module mac_multigrid_module

  use bl_types
  use ml_layout_module
  use define_bc_module
  use multifab_module
  use bndry_reg_module
  use bl_constants_module
  use bl_error_module
  use sparse_solve_module

  implicit none

  private

  public :: mac_multigrid

contains

  subroutine mac_multigrid(mla,rh,phi,fine_flx,alpha,beta,dx,&
       the_bc_tower,bc_comp,stencil_order,ref_ratio,umac_norm)

     use stencil_fill_module, only: stencil_fill_cc_all_mglevels
     use mg_module, only: mg_tower, mg_tower_build, mg_tower_destroy
     use ml_solve_module
     use probin_module, only: mg_verbose, cg_verbose

    type(ml_layout),intent(in   ) :: mla
    integer        ,intent(in   ) :: stencil_order
    integer        ,intent(in   ) :: ref_ratio(:,:)

    real(dp_t), intent(in) :: dx(:,:)
    type(bc_tower), intent(in) :: the_bc_tower
    integer     ,intent(in   ) :: bc_comp
    real(dp_t), intent(in), optional :: umac_norm(:)

    type(layout  ) :: la
    type(box     ) :: pd

    type(multifab), allocatable :: coeffs(:)

    type(multifab) , intent(in   ) :: alpha(:), beta(:,:)
    type(multifab) , intent(inout) ::    rh(:),  phi(:)
    type(bndry_reg), intent(inout) :: fine_flx(2:)

    type( multifab) :: ss
    type(imultifab) :: mm
    type(sparse)    :: sparse_object
    type(mg_tower)  :: mgt(mla%nlevel)
    integer         :: dm, ns, nlevs, test

    ! MG solver defaults
    integer :: bottom_solver, bottom_max_iter
    integer    :: max_iter
    integer    :: min_width
    integer    :: max_nlevel
    integer    :: n, nu1, nu2, gamma, cycle_type, smoother
    integer    :: max_nlevel_in,do_diagnostics
    real(dp_t) :: rel_eps,abs_eps,omega,bottom_solver_eps
    real(dp_t) ::  xa(mla%dim),  xb(mla%dim)
    real(dp_t) :: pxa(mla%dim), pxb(mla%dim)

    !! Defaults:

    nlevs = mla%nlevel
    dm    = mla%dim

    test           = 0

    max_nlevel        = mgt(nlevs)%max_nlevel
    max_iter          = mgt(nlevs)%max_iter
    rel_eps           = mgt(nlevs)%eps
    abs_eps           = mgt(nlevs)%abs_eps
    smoother          = mgt(nlevs)%smoother
    nu1               = mgt(nlevs)%nu1
    nu2               = mgt(nlevs)%nu2
    gamma             = mgt(nlevs)%gamma
    omega             = mgt(nlevs)%omega
    cycle_type        = mgt(nlevs)%cycle_type
    bottom_solver     = mgt(nlevs)%bottom_solver
    bottom_solver_eps = mgt(nlevs)%bottom_solver_eps
    bottom_max_iter   = mgt(nlevs)%bottom_max_iter
    min_width         = mgt(nlevs)%min_width

    ! Note: put this here to minimize asymmetries - ASA
    if (nlevs .eq. 1) then
       rel_eps = 1.d-12
    else if (nlevs .eq. 2) then
       rel_eps = 1.d-11
    else
       rel_eps = 1.d-10
    endif

    abs_eps = -1.0_dp_t
    if (present(umac_norm)) then
       do n = 1,nlevs
          abs_eps = max(abs_eps, umac_norm(n) / dx(n,1))
       end do
       abs_eps = rel_eps * abs_eps
    end if

    bottom_solver = 1
    bottom_solver_eps = 1.d-3

    if ( test /= 0 .AND. max_iter == mgt(nlevs)%max_iter ) &
         max_iter = 1000

    ns = 1 + dm*3

    do n = nlevs, 1, -1

       if (n == 1) then
          max_nlevel_in = max_nlevel
       else
          if ( all(ref_ratio(n-1,:) == 2) ) then
             max_nlevel_in = 1
          else if ( all(ref_ratio(n-1,:) == 4) ) then
             max_nlevel_in = 2
          else
             call bl_error("MAC_MULTIGRID: confused about ref_ratio")
          end if
       end if

       pd = layout_get_pd(mla%la(n))

       call mg_tower_build(mgt(n), mla%la(n), pd, &
                           the_bc_tower%bc_tower_array(n)%ell_bc_level_array(0,:,:,bc_comp),&
                           dh = dx(n,:), &
                           ns = ns, &
                           smoother = smoother, &
                           nu1 = nu1, &
                           nu2 = nu2, &
                           gamma = gamma, &
                           cycle_type = cycle_type, &
                           omega = omega, &
                           bottom_solver = bottom_solver, &
                           bottom_max_iter = bottom_max_iter, &
                           bottom_solver_eps = bottom_solver_eps, &
                           max_iter = max_iter, &
                           max_nlevel = max_nlevel_in, &
                           min_width = min_width, &
                           eps = rel_eps, &
                           abs_eps = abs_eps, &
                           verbose = mg_verbose, &
                           cg_verbose = cg_verbose, &
                           nodal = rh(nlevs)%nodal)

    end do

    !! Fill coefficient array

    do n = nlevs,1,-1

       allocate(coeffs(mgt(n)%nlevels))

       la = mla%la(n)
       pd = layout_get_pd(la)

       call multifab_build(coeffs(mgt(n)%nlevels), la, 1+dm, 1)
       call multifab_copy_c(coeffs(mgt(n)%nlevels),1,alpha(n),1, 1,ng=alpha(n)%ng)
       call multifab_copy_c(coeffs(mgt(n)%nlevels),2, beta(n,1),1,1,ng=beta(n,1)%ng)
       call multifab_copy_c(coeffs(mgt(n)%nlevels),3, beta(n,2),1,1,ng=beta(n,2)%ng)
       if (dm > 2) &
          call multifab_copy_c(coeffs(mgt(n)%nlevels),4, beta(n,3),1,1,ng=beta(n,3)%ng)

       if (n > 1) then
          xa = HALF*ref_ratio(n-1,:)*mgt(n)%dh(:,mgt(n)%nlevels)
          xb = HALF*ref_ratio(n-1,:)*mgt(n)%dh(:,mgt(n)%nlevels)
       else
          xa = ZERO
          xb = ZERO
       end if

       pxa = ZERO
       pxb = ZERO

       call stencil_fill_cc_all_mglevels(mgt(n), coeffs, xa, xb, pxa, pxb, pd, stencil_order, &
                                         the_bc_tower%bc_tower_array(n)%ell_bc_level_array(0,:,:,bc_comp))
 
       call destroy(coeffs(mgt(n)%nlevels))
       deallocate(coeffs)

    end do

    if (mg_verbose >= 3) then
       do_diagnostics = 1
    else
       do_diagnostics = 0
    end if

    call ml_cc_solve(mla, mgt, rh, phi, fine_flx, ref_ratio,do_diagnostics)

    do n = 1,nlevs
       call multifab_fill_boundary(phi(n))
    end do

    if ( test == 3 ) then
       call sparse_destroy(sparse_object)
    end if
    if ( test > 0 ) then
       call destroy(ss)
       call destroy(mm)
    end if

    do n = 1, nlevs
       call mg_tower_destroy(mgt(n))
    end do

  end subroutine mac_multigrid

end module mac_multigrid_module
