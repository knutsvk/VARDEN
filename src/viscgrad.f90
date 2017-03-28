module viscgrad_module
  
  use bl_types
  use multifab_module
  use define_bc_module
  use ml_layout_module

  implicit none

  private

  public :: update_viscgrad

contains

   subroutine update_viscgrad(mla, visc_grad_term, viscosity, u, dx, the_bc_level)

    use bl_constants_module
    use ml_restrict_fill_module
    use probin_module, only: extrap_comp, prob_lo, prob_hi

    type(ml_layout)   , intent(in   ) :: mla
    type(multifab)    , intent(inout) :: visc_grad_term(:)
    type(multifab)    , intent(in   ) ::      viscosity(:)
    type(multifab)    , intent(in   ) ::              u(:)
    real(kind = dp_t) , intent(in   ) :: dx(:,:)
    type(bc_level)    , intent(in   ) :: the_bc_level(:)

    ! local
    real(kind=dp_t), pointer :: vgp(:,:,:,:)
    real(kind=dp_t), pointer ::  vp(:,:,:,:)
    real(kind=dp_t), pointer ::  up(:,:,:,:)

    integer :: lo(mla%dim), hi(mla%dim)
    integer :: i, dm, n, nlevs
    integer :: ng_u, ng_v
    
    type(bl_prof_timer), save :: bpt

    call build(bpt,"viscgrad")

    nlevs = mla%nlevel
    dm    = mla%dim

    ng_u = u(1)%ng
    ng_v = viscosity(1)%ng

    do n=1, nlevs

       do i = 1, nfabs(visc_grad_term(n))
          vgp => dataptr(visc_grad_term(n),i)
          vp =>  dataptr(     viscosity(n),i)
          up =>  dataptr(             u(n),i)
          lo = lwb(get_box(visc_grad_term(n),i))
          hi = upb(get_box(visc_grad_term(n),i))
          select case (dm)
          case (2)
             call update_viscgrad_2d(vgp(:,:,1,:), vp(:,:,1,1), up(:,:,1,:), lo, hi, &
                                       ng_u, ng_v, dx(n,:), &
                                       the_bc_level(n)%phys_bc_level_array(i,:,:))
          case (3) !TODO
             !call update_strainrate_3d(dgp(:,:,:,:), up(:,:,:,:), lo, hi, ng_dg, ng_u, dx(n,:), &
             !                          the_bc_level(n)%phys_bc_level_array(i,:,:))
          end select
       end do

    enddo ! end loop over levels

    ! restrict cell-centered multifab data, fill all boundaries
    ! TODO: Maybe need to change since visc_grad_term has dm comps? 
    call ml_restrict_and_fill(nlevs, visc_grad_term, mla%mba%rr, the_bc_level, bcomp=extrap_comp, &
                              dx_in=dx(1,:), prob_lo_in=prob_lo, prob_hi_in=prob_hi)

    call destroy(bpt)

  end subroutine update_viscgrad

  subroutine update_viscgrad_2d(visc_grad_term, viscosity, u, lo, hi, ng_u, ng_v, dx, bc)

    use bl_constants_module
    use bc_module
    use probin_module, only: u_bcx_lo, v_bcx_lo, u_bcx_hi, v_bcx_hi, &
                             u_bcy_lo, v_bcy_lo, u_bcy_hi, v_bcy_hi

    integer           , intent(in   ) :: lo(:), hi(:), ng_u, ng_v
    real (kind = dp_t), intent(inout) :: visc_grad_term(lo(1):     , lo(2):     ,:)  
    real (kind = dp_t), intent(in   ) ::      viscosity(lo(1)-ng_v:, lo(2)-ng_v:)  
    real (kind = dp_t), intent(in   ) ::              u(lo(1)-ng_u:, lo(2)-ng_u:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    integer           , intent(in   ) :: bc(:,:)

    integer :: i, j
    real (kind = dp_t) ux, uy, vx, vy, viscx, viscy

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          ux = (u(i+1,j,1) - u(i-1,j,1)) / (TWO * dx(1))
          vx = (u(i+1,j,2) - u(i-1,j,2)) / (TWO * dx(1))
          uy = (u(i,j+1,1) - u(i,j-1,1)) / (TWO * dx(2))
          vy = (u(i,j+1,2) - u(i,j-1,2)) / (TWO * dx(2))
          viscx = (viscosity(i+1,j) - viscosity(i-1,j)) / (TWO * dx(1))
          viscy = (viscosity(i,j+1) - viscosity(i,j-1)) / (TWO * dx(2))
          visc_grad_term(i,j,1) = TWO * viscx * ux + viscy * (uy + vx)
          visc_grad_term(i,j,2) = TWO * viscy * vy + viscx * (uy + vx)
       enddo
    enddo

    if ( bc(1,1) == NO_SLIP_WALL ) then
       i = lo(1)
       do j = lo(2), hi(2)
          ux = TWO * (u(i,j,1) - u_bcx_lo) / dx(1)
          vx = TWO * (u(i,j,2) - v_bcx_lo) / dx(1)
          uy = (u(i,j+1,1) - u(i,j-1,1)) / (TWO * dx(2))
          vy = (u(i,j+1,2) - u(i,j-1,2)) / (TWO * dx(2))
          viscx = (viscosity(i+1,j) - viscosity(i-1,j)) / (TWO * dx(1))
          viscy = (viscosity(i,j+1) - viscosity(i,j-1)) / (TWO * dx(2))
          visc_grad_term(i,j,1) = TWO * viscx * ux + viscy * (uy + vx)
          visc_grad_term(i,j,2) = TWO * viscy * vy + viscx * (uy + vx)
       enddo
    endif

    if ( bc(1,2) == NO_SLIP_WALL ) then
       i = hi(1)
       do j = lo(2), hi(2)
          ux = -TWO * (u(i,j,1) - u_bcx_hi) / dx(1)
          vx = -TWO * (u(i,j,2) - v_bcx_hi) / dx(1)
          uy = (u(i,j+1,1) - u(i,j-1,1)) / (TWO * dx(2))
          vy = (u(i,j+1,2) - u(i,j-1,2)) / (TWO * dx(2))
          viscx = (viscosity(i+1,j) - viscosity(i-1,j)) / (TWO * dx(1))
          viscy = (viscosity(i,j+1) - viscosity(i,j-1)) / (TWO * dx(2))
          visc_grad_term(i,j,1) = TWO * viscx * ux + viscy * (uy + vx)
          visc_grad_term(i,j,2) = TWO * viscy * vy + viscx * (uy + vx)
       enddo
    endif

    if ( bc(2,1) == NO_SLIP_WALL ) then
       j = lo(2)
       do i = lo(1), hi(1)
          ux = (u(i+1,j,1) - u(i-1,j,1)) / (TWO * dx(1))
          vx = (u(i+1,j,2) - u(i-1,j,2)) / (TWO * dx(1))
          uy = TWO * (u(i,j,1) - u_bcy_lo) / dx(2)
          vy = TWO * (u(i,j,2) - v_bcy_lo) / dx(2)
          viscx = (viscosity(i+1,j) - viscosity(i-1,j)) / (TWO * dx(1))
          viscy = (viscosity(i,j+1) - viscosity(i,j-1)) / (TWO * dx(2))
          visc_grad_term(i,j,1) = TWO * viscx * ux + viscy * (uy + vx)
          visc_grad_term(i,j,2) = TWO * viscy * vy + viscx * (uy + vx)
       enddo
    endif

    if ( bc(2,2) == NO_SLIP_WALL ) then
       j = hi(2)
       do i = lo(1), hi(1)
          ux = (u(i+1,j,1) - u(i-1,j,1)) / (TWO * dx(1))
          vx = (u(i+1,j,2) - u(i-1,j,2)) / (TWO * dx(1))
          uy = -TWO * (u(i,j,1) - u_bcy_hi) / dx(2)
          vy = -TWO * (u(i,j,2) - v_bcy_hi) / dx(2)
          viscx = (viscosity(i+1,j) - viscosity(i-1,j)) / (TWO * dx(1))
          viscy = (viscosity(i,j+1) - viscosity(i,j-1)) / (TWO * dx(2))
          visc_grad_term(i,j,1) = TWO * viscx * ux + viscy * (uy + vx)
          visc_grad_term(i,j,2) = TWO * viscy * vy + viscx * (uy + vx)
       enddo
    endif

  end subroutine update_viscgrad_2d

  !TODO: 3D

end module viscgrad_module
