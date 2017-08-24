module viscosity_module
  
  use bl_types
  use multifab_module
  use define_bc_module
  use ml_layout_module

  implicit none

  private

  public :: update_viscosity

contains

  subroutine update_viscosity(mla, viscosity, phi, dx, the_bc_level)

    use bl_constants_module
    use ml_restrict_fill_module
    use probin_module, only: extrap_comp, prob_lo, prob_hi

    type(ml_layout)   , intent(in   ) :: mla
    type(multifab)    , intent(inout) :: viscosity(:)
    type(multifab)    , intent(in   ) :: phi(:)
    real(kind = dp_t) , intent(in   ) :: dx(:,:)
    type(bc_level)    , intent(in   ) :: the_bc_level(:)

    ! local
    real(kind=dp_t), pointer :: vp(:,:,:,:)
    real(kind=dp_t), pointer :: pp(:,:,:,:)

    integer :: lo(mla%dim), hi(mla%dim)
    integer :: i, dm, n, nlevs, ng
    
    type(bl_prof_timer), save :: bpt

    call build(bpt,"viscosity")

    nlevs = mla%nlevel
    dm    = mla%dim
    ng = nghost(viscosity(1))

    do n=1, nlevs

       do i = 1, nfabs(viscosity(n))
          vp => dataptr(  viscosity(n),i)
          pp => dataptr(        phi(n),i)
          lo = lwb(get_box(viscosity(n),i))
          hi = upb(get_box(viscosity(n),i))
          select case (dm)
          case (2)
             call update_viscosity_2d(vp(:,:,1,1), pp(:,:,1,1), lo, hi, ng, dx(n,:))
          case (3)
          !   call update_viscosity_3d(viscp(:,:,:,:), dgp(:,:,:,:), lo, hi, ng_visc, ng_dg, dx(n,:))
          end select
       end do

    enddo ! end loop over levels

    ! restrict cell-centered multifab data, fill all boundaries
    call ml_restrict_and_fill(nlevs, viscosity, mla%mba%rr, the_bc_level, bcomp=extrap_comp, &
                              dx_in=dx(1,:), prob_lo_in=prob_lo, prob_hi_in=prob_hi)

    call destroy(bpt)

  end subroutine update_viscosity

  subroutine update_viscosity_2d(viscosity, phi, lo, hi, ng, dx)

    use bl_constants_module
    use probin_module, only: visc_coef

    integer           , intent(in   ) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(inout) ::   viscosity(lo(1)-ng:,lo(2)-ng:)  
    real (kind = dp_t), intent(in   ) ::         phi(lo(1)   :,lo(2)   :)  
    real (kind = dp_t), intent(in   ) :: dx(:)

    integer :: i, j
    real (kind = dp_t) frac

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
       ! TODO: make level-set dependent
          viscosity(i,j) = visc_coef
       enddo
    enddo

  end subroutine update_viscosity_2d

  ! TODO: 
  subroutine update_viscosity_3d(visc, phi, lo, hi, ng, dx)

    use bl_constants_module
    use probin_module, only: visc_coef

    integer           , intent(in   ) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(inout) ::    visc(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)  
    real (kind = dp_t), intent(in   ) ::     phi(lo(1)   :,lo(2)   :,lo(3)   :)  
    real (kind = dp_t), intent(in   ) :: dx(:)

    integer :: i, j, k

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
        ! TODO: make level-set dependent
            visc(i,j,k) = visc_coef
          enddo
       enddo
    enddo

  end subroutine update_viscosity_3d

end module viscosity_module
