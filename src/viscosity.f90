module viscosity_module
  
  use bl_types
  use multifab_module
  use define_bc_module
  use ml_layout_module

  implicit none

  private

  public :: update_viscosity

contains

  subroutine update_viscosity(mla, visc, dotgam, dx, the_bc_level)

    use bl_constants_module
    use ml_restrict_fill_module
    use probin_module, only: extrap_comp

    type(ml_layout)   , intent(in   ) :: mla
    type(multifab)    , intent(inout) :: visc(:)
    type(multifab)    , intent(in   ) :: dotgam(:)
    real(kind = dp_t) , intent(in   ) :: dx(:,:)
    type(bc_level)    , intent(in   ) :: the_bc_level(:)

    ! local
    real(kind=dp_t), pointer :: viscp(:,:,:,:)
    real(kind=dp_t), pointer ::   dgp(:,:,:,:)

    integer :: lo(get_dim(visc(1))), hi(get_dim(visc(1)))
    integer :: i, dm, nscal, n, nlevs
    integer :: ng_visc, ng_dg
    
    type(bl_prof_timer), save :: bpt

    call build(bpt,"viscosity")

    nlevs = mla%nlevel
    dm    = mla%dim

    ng_visc = visc(1)%ng
    ng_dg = dotgam(1)%ng

    do n=1, nlevs

       do i = 1, nfabs(visc(n))
          viscp => dataptr(  visc(n),i)
          dgp   => dataptr(dotgam(n),i)
          lo = lwb(get_box(visc(n),i))
          hi = upb(get_box(visc(n),i))
          select case (dm)
          case (2)
             call update_viscosity_2d(viscp(:,:,1,:), dgp(:,:,1,:), lo, hi, ng_visc, ng_dg, dx(n,:))
          case (3)
             call update_viscosity_3d(viscp(:,:,:,:), dgp(:,:,:,:), lo, hi, ng_visc, ng_dg, dx(n,:))
          end select
       end do

    enddo ! end loop over levels

    ! restrict cell-centered multifab data, fill all boundaries
    call ml_restrict_and_fill(nlevs, visc, mla%mba%rr, the_bc_level, bcomp=extrap_comp)

    call destroy(bpt)

  end subroutine update_viscosity

  subroutine update_viscosity_2d(visc, dotgam, lo, hi, ng_v, ng_dg, dx)

    use bl_constants_module
    use probin_module, only: visc_coef, yield_stress, papa_reg

    integer           , intent(in   ) :: lo(:), hi(:), ng_v, ng_dg
    real (kind = dp_t), intent(inout) ::    visc(lo(1)-ng_v: ,lo(2)-ng_v:,:)  
    real (kind = dp_t), intent(in   ) ::  dotgam(lo(1)-ng_dg:,lo(2)-ng_dg:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)

    integer :: i, j
    real (kind = dp_t) frac

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          frac = dotgam(i,j,1) / papa_reg
          if (frac < 1e-9) then
             ! Taylor expand exp() to avoid dividing by small number
             visc(i,j,1) = visc_coef + yield_stress / papa_reg * &
                (ONE - HALF * frac + SIXTH * frac**2 - SIXTH * FOURTH * frac**3)
          else 
             visc(i,j,1) = visc_coef + yield_stress / papa_reg * &
                (ONE - exp(-frac)) / frac
          endif
       enddo
    enddo

  end subroutine update_viscosity_2d

  ! TODO: 
  subroutine update_viscosity_3d(visc, u, lo, hi, ng_v, ng_u, dx)

    use bl_constants_module
    use probin_module, only: visc_coef, yield_stress, papa_reg

    integer           , intent(in   ) :: lo(:), hi(:), ng_v, ng_u
    real (kind = dp_t), intent(inout) ::    visc(lo(1)-ng_v:,lo(2)-ng_v:,lo(3)-ng_v:,:)  
    real (kind = dp_t), intent(in   ) ::       u(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)

    integer :: i, j, k
    real (kind = dp_t) ux, uy, uz, vx, vy, vz, wx, wy, wz, dotgamma, frac

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             ux = (u(i+1,j,k,1) - u(i-1,j,k,1)) / (TWO * dx(1))
             uy = (u(i,j+1,k,1) - u(i,j-1,k,1)) / (TWO * dx(2))
             uz = (u(i,j,k+1,1) - u(i,j,k-1,1)) / (TWO * dx(3))
             vx = (u(i+1,j,k,2) - u(i-1,j,k,2)) / (TWO * dx(1))
             vy = (u(i,j+1,k,2) - u(i,j-1,k,2)) / (TWO * dx(2))
             vz = (u(i,j,k+1,2) - u(i,j,k-1,2)) / (TWO * dx(3))
             wx = (u(i+1,j,k,3) - u(i-1,j,k,3)) / (TWO * dx(1))
             wy = (u(i,j+1,k,3) - u(i,j-1,k,3)) / (TWO * dx(2))
             wz = (u(i,j,k+1,3) - u(i,j,k-1,3)) / (TWO * dx(3))
             dotgamma = sqrt(TWO * (ux**2 + vy**2 + wz**2) + (uy + vx)**2 + (vz + wy)**2 + (wx + uz)**2)
             frac = dotgamma / papa_reg
             if (frac < 1e-5) then
                visc(i,j,k,1) = visc_coef + yield_stress / papa_reg * (ONE - HALF * frac + SIXTH * frac**2)
             else 
                visc(i,j,k,1) = visc_coef + yield_stress / dotgamma * (ONE - exp(-frac))
             endif
          enddo
       enddo
    enddo

  end subroutine update_viscosity_3d

end module viscosity_module
