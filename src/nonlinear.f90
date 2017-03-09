module nonlinear_module
  
  use bl_types
  use multifab_module
  use define_bc_module
  use ml_layout_module

  implicit none

  private

  public :: update_nonlinear

contains

  subroutine update_nonlinear(mla, nonlinear_term, u, viscosity, dx, the_bc_level)

    use bl_constants_module
    use ml_restrict_fill_module
    use probin_module, only: extrap_comp

    type(ml_layout)   , intent(in   ) :: mla
    type(multifab)    , intent(inout) :: nonlinear_term(:)
    type(multifab)    , intent(in   ) :: u(:)
    type(multifab)    , intent(in   ) :: viscosity(:)
    real(kind = dp_t) , intent(in   ) :: dx(:,:)
    type(bc_level)    , intent(in   ) :: the_bc_level(:)

    ! local
    real(kind=dp_t), pointer :: np(:,:,:,:)
    real(kind=dp_t), pointer :: up(:,:,:,:)
    real(kind=dp_t), pointer :: vp(:,:,:,:)

    integer :: lo(mla%dim), hi(mla%dim)
    integer :: i, dm, n, nlevs
    integer :: ng_u, ng_v
    
    type(bl_prof_timer), save :: bpt

    call build(bpt,"nonlinear")

    nlevs = mla%nlevel
    dm    = mla%dim

    ng_u = u(1)%ng
    ng_v = viscosity(1)%ng

    do n=1, nlevs

       do i = 1, nfabs(nonlinear_term(n))
          np => dataptr(nonlinear_term(n),i)
          up => dataptr(             u(n),i)
          vp => dataptr(     viscosity(n),i)
          lo = lwb(get_box(nonlinear_term(n),i))
          hi = upb(get_box(nonlinear_term(n),i))
          select case (dm)
          case (2)
             call update_nonlinear_2d(np(:,:,1,:), up(:,:,1,:), vp(:,:,1,1), &
                                      lo, hi, ng_u, ng_v, dx(n,:))
          case (3)
             !call update_nonlinear_3d(nlp(:,:,:,:), up(:,:,:,:), viscp(:,:,:,:), &
             !                         lo, hi, ng_nl, ng_u, ng_v, dx(n,:))
          end select
       end do

    enddo ! end loop over levels

    ! restrict cell-centered multifab data, fill all boundaries
    call ml_restrict_and_fill(nlevs, nonlinear_term, mla%mba%rr, the_bc_level, bcomp=extrap_comp, &
                              same_boundary=.true.)

    call destroy(bpt)

  end subroutine update_nonlinear

  subroutine update_nonlinear_2d(nonlinear_term, u, viscosity, lo, hi, ng_u, ng_v, dx)

    use bl_constants_module

    integer           , intent(in   ) :: lo(:), hi(:), ng_u, ng_v
    real (kind = dp_t), intent(inout) :: nonlinear_term(lo(1)   :, lo(2)   :,:)  
    real (kind = dp_t), intent(in   ) ::              u(lo(1)-ng_u:, lo(2)-ng_u:,:)  
    real (kind = dp_t), intent(in   ) ::      viscosity(lo(1)-ng_v:, lo(2)-ng_v:)  
    real (kind = dp_t), intent(in   ) :: dx(:)

    integer :: i, j
    real (kind = dp_t) ux, uy, vx, vy, viscx, viscy

    do j = lo(2)+1, hi(2)-1
       do i = lo(1), hi(1)
          ux = (u(i+1,j,1) - u(i-1,j,1)) / (TWO * dx(1))
          vx = (u(i+1,j,2) - u(i-1,j,2)) / (TWO * dx(1))
          uy = (u(i,j+1,1) - u(i,j-1,1)) / (TWO * dx(2))
          vy = (u(i,j+1,2) - u(i,j-1,2)) / (TWO * dx(2))
          viscx = (viscosity(i+1,j) - viscosity(i-1,j)) / (TWO * dx(1))
          viscy = (viscosity(i,j+1) - viscosity(i,j-1)) / (TWO * dx(2))
          nonlinear_term(i,j,1) = TWO * viscx * ux + viscy * (uy + vx)
          nonlinear_term(i,j,2) = TWO * viscy * vy + viscx * (uy + vx)
       enddo
    enddo

    ! TODO: Also for boundaries in x-direction
    j = lo(2)
    do i = lo(1), hi(1)
       ux = (u(i+1,j,1) - u(i-1,j,1)) / (TWO * dx(1))
       vx = (u(i+1,j,2) - u(i-1,j,2)) / (TWO * dx(1))
       uy = -(THREE * u(i,j,1) - FOUR * u(i,j+1,1) + u(i,j+2,1)) / (TWO * dx(2))
       vy = -(THREE * u(i,j,2) - FOUR * u(i,j+1,2) + u(i,j+2,2)) / (TWO * dx(2))
       viscx = (viscosity(i+1,j) - viscosity(i-1,j)) / (TWO * dx(1))
       viscy = -(THREE * viscosity(i,j) - FOUR * viscosity(i,j+1) + viscosity(i,j+2)) / (TWO * dx(2))
       nonlinear_term(i,j,1) = TWO * viscx * ux + viscy * (uy + vx)
       nonlinear_term(i,j,2) = TWO * viscy * vy + viscx * (uy + vx)
    enddo

    j = hi(2)
    do i = lo(1), hi(1)
       ux = (u(i+1,j,1) - u(i-1,j,1)) / (TWO * dx(1))
       vx = (u(i+1,j,2) - u(i-1,j,2)) / (TWO * dx(1))
       uy = (THREE * u(i,j,1) - FOUR * u(i,j-1,1) + u(i,j-2,1)) / (TWO * dx(2))
       vy = (THREE * u(i,j,2) - FOUR * u(i,j-1,2) + u(i,j-2,2)) / (TWO * dx(2))
       viscx = (viscosity(i+1,j) - viscosity(i-1,j)) / (TWO * dx(1))
       viscy = (THREE * viscosity(i,j) - FOUR * viscosity(i,j-1) + viscosity(i,j-2)) / (TWO * dx(2))
       nonlinear_term(i,j,1) = TWO * viscx * ux + viscy * (uy + vx)
       nonlinear_term(i,j,2) = TWO * viscy * vy + viscx * (uy + vx)
    enddo

  end subroutine update_nonlinear_2d

  subroutine update_nonlinear_3d(nonlin, u, visc, lo, hi, ng_nl, ng_u, ng_v, dx)

    use bl_constants_module

    integer           , intent(in   ) :: lo(:), hi(:), ng_nl, ng_v, ng_u
    real (kind = dp_t), intent(inout) :: nonlin(lo(1)-ng_nl:,lo(2)-ng_nl:,lo(3)-ng_nl:,:)  
    real (kind = dp_t), intent(in   ) ::      u(lo(1)-ng_u: ,lo(2)-ng_u: ,lo(3)-ng_u: ,:)  
    real (kind = dp_t), intent(in   ) ::   visc(lo(1)-ng_v: ,lo(2)-ng_v: ,lo(3)-ng_v: ,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)

    integer :: i, j, k
    real (kind = dp_t) ux, uy, uz, vx, vy, vz, wx, wy, wz, viscx, viscy, viscz

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
             viscx = (visc(i+1,j,k,1) - visc(i-1,j,k,1)) / (TWO * dx(1))
             viscy = (visc(i,j+1,k,1) - visc(i,j-1,k,1)) / (TWO * dx(2))
             viscz = (visc(i,j,k+1,1) - visc(i,j,k-1,1)) / (TWO * dx(3))
             nonlin(i,j,k,1) = TWO * viscx * ux + viscy * (uy + vx) + viscz * (uz + wx)
             nonlin(i,j,k,2) = TWO * viscy * vy + viscz * (vz + wy) + viscx * (vx + uy)
             nonlin(i,j,k,3) = TWO * viscz * wz + viscx * (wx + uz) + viscy * (wy + vz)
          enddo
       enddo
    enddo

  end subroutine update_nonlinear_3d

end module nonlinear_module
