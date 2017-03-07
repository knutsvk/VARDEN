module nonlinear_module
  
  use bl_types
  use multifab_module
  use define_bc_module
  use ml_layout_module

  implicit none

  private

  public :: update_nonlinear

contains

  subroutine update_nonlinear(mla, nonlin, u, visc, dx, the_bc_level)

    use bl_constants_module
    use ml_restrict_fill_module
    use probin_module, only: extrap_comp

    type(ml_layout)   , intent(in   ) :: mla
    type(multifab)    , intent(inout) :: nonlin(:)
    type(multifab)    , intent(in   ) :: u(:)
    type(multifab)    , intent(in   ) :: visc(:)
    real(kind = dp_t) , intent(in   ) :: dx(:,:)
    type(bc_level)    , intent(in   ) :: the_bc_level(:)

    ! local
    real(kind=dp_t), pointer ::   nlp(:,:,:,:)
    real(kind=dp_t), pointer ::    up(:,:,:,:)
    real(kind=dp_t), pointer :: viscp(:,:,:,:)

    integer :: lo(get_dim(nonlin(1))), hi(get_dim(nonlin(1)))
    integer :: i, dm, nscal, n, nlevs
    integer :: ng_nl, ng_u, ng_v
    
    type(bl_prof_timer), save :: bpt

    call build(bpt,"nonlinear")

    nlevs = mla%nlevel
    dm    = mla%dim

    ng_nl = nonlin(1)%ng
    ng_u  =      u(1)%ng
    ng_v  =   visc(1)%ng

    do n=1, nlevs

       do i = 1, nfabs(visc(n))
            nlp => dataptr(nonlin(n),i)
          viscp => dataptr(  visc(n),i)
             up => dataptr(     u(n),i)
          lo = lwb(get_box(nonlin(n),i))
          hi = upb(get_box(nonlin(n),i))
          select case (dm)
          case (2)
             call update_nonlinear_2d(nlp(:,:,1,:), up(:,:,1,:), viscp(:,:,1,:), &
                                      lo, hi, ng_nl, ng_u, ng_v, dx(n,:))
          case (3)
             call update_nonlinear_3d(nlp(:,:,:,:), up(:,:,:,:), viscp(:,:,:,:), &
                                      lo, hi, ng_nl, ng_u, ng_v, dx(n,:))
          end select
       end do

    enddo ! end loop over levels

    ! restrict cell-centered multifab data, fill all boundaries
    call ml_restrict_and_fill(nlevs, nonlin, mla%mba%rr, the_bc_level, &
                              bcomp=extrap_comp, same_boundary=.true.)

    call destroy(bpt)

  end subroutine update_nonlinear

  subroutine update_nonlinear_2d(nonlin, u, visc, lo, hi, ng_nl, ng_u, ng_v, dx)

    use bl_constants_module

    integer           , intent(in   ) :: lo(:), hi(:), ng_nl, ng_v, ng_u
    real (kind = dp_t), intent(inout) :: nonlin(lo(1)-ng_nl:, lo(2)-ng_nl:,:)  
    real (kind = dp_t), intent(in   ) ::      u(lo(1)-ng_u: , lo(2)-ng_u: ,:)  
    real (kind = dp_t), intent(in   ) ::   visc(lo(1)-ng_v: , lo(2)-ng_v: ,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)

    integer :: i, j
    real (kind = dp_t) ux, uy, vx, vy, viscx, viscy

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          ux = (u(i+1,j,1) - u(i-1,j,1)) / (TWO * dx(1))
          uy = (u(i,j+1,1) - u(i,j-1,1)) / (TWO * dx(2))
          vx = (u(i+1,j,2) - u(i-1,j,2)) / (TWO * dx(1))
          vy = (u(i,j+1,2) - u(i,j-1,2)) / (TWO * dx(2))
          viscx = (visc(i+1,j,1) - visc(i-1,j,1)) / (TWO * dx(1))
          viscy = (visc(i,j+1,1) - visc(i,j-1,1)) / (TWO * dx(2))
          nonlin(i,j,1) = TWO * viscx * ux + viscy * (uy + vx)
          nonlin(i,j,2) = TWO * viscy * vy + viscx * (uy + vx)
       enddo
    enddo

    ! FILL GHOST CELLS AS WELL

    ! Left
    i = lo(1)
    do j = lo(2), hi(2)
       ux = -(THREE * u(i,j,1) - FOUR * u(i+1,j,1) + u(i+2,j,1)) / (TWO * dx(1))
       uy = (u(i,j+1,1) - u(i,j-1,1)) / (TWO * dx(2))
       vx = -(THREE * u(i,j,2) - FOUR * u(i+1,j,2) + u(i+2,j,2)) / (TWO * dx(1))
       vy = (u(i,j+1,2) - u(i,j-1,2)) / (TWO * dx(2))
       viscx = -(THREE * visc(i,j,1) - FOUR * visc(i+1,j,1) + visc(i+2,j,1)) / (TWO * dx(1))
       viscy = (visc(i,j+1,1) - visc(i,j-1,1)) / (TWO * dx(2))
       nonlin(i,j,1) = TWO * viscx * ux + viscy * (uy + vx)
       nonlin(i,j,2) = TWO * viscy * vy + viscx * (uy + vx)
    enddo

    ! Right
    i = hi(1)
    do j = lo(2), hi(2)
       ux = (THREE * u(i,j,1) - FOUR * u(i-1,j,1) + u(i-2,j,1)) / (TWO * dx(1))
       uy = (u(i,j+1,1) - u(i,j-1,1)) / (TWO * dx(2))
       vx = (THREE * u(i,j,2) - FOUR * u(i-1,j,2) + u(i-2,j,2)) / (TWO * dx(1))
       vy = (u(i,j+1,2) - u(i,j-1,2)) / (TWO * dx(2))
       viscx = (THREE * visc(i,j,1) - FOUR * visc(i-1,j,1) + visc(i-2,j,1)) / (TWO * dx(1))
       viscy = (visc(i,j+1,1) - visc(i,j-1,1)) / (TWO * dx(2))
       nonlin(i,j,1) = TWO * viscx * ux + viscy * (uy + vx)
       nonlin(i,j,2) = TWO * viscy * vy + viscx * (uy + vx)
    enddo

    ! Bottom
    j = lo(2)
    do i = lo(1), hi(1)
       ux = (u(i+1,j,1) - u(i-1,j,1)) / (TWO * dx(1))
       uy = -(THREE * u(i,j,1) - FOUR * u(i,j+1,1) + u(i,j+2,1)) / (TWO * dx(2))
       vx = (u(i+1,j,2) - u(i-1,j,2)) / (TWO * dx(1))
       vy = -(THREE * u(i,j,2) - FOUR * u(i,j+1,2) + u(i,j+2,2)) / (TWO * dx(2))
       viscx = (visc(i+1,j,1) - visc(i-1,j,1)) / (TWO * dx(1))
       viscy = -(THREE * visc(i,j,1) - FOUR * visc(i,j+1,1) + visc(i,j+2,1)) / (TWO * dx(2))
       nonlin(i,j,1) = TWO * viscx * ux + viscy * (uy + vx)
       nonlin(i,j,2) = TWO * viscy * vy + viscx * (uy + vx)
    enddo

    ! Top
    j = hi(2)
    do i = lo(1), hi(1)
       ux = (u(i+1,j,1) - u(i-1,j,1)) / (TWO * dx(1))
       uy = (THREE * u(i,j,1) - FOUR * u(i,j-1,1) + u(i,j-2,1)) / (TWO * dx(2))
       vx = (u(i+1,j,2) - u(i-1,j,2)) / (TWO * dx(1))
       vy = (THREE * u(i,j,2) - FOUR * u(i,j-1,2) + u(i,j-2,2)) / (TWO * dx(2))
       viscx = (visc(i+1,j,1) - visc(i-1,j,1)) / (TWO * dx(1))
       viscy = (THREE * visc(i,j,1) - FOUR * visc(i,j-1,1) + visc(i,j-2,1)) / (TWO * dx(2))
       nonlin(i,j,1) = TWO * viscx * ux + viscy * (uy + vx)
       nonlin(i,j,2) = TWO * viscy * vy + viscx * (uy + vx)
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
