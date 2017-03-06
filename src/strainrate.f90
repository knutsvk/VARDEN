module strainrate_module
  
  use bl_types
  use multifab_module
  use define_bc_module
  use ml_layout_module

  implicit none

  private

  public :: update_strainrate

contains

   subroutine update_strainrate(mla, dotgam, u, dx, the_bc_level)

    use bl_constants_module
    use ml_restrict_fill_module
    use probin_module, only: extrap_comp

    type(ml_layout)   , intent(in   ) :: mla
    type(multifab)    , intent(inout) :: dotgam(:)
    type(multifab)    , intent(in   ) :: u(:)
    real(kind = dp_t) , intent(in   ) :: dx(:,:)
    type(bc_level)    , intent(in   ) :: the_bc_level(:)

    ! local
    real(kind=dp_t), pointer :: dgp(:,:,:,:)
    real(kind=dp_t), pointer ::  up(:,:,:,:)

    integer :: lo(get_dim(dotgam(1))), hi(get_dim(dotgam(1)))
    integer :: i, dm, nscal, n, nlevs
    integer :: ng_dg, ng_u
    
    type(bl_prof_timer), save :: bpt

    call build(bpt,"strainrate")

    nlevs = mla%nlevel
    dm    = mla%dim

    ng_dg = dotgam(1)%ng
    ng_u = u(1)%ng

    do n=1, nlevs

       do i = 1, nfabs(dotgam(n))
          dgp => dataptr(dotgam(n),i)
          up  => dataptr( u(n),i)
          lo = lwb(get_box(dotgam(n),i))
          hi = upb(get_box(dotgam(n),i))
          select case (dm)
          case (2)
             call update_strainrate_2d(dgp(:,:,1,:), up(:,:,1,:), lo, hi, ng_dg, ng_u, dx(n,:))
          case (3)
             call update_strainrate_3d(dgp(:,:,:,:), up(:,:,:,:), lo, hi, ng_dg, ng_u, dx(n,:))
          end select
       end do

    enddo ! end loop over levels

    ! restrict cell-centered multifab data, fill all boundaries
    call ml_restrict_and_fill(nlevs, dotgam, mla%mba%rr, the_bc_level, bcomp=extrap_comp)

    call destroy(bpt)

  end subroutine update_strainrate

  subroutine update_strainrate_2d(dotgam, u, lo, hi, ng_dg, ng_u, dx)

    use bl_constants_module

    integer           , intent(in   ) :: lo(:), hi(:), ng_dg, ng_u
    real (kind = dp_t), intent(inout) :: dotgam(lo(1)-ng_dg:,lo(2)-ng_dg:,:)  
    real (kind = dp_t), intent(in   ) ::      u(lo(1)-ng_u:,lo(2)-ng_u:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)

    integer :: i, j
    real (kind = dp_t) ux, uy, vx, vy

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          ux = (u(i+1,j,1) - u(i-1,j,1)) / (TWO * dx(1))
          uy = (u(i,j+1,1) - u(i,j-1,1)) / (TWO * dx(2))
          vx = (u(i+1,j,2) - u(i-1,j,2)) / (TWO * dx(1))
          vy = (u(i,j+1,2) - u(i,j-1,2)) / (TWO * dx(2))
          dotgam(i,j,1) = sqrt(TWO * (ux**2 + vy**2) + (uy + vx)**2)
       enddo
    enddo

    if (bc(1,1) .eq. INLET .or. bc(1,1) .eq. SLIP_WALL .or. &
         bc(1,1) .eq. NO_SLIP_WALL) then
       i = lo(1)
       do j = lo(2), hi(2)
          ux = -(THREE * u(i,j,1) - FOUR * u(i+1,j,1) + u(i+2,j,1)) / (TWO * dx(1))
          uy = (u(i,j+1,1) - u(i,j-1,1)) / (TWO * dx(2))
          vx = -(THREE * u(i,j,2) - FOUR * u(i+1,j,2) + u(i+2,j,2)) / (TWO * dx(1))
          vy = (u(i,j+1,2) - u(i,j-1,2)) / (TWO * dx(2))
          dotgam(i,j,1) = sqrt(TWO*(ux**2 + vy**2) + (uy + vx)**2)
       end do
    end if

    if (bc(1,2) .eq. INLET .or. bc(1,2) .eq. SLIP_WALL .or. &
         bc(1,2) .eq. NO_SLIP_WALL) then
       i = hi(1)
       do j = lo(2), hi(2)
          ux = (THREE * u(i,j,1) - FOUR * u(i-1,j,1) + u(i-2,j,1)) / (TWO * dx(1))
          uy = (u(i,j+1,1) - u(i,j-1,1)) / (TWO * dx(2))
          vx = (THREE * u(i,j,2) - FOUR * u(i-1,j,2) + u(i-2,j,2)) / (TWO * dx(1))
          vy = (u(i,j+1,2) - u(i,j-1,2)) / (TWO * dx(2))
          dotgam(i,j,1) = sqrt(TWO*(ux**2 + vy**2) + (uy + vx)**2)
       end do
    end if

    if (bc(2,1) .eq. INLET .or. bc(2,1) .eq. SLIP_WALL .or. &
         bc(2,1) .eq. NO_SLIP_WALL) then
       j = lo(2)
       do i = lo(1), hi(1)
          ux = (u(i+1,j,1) - u(i-1,j,1)) / (TWO * dx(1)) 
          uy = -(THREE * u(i,j,1) - FOUR * u(i,j+1,1) + u(i,j+2,1)) / (TWO * dx(2))
          vx = (u(i+1,j,2) - u(i-1,j,2)) / (TWO * dx(1)) 
          vy = -(THREE * u(i,j,2) - FOUR * u(i,j+1,2) + u(i,j+2,2)) / (TWO * dx(2))
          dotgam(i,j,1) = sqrt(TWO*(ux**2 + vy**2) + (uy + vx)**2)
       end do
    end if

    if (bc(2,2) .eq. INLET .or. bc(2,2) .eq. SLIP_WALL .or. &
         bc(2,2) .eq. NO_SLIP_WALL) then
       j = hi(2)
       do i = lo(1), hi(1)
          ux = (u(i+1,j,1) - u(i-1,j,1)) / (TWO * dx(1)) 
          uy = (THREE * u(i,j,1) - FOUR * u(i,j-1,1) + u(i,j-2,1)) / (TWO * dx(2))
          vx = (u(i+1,j,2) - u(i-1,j,2)) / (TWO * dx(1)) 
          vy = (THREE * u(i,j,2) - FOUR * u(i,j-1,2) + u(i,j-2,2)) / (TWO * dx(2))
          dotgam(i,j,1) = sqrt(TWO*(ux**2 + vy**2) + (uy + vx)**2)
       end do
    end if

  end subroutine update_strainrate_2d

  subroutine update_strainrate_3d(dotgam, u, lo, hi, ng_dg, ng_u, dx)

    use bl_constants_module

    integer           , intent(in   ) :: lo(:), hi(:), ng_dg, ng_u
    real (kind = dp_t), intent(inout) :: dotgam(lo(1)-ng_dg:,lo(2)-ng_dg:,lo(3)-ng_dg:,:)  
    real (kind = dp_t), intent(in   ) ::      u(lo(1)-ng_u: ,lo(2)-ng_u: ,lo(3)-ng_u:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)

    integer :: i, j, k
    real (kind = dp_t) ux, uy, uz, vx, vy, vz, wx, wy, wz

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
             dotgam(i,j,k,1) = sqrt(TWO * (ux**2 + vy**2 + wz**2) + (uy + vx)**2 + (vz + wy)**2 + (wx + uz)**2)
          enddo
       enddo
    enddo
     ! TODO: BCs for 3D

  end subroutine update_strainrate_3d

end module strainrate_module
