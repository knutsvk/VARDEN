module init_module

  use bl_types
  use bl_constants_module
  use bc_module
  use inlet_bc
  use define_bc_module
  use multifab_module
  use ml_restrict_fill_module
  use ml_layout_module

  implicit none

  private
  public :: initdata, initdata_on_level

contains

  subroutine initdata_on_level(u,s,dx,bc)

    use probin_module, only : nscal

    type(multifab) , intent(inout) :: u,s
    real(kind=dp_t), intent(in   ) :: dx(:)
    type(bc_level) , intent(in   ) :: bc

    real(kind=dp_t), pointer:: uop(:,:,:,:), sop(:,:,:,:)
    integer :: lo(u%dim),hi(u%dim)
    integer :: i,ng,dm

    type(bl_prof_timer), save :: bpt

    call build(bpt,"initdata_on_level")

    ng = u%ng
    dm = u%dim

    do i = 1, nfabs(u)
       uop => dataptr(u,i)
       sop => dataptr(s,i)
       lo =  lwb(get_box(u,i))
       hi =  upb(get_box(u,i))
       select case (dm)
       case (2)
          call initdata_2d(uop(:,:,1,:), sop(:,:,1,:), lo, hi, ng, dx)
       case (3)
          call initdata_3d(uop(:,:,:,:), sop(:,:,:,:), lo, hi, ng, dx)
       end select
    end do

    call multifab_fill_boundary(u)
    call multifab_fill_boundary(s)

    call multifab_physbc(u,1,1,   dm,   bc)
    call multifab_physbc(s,1,dm+1,nscal,bc)

    call destroy(bpt)

  end subroutine initdata_on_level

  subroutine initdata(nlevs,u,s,dx,bc,mla)

    use multifab_physbc_module
    use probin_module, only: nscal

    integer        , intent(in   ) :: nlevs
    type(multifab) , intent(inout) :: u(:),s(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: bc(:)
    type(ml_layout), intent(inout) :: mla


    real(kind=dp_t), pointer:: uop(:,:,:,:), sop(:,:,:,:)
    integer :: lo(u(1)%dim),hi(u(1)%dim)
    integer :: ng,dm,i,n

    type(bl_prof_timer), save :: bpt

    call build(bpt,"initdata")

    ng = u(1)%ng
    dm = u(1)%dim

    do n=1,nlevs

       do i = 1, nfabs(u(n))
          uop => dataptr(u(n),i)
          sop => dataptr(s(n),i)
          lo =  lwb(get_box(u(n),i))
          hi =  upb(get_box(u(n),i))
          select case (dm)
          case (2)
             call initdata_2d(uop(:,:,1,:), sop(:,:,1,:), lo, hi, ng, dx(n,:))
          case (3)
             call initdata_3d(uop(:,:,:,:), sop(:,:,:,:), lo, hi, ng, dx(n,:))
          end select
       end do

    enddo

    call ml_restrict_and_fill(nlevs, u, mla%mba%rr, bc, bcomp=1)
    call ml_restrict_and_fill(nlevs, s, mla%mba%rr, bc, bcomp=dm+1)

    call destroy(bpt)

  end subroutine initdata

  subroutine initdata_2d(u,s,lo,hi,ng,dx)

    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(out) :: u(lo(1)-ng:,lo(2)-ng:,:)  
    real (kind = dp_t), intent(out) :: s(lo(1)-ng:,lo(2)-ng:,:)  
    real (kind = dp_t), intent(in ) :: dx(:)

    !     Local variables
    integer :: i, j
    real (kind = dp_t)  :: x,y,dist
    real (kind = dp_t)  :: xblob = 0.5d0, yblob = 0.5d0, densfact = 2.0d0
    real (kind = dp_t)  :: blobrad = 0.1d0
   
    u = 0.d0

    s(:,:,1) = ONE
    s(:,:,2) = ZERO

    ! add a density and tracer perturbation
    do j=lo(2),hi(2)
       y = dx(2)*(j + HALF)
       do i=lo(1),hi(1)
          x = dx(1)*((i) + HALF)
          dist = SQRT((x-xblob)**2 + (y-yblob)**2)
          s(i,j,1) = ONE + HALF*(densfact-ONE)*(ONE-TANH(30.*(dist-blobrad))) 
          s(i,j,2) = s(i,j,1)
       enddo
    enddo

  end subroutine initdata_2d

  subroutine initdata_3d(u,s,lo,hi,ng,dx)

    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(out) :: u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
    real (kind = dp_t), intent(out) :: s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
    real (kind = dp_t), intent(in ) :: dx(:)

    !     Local variables
    integer :: i, j, k
    real (kind = dp_t)  :: x,y,z,dist
    real (kind = dp_t)  :: xblob = 0.5d0, yblob = 0.5d0, zblob = 0.5d0, densfact = 10.0d0
    real (kind = dp_t)  :: blobrad = 0.1d0

    u = 0.d0

    s(:,:,:,1) = ONE
    s(:,:,:,2) = ZERO

    ! add a density and tracer perturbation
    do k=lo(3),hi(3)
       z = dx(3)*(k + HALF)
       do j=lo(2),hi(2)
          y = dx(2)*(j + HALF)
          do i=lo(1),hi(1)
             x = dx(1)*((i) + HALF)
             dist = SQRT((x-xblob)**2 + (y-yblob)**2 + (z-zblob)**2)
             s(i,j,k,1) = ONE + HALF*(densfact-ONE)*(ONE-TANH(30.*(dist-blobrad))) 
             s(i,j,k,2) = s(i,j,k,1)
          enddo
       enddo
    enddo

  end subroutine initdata_3d

end module init_module
