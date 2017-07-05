module tag_boxes_module

  use multifab_module
  use bl_error_module
  use probin_module, only: prob_type

  implicit none 

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! MUST set this to .true. if tagging uses ghost cells (e.g., tagging on gradient). !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  logical, save :: tagging_needs_ghost_cells = .false.

contains

  subroutine tag_boxes(mf,tagboxes,dx,lev,aux_tag_mf)

    type( multifab)          , intent(in   ) :: mf
    type(lmultifab)          , intent(inout) :: tagboxes
    real(dp_t)               , intent(in   ) :: dx
    integer                  , intent(in   ) :: lev
    type( multifab), optional, intent(in   ) :: aux_tag_mf
    ! aux_tag_mf allows user to pass in additional multifabs for tagging logic
    
    ! local variables
    real(kind = dp_t), pointer :: mfp(:,:,:,:)
    real(kind = dp_t), pointer ::  up(:,:,:,:)
    logical          , pointer ::  tp(:,:,:,:)
    integer           :: i, lo(get_dim(mf)), hi(get_dim(mf)), ng

    ng = nghost(mf)

    do i = 1, nfabs(mf)
       mfp => dataptr(mf, i)
       tp  => dataptr(tagboxes, i)
       up  => dataptr(aux_tag_mf, i)
       lo =  lwb(get_box(tagboxes, i))
       hi =  upb(get_box(tagboxes, i))
       select case (get_dim(mf))
       case (2)
          call tag_boxes_2d(tp(:,:,1,1),mfp(:,:,1,1),up(:,:,1,:),lo,hi,ng,dx,lev)
       case  (3)
          call tag_boxes_3d(tp(:,:,:,1),mfp(:,:,:,1),lo,hi,ng,dx,lev)
       end select
    end do

  end subroutine tag_boxes

  subroutine tag_boxes_2d(tagbox,mf,u,lo,hi,ng,dx,lev)

    integer          , intent(in   ) :: lo(:),hi(:),ng
    logical          , intent(  out) :: tagbox(lo(1)   :,lo(2)   :)
    real(kind = dp_t), intent(in   ) ::     mf(lo(1)-ng:,lo(2)-ng:)
    real(kind = dp_t), intent(in   ) ::      u(lo(1)-ng:,lo(2)-ng:,:)
    real(dp_t)       , intent(in   ) :: dx
    integer          , intent(in   ) :: lev

    ! local variables
    integer :: i,j

    ! initially say that we do not want to tag any cells for refinement
    tagbox = .false.

    if (prob_type .eq. 1 .or. prob_type .eq. 2) then
       select case(lev)
       case (1)
          ! tag all boxes where the first component of mf >= 1.01
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                if (mf(i,j) .gt. 1.01_dp_t) then
                   tagbox(i,j) = .true.
                end if
             end do
          enddo
       case (2)
          ! for level 2 tag all boxes where the first component of mf >= 1.1
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                if (mf(i,j) .gt. 1.1_dp_t) then
                   tagbox(i,j) = .true.
                end if
             end do
          end do
       case default
          ! for level 3 or greater tag all boxes where the first component of mf >= 1.5
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                if (mf(i,j) .gt. 1.5_dp_t) then
                   tagbox(i,j) = .true.
                end if
             end do
          end do
       end select
    else if (prob_type .eq. 8) then
       ! two-phase flow test case where density = 2 to the left of interface and density = 1 to its right. 
       ! cells near the interface are tagged for refinement
       select case(lev)
       case (1)
          ! level 1: tag all boxes where density is 1% different from 1 or 2
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                if (abs(mf(i,j)-1.5_dp_t) .lt. 0.49_dp_t) then
                   tagbox(i,j) = .true.
                end if
             end do
          end do
       case (2)
          ! level 2: tag all boxes where density is 5% different from 1 or 2
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                if (abs(mf(i,j)-1.5_dp_t) .lt. 0.45_dp_t) then
                   tagbox(i,j) = .true.
                end if
             end do
          end do
       case default
          ! level 3++: tag all boxes where density is 10% different from 1 or 2
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                if (abs(mf(i,j)-1.5_dp_t) .lt. 0.4_dp_t) then
                   tagbox(i,j) = .true.
                end if
             end do
          end do
       end select
    else if (prob_type .eq. 9) then
       select case(lev)
       case (1)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                ! vx = (mf(i+1,j,2) - mf(i-1,j,2)) / (2.0d0 * dx)
                ! uy = (mf(i,j+1,1) - mf(i,j-1,1)) / (2.0d0 * dx) 
                ! if (abs(vx-uy) > 1.0d0)  then
                !    tagbox(i,j) = .true.
                ! end if
                if(abs(u(i,j,1)) < 1.0d0) then 
                   tagbox(i,j) = .true. 
                end if
             end do
          enddo
       case (2)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                ! vx = (mf(i+1,j,2) - mf(i-1,j,2)) / (2.0d0 * dx)
                ! uy = (mf(i,j+1,1) - mf(i,j-1,1)) / (2.0d0 * dx) 
                ! if (abs(vx-uy) > 10.0d0)  then
                !    tagbox(i,j) = .true.
                ! end if
                if(abs(u(i,j,1)) < 0.1d0) then 
                   tagbox(i,j) = .true. 
                end if
             end do
          end do
       case default
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                ! vx = (mf(i+1,j,2) - mf(i-1,j,2)) / (2.0d0 * dx)
                ! uy = (mf(i,j+1,1) - mf(i,j-1,1)) / (2.0d0 * dx) 
                ! if (abs(vx-uy) > 20.0d0)  then
                !    tagbox(i,j) = .true.
                ! end if
                if(abs(u(i,j,1)) < 0.01d0) then 
                   tagbox(i,j) = .true. 
                end if
             end do
          end do
       end select
    else
       call bl_error('Unsupported prob_type')
    end if

  end subroutine tag_boxes_2d

  subroutine tag_boxes_3d(tagbox,mf,lo,hi,ng,dx,lev)

    integer          , intent(in   ) :: lo(:),hi(:),ng
    logical          , intent(  out) :: tagbox(lo(1)   :,lo(2)   :,lo(3)   :)
    real(kind = dp_t), intent(in   ) ::     mf(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    real(dp_t)       , intent(in   ) :: dx
    integer          , intent(in   ) :: lev

    ! local variables
    integer :: i,j,k

    ! initially say that we do not want to tag any cells for refinement
    tagbox = .false.

    if (prob_type .eq. 1 .or. prob_type .eq. 2) then
    select case(lev)
    case (1)
       ! tag all boxes where the first component of mf >= 1.01
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                if (mf(i,j,k) .gt. 1.01_dp_t) then
                   tagbox(i,j,k) = .true.
                end if
             end do
          enddo
       end do
    case (2)
       ! for level 2 tag all boxes where the first component of mf >= 1.1
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                if (mf(i,j,k) .gt. 1.1_dp_t) then
                   tagbox(i,j,k) = .true.
                end if
             end do
          end do
       end do
    case default
       ! for level 3 or greater tag all boxes where the first component of mf >= 1.5
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                if (mf(i,j,k) .gt. 1.5_dp_t) then
                   tagbox(i,j,k) = .true.
                end if
             end do
          end do
       end do
    end select
    else if (prob_type .eq. 3) then
    select case(lev)
    case (1)
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                if (mf(i,j,k) .gt. 1.2d0 .and. mf(i,j,k) .lt. 1.8d0) then
                   tagbox(i,j,k) = .true.
                end if
             end do
          enddo
       end do
    case (2)
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                if (mf(i,j,k) .gt. 1.2d0 .and. mf(i,j,k) .lt. 1.8d0) then
                   tagbox(i,j,k) = .true.
                end if
             end do
          end do
       end do
    case default
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                if (mf(i,j,k) .gt. 1.2d0 .and. mf(i,j,k) .lt. 1.8d0) then
                   tagbox(i,j,k) = .true.
                end if
             end do
          end do
       end do
    end select
    else
       call bl_error('Unsupported prob_type')
    end if


  end subroutine tag_boxes_3d

end module tag_boxes_module
