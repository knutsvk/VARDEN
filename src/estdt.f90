module estdt_module 

  use bl_types
  use multifab_module
  use probin_module, only : verbose

  implicit none

  private

  public :: estdt

contains

  subroutine estdt (lev, u, s, gp, ext_vel_force, lapu, viscosity, visc_grad_term, dx, dtold, dt)

    use probin_module, only: max_dt_growth, cflfac

    integer        , intent( in) :: lev
    type(multifab) , intent( in) :: u
    type(multifab) , intent( in) :: s
    type(multifab) , intent( in) :: gp
    type(multifab) , intent( in) :: ext_vel_force
    type(multifab) , intent( in) :: lapu
    type(multifab) , intent( in) :: viscosity
    type(multifab) , intent( in) :: visc_grad_term
    real(kind=dp_t), intent( in) :: dx(:)
    real(kind=dp_t), intent( in) :: dtold
    real(kind=dp_t), intent(out) :: dt

    real(kind=dp_t), pointer::  up(:,:,:,:), sp(:,:,:,:)
    real(kind=dp_t), pointer:: gpp(:,:,:,:), fp(:,:,:,:)
    real(kind=dp_t), pointer::  lp(:,:,:,:), vp(:,:,:,:)
    real(kind=dp_t), pointer:: vgp(:,:,:,:)
    integer :: lo(get_dim(u)), hi(get_dim(u)), dm
    integer :: ng_u, ng_s, ng_g, ng_f, ng_l, ng_v, ng_vg
    integer :: i
    real(kind=dp_t) :: dt_proc, dt_grid, dt_start

    type(bl_prof_timer), save :: bpt

    call build(bpt,"estdt")

    dm = get_dim(u)

    ng_u  =              u%ng
    ng_s  =              s%ng
    ng_g  =             gp%ng
    ng_f  =  ext_vel_force%ng
    ng_l  =           lapu%ng
    ng_v  =      viscosity%ng
    ng_vg = visc_grad_term%ng

    dt_proc  = 1.d20
    dt_start = 1.d20

    do i = 1, nfabs(u)
       up  => dataptr(u, i)
       sp  => dataptr(s, i)
       gpp => dataptr(gp, i)
       fp  => dataptr(ext_vel_force, i)
       lp  => dataptr(lapu, i)
       vp  => dataptr(viscosity, i)
       vgp => dataptr(visc_grad_term, i)
       lo = lwb(get_box(u, i))
       hi = upb(get_box(u, i))

       dt_grid = 1.d20
       select case (dm)
       case (2)
          call estdt_2d(up(:,:,1,:), ng_u, sp(:,:,1,1), ng_s, gpp(:,:,1,:), ng_g, fp(:,:,1,:), &
                        ng_f, lp(:,:,1,:), ng_l, vp(:,:,1,1), ng_v, vgp(:,:,1,:), ng_vg, &
                        lo, hi, dx, dt_grid)
       case (3)
          call estdt_3d(up(:,:,:,:), ng_u, sp(:,:,:,1), ng_s, gpp(:,:,:,:), ng_g, fp(:,:,:,:), &
                        ng_f, lp(:,:,:,:), ng_l, vp(:,:,:,1), ng_v, lo, hi, dx, dt_grid)
       end select

       dt_proc = min(dt_grid, dt_proc)
    end do

    ! This sets dt to be the min of dt_proc over all processors.
    call parallel_reduce(dt, dt_proc, MPI_MIN)

    if ( dt .eq. dt_start ) then
       dt = min(dx(1), dx(2))
       if ( dm .eq. 3 ) dt = min(dt, dx(3))
    end if

    dt = dt * cflfac

    if ( dtold .gt. 0.0d0 ) dt = min(dt, max_dt_growth * dtold)

    if ( parallel_IOProcessor() .and. verbose .ge. 1 ) then
       write(6,1000) lev, dt
    end if
1000 format("Computing dt at level ", i2," to be ... ", e15.8)

    call destroy(bpt)

  end subroutine estdt

  subroutine estdt_2d(u, ng_u, s, ng_s, gp, ng_g, ext_vel_force, ng_f, lapu, ng_l, &
                      viscosity, ng_v, visc_grad_term, ng_vg, lo, hi, dx, dt)

    integer, intent(in) :: lo(:), hi(:), ng_u, ng_s, ng_g, ng_f, ng_l, ng_v, ng_vg
    real (kind = dp_t), intent(in   ) ::              u(lo(1)-ng_u: ,lo(2)-ng_u: ,:)  
    real (kind = dp_t), intent(in   ) ::              s(lo(1)-ng_s: ,lo(2)-ng_s:   )  
    real (kind = dp_t), intent(in   ) ::             gp(lo(1)-ng_g: ,lo(2)-ng_g: ,:)  
    real (kind = dp_t), intent(in   ) ::  ext_vel_force(lo(1)-ng_f: ,lo(2)-ng_f: ,:)  
    real (kind = dp_t), intent(in   ) ::           lapu(lo(1)-ng_l: ,lo(2)-ng_l: ,:)  
    real (kind = dp_t), intent(in   ) ::      viscosity(lo(1)-ng_v: ,lo(2)-ng_v:   )  
    real (kind = dp_t), intent(in   ) :: visc_grad_term(lo(1)-ng_vg:,lo(2)-ng_vg:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    real (kind = dp_t), intent(inout) :: dt

    !     Local variables
    real (kind = dp_t) u_max, v_max, fx, fy, fx_max, fy_max
    real (kind = dp_t) eps
    integer :: i, j

    eps = 1.0e-8

    u_max  = 0.0d0 
    v_max  = 0.0d0 
    fx_max = 0.0d0 
    fy_max = 0.0d0 

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          u_max  = max(u_max , abs(u(i,j,1)))
          v_max  = max(v_max , abs(u(i,j,2)))
          fx = (viscosity(i,j) * lapu(i,j,1) + visc_grad_term(i,j,1) - gp(i,j,1)) / s(i,j) &
             + ext_vel_force(i,j,1)
          fy = (viscosity(i,j) * lapu(i,j,2) + visc_grad_term(i,j,2) - gp(i,j,2)) / s(i,j) &
             + ext_vel_force(i,j,2)
          fx_max = max(fx_max, abs(fx))
          fy_max=  max(fy_max, abs(fy))
       enddo  
    enddo

    if ( u_max .gt. eps ) dt = min(dt, dx(1) / u_max)
    if ( v_max .gt. eps ) dt = min(dt, dx(2) / v_max)

    if ( fx_max .gt. eps ) dt = min(dt, sqrt(2.0d0 * dx(1) / fx_max))
    if ( fy_max .gt. eps ) dt = min(dt, sqrt(2.0d0 * dx(2) / fy_max))

  end subroutine estdt_2d

  subroutine estdt_3d(u, ng_u, s, ng_s, gp, ng_g, ext_vel_force, ng_f, lapu, ng_l, &
                      viscosity, ng_v, lo, hi, dx, dt)

    integer, intent(in) :: lo(:), hi(:), ng_u, ng_s, ng_g, ng_f, ng_l, ng_v
    real (kind = dp_t), intent(in   ) ::             u(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:,:)  
    real (kind = dp_t), intent(in   ) ::             s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:)  
    real (kind = dp_t), intent(in   ) ::            gp(lo(1)-ng_g:,lo(2)-ng_g:,lo(3)-ng_g:,:)  
    real (kind = dp_t), intent(in   ) :: ext_vel_force(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:,:)  
    real (kind = dp_t), intent(in   ) ::          lapu(lo(1)-ng_l:,lo(2)-ng_l:,lo(3)-ng_l:,:)  
    real (kind = dp_t), intent(in   ) ::     viscosity(lo(1)-ng_v:,lo(2)-ng_v:,lo(3)-ng_v:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    real (kind = dp_t), intent(inout) :: dt

    !     Local variables
    real (kind = dp_t)  u_max, v_max, w_max, fx, fy, fz, fx_max, fy_max, fz_max
    real (kind = dp_t)  eps
    integer :: i, j, k

    eps = 1.0e-8

    u_max  = 0.0d0 
    v_max  = 0.0d0 
    w_max  = 0.0d0
    fx_max = 0.0d0 
    fy_max = 0.0d0 
    fz_max = 0.0d0 

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             u_max  = max(u_max, abs(u(i,j,k,1)))
             v_max  = max(v_max, abs(u(i,j,k,2)))
             w_max  = max(w_max, abs(u(i,j,k,3)))
             fx = (viscosity(i,j,k)*lapu(i,j,k,1) - gp(i,j,k,1)) / s(i,j,k) + ext_vel_force(i,j,k,1)
             fy = (viscosity(i,j,k)*lapu(i,j,k,2) - gp(i,j,k,2)) / s(i,j,k) + ext_vel_force(i,j,k,2)
             fz = (viscosity(i,j,k)*lapu(i,j,k,3) - gp(i,j,k,3)) / s(i,j,k) + ext_vel_force(i,j,k,3)
             fx_max = max(fx_max, abs(fx))
             fy_max = max(fy_max, abs(fy))
             fz_max = max(fz_max, abs(fz))
          enddo
       enddo
    enddo

    if ( u_max .gt. eps ) dt = min(dt, dx(1) / u_max)
    if ( v_max .gt. eps ) dt = min(dt, dx(2) / v_max)
    if ( w_max .gt. eps ) dt = min(dt, dx(3) / w_max)

    if ( fx_max .gt. eps ) dt = min(dt, sqrt(2.0d0 * dx(1) / fx_max))
    if ( fy_max .gt. eps ) dt = min(dt, sqrt(2.0d0 * dx(2) / fy_max))
    if ( fz_max .gt. eps ) dt = min(dt, sqrt(2.0d0 * dx(3) / fz_max))

  end subroutine estdt_3d

end module estdt_module
