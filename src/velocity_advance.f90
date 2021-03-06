module velocity_advance_module

  use bl_types
  use bl_constants_module
  use define_bc_module
  use multifab_module
  use ml_layout_module

  implicit none

  private

  public :: velocity_advance

contains

  subroutine velocity_advance(mla, uold, unew, sold, lapu, rhohalf, umac, gp, ext_vel_force, &
                              viscosity, dx, dt, the_bc_tower)

    use viscous_module, only : visc_solve
    use mkflux_module , only : mkflux
    use mkforce_module, only : mkvelforce
    use update_module , only : update
    use probin_module , only : verbose, diffusion_type

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) ::           uold(:)
    type(multifab) , intent(inout) ::           unew(:)
    type(multifab) , intent(in   ) ::           sold(:)
    type(multifab) , intent(inout) ::           lapu(:)
    type(multifab) , intent(in   ) ::        rhohalf(:)
    type(multifab) , intent(in   ) ::           umac(:,:)
    type(multifab) , intent(in   ) ::             gp(:)
    type(multifab) , intent(in   ) ::  ext_vel_force(:)
    type(multifab) , intent(inout) ::      viscosity(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:), dt
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local
    type(multifab)  :: vel_force(mla%nlevel)
    type(multifab)  :: uflux(mla%nlevel, mla%dim)
    type(multifab)  :: uedge(mla%nlevel, mla%dim)
    integer         :: i, n, dm, comp, nlevs
    logical         :: is_vel, is_conservative(mla%dim)
    real(kind=dp_t) :: visc_fac, visc_dt
    real(kind=dp_t) :: umin, umax

    nlevs = mla%nlevel
    dm    = mla%dim

    is_conservative = .false.
    is_vel = .true.

    do n = 1, nlevs
       call multifab_build(vel_force(n), get_layout(ext_vel_force(n)), dm, 1)

       do i = 1,dm
         call multifab_build_edge(uflux(n,i), mla%la(n), dm, 0, i)
         call multifab_build_edge(uedge(n,i), mla%la(n), dm, 0, i)
         call setval(uflux(n,i), ZERO, all=.true.)
         call setval(uedge(n,i), ZERO, all=.true.)
       end do
    enddo

    !********************************************************
    ! Create the velocity forcing term at time n using rho and the full viscous term.
    !********************************************************

    visc_fac = ONE
    call mkvelforce(mla, vel_force, ext_vel_force, sold, gp, lapu, viscosity, visc_fac, dx, &
                    the_bc_tower)

    !********************************************************
    ! Create the edge state velocities
    !********************************************************

    call mkflux(mla, uold, uedge, uflux, umac, vel_force, dx, dt, the_bc_tower%bc_tower_array, &
                is_vel, is_conservative)

    !********************************************************
    ! Now create vel_force at half-time using rhohalf and half the viscous term.
    !********************************************************

    ! The lapu term will be added to the rhs in visc_solve
    ! for Crank-Nicolson
    visc_fac = ZERO
    call mkvelforce(mla, vel_force, ext_vel_force, rhohalf, gp, lapu, viscosity, visc_fac, dx, &
                    the_bc_tower)

    !********************************************************
    ! Update the velocity with convective differencing
    !********************************************************

    call update(mla, uold, umac, uedge, uflux, vel_force, unew, dx, dt, is_vel, is_conservative, &
                the_bc_tower%bc_tower_array)

    do n = 1, nlevs
       call multifab_destroy(vel_force(n))
       do i = 1,dm
          call multifab_destroy(uflux(n,i))
          call multifab_destroy(uedge(n,i))
       end do
    enddo

    ! Crank-Nicolson
    if (diffusion_type .eq. 1) then
       visc_dt = HALF*dt
       
    ! backward Euler
    else if (diffusion_type .eq. 2) then
       visc_dt = dt
       
    else
       call bl_error('BAD DIFFUSION TYPE ')
    end if
    
    if (verbose .ge. 1) then
       if (parallel_IOProcessor()) write(6,*) 'Doing the viscous solve ...'
    end if
    call visc_solve(mla, unew, lapu, rhohalf, viscosity, dx, visc_dt, the_bc_tower)

    if (verbose .ge. 1) then
       do n = 1, nlevs
          do comp = 1, dm
             umin = multifab_min_c(unew(n),comp) 
             umax = multifab_max_c(unew(n),comp)
             if (comp .eq. 1) then
                if (parallel_IOProcessor()) write(6,2001) n,umin,umax
             else if (comp .eq. 2) then
                if (parallel_IOProcessor()) write(6,2002) n,umin,umax
             else if (comp .eq. 3) then
                if (parallel_IOProcessor()) write(6,2003) n,umin,umax
             end if
          end do
       end do
    end if

2001 format('... level ', i2,' new min/max : x-vel           ',e17.10,2x,e17.10)
2002 format('... level ', i2,' new min/max : y-vel           ',e17.10,2x,e17.10)
2003 format('... level ', i2,' new min/max : z-vel           ',e17.10,2x,e17.10)

  end subroutine velocity_advance

end module velocity_advance_module
