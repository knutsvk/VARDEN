subroutine varden()

  use BoxLib
  use bl_constants_module
  use list_box_module
  use ml_boxarray_module
  use layout_module
  use multifab_module
  use initialize_module
  use estdt_module
  use vort_module
  use proj_parameters
  use ml_restrict_fill_module
  use bc_module
  use define_bc_module
  use bl_mem_stat_module
  use bl_timer_module
  use box_util_module
  use bl_IO_module
  use fabio_module
  use plotfile_module
  use checkpoint_module
  use advance_module
  use regrid_module
  use explicit_diffusive_module, only : get_explicit_diffusive_term
  use strainrate_module
  use viscosity_module
  use stress_module

  use probin_module, only : dim_in, max_levs, nlevs, ng_cell, ng_grow, pmask, init_iter, max_step, &
                            stop_time, restart, chk_int, plot_int, regrid_int, init_shrink, & 
                            fixed_dt, nodal, ref_ratio, fixed_grids, grids_file_name, & 
                            do_initial_projection, grav, probin_init, probin_close, dpdx, &
                            visc_coef, yield_stress, prob_lo, prob_hi

  implicit none

  integer    :: dm, comp
  real(dp_t) :: time, dt, dtold, dt_lev, dt_temp
  integer    :: n, istep
  integer    :: n_chk_comps
  integer    :: last_plt_written, last_chk_written
  integer    :: init_step
  integer    :: press_comp

  real(dp_t)  , pointer     :: dx(:,:)
  type(ml_layout)           :: mla

  ! Cell-based quantities
  type(multifab), pointer     ::            uold(:)
  type(multifab), pointer     ::            sold(:)
  type(multifab), pointer     ::              gp(:)
  type(multifab), pointer     ::               p(:)
  type(multifab), allocatable ::            unew(:)
  type(multifab), allocatable ::            snew(:)
  type(multifab), allocatable ::         rhohalf(:)
  type(multifab), allocatable ::   ext_vel_force(:)
  type(multifab), allocatable ::  ext_scal_force(:)
  type(multifab), allocatable ::            lapu(:)
  type(multifab), allocatable ::     strain_rate(:)
  type(multifab), allocatable ::       viscosity(:)
  type(multifab), allocatable ::      stress_old(:)
  type(multifab), allocatable ::      stress_new(:)
  type(multifab), allocatable ::        plotdata(:)

  character(len=5)               :: plot_index, check_index
  character(len=256)             :: plot_file_name, check_file_name
  character(len=20), allocatable :: plot_names(:)

  type(bc_tower) ::  the_bc_tower
  type(bc_level) ::  bc

  last_plt_written = -1
  last_chk_written = -1

  call probin_init()

  dm = dim_in
  press_comp = dm + nscal + 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Set up plot_names for writing plot files.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  allocate(plot_names(2*dm+nscal+5))

  plot_names(1) = "x_vel"
  plot_names(2) = "y_vel"
  if (dm > 2) plot_names(3) = "z_vel"
  plot_names(dm+1) = "density"
  if (nscal > 1) plot_names(dm+2) = "tracer"
  plot_names(dm+nscal+1) = "magvel"
  plot_names(dm+nscal+2) = "vort"
  plot_names(dm+nscal+3) = "strainrate"
  plot_names(dm+nscal+4) = "viscosity"
  plot_names(dm+nscal+5) = "stress"
  plot_names(dm+nscal+6) = "gpx"
  plot_names(dm+nscal+7) = "gpy"
  if (dm > 2) plot_names(dm+nscal+8) = "gpz"


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Initialize the grids and the data.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (restart >= 0) then

     call initialize_from_restart(mla, restart, time, dt, dx, pmask, uold, sold, gp, p, &
                                  the_bc_tower)

  else if (fixed_grids /= '') then

     call initialize_with_fixed_grids(mla, pmask, dx, uold, sold, gp, p, the_bc_tower)

  else  ! Adaptive gridding

     call initialize_with_adaptive_grids(mla, pmask, dx, uold, sold, gp, p, the_bc_tower)

  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Allocate new-time state and temp variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  allocate(unew(nlevs), snew(nlevs))
  allocate(ext_vel_force(nlevs), ext_scal_force(nlevs))
  allocate(lapu(nlevs))
  allocate(strain_rate(nlevs), viscosity(nlevs))
  allocate(stress_old(nlevs), stress_new(nlevs))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Initial projection if not restart
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (restart < 0) then

     time    = ZERO
     dt_temp = ONE

     ! Note that we use rhohalf, filled with 1 at this point, as a temporary
     ! in order to do a constant-density initial projection.
     if (do_initial_projection > 0) then
        allocate(rhohalf(nlevs))
        do n = 1,nlevs
          call multifab_build(rhohalf(n), mla%la(n),1, 1)
          call setval(rhohalf(n),ONE, all=.true.)
        end do
        call hgproject(initial_projection,mla,uold,uold,rhohalf,p,gp,dx,dt_temp, &
                       the_bc_tower,press_comp)
         do n = 1,nlevs
           call multifab_destroy(rhohalf(n))
        end do
        deallocate(rhohalf)
     end if

     do n = 1,nlevs
        call setval( p(n)  ,0.0_dp_t, all=.true.)
        call setval(gp(n)  ,0.0_dp_t, all=.true.)
     end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Write the grids into a "grdlog" file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     if (grids_file_name /= '') &
        call write_grids(grids_file_name,mla,0)

  else

     if (grids_file_name /= '') &
        call write_grids(grids_file_name,mla,restart)

  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Create "new" unew, snew, ext_vel_force, ext_scal_force, lapu, strain_rate, viscosity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call make_temps(mla)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Impose bc's on uold and copy to unew
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do n = 1,nlevs

     call multifab_fill_boundary(uold(n))
     call multifab_fill_boundary(sold(n))

     bc = the_bc_tower%bc_tower_array(n)
     call multifab_physbc(uold(n), 1, 1,    dm,    bc, &
                          dx_in=dx(n,:), prob_lo_in=prob_lo, prob_hi_in=prob_hi)
     call multifab_physbc(sold(n), 1, dm+1, nscal, bc, &
                          dx_in=dx(n,:), prob_lo_in=prob_lo, prob_hi_in=prob_hi)

     !    This is done to impose any Dirichlet bc's on unew or snew.
     call multifab_copy_c(unew(n),1,uold(n),1,dm   ,ng=unew(n)%ng)
     call multifab_copy_c(snew(n),1,sold(n),1,nscal,ng=snew(n)%ng)
  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Update lapu, strain rate and viscosity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 ! compute lapu
  do comp = 1, dm
    call get_explicit_diffusive_term(mla, lapu, uold, comp, comp, dx, the_bc_tower)
  end do

  if ( yield_stress > 0.0d0 ) then
     ! compute rate-of-strain magnitude
     call update_strainrate(mla, strain_rate, uold, dx, the_bc_tower%bc_tower_array)

     ! compute viscosity
     call update_viscosity(mla, viscosity, strain_rate, dx, the_bc_tower%bc_tower_array)

     ! compute stress magnitude
     call update_stress(mla, stress_new, viscosity, strain_rate, dx, the_bc_tower%bc_tower_array)
  endif

  if (restart < 0) then

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Compute the time step.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     dt = 1.d20
     dtold = dt
     do n = 1, nlevs
        call estdt(n, uold(n), sold(n), gp(n), ext_vel_force(n), lapu(n), viscosity(n), dx(n,:), &
                   dtold, dt_lev)
        dt = min(dt, dt_lev)
     end do

     dt = dt * init_shrink

     if ( fixed_dt > 0.d0 ) dt = fixed_dt
     if ( stop_time >= 0.d0 ) then
        if (time+dt > stop_time) dt = min(dt, stop_time - time)
     end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Begin the initial iterations to define an initial pressure field.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     if (init_iter > 0) call initial_iters()

     istep = 0

     if ( plot_int > 0 ) then
        if ( mod(istep,plot_int) .eq. 0 ) then
           call write_plotfile(istep)
           last_plt_written = istep
        end if
     end if

     if ( chk_int > 0 ) then
        if ( mod(istep,chk_int) .eq. 0 ) then
           call write_checkfile(istep)
           last_chk_written = istep
        end if
     end if

     init_step = 1

  else 

     init_step = restart+1

  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Begin the real integration.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if ((max_step >= init_step) .and. (time < stop_time .or. stop_time < 0.d0)) then

     do istep = init_step, max_step

        if (max_levs > 1 .and. regrid_int > 0 .and. &
            (mod(istep-1,regrid_int) .eq. 0)) then

           ! Delete everything defined on the old mla.
           call delete_temps()

           ! Keep the state on the previous grid in unew,snew
           ! Create the state on the new grid in uold,sold
           call regrid(mla, uold, sold, gp, p, dx, the_bc_tower)
           if (grids_file_name /= '') &
              call write_grids(grids_file_name, mla, istep)

           ! Create "new" unew, snew, ext_vel_force, ext_scal_force, lapu, strain_rate, viscosity
           call make_temps(mla)

        end if  

        do n = 2, nlevs
           call multifab_fill_ghost_cells(uold(n),uold(n-1), &
                                          ng_cell,mla%mba%rr(n-1,:), &
                                          the_bc_tower%bc_tower_array(n-1), &
                                          the_bc_tower%bc_tower_array(n  ), &
                                          1,1,dm)
           call multifab_fill_ghost_cells(sold(n),sold(n-1), &
                                          ng_cell,mla%mba%rr(n-1,:), &
                                          the_bc_tower%bc_tower_array(n-1), &
                                          the_bc_tower%bc_tower_array(n  ), &
                                          1,dm+1,nscal)
           call multifab_fill_ghost_cells(gp(n),gp(n-1), &
                                          ng_grow,mla%mba%rr(n-1,:), &
                                          the_bc_tower%bc_tower_array(n-1), &
                                          the_bc_tower%bc_tower_array(n  ), &
                                          1,1,dm)
        end do

        do n = 1,nlevs
           call multifab_fill_boundary(uold(n))
           call multifab_fill_boundary(sold(n))
           call multifab_fill_boundary(gp(n))

           bc = the_bc_tower%bc_tower_array(n)
           call multifab_physbc(uold(n), 1,    1,    dm, bc, &
                                dx_in=dx(n,:), prob_lo_in=prob_lo, prob_hi_in=prob_hi)
           call multifab_physbc(sold(n), 1, dm+1, nscal, bc, &
                                dx_in=dx(n,:), prob_lo_in=prob_lo, prob_hi_in=prob_hi)
        end do

        ! compute lapu
        do comp = 1, dm
          call get_explicit_diffusive_term(mla, lapu, uold, comp, comp, dx, the_bc_tower)
        end do

        if ( yield_stress > 0.0d0 ) then
           ! compute rate-of-strain magnitude
           call update_strainrate(mla, strain_rate, uold, dx, the_bc_tower%bc_tower_array)

           ! compute viscosity
           call update_viscosity(mla, viscosity, strain_rate, dx, the_bc_tower%bc_tower_array)

           ! compute stress magnitude
           call update_stress(mla, stress_new, viscosity, strain_rate, dx, &
                              the_bc_tower%bc_tower_array)
        endif

        if (istep > 1) then
           dtold = dt
           dt = 1.d20
           do n = 1,nlevs
              call estdt(n, uold(n), sold(n), gp(n), ext_vel_force(n), lapu(n), viscosity(n), &
                         dx(n,:), dtold, dt_lev)
              dt = min(dt, dt_lev)
           end do
           if (fixed_dt > 0.d0) dt = fixed_dt
           if (stop_time >= 0.d0) then
              if (time+dt > stop_time) then
                 dt = stop_time - time
                 if (parallel_IOProcessor()) &
                    print*, "Stop time limits dt =",dt
              end if
           end if
        end if

        call advance_timestep(istep, mla, sold, uold, snew, unew, gp, p, ext_vel_force, &
                              ext_scal_force, lapu, viscosity, the_bc_tower, dt, time, dx, &
                              press_comp, regular_timestep)

        time = time + dt

         if (parallel_IOProcessor()) then
            write(6,1000) istep,time,dt
         end if

         if ( plot_int > 0 ) then
           if ( mod(istep,plot_int) .eq. 0 ) then
              call write_plotfile(istep)
              last_plt_written = istep
           end if
         end if

         if ( chk_int > 0 ) then
           if ( mod(istep,chk_int) .eq. 0 ) then
              call write_checkfile(istep)
              last_chk_written = istep
           end if
         end if
 
         if ( verbose > 0 ) call print_and_reset_fab_byte_spread()

         if ( stop_time >= 0.d0 ) then
            if ( time >= stop_time ) goto 999
         else if ( istep > init_step+1 ) then
            if ( yield_stress > 0.d0 ) then
               if ( steady_state(mla, stress_old, stress_new) ) goto 999
            else
               if ( steady_state(mla, uold, unew) ) goto 999
            end if
         end if

         do n = 1, nlevs
            call multifab_copy_c(uold(n),1,unew(n),1,dm)
            call multifab_copy_c(sold(n),1,snew(n),1,nscal)
            call multifab_copy(stress_old(n), stress_new(n))
         end do

     end do ! istep loop

999  continue

     if (istep > max_step) istep = max_step

     if (last_plt_written .ne. istep .and. plot_int > 0) call write_plotfile(istep)
     if (last_chk_written .ne. istep .and. chk_int  > 0) call write_checkfile(istep)

  end if
  
  call delete_state(uold,sold,gp,p)

  deallocate(uold,sold,p,gp)
  deallocate(dx)

  call delete_temps()

  call bc_tower_destroy(the_bc_tower)

  call destroy(mla)

  call probin_close()

  if ( verbose > 0 ) then
     if ( parallel_IOProcessor() ) then
        print *, 'MEMORY STATS AT END OF RUN '
        print*, ' '
     end if
     call print(multifab_mem_stats(),    "    multifab")
     call print(fab_mem_stats(),         "         fab")
     call print(boxarray_mem_stats(),    "    boxarray")
     call print(layout_mem_stats(),      "      layout")
     call print(boxassoc_mem_stats(),    "    boxassoc")
     call print(fgassoc_mem_stats(),     "     fgassoc")
     call print(syncassoc_mem_stats(),   "   syncassoc")
     call print(copyassoc_mem_stats(),   "   copyassoc")
     call print(fluxassoc_mem_stats(),   "   fluxassoc")
     if ( parallel_IOProcessor() ) print*, ''
  end if

1000 format('STEP = ',i4,1x,' TIME = ',f14.10,1x,'DT = ',f14.9)

contains

  subroutine make_temps(mla_loc)

    type(ml_layout),intent(in   ) :: mla_loc

    do n = nlevs,1,-1
       call multifab_build(          unew(n), mla_loc%la(n),    dm, ng_cell)
       call multifab_build(          snew(n), mla_loc%la(n), nscal, ng_cell)
       call multifab_build( ext_vel_force(n), mla_loc%la(n),    dm, 1)
       call multifab_build(ext_scal_force(n), mla_loc%la(n), nscal, 1)
       call multifab_build(          lapu(n), mla_loc%la(n),    dm, 0)
       call multifab_build(   strain_rate(n), mla_loc%la(n),     1, 0)
       call multifab_build(     viscosity(n), mla_loc%la(n),     1, 1)
       call multifab_build(    stress_old(n), mla_loc%la(n),     1, 0)
       call multifab_build(    stress_new(n), mla_loc%la(n),     1, 0)

       call setval(          unew(n),      ZERO,           all=.true.)
       call setval(          snew(n),      ZERO,           all=.true.)
       call setval( ext_vel_force(n),      dpdx,  1, dm-1, all=.true.)
       call setval( ext_vel_force(n),      grav, dm,    1, all=.true.)
       call setval(ext_scal_force(n),      ZERO,           all=.true.)
       call setval(          lapu(n),      ZERO,           all=.true.)
       call setval(   strain_rate(n),      ZERO,           all=.true.)
       call setval(      viscosity(n), visc_coef,           all=.true.)
       call setval(    stress_old(n),      ZERO,           all=.true.)
       call setval(    stress_new(n),      ZERO,           all=.true.)
    end do

  end subroutine make_temps

  subroutine delete_temps()

    do n = 1,nlevs
       call multifab_destroy(unew(n))
       call multifab_destroy(snew(n))
       call multifab_destroy( ext_vel_force(n))
       call multifab_destroy(ext_scal_force(n))
       call multifab_destroy(          lapu(n))
       call multifab_destroy(   strain_rate(n))
       call multifab_destroy(     viscosity(n))
       call multifab_destroy(    stress_old(n))
       call multifab_destroy(    stress_new(n))
    end do

  end subroutine delete_temps

  subroutine delete_state(u,s,gp,p)

    type(multifab), intent(inout) :: u(:),s(:),p(:),gp(:)

    do n = 1,nlevs
       call multifab_destroy( u(n))
       call multifab_destroy( s(n))
       call multifab_destroy( p(n))
       call multifab_destroy(gp(n))
    end do

  end subroutine delete_state

  subroutine initial_iters()

    if (parallel_IOProcessor() .and. verbose .ge. 1) &
         print *,'DOING ',init_iter,' INITIAL ITERATIONS ' 

    do istep = 1,init_iter

       do n = 2, nlevs
          call multifab_fill_ghost_cells(uold(n),uold(n-1), &
                                         ng_cell,mla%mba%rr(n-1,:), &
                                         the_bc_tower%bc_tower_array(n-1), &
                                         the_bc_tower%bc_tower_array(n  ), &
                                         1,1,dm)
          call multifab_fill_ghost_cells(sold(n),sold(n-1), &
                                         ng_cell,mla%mba%rr(n-1,:), &
                                         the_bc_tower%bc_tower_array(n-1), &
                                         the_bc_tower%bc_tower_array(n  ), &
                                         1,dm+1,nscal)
          call multifab_fill_ghost_cells(gp(n),gp(n-1), &
                                         ng_grow,mla%mba%rr(n-1,:), &
                                         the_bc_tower%bc_tower_array(n-1), &
                                         the_bc_tower%bc_tower_array(n  ), &
                                         1,1,dm)
       end do

       call advance_timestep(istep, mla, sold, uold, snew, unew, gp, p, ext_vel_force, & 
                             ext_scal_force, lapu, viscosity, the_bc_tower, dt, time, dx, &
                             press_comp, pressure_iters)

    end do

  end subroutine initial_iters

  subroutine write_plotfile(istep_to_write)

    use probin_module, only : prob_lo, prob_hi, plot_base_name
    use bc_module
    use define_bc_module
    use ml_boxarray_module

    integer, intent(in   )  :: istep_to_write

    integer                 :: n, n_plot_comps
    integer                 :: mvel_comp, vort_comp, gpx_comp
    integer                 :: strainrate_comp, viscosity_comp, stress_comp
    logical                 :: coarsen_plot_data
 
    ! These are only used if you want to coarsen your plotdata before writing
    ! Start crse
    type(box)                   :: pd_crse
    type(boxarray)              :: ba_crse
    type(layout)                :: la_crse
    type(multifab), allocatable :: mf_crse(:)
    real(dp_t),     allocatable :: dx_crse(:)
    integer,        allocatable :: ref_ratio(:)
    integer                     :: rr_for_write(nlevs-1)
    integer                     :: coarsening_factor

    allocate(ref_ratio(dm))
    allocate(dx_crse(dm))
    allocate(mf_crse(nlevs))
    ! End crse
  
    coarsen_plot_data = .false.
    coarsening_factor = 2

    allocate(plotdata(nlevs))

    n_plot_comps = 2 * dm + nscal + 5

    do n = 1, nlevs
       call multifab_build(plotdata(n), mla%la(n), n_plot_comps, 0)
       call multifab_copy_c(plotdata(n),    1, uold(n), 1, dm)
       call multifab_copy_c(plotdata(n), 1+dm, sold(n), 1, nscal)

       mvel_comp = 1 + dm + nscal
       call make_magvel(plotdata(n), mvel_comp, uold(n))

       vort_comp = mvel_comp + 1
       call make_vorticity(plotdata(n), vort_comp, uold(n), dx(n,:), &
                           the_bc_tower%bc_tower_array(n))

       strainrate_comp = vort_comp + 1
       call multifab_copy_c(plotdata(n), strainrate_comp, strain_rate(n), 1, 1)

       viscosity_comp = strainrate_comp + 1
       call multifab_copy_c(plotdata(n), viscosity_comp, viscosity(n), 1, 1)

       stress_comp = viscosity_comp + 1
       call multifab_copy_c(plotdata(n), stress_comp, stress_old(n), 1, 1)

       gpx_comp = stress_comp + 1
       call multifab_copy_c(plotdata(n), gpx_comp, gp(n), 1, dm)
    end do

    write(unit=plot_index,fmt='(i5.5)') istep_to_write
    plot_file_name = trim(plot_base_name) // plot_index

    if (coarsen_plot_data) then

       ! We have only implemented this for nlevs = 1 right now
       ref_ratio(1:dm) = coarsening_factor
       rr_for_write(:) = coarsening_factor
       dx_crse(:) = dx(1,:) / ref_ratio

       pd_crse = coarsen(mla%mba%pd(1),ref_ratio)

       do n = 1, nlevs
          call boxarray_build_copy(ba_crse,get_boxarray(mla%la(n)))
          call boxarray_coarsen(ba_crse,ref_ratio)
          call layout_build_ba(la_crse,ba_crse,coarsen(mla%mba%pd(n),ref_ratio))
          call print(mla%la(n),'LA FINE')
          call print(la_crse,'LA CRSE')
          call multifab_build(mf_crse(n), la_crse, n_plot_comps, 0)
          call destroy(ba_crse)

          call ml_cc_restriction(mf_crse(n),plotdata(n),ref_ratio)
       end do

       call fabio_ml_multifab_write_d(mf_crse, rr_for_write, plot_file_name, plot_names, &
                                      pd_crse, prob_lo, prob_hi, time, dx_crse)
    else

       call fabio_ml_multifab_write_d(plotdata, mla%mba%rr(:,1), plot_file_name, plot_names, &
                                      mla%mba%pd(1), prob_lo, prob_hi, time, dx(1,:))
    end if

    do n = 1,nlevs
      call multifab_destroy(plotdata(n))
    end do
    deallocate(plotdata)

    if (coarsen_plot_data) then
       do n = 1, nlevs
          call destroy(mf_crse(n))
       end do
       deallocate(mf_crse)
       deallocate(ref_ratio)
       deallocate(dx_crse)
    end if

    call write_job_info(plot_file_name, mla%mba, the_bc_tower)

  end subroutine write_plotfile

  subroutine write_checkfile(istep_to_write)

    use probin_module, only : check_base_name

    integer, intent(in   ) :: istep_to_write

    type(multifab), pointer     ::  chkdata(:)

    allocate(chkdata(nlevs))
    n_chk_comps = 2*dm + nscal
    do n = 1,nlevs
       call multifab_build(chkdata(n), mla%la(n), n_chk_comps, 0)
       call multifab_copy_c(chkdata(n),1         ,uold(n),1,dm)
       call multifab_copy_c(chkdata(n),1+dm      ,sold(n),1,nscal)
       call multifab_copy_c(chkdata(n),1+dm+nscal,  gp(n),1,dm)
    end do
    write(unit=check_index,fmt='(i5.5)') istep_to_write
    check_file_name = trim(check_base_name) // check_index

    call checkpoint_write(nlevs, check_file_name, chkdata, p, mla%mba%rr, time, dt, verbose)

    do n = 1,nlevs
       call multifab_destroy(chkdata(n))
    end do
    deallocate(chkdata)

  end subroutine write_checkfile

  subroutine write_grids(grids_file_name,mla,nstep)

    character(len=128), intent(inout) :: grids_file_name
    type(ml_layout)   , intent(in   ) :: mla
    integer           , intent(in   ) :: nstep

    integer        :: i,un,nb,tp(mla%dim)
    type(box)      :: bx

    un = 11
    tp = 0
   
    if ( parallel_IOProcessor() ) then
       open(un,file=grids_file_name, position='append')
       write(unit=un, fmt='("At step ",i5,":")') nstep
!      write(unit=un, fmt='(i1," levels ")') mla%mba%nlevel
       write(unit=un, fmt='(i2)') mla%mba%nlevel
       do n = 1, mla%mba%nlevel
          nb = nboxes(mla%mba%bas(n))
          bx = ml_layout_get_pd(mla,n)
          write(unit=un, fmt='("   (")', advance = 'no') 
          write(unit=un, fmt='("(", 3(I0,:,", "))', advance = 'no') bx%lo(1:bx%dim)
          write(unit=un, fmt='(") (", 3(I0,:,", "))', advance = 'no') bx%hi(1:bx%dim)
          write(unit=un, fmt='(") (", 3(I0,:,","))', advance = 'no') tp(1:bx%dim)
          write(unit=un, fmt='("))")', advance = 'no' )
          write(unit=un, fmt='(" ",i4)', advance = 'yes') nb
          do i = 1, nb
             bx = get_box(mla%mba%bas(n),i)
             tp = 0
             write(unit=un, fmt='("      (")', advance = 'no') 
             write(unit=un, fmt='("(", 3(I0,:,", "))', advance = 'no') bx%lo(1:bx%dim)
             write(unit=un, fmt='(") (", 3(I0,:,", "))', advance = 'no') bx%hi(1:bx%dim)
             write(unit=un, fmt='(") (", 3(I0,:,","))', advance = 'no') tp(1:bx%dim)
             write(unit=un, fmt='("))")', advance = 'no' )
             write(unit=un, fmt='(" ")')
          end do
       end do
       write(unit=un, fmt='(" ")')
       close(un)
    end if

  end subroutine write_grids

  function steady_state(mla, old, new) result(r)

      type(ml_layout), intent(in) :: mla
      type(multifab ), intent(in) :: old(mla%nlevel)
      type(multifab ), intent(in) :: new(mla%nlevel)

      type(multifab) :: diff(mla%nlevel)
      integer :: nlevs, nc, ng
      real(dp_t) :: max_diff, tol
      logical :: r

      nlevs = mla%nlevel
      nc = ncomp(old(1))
      ng = nghost(old(1))
      max_diff = 0.0d0
      tol = 1.0e-5
      r = .false.

      do n = 1, nlevs
         call multifab_build(diff(n), mla%la(n), nc, ng)
         call multifab_copy(diff(n), old(n))
         call multifab_sub_sub(diff(n), new(n))
         max_diff = max(max_diff, multifab_norm_inf(diff(n)))
      end do

      if ( max_diff <= tol ) r = .true.

  end function steady_state

end subroutine varden
