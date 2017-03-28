module pre_advance_module

  use bl_types

  use define_bc_module
  use multifab_module
  use ml_layout_module

  implicit none

  private

  public :: advance_premac

contains

  subroutine advance_premac(mla, uold, sold, lapu, umac, gp, ext_vel_force, viscosity, &
                            visc_grad_term, dx, dt, the_bc_tower)

    use velpred_module
    use mkforce_module

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) ::           uold(:)
    type(multifab) , intent(inout) ::           sold(:)
    type(multifab) , intent(in   ) ::           lapu(:)
    type(multifab) , intent(inout) ::           umac(:,:)
    type(multifab) , intent(in   ) ::             gp(:)
    type(multifab) , intent(in   ) ::  ext_vel_force(:)
    type(multifab) , intent(in   ) ::      viscosity(:)
    type(multifab) , intent(in   ) :: visc_grad_term(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:), dt
    type(bc_tower) , intent(in   ) :: the_bc_tower

    type(multifab)  :: vel_force(mla%nlevel)
    integer         :: dm, n, nlevs
    real(kind=dp_t) :: visc_fac

    type(bl_prof_timer), save :: bpt

    call build(bpt, "advance_premac")

    dm    = mla%dim
    nlevs = mla%nlevel

    do n = 1, nlevs
       call multifab_build(vel_force(n), get_layout(ext_vel_force(n)), dm, 1)
       call setval(vel_force(n), 0.0_dp_t, all=.true.)
    enddo

    visc_fac = 1.0d0
    call mkvelforce(mla, vel_force, ext_vel_force, sold, gp, lapu, viscosity, visc_grad_term, &
                    visc_fac, dx, the_bc_tower)

    call velpred(nlevs, uold, umac, vel_force, dx, dt, the_bc_tower%bc_tower_array, mla)

    do n = 1, nlevs
       call multifab_destroy(vel_force(n))
    enddo

    call destroy(bpt)

  end subroutine advance_premac

end module pre_advance_module
