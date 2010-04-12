! Demonstration of the use of CallbackPot Potential type.

! my_callback() is an example callback routine which calculates
! energies, forces and virial given an atoms object.  After the
! potential has been constructed with args_str "CallbackPot" and then
! callback routine has been assigned with set_callback(), it can be
! used like any other Potential object.

! James Kermode <james.kermode@kcl.ac.uk>
! 12 April 2010.

module callback_module

  use libAtoms_module
  use QUIP_module
  implicit None

contains

  subroutine my_callback(at)
    type(Atoms), intent(inout) :: at

    logical :: calc_energy, calc_local_e, calc_force, calc_virial, dummy
    real(dp), pointer, dimension(:,:) :: force
    real(dp) :: virial(3,3)

    calc_energy = .false.
    calc_local_e = .false.
    calc_force = .false.
    calc_virial = .false.
    dummy = get_value(at%params, 'calc_energy', calc_energy)
    dummy = get_value(at%params, 'calc_local_e', calc_local_e)
    dummy = get_value(at%params, 'calc_force', calc_force)
    dummy = get_value(at%params, 'calc_virial', calc_virial)

    if (calc_energy) then
       call set_value(at%params, 'energy', 1.0_dp)
    end if

    if (calc_force) then
       call add_property(at, 'force', 0.0_dp, n_cols=3)
       dummy = assign_pointer(at, 'force', force)
       force = 1.0_dp
    end if

    if (calc_virial) then
       virial = 1.0_dp
       call set_value(at%params, 'virial', virial)
    end if

  end subroutine my_callback

end module callback_module

program callback_test

  use libAtoms_module
  use QUIP_module
  use callback_module
  
  implicit none

  type(Atoms) :: at
  type(Potential) :: pot
  real(dp) :: energy, virial(3,3)
  real(dp), allocatable, dimension(:,:) :: force
  real(dp), pointer, dimension(:,:) :: force_ptr  
  logical :: dummy
  
  call system_initialise

  call diamond(at, 5.44_dp, (/14/))
  call calc_connect(at)
  call initialise(pot, 'CallbackPot')
  call set_callback(pot, my_callback)

  call calc(pot, at, e=energy)
  call print('energy = '//energy)
  
  allocate(force(3,at%N))
  call calc(pot, at, f=force)
  call print('force =')
  call print(force)
  deallocate(force)
  
  call calc(pot, at, args_str="calc_force=T calc_virial=T")
  dummy = assign_pointer(at, 'force', force_ptr)
  call print('virial=')
  call print(virial)
  call print('force_ptr=')
  call print(force_ptr)

  call system_finalise

end program callback_test
