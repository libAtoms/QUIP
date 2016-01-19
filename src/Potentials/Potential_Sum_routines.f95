
  !*************************************************************************
  !*
  !*  Potential_Sum routines
  !*
  !*************************************************************************

  recursive subroutine Potential_Sum_Initialise(this, args_str, pot1, pot2, mpi, error)
    type(Potential_Sum), intent(inout) :: this
    character(len=*), intent(in) :: args_str
    type(Potential), intent(in), target :: pot1, pot2
    type(MPI_Context), intent(in), optional :: mpi
    integer, intent(out), optional :: error

    INIT_ERROR(error)

    call finalise(this)

    this%pot1 => pot1
    this%pot2 => pot2

    if (present(mpi)) this%mpi = mpi

  end subroutine Potential_Sum_Initialise

  recursive subroutine Potential_Sum_Finalise(this)
    type(Potential_Sum), intent(inout) :: this
    
    nullify(this%pot1)
    nullify(this%pot2)

  end subroutine Potential_Sum_Finalise

  recursive subroutine Potential_Sum_Print(this, file)
    type(Potential_Sum), intent(inout) :: this
    type(Inoutput), intent(inout), optional :: file

    call print('Potential_Sum:', file=file)
    call print('', file=file)
    if (associated(this%pot1)) then
       call print('Potential 1:', file=file)
       call print(this%pot1, file=file)
       call print('', file=file)
    else
       call print('Potential 1 not initialised', file=file)
       call print('', file=file)
    end if
    if (associated(this%pot2)) then
       call print('Potential 2:', file=file)
       call Print(this%pot2, file=file)
       call print('', file=file)
    else
       call print('Potential 2 not initialised', file=file)
       call print('', file=file)
    end if

  end subroutine Potential_Sum_Print

  recursive subroutine Potential_Sum_Calc(this, at, args_str, error)
    type(Potential_Sum), intent(inout) :: this
    type(Atoms), intent(inout) :: at
    character(*), intent(in), optional :: args_str
    integer, intent(out), optional :: error

    real(dp) :: energy, virial(3,3)
    real(dp), pointer :: at_force_ptr(:,:), at_local_energy_ptr(:), at_local_virial_ptr(:,:)

    real(dp) :: my_e_1, my_e_2
    real(dp), allocatable :: my_local_e_1(:)
    real(dp), allocatable :: my_f_1(:,:), my_local_virial_1(:,:)
    real(dp) :: my_virial_1(3,3)
    type(Dictionary) :: params
    character(STRING_LENGTH) :: calc_energy, calc_force, calc_local_energy, calc_virial, calc_local_virial, calc_args_pot1, calc_args_pot2, my_args_str
    logical :: store_contributions

    INIT_ERROR(error)

    call initialise(params)
    call param_register(params,"energy", "", calc_energy, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params,"force", "", calc_force, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params,"virial", "", calc_virial, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params,"local_energy", "", calc_local_energy, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params,"local_virial", "", calc_local_virial, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params,"calc_args_pot1", "", calc_args_pot1, help_string="additional args_str to pass along to pot1")
    call param_register(params,"calc_args_pot2", "", calc_args_pot2, help_string="additional args_str to pass along to pot2")
    call param_register(params,"store_contributions", "F", store_contributions, help_string="if true, store contributions to sum with _pot1 and _pot2 suffixes")
    if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='Potential_Sum_calc args_str')) then
       RAISE_ERROR('Potential_Sum_calc failed to parse args_str="'//trim(args_str)//'"', error)
    endif
    call finalise(params)

    my_args_str = optional_default("", args_str)

    call calc(this%pot1, at, args_str=trim(my_args_str)//" "//calc_args_pot1, error=error)
    PASS_ERROR(error)
    if (len_trim(calc_energy) > 0) then
       call get_param_value(at, trim(calc_energy), my_e_1)
       call print("Potential_sum my_e_1 " // my_e_1, PRINT_VERBOSE)
       if (store_contributions) call set_param_value(at, trim(calc_energy)//"_pot1", my_e_1)
     endif
    if (len_trim(calc_virial) > 0) then
       call get_param_value(at, trim(calc_virial), my_virial_1)
       if (store_contributions) call set_param_value(at, trim(calc_virial)//"_pot1", my_virial_1)
    endif
    if (len_trim(calc_local_energy) > 0) then
       call assign_property_pointer(at, trim(calc_local_energy), at_local_energy_ptr, error=error)
       PASS_ERROR(error)
       allocate(my_local_e_1(at%N))
       my_local_e_1 = at_local_energy_ptr
       if (store_contributions) call add_property(at, trim(calc_local_energy)//"_pot1", at_local_energy_ptr, overwrite=.true.)
    endif
    if (len_trim(calc_force) > 0) then
       call assign_property_pointer(at, trim(calc_force), at_force_ptr, error=error)
       PASS_ERROR(error)
       allocate(my_f_1(3, at%N))
       my_f_1 = at_force_ptr
       if (store_contributions) call add_property(at, trim(calc_force)//"_pot1", at_force_ptr, overwrite=.true.)
    endif
    if (len_trim(calc_local_virial) > 0) then
       call assign_property_pointer(at, trim(calc_local_virial), at_local_virial_ptr, error=error)
       PASS_ERROR(error)
       allocate(my_local_virial_1(9, at%N))
       my_local_virial_1 = at_local_virial_ptr
       if (store_contributions) call add_property(at, trim(calc_local_virial)//"_pot1", at_local_virial_ptr, overwrite=.true.)
    endif


    call calc(this%pot2, at, args_str=trim(my_args_str)//" "//calc_args_pot2, error=error)
    PASS_ERROR(error)
    if (len_trim(calc_energy) > 0) then
       call get_param_value(at, trim(calc_energy), energy)
       call print("Potential_sum my_e_2 " // energy, PRINT_VERBOSE)
       if (store_contributions) call set_param_value(at, trim(calc_energy)//"_pot2", energy)
       energy = my_e_1 + energy
       call set_param_value(at, trim(calc_energy), energy)
    endif
    if (len_trim(calc_virial) > 0) then
       call get_param_value(at, trim(calc_virial), virial)
       if (store_contributions) call set_param_value(at, trim(calc_virial)//"_pot2", virial)
       virial = my_virial_1 + virial
       call set_param_value(at, trim(calc_virial), virial)
    endif
    if (len_trim(calc_local_energy) > 0) then
       if (store_contributions) call add_property(at, trim(calc_local_energy)//"_pot2", at_local_energy_ptr, overwrite=.true.)
       at_local_energy_ptr = my_local_e_1 + at_local_energy_ptr
    endif
    if (len_trim(calc_force) > 0) then
       if (store_contributions) call add_property(at, trim(calc_force)//"_pot2", at_force_ptr, overwrite=.true.)
       at_force_ptr = my_f_1 + at_force_ptr
    end if
    if (len_trim(calc_local_virial) > 0) then
       if (store_contributions) call add_property(at, trim(calc_local_virial)//"_pot2", at_local_virial_ptr, overwrite=.true.)
       at_local_virial_ptr = my_local_virial_1 + at_local_virial_ptr
    end if

    if (allocated(my_local_e_1)) deallocate(my_local_e_1)
    if (allocated(my_f_1)) deallocate(my_f_1)
    if (allocated(my_local_virial_1)) deallocate(my_local_virial_1)

  end subroutine Potential_Sum_Calc

  recursive function Potential_Sum_Cutoff(this)
    type(Potential_Sum), intent(in) :: this
    real(dp) :: potential_sum_cutoff

    if(associated(this%pot1) .and. associated(this%pot2)) then
       potential_sum_cutoff = max(cutoff(this%pot1), cutoff(this%pot2))
    else
       potential_sum_cutoff = 0.0_dp
    endif

  end function Potential_Sum_Cutoff

