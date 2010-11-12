
  !*************************************************************************
  !*
  !*  Potential_EVB routines
  !*
  !*************************************************************************

  recursive subroutine Potential_EVB_Initialise(this, args_str, pot1, mpi, error)
    type(Potential_EVB), intent(inout)         :: this
    character(len=*),    intent(in)            :: args_str
    type(Potential),     intent(in),  target   :: pot1
    type(MPI_Context),   intent(in),  optional :: mpi
    integer,             intent(out), optional :: error

    type(Dictionary) :: params

    INIT_ERROR(error)

    call finalise(this)

    !Only 1 potential for the 2 MM calculations
    this%pot1 => pot1

    call initialise(params)
    call param_register(params, 'mm_args_str', '', this%mm_args_str, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'topology_suffix1', '_EVB1', this%topology_suffix1, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'topology_suffix2', '_EVB2', this%topology_suffix2, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'form_bond', '0 0', this%form_bond, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'break_bond', '0 0', this%break_bond, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'diagonal_dE2', '0.0', this%diagonal_dE2, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'offdiagonal_A12', '0.0', this%offdiagonal_A12, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'offdiagonal_mu12', '0.0', this%offdiagonal_mu12, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'offdiagonal_mu12_square', '0.0', this%offdiagonal_mu12_square, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'offdiagonal_r0', '0.0', this%offdiagonal_r0, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'save_forces', 'T', this%save_forces, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'save_energies', 'T', this%save_energies, help_string="No help yet.  This source file was $LastChangedBy$")
    if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='Potential_EVB_initialise args_str')) then
       RAISE_ERROR('Potential_EVB_initialise failed to parse args_str="'//trim(args_str)//'"', error)
    endif
    call finalise(params)

    if (present(mpi)) this%mpi = mpi

  end subroutine Potential_EVB_Initialise

  recursive subroutine Potential_EVB_Finalise(this)
    type(Potential_EVB), intent(inout) :: this
    
    nullify(this%pot1)

    this%mm_args_str = ""
    this%topology_suffix1 = ""
    this%topology_suffix2 = ""
    this%form_bond(1:2) = 0
    this%break_bond(1:2) = 0
    this%diagonal_dE2 = 0._dp
    this%offdiagonal_A12 = 0._dp
    this%offdiagonal_mu12 = 0._dp
    this%offdiagonal_mu12_square = 0._dp
    this%offdiagonal_r0 = 0._dp
    this%save_forces = .false.
    this%save_energies = .false.

  end subroutine Potential_EVB_Finalise

  recursive subroutine Potential_EVB_Print(this, file)
    type(Potential_EVB), intent(inout)           :: this
    type(Inoutput),      intent(inout), optional :: file

    call print('Potential_EVB:', file=file)
    call print('  mm_args_str='//trim(this%mm_args_str), file=file)
    call print('  topology_suffix1='//trim(this%topology_suffix1), file=file)
    call print('  topology_suffix2='//trim(this%topology_suffix2), file=file)
    call print('  evb1-form and evb2-break bond: '//this%form_bond(1:2), file=file)
    call print('  evb1-break and evb2-form bond: '//this%break_bond(1:2), file=file)
    call print('  diagonal E2(shift to E2): '//this%diagonal_dE2, file=file)
    call print('  offdiagonal A12(pre-exponent factor): '//this%offdiagonal_A12, file=file)
    call print('  offdiagonal mu12(exponent factor): '//this%offdiagonal_mu12, file=file)
    call print('  offdiagonal mu12_square(exponent factor): '//this%offdiagonal_mu12_square, file=file)
    call print('  offdiagonal r0(exponent): '//this%offdiagonal_r0, file=file)
    call print('  save_forces: '//this%save_forces, file=file)
    call print('  save_energies: '//this%save_energies, file=file)
    call print('', file=file)
    if (associated(this%pot1)) then
       call print('Potential 1:', file=file)
       call print(this%pot1, file=file)
       call print('', file=file)
    else
       call print('Potential 1 not initialised', file=file)
       call print('', file=file)
    end if
    call print('Potential 2: same as Potential 1', file=file)
    call print('', file=file)

  end subroutine Potential_EVB_Print

  recursive subroutine Potential_EVB_Calc(this, at, args_str, error)
    type(Potential_EVB), intent(inout)          :: this
    type(Atoms),         intent(inout)          :: at
    character(*),        intent(in),  optional  :: args_str
    integer,             intent(out), optional  :: error

    real(dp) :: e, virial
    real(dp), pointer       :: at_force_ptr(:,:), dgap_dr_ptr(:,:)

    real(dp) :: gap, term1
    type(Dictionary)        :: params
    character(FIELD_LENGTH) :: mm_args_str
    character(FIELD_LENGTH) :: topology_suffix1, topology_suffix2
    integer                 :: form_bond(2), break_bond(2), atom1, atom3
    logical                 :: have_form_bond, have_break_bond

    logical                 :: save_energies, save_forces
    real(dp)                :: my_e_1, my_e_2, e_offdiag
    real(dp), allocatable   :: my_f_1(:,:), my_f_2(:,:), de_offdiag_dr(:,:)
    real(dp)                :: offdiagonal_A12, offdiagonal_r0, offdiagonal_mu12, offdiagonal_mu12_square, diagonal_dE2, &
                               rab, d_rab_dx(3)
    logical                 :: no_coupling, dummy
    character(STRING_LENGTH) :: extra_calc_args

    character(FIELD_LENGTH) :: psf_print
    character(STRING_LENGTH) :: calc_energy, calc_force, calc_virial, calc_local_energy, calc_local_virial, calc_EVB_gap
    character(STRING_LENGTH) :: use_calc_energy
    character(10240) :: new_args_str

    INIT_ERROR(error)

    !read args_str
    call initialise(params)
    call param_register(params, 'mm_args_str', ''//this%mm_args_str, mm_args_str, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'topology_suffix1', ''//this%topology_suffix1, topology_suffix1, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'topology_suffix2', ''//this%topology_suffix2, topology_suffix2, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'form_bond', ''//this%form_bond, form_bond, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'break_bond', ''//this%break_bond, break_bond, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'diagonal_dE2', ''//this%diagonal_dE2, diagonal_dE2, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'offdiagonal_A12', ''//this%offdiagonal_A12, offdiagonal_A12, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'offdiagonal_mu12', ''//this%offdiagonal_mu12, offdiagonal_mu12, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'offdiagonal_mu12_square', ''//this%offdiagonal_mu12_square, offdiagonal_mu12_square, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'offdiagonal_r0', ''//this%offdiagonal_r0, offdiagonal_r0, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'save_forces', ''//this%save_forces, save_forces, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'save_energies', ''//this%save_energies, save_energies, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'energy', '', calc_energy, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'force', '', calc_force, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'virial', '', calc_virial, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'local_energy', '', calc_local_energy, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'local_virial', '', calc_local_virial, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'EVB_gap', '', calc_EVB_gap, help_string="No help yet.  This source file was $LastChangedBy$")
    if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='Potential_EVB_calc args_str')) then
       RAISE_ERROR('Potential_EVB_calc failed to parse args_str="'//trim(args_str)//'"', error)
    endif
    call finalise(params)

    !CHECK ARGUMENTS

    if (len_trim(calc_virial) > 0 .or. len_trim(calc_local_energy) > 0 .or. len_trim(calc_local_virial) > 0) then
       RAISE_ERROR('Potential_EVB_calc: supports only energy and forces, not virial, local_energy or local_virial', error)
    endif

    if (len_trim(calc_force) > 0) then
       call assign_property_pointer(at, trim(calc_force), at_force_ptr, error=error)
       PASS_ERROR_WITH_INFO("Potential_EVB_Calc assigning pointer for force property '"//trim(calc_force)//"'", error)
    endif

    !coupling parameters
    if (offdiagonal_A12 .feq. 0._dp) then
       call print('WARNING! Offdiagonal A12 is set to 0. No coupling between resonance states.')
       no_coupling = .true.
    else
       if (offdiagonal_A12 < 0._dp .or. offdiagonal_mu12 < 0._dp .or. offdiagonal_mu12_square < 0._dp .or. offdiagonal_r0 < 0._dp ) then
          RAISE_ERROR('Potential_EVB_calc offdiagonal parameters must be positive or 0 for no coupling. Got offdiagonal_A12: '//offdiagonal_A12 //' and offdiagonal_mu12: '//offdiagonal_mu12//' and offdiagonal_mu12_square: '//offdiagonal_mu12_square//' and offdiagonal_r0: '//offdiagonal_r0, error)
       endif
       no_coupling = .false.
    endif

    !form_bond
    have_form_bond = .true.
    if(any(form_bond<1) .or. any(form_bond>at%N)) then
       !check whether all 0 (skip)
       if (all(form_bond==0)) then
          have_form_bond = .false.
       else
          RAISE_ERROR('Potential_EVB_calc form_bond is out of range 1--'//at%N//': '//form_bond, error)
       endif
    endif

    !break_bond
    have_break_bond = .true.
    if(any(break_bond<1) .or. any(break_bond>at%N)) then
       !check whether all 0 (skip)
       if (all(break_bond==0)) then
          have_break_bond = .false.
       else
          RAISE_ERROR('Potential_EVB_calc break_bond is out of range 1--'//at%N//': '//break_bond, error)
       endif
    endif

    if (.not.(have_form_bond .or. have_break_bond)) then
       RAISE_ERROR('Potential_EVB_calc no bonds to form neither to break. ', error)
    endif

    !CALCULATE E,F WITH DIFFERENT TOPOLOGIES

    ! allocate local arrays
    if (len_trim(calc_force) > 0) then
      allocate(my_f_1(3,at%N))
      allocate(my_f_2(3,at%N))
      allocate(de_offdiag_dr(3,at%N))
    endif

    ! SETUP CALC_ARGS, AND CALL CALC FOR TOPOLOGY 1
    ! topology suffix
    extra_calc_args="topology_suffix="//trim(topology_suffix1)
    ! add energy= arg if needed
    use_calc_energy=trim(calc_energy)
    if (len_trim(calc_energy) == 0) then
      use_calc_energy="energy"
      do while (has_key(at%params, trim(use_calc_energy)))
	 use_calc_energy = "T"//trim(use_calc_energy)
      end do
    endif
    extra_calc_args=trim(extra_calc_args)//" energy="//trim(use_calc_energy)
    if (len_trim(calc_force) > 0) extra_calc_args=trim(extra_calc_args)//" force="//trim(calc_force)
    ! add args to form/break bonds for topology 1
    if (have_form_bond) extra_calc_args=trim(extra_calc_args)//" form_bond={"//form_bond(1:2)//"}"
    if (have_break_bond) extra_calc_args=trim(extra_calc_args)//" break_bond={"//break_bond(1:2)//"}"
    !calc with topology1
call print("EVB1 ARGS_STR "//trim(mm_args_str)//" "//trim(extra_calc_args))
    call calc(this%pot1, at, args_str=trim(mm_args_str)//" "//trim(extra_calc_args), error=error)
    PASS_ERROR(error)

    call get_param_value(at, trim(use_calc_energy), my_e_1, error=error)
    PASS_ERROR_WITH_INFO("getting energy parameter '"//trim(use_calc_energy)//"' for topology 1", error)
    if (len_trim(calc_energy) == 0) call remove_value(at%params, trim(use_calc_energy))
    if (len_trim(calc_force) > 0) my_f_1 = at_force_ptr

    ! SETUP CALC_ARGS, AND CALL CALC FOR TOPOLOGY 2
    ! topology suffix
    extra_calc_args="topology_suffix="//trim(topology_suffix2)
    ! add energy and force args if needed
    use_calc_energy=trim(calc_energy)
    if (len_trim(calc_energy) == 0) then
      use_calc_energy="energy"
      do while (has_key(at%params, trim(use_calc_energy)))
	 use_calc_energy = "T"//trim(use_calc_energy)
      end do
    endif
    extra_calc_args=trim(extra_calc_args)//" energy="//trim(use_calc_energy)
    if (len_trim(calc_force) > 0) extra_calc_args=trim(extra_calc_args)//" force="//trim(calc_force)
    ! form/break bonds for topology 2
    if (have_form_bond) extra_calc_args=trim(extra_calc_args)//" break_bond={"//form_bond(1:2)//"}"
    if (have_break_bond) extra_calc_args=trim(extra_calc_args)//" form_bond={"//break_bond(1:2)//"}"
    !calc with topology2
call print("EVB2 ARGS_STR "//trim(mm_args_str)//" "//trim(extra_calc_args))
    call calc(this%pot1, at, args_str=trim(mm_args_str)//" "//trim(extra_calc_args), error=error)
    PASS_ERROR(error)

    call get_param_value(at, trim(use_calc_energy), my_e_2, error=error)
    PASS_ERROR_WITH_INFO("getting energy parameter '"//trim(use_calc_energy)//"' for topology 2", error)
    my_e_2 = my_e_2 + diagonal_dE2
    if (len_trim(calc_energy) == 0) call remove_value(at%params, trim(use_calc_energy))
    if (len_trim(calc_force) > 0) my_f_2 = at_force_ptr

    call print("EVB energies my_e_1 " // my_e_1 // " " // my_e_2, PRINT_VERBOSE)

    !CALCULATE EVB ENERGY AND FORCES

    !distance or distance difference
    if (no_coupling) then
       call print("EVB no coupling", PRINT_VERBOSE)
       !take the E, F of the resonance state with the smaller E
       if (my_e_1 < my_e_2) then
	  if (len_trim(calc_energy) > 0) call set_param_value(at, trim(calc_energy), my_e_1)
	  if (len_trim(calc_force) > 0) at_force_ptr = my_f_1
       else
	  if (len_trim(calc_energy) > 0) call set_param_value(at, trim(calc_energy), my_e_2)
	  if (len_trim(calc_force) > 0) at_force_ptr = my_f_2
       endif
       if (len_trim(calc_EVB_gap) > 0) then
          gap=my_e_1-my_e_2
          call set_param_value(at, trim(calc_EVB_gap), gap)
          call print("EVB gap " // gap, PRINT_VERBOSE)
          if (len_trim(calc_force) > 0) then
             call add_property(at, trim(calc_EVB_gap)//"_force", 0.0_dp, n_cols=3, ptr2=dgap_dr_ptr, error=error)
             PASS_ERROR(error)
             dgap_dr_ptr = my_f_1 - my_f_2
             !dgap_dr_ptr = ((my_e_1 - my_e_2)*(my_f_1 - my_f_2) - 4.0_dp*e_offdiag*de_offdiag_dr) / sqrt((my_e_1 - my_e_2)**2.0_dp + 4._dp*e_offdiag**2)
          endif
       endif
    else
       !calculate coupling terms -- rab is the bondlength or the sum of the bond length
       rab = 0._dp
       if (have_form_bond.and.have_break_bond) then
          !the distance of the two atoms around the central one
          atom1=form_bond(1)
          atom3=form_bond(2)
          if (break_bond(1)==atom1) then
             atom1=break_bond(2)
          elseif (break_bond(1)==atom3) then
             atom3=break_bond(2)
          elseif (break_bond(2)==atom1) then
             atom1=break_bond(1)
          elseif (break_bond(2)==atom3) then
             atom3=break_bond(1)
          else
             call system_abort("The forming and breaking bonds are not around a common atom.  This case has not yet been implemented.")
          endif
          rab = distance_min_image(at,atom1,atom3)
       else
          !use the only bond that breaks/forms
          if (have_form_bond) then
             rab = distance_min_image(at,form_bond(1),form_bond(2))
          elseif (have_break_bond) then
             rab = distance_min_image(at,break_bond(1),break_bond(2))
          endif
       endif


       if (len_trim(calc_energy) > 0 .or. len_trim(calc_force) > 0) then
          !original by Warshel
          !e_offdiag = offdiagonal_A12 * exp(-offdiagonal_mu12 * abs(rab))
          !Letif's modification (the original is too strong and sucks the two Cl atoms together)
          e_offdiag = offdiagonal_A12 * exp(-offdiagonal_mu12 * (rab - offdiagonal_r0) - offdiagonal_mu12_square * (rab - offdiagonal_r0)**2._dp)
	  call print("EVB e_offidag " // e_offdiag, PRINT_VERBOSE)
       endif

       !energy
       term1 =  sqrt((my_e_1 - my_e_2)**2._dp + 4._dp*e_offdiag**2)
       if (len_trim(calc_energy) > 0) then
          e = 0.5_dp * (my_e_1 + my_e_2) - 0.5_dp * term1
	  call set_param_value(at, trim(calc_energy), e)
	  call print("EVB coupled energy " // e, PRINT_VERBOSE)
       endif
       if (len_trim(calc_EVB_gap) > 0) then
          gap=my_e_1-my_e_2
	  call set_param_value(at, trim(calc_EVB_gap), gap)
	  call print("EVB gap " // gap, PRINT_VERBOSE)
       endif

       !force
       if (len_trim(calc_force) > 0) then
          !force coupling term
          de_offdiag_dr = 0._dp
          if (have_form_bond.and.have_break_bond) then !breaking and forming bonds around the same atom
             d_rab_dx = diff_min_image(at,atom1,atom3)/distance_min_image(at,atom1,atom3)
             de_offdiag_dr(1:3,atom1) = de_offdiag_dr(1:3,atom1) + e_offdiag * ( offdiagonal_mu12 + 2._dp * offdiagonal_mu12_square * (rab-offdiagonal_r0) ) * d_rab_dx(1:3)
             de_offdiag_dr(1:3,atom3) = de_offdiag_dr(1:3,atom3) - e_offdiag * ( offdiagonal_mu12 + 2._dp * offdiagonal_mu12_square * (rab-offdiagonal_r0) ) * d_rab_dx(1:3)
          elseif (have_form_bond) then !only 1 forming bond
             d_rab_dx = diff_min_image(at,form_bond(1),form_bond(2))/distance_min_image(at,form_bond(1),form_bond(2))
             de_offdiag_dr(1:3,form_bond(1)) = de_offdiag_dr(1:3,form_bond(1)) + e_offdiag * ( offdiagonal_mu12 + 2._dp * offdiagonal_mu12_square * (rab-offdiagonal_r0) ) * d_rab_dx(1:3)
             de_offdiag_dr(1:3,form_bond(2)) = de_offdiag_dr(1:3,form_bond(2)) - e_offdiag * ( offdiagonal_mu12 + 2._dp * offdiagonal_mu12_square * (rab-offdiagonal_r0) ) * d_rab_dx(1:3)
          elseif (have_break_bond) then !only 1 breaking bond
             d_rab_dx = diff_min_image(at,break_bond(1),break_bond(2))/distance_min_image(at,break_bond(1),break_bond(2))
             de_offdiag_dr(1:3,break_bond(1)) = de_offdiag_dr(1:3,break_bond(1)) + e_offdiag * ( offdiagonal_mu12 + 2._dp * offdiagonal_mu12_square * (rab-offdiagonal_r0) ) * d_rab_dx(1:3)
             de_offdiag_dr(1:3,break_bond(2)) = de_offdiag_dr(1:3,break_bond(2)) - e_offdiag * ( offdiagonal_mu12 + 2._dp * offdiagonal_mu12_square * (rab-offdiagonal_r0) ) * d_rab_dx(1:3)
          endif
          !force
          at_force_ptr = 0.5_dp * (my_f_1 + my_f_2) - &
              0.5_dp * ((my_e_1 - my_e_2)*(my_f_1 - my_f_2) - 4.0_dp*e_offdiag*de_offdiag_dr) / sqrt((my_e_1 - my_e_2)**2.0_dp + 4._dp*e_offdiag**2)
	  if (len_trim(calc_EVB_gap) > 0) then
	    call add_property(at, trim(calc_EVB_gap)//"_force", 0.0_dp, n_cols=3, ptr2=dgap_dr_ptr, error=error)
	    PASS_ERROR(error)
            dgap_dr_ptr = my_f_1 - my_f_2
	    !dgap_dr_ptr = ((my_e_1 - my_e_2)*(my_f_1 - my_f_2) - 4.0_dp*e_offdiag*de_offdiag_dr) / sqrt((my_e_1 - my_e_2)**2.0_dp + 4._dp*e_offdiag**2)
	  endif
       endif
    endif

    !SAVE E,F IF NEEDED

if (len_trim(calc_energy) > 0) then
call print("DEBUG EVB E1 " // my_e_1 // " E2 " // my_e_2 // " e " // e, PRINT_ALWAYS)
endif
    if (save_energies .and. len_trim(calc_energy) > 0) then
       call set_param_value(at,'EVB1_'//trim(calc_energy),my_e_1)
       call set_param_value(at,'EVB2_'//trim(calc_energy),my_e_2)
    endif

    if (save_forces .and. len_trim(calc_force) > 0) then
      call add_property(at, 'EVB1_'//trim(calc_force), my_f_1)
      call add_property(at, 'EVB2_'//trim(calc_force), my_f_2)
    endif

    if (allocated(my_f_1)) deallocate(my_f_1)
    if (allocated(my_f_2)) deallocate(my_f_2)
    if (allocated(de_offdiag_dr)) deallocate(de_offdiag_dr)

  end subroutine Potential_EVB_Calc

  recursive function Potential_EVB_Cutoff(this)
    type(Potential_EVB), intent(in) :: this
    real(dp) :: potential_EVB_cutoff

    if(associated(this%pot1)) then
       potential_EVB_cutoff = cutoff(this%pot1)
    else
       potential_EVB_cutoff = 0.0_dp
    endif

  end function Potential_EVB_Cutoff

