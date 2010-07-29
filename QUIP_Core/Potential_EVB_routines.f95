
  !*************************************************************************
  !*
  !*  Potential_EVB routines
  !*
  !*************************************************************************

  subroutine Potential_EVB_Initialise(this, args_str, pot1, mpi, error)
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
    call param_register(params, 'mm_args_str', '', this%mm_args_str)
    call param_register(params, 'topology_suffix1', '', this%topology_suffix1)
    call param_register(params, 'topology_suffix2', '', this%topology_suffix2)
    call param_register(params, 'form_bond', '0 0', this%form_bond)
    call param_register(params, 'break_bond', '0 0', this%break_bond)
    !call param_register(params, 'energy_offset', '0.0' , this%energy_offset)
    call param_register(params, 'offdiagonal_A12', '0.0', this%offdiagonal_A12)
    call param_register(params, 'offdiagonal_mu12', '0.0', this%offdiagonal_mu12)
    call param_register(params, 'save_forces', 'T', this%save_forces)
    call param_register(params, 'save_energies', 'T', this%save_energies)
    if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='Potential_EVB_initialise args_str')) then
       RAISE_ERROR('Potential_EVB_initialise failed to parse args_str="'//trim(args_str)//'"', error)
    endif
    call finalise(params)

    if (present(mpi)) this%mpi = mpi

  end subroutine Potential_EVB_Initialise

  subroutine Potential_EVB_Finalise(this)
    type(Potential_EVB), intent(inout) :: this
    
    nullify(this%pot1)

    this%mm_args_str = ""
    this%topology_suffix1 = ""
    this%topology_suffix2 = ""
    this%form_bond(1:2) = 0
    this%break_bond(1:2) = 0
    !this%energy_offset = 0._dp
    this%offdiagonal_A12 = 0._dp
    this%offdiagonal_mu12 = 0._dp
    this%save_forces = .false.
    this%save_energies = .false.

  end subroutine Potential_EVB_Finalise

  subroutine Potential_EVB_Print(this, file)
    type(Potential_EVB), intent(inout)           :: this
    type(Inoutput),      intent(inout), optional :: file

    call print('Potential_EVB:', file=file)
    call print('  mm_args_str'//trim(this%mm_args_str), file=file)
    call print('  topology_suffix1'//trim(this%topology_suffix1), file=file)
    call print('  topology_suffix2'//trim(this%topology_suffix2), file=file)
    call print('  evb1-form and evb2-break bond:  '//this%form_bond(1:2), file=file)
    call print('  evb1-break and evb2-form bond: '//this%break_bond(1:2), file=file)
    !call print('  energy offset: '//this%energy_offset, file=file)
    call print('  offdiagonal A12(pre-exponent factor): '//this%offdiagonal_A12, file=file)
    call print('  offdiagonal mu12(exponent factor): '//this%offdiagonal_mu12, file=file)
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

  subroutine Potential_EVB_Calc(this, at, e, local_e, f, df, virial, args_str, error)
    type(Potential_EVB), intent(inout)          :: this
    type(Atoms),         intent(inout)          :: at
    real(dp),            intent(out), optional  :: e
    real(dp),            intent(out), optional  :: local_e(:)
    real(dp),            intent(out), optional  :: f(:,:)
    real(dp),            intent(out), optional  :: df(:,:)
    real(dp),            intent(out), optional  :: virial(3,3)
    character(*),        intent(in),  optional  :: args_str
    integer,             intent(out), optional  :: error

    type(Dictionary)        :: params
    character(FIELD_LENGTH) :: mm_args_str
    character(FIELD_LENGTH) :: topology_suffix1, topology_suffix2
    !real(dp)                :: energy_offset
    integer                 :: form_bond(2), break_bond(2)
    logical                 :: skip_form_bond, skip_break_bond

    logical                 :: save_energies, save_forces
    real(dp)                :: my_e_1, my_e_2, e_offdiag
    real(dp), allocatable   :: my_f_1(:,:), my_f_2(:,:), delta_eoffdiag(:,:)
    real(dp)                :: offdiagonal_A12, offdiagonal_mu12, &
                               rab, d_rab_dx(3)
    real(dp), pointer       :: force_ptr(:,:)
    logical                 :: no_coupling, dummy

    character(FIELD_LENGTH) :: psf_print

    INIT_ERROR(error)

    !read args_str
    call initialise(params)
    call param_register(params, 'mm_args_str', ''//this%mm_args_str, mm_args_str)
    call param_register(params, 'topology_suffix1', ''//this%topology_suffix1, topology_suffix1)
    call param_register(params, 'topology_suffix2', ''//this%topology_suffix2, topology_suffix2)
    call param_register(params, 'form_bond', ''//this%form_bond, form_bond)
    call param_register(params, 'break_bond', ''//this%break_bond, break_bond)
    !call param_register(params, 'energy_offset', ''//this%energy_offset, energy_offset)
    call param_register(params, 'offdiagonal_A12', ''//this%offdiagonal_A12, offdiagonal_A12)
    call param_register(params, 'offdiagonal_mu12', ''//this%offdiagonal_mu12, offdiagonal_mu12)
    call param_register(params, 'save_forces', ''//this%save_forces, save_forces)
    call param_register(params, 'save_energies', ''//this%save_energies, save_energies)
    if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='Potential_EVB_initialise args_str')) then
       RAISE_ERROR('Potential_EVB_calc failed to parse args_str="'//trim(args_str)//'"', error)
    endif
    call finalise(params)

    !!!!!
    !at the moment it only works with existing topologies
    call initialise(params)
      psf_print = ''
      call param_register(params, 'PSF_print', 'NO_PSF', psf_print)
    if (.not. param_read_line(params, mm_args_str, ignore_unknown=.true.,task='Potential_EVB_calc mm_args_str')) then
       RAISE_ERROR('Potential_EVB_calc failed to parse mm_args_str="'//trim(mm_args_str)//'"', error)
    endif
    call finalise(params)
    if ( (trim(psf_print) /= 'USE_EXISTING_PSF')) &
      call system_abort("EVB calculation has been only implemented for using existing topology files yet.  Generate 2 topology files with the resonance structures, and try with them.")
    !!!!!


    !CHECK ARGUMENTS

    if (present(local_e) .or. present(virial) .or. present(df)) then
         RAISE_ERROR('Potential_EVB_calc: supports only energy and forces, not virial, energy or local_e', error)
    endif

    !coupling parameters
    if (offdiagonal_A12 .feq. 0._dp) then
       call print('WARNING! Offdiagonal A12 is set to 0. No coupling between resonance states.')
       no_coupling = .true.
    else
       if (offdiagonal_A12 < 0._dp .or. offdiagonal_mu12 < 0._dp) &
          RAISE_ERROR('Potential_EVB_calc offdiagonal parameters must be positive or 0 for no coupling. '// &
                      'Got offdiagonal_A12: '//offdiagonal_A12//' and offdiagonal_mu12: '//offdiagonal_mu12, error)

       no_coupling = .false.
    endif

    !form_bond
    skip_form_bond = .false.
    if(any(form_bond<1) .or. any(form_bond>at%N)) then
       !check whether all 0 (skip)
       if (all(form_bond==0)) then
          skip_form_bond = .true.
       else
          RAISE_ERROR('Potential_EVB_calc form_bond is out of range 1--'//at%N//': '//form_bond, error)
       endif
    endif

    !break_bond
    skip_break_bond = .false.
    if(any(break_bond<1) .or. any(break_bond>at%N)) then
       !check whether all 0 (skip)
       if (all(break_bond==0)) then
          skip_break_bond = .true.
       else
          RAISE_ERROR('Potential_EVB_calc break_bond is out of range 1--'//at%N//': '//break_bond, error)
       endif
    endif

    if (skip_form_bond .and. skip_break_bond) then
       RAISE_ERROR('Potential_EVB_calc no bonds to form neither to break. ', error)
    endif


    !CALCULATE E,F WITH DIFFERENT TOPOLOGIES

    if (present(f)) allocate(my_f_1(size(f, 1),size(f, 2)))
    if (present(f)) allocate(my_f_2(size(f, 1),size(f, 2)))
    if (present(f)) allocate(delta_eoffdiag(size(f, 1),size(f, 2)))

    !set breaking and forming bonds for topology1
    call remove_value(at%params, 'form_bond')
    call remove_value(at%params, 'break_bond')
    if (.not.skip_form_bond) call set_value(at%params,'form_bond',form_bond(1:2))
    if (.not.skip_break_bond) call set_value(at%params,'break_bond',break_bond(1:2))
    !calc with topology1
    call calc(this%pot1, at, e, local_e, f, df, virial, trim(mm_args_str)//" topology_suffix="//trim(topology_suffix1), error)
    PASS_ERROR(error)

    if (present(e)) my_e_1 = e
    if (present(f)) my_f_1 = f

    !set breaking and forming bonds for topology2
    call remove_value(at%params, 'form_bond')
    call remove_value(at%params, 'break_bond')
    if (.not.skip_break_bond) call set_value(at%params,'form_bond',break_bond(1:2))
    if (.not.skip_form_bond) call set_value(at%params,'break_bond',form_bond(1:2))
    !calc with topology2
    call calc(this%pot1, at, e, local_e, f, df, virial, trim(mm_args_str)//" topology_suffix="//trim(topology_suffix2), error)
    PASS_ERROR(error)

    if (present(e)) my_e_2 = e
    if (present(f)) my_f_2 = f


    !CALCULATE EVB ENERGY AND FORCES

    !distance or distance difference
    if (no_coupling) then
       !take the E, F of the resonance state with the smaller E
       if (my_e_1 < my_e_2) then
          if (present(e)) e = my_e_1
          if (present(f)) f = my_f_1
       else
          if (present(e)) e = my_e_2
          if (present(f)) f = my_f_2
       endif
    else
       !calculate coupling terms
       rab = 0._dp
       if (.not. skip_form_bond) rab = rab + distance_min_image(at,form_bond(1),form_bond(2))
       if (.not. skip_form_bond) rab = rab - distance_min_image(at,break_bond(1),break_bond(2))
   
       if (present(e) .or. present(f)) then
          e_offdiag = offdiagonal_A12 * exp(-offdiagonal_mu12 * abs(rab))
       endif

       !energy
       if (present(e)) then
          e = 0.5_dp * (my_e_1 + my_e_2) - 0.5_dp * sqrt((my_e_1 - my_e_2)**2._dp + 4._dp*e_offdiag)
       endif

       !force
       if (present(f)) then
          !force coupling term
          delta_eoffdiag(1:size(f, 1),1:size(f, 2)) = 0._dp
          if (.not. skip_form_bond) then
             d_rab_dx = diff_min_image(at,form_bond(1),form_bond(2))/distance_min_image(at,form_bond(1),form_bond(2))
             if (rab < 0) d_rab_dx = - d_rab_dx
             delta_eoffdiag(1:size(f, 1),form_bond(1)) = delta_eoffdiag(1:size(f, 1),form_bond(1)) + 4._dp * e_offdiag * (offdiagonal_mu12) * d_rab_dx(1:3)
             delta_eoffdiag(1:size(f, 1),form_bond(2)) = delta_eoffdiag(1:size(f, 1),form_bond(2)) - 4._dp * e_offdiag * (offdiagonal_mu12) * d_rab_dx(1:3)
          endif
          if (.not. skip_break_bond) then
             d_rab_dx = diff_min_image(at,break_bond(1),break_bond(2))/distance_min_image(at,break_bond(1),break_bond(2))
             if (rab > 0) d_rab_dx = - d_rab_dx
             delta_eoffdiag(1:size(f, 1),break_bond(1)) = delta_eoffdiag(1:size(f, 1),break_bond(1)) + 4._dp * e_offdiag * (offdiagonal_mu12) * d_rab_dx(1:3)
             delta_eoffdiag(1:size(f, 1),break_bond(2)) = delta_eoffdiag(1:size(f, 1),break_bond(2)) - 4._dp * e_offdiag * (offdiagonal_mu12) * d_rab_dx(1:3)
          endif
          !force
          f = 0.5_dp * (my_f_1(1:size(f, 1),1:size(f, 2)) + my_f_2(1:size(f, 1),1:size(f, 2))) - &
              0.5_dp * (my_e_1 - my_e_2)*(my_f_1(1:size(f, 1),1:size(f, 2)) - my_f_2(1:size(f, 1),1:size(f, 2)) - 4._dp * delta_eoffdiag) / sqrt((my_e_1 - my_e_2)**2.0_dp + 4._dp*e_offdiag**2._dp)
       endif
    endif

    !SAVE E,F IF NEEDED

    if (save_energies) then
       call set_value(at%params,'EVB1_energy',my_e_1)
       call set_value(at%params,'EVB2_energy',my_e_2)
    endif

    if (save_forces) then
       !only save forces, if requested and calculated
       if (present(f)) then
          if (.not. has_property(at, 'EVB1_force')) &
               call add_property(at, 'EVB1_force', 0.0_dp, n_cols=3)
          dummy = assign_pointer(at, 'EVB1_force', force_ptr)
          call print(size(force_ptr,1)//" "//size(force_ptr,2))
          call print(size(my_f_1,1)//" "//size(my_f_1,2))
          force_ptr = my_f_1
   
          if (.not. has_property(at, 'EVB2_force')) &
               call add_property(at, 'EVB2_force', 0.0_dp, n_cols=3)
          dummy = assign_pointer(at, 'EVB2_force', force_ptr)
          force_ptr = my_f_2

          if (.not. has_property(at, 'force')) &
               call add_property(at, 'force', 0.0_dp, n_cols=3)
          dummy = assign_pointer(at, 'force', force_ptr)
          force_ptr = f
       else
          call print('WARNING! Potential_EVB_calc did not calculate EVB forces, so it cannot and will not save them.',PRINT_ALWAYS)
       endif
    endif


    if (present(f)) deallocate(my_f_1)
    if (present(f)) deallocate(my_f_2)
    if (present(f)) deallocate(delta_eoffdiag)
  
  end subroutine Potential_EVB_Calc

  function Potential_EVB_Cutoff(this)
    type(Potential_EVB), intent(in) :: this
    real(dp) :: potential_EVB_cutoff

    if(associated(this%pot1)) then
       potential_EVB_cutoff = cutoff(this%pot1)
    else
       potential_EVB_cutoff = 0.0_dp
    endif

  end function Potential_EVB_Cutoff

