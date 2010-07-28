
  !*************************************************************************
  !*
  !*  Potential_EVB routines
  !*
  !*************************************************************************

  subroutine Potential_EVB_Initialise(this, args_str, pot1, mpi, error)
    type(Potential_EVB), intent(inout) :: this
    character(len=*), intent(in) :: args_str
    type(Potential), intent(in), target :: pot1
    type(MPI_Context), intent(in), optional :: mpi
    integer, intent(out), optional :: error

    type(Dictionary) :: params

    INIT_ERROR(error)

    call finalise(this)

    !the two potentials are the same
    this%pot1 => pot1
    this%pot2 => pot1

    call initialise(params)
    call param_register(params, 'evb_step', '-1', this%evb_step)
    call param_register(params, 'mm_args_str', '', this%mm_args_str)
    call param_register(params, 'topology_suffix1', '', this%topology_suffix1)
    call param_register(params, 'topology_suffix2', '', this%topology_suffix2)
    call param_register(params, 'form_bond', '0 0', this%form_bond)
    call param_register(params, 'break_bond', '0 0', this%break_bond)
    !call param_register(params, 'energy_offset', '0.0',this%energy_offset)
    if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='Potential_EVB_initialise args_str')) then
       RAISE_ERROR('Potential_FM_initialise failed to parse args_str="'//trim(args_str)//'"', error)
    endif
    call finalise(params)

    if (present(mpi)) this%mpi = mpi

  end subroutine Potential_EVB_Initialise

  subroutine Potential_EVB_Finalise(this)
    type(Potential_EVB), intent(inout) :: this
    
    nullify(this%pot1)
    nullify(this%pot2)

    this%mm_args_str = ""
    this%topology_suffix1 = ""
    this%topology_suffix2 = ""
    this%form_bond(1:2) = 0
    this%break_bond(1:2) = 0

  end subroutine Potential_EVB_Finalise

  subroutine Potential_EVB_Print(this, file)
    type(Potential_EVB), intent(inout) :: this
    type(Inoutput), intent(inout), optional :: file

    call print('Potential_EVB:', file=file)
    call print('  evb_step'//trim(this%evb_step))
    call print('  mm_args_str'//trim(this%mm_args_str))
    call print('  topology_suffix1'//trim(this%topology_suffix1))
    call print('  topology_suffix2'//trim(this%topology_suffix2))
    call print('  pot1 - form_bond:  '//this%form_bond(1:2))
    call print('  pot1 - break_bond: '//this%break_bond(1:2))
    call print('  pot2 - form_bond:  '//this%break_bond(1:2))
    call print('  pot2 - break_bond: '//this%form_bond(1:2))
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

  end subroutine Potential_EVB_Print

  subroutine Potential_EVB_Calc(this, at, e, local_e, f, df, virial, args_str, error)
    type(Potential_EVB), intent(inout) :: this
    type(Atoms), intent(inout) :: at
    real(dp), intent(out), optional :: e
    real(dp), intent(out), optional :: local_e(:)
    real(dp), intent(out), optional :: f(:,:)
    real(dp), intent(out), optional :: df(:,:)
    real(dp), intent(out), optional :: virial(3,3)
    character(*), intent(in), optional :: args_str
    integer, intent(out), optional :: error

    type(Dictionary) :: params
    character(FIELD_LENGTH) :: mm_args_str
    character(FIELD_LENGTH) :: topology_suffix1, topology_suffix2
    character(FIELD_LENGTH) :: evb_step
    real(dp) :: my_e_1
    real(dp), allocatable :: my_f_1(:,:)
    !real(dp) :: energy_offset
    logical :: skip_form_bond
    logical :: skip_break_bond
    integer :: form_bond(2)
    integer :: break_bond(2)
    logical, save :: first_call = .true.
    character(FIELD_LENGTH) :: psf_print

    INIT_ERROR(error)

    !read args_str
    call initialise(params)
    call param_register(params, 'evb_step', ''//this%evb_step, evb_step)
    call param_register(params, 'mm_args_str', ''//this%mm_args_str, mm_args_str)
    call param_register(params, 'topology_suffix1', ''//this%topology_suffix1, topology_suffix1)
    call param_register(params, 'topology_suffix2', ''//this%topology_suffix2, topology_suffix2)
    call param_register(params, 'form_bond', ''//this%form_bond, form_bond)
    call param_register(params, 'break_bond', ''//this%break_bond, break_bond)
    !call param_register(params, 'energy_offset', ''//this%energy_offset, energy_offset)
    if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='Potential_EVB_initialise args_str')) then
       RAISE_ERROR('Potential_FM_initialise failed to parse args_str="'//trim(args_str)//'"', error)
    endif
    call finalise(params)

    !!!!!
    !at the moment it only works with existing topologies
    call initialise(params)
      psf_print = ''
      call param_register(params, 'PSF_print', 'NO_PSF', psf_print)
    if (.not. param_read_line(params, mm_args_str, ignore_unknown=.true.,task='Potential_EVB_initialise mm1_args_str')) then
       RAISE_ERROR('Potential_FM_initialise failed to parse mm1_args_str="'//trim(mm_args_str)//'"', error)
    endif
    call finalise(params)
    if ( (trim(psf_print) /= 'USE_EXISTING_PSF')) &
      call system_abort("EVB calculation has been only implemented for using existing topology files yet.  Generate 2 topology files with the resonance structures, and try with them.")
    !!!!!

    !check arguments
    if (present(local_e) .or. present(virial) .or. present(df)) then
         RAISE_ERROR('Potential_EVB_calc: supports only energy and forces, not virial, energy or local_e', error)
    endif

    !form_bond
    skip_form_bond = .false.
    if(any(form_bond<1) .or. any(form_bond>at%N)) then
       !check whether all 0 (skip)
       if (all(form_bond==0)) then
          skip_form_bond = .true.
       else
          RAISE_ERROR('Potential_EVB_initialise form_bond is out of range 1--'//at%N//': '//form_bond, error)
       endif
    endif

    !break_bond
    skip_break_bond = .false.
    if(any(break_bond<1) .or. any(break_bond>at%N)) then
       !check whether all 0 (skip)
       if (all(break_bond==0)) then
          skip_break_bond = .true.
       else
          RAISE_ERROR('Potential_EVB_initialise break_bond is out of range 1--'//at%N//': '//break_bond, error)
       endif
    endif


    if (present(f)) allocate(my_f_1(size(f, 1),size(f, 2)))

    !set breaking and forming bonds for topology1
    call remove_value(at%params, 'form_bond')
    call remove_value(at%params, 'break_bond')
    if (.not.skip_form_bond) call set_value(at%params,'form_bond',form_bond(1:2))
    if (.not.skip_break_bond) call set_value(at%params,'break_bond',break_bond(1:2))
    call calc(this%pot1, at, e, local_e, f, df, virial, trim(mm_args_str)//" topology_suffix="//trim(topology_suffix1), error)
    PASS_ERROR(error)

    if (present(e)) my_e_1 = e
    if (present(f)) my_f_1 = f

    !set breaking and forming bonds for topology2
    call remove_value(at%params, 'form_bond')
    call remove_value(at%params, 'break_bond')
    if (.not.skip_break_bond) call set_value(at%params,'form_bond',break_bond(1:2))
    if (.not.skip_form_bond) call set_value(at%params,'break_bond',form_bond(1:2))
    call calc(this%pot2, at, e, local_e, f, df, virial, trim(mm_args_str)//" topology_suffix="//trim(topology_suffix2), error)
    PASS_ERROR(error)

    !print energies here, cannot pass them back to main
    if (first_call) then
       if (len_trim(evb_step)/=0) then
          call print('EVB | evb_step   energy1(eV)        energy2(eV)        E_GAP(eV)',PRINT_ALWAYS)
          call print('EVB | '//trim(evb_step)//' '//my_e_1//' '//e//' '//(my_e_1-e),PRINT_ALWAYS)
       endif
       first_call = .false.
    else
       if (len_trim(evb_step)/=0) then
          call print('EVB | '//trim(evb_step)//' '//my_e_1//' '//e//' '//(my_e_1-e),PRINT_ALWAYS)
       endif
    endif

    !E_GAP energy
    if (present(e)) e = my_e_1 - e
    !EVB force
    if (present(f)) f = my_f_1 - f

    if (present(f)) deallocate(my_f_1)
  
  end subroutine Potential_EVB_Calc

  function Potential_EVB_Cutoff(this)
    type(Potential_EVB), intent(in) :: this
    real(dp) :: potential_EVB_cutoff

    if(associated(this%pot1) .and. associated(this%pot2)) then
       potential_EVB_cutoff = max(cutoff(this%pot1), cutoff(this%pot2))
    else
       potential_EVB_cutoff = 0.0_dp
    endif

  end function Potential_EVB_Cutoff

