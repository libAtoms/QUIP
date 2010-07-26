
  !*************************************************************************
  !*
  !*  Potential_Sum routines
  !*
  !*************************************************************************

  subroutine Potential_Sum_Initialise(this, args_str, pot1, pot2, mpi, error)
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

  subroutine Potential_Sum_Finalise(this)
    type(Potential_Sum), intent(inout) :: this
    
    nullify(this%pot1)
    nullify(this%pot2)

  end subroutine Potential_Sum_Finalise

  subroutine Potential_Sum_Print(this, file)
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

  subroutine Potential_Sum_Calc(this, at, e, local_e, f, df, virial, args_str, error)
    type(Potential_Sum), intent(inout) :: this
    type(Atoms), intent(inout) :: at
    real(dp), intent(out), optional :: e
    real(dp), intent(out), optional :: local_e(:)
    real(dp), intent(out), optional :: f(:,:)
    real(dp), intent(out), optional :: df(:,:)
    real(dp), intent(out), optional :: virial(3,3)
    character(*), intent(in), optional :: args_str
    integer, intent(out), optional :: error

    real(dp) :: my_e_1, my_e_2
    real(dp), allocatable :: my_local_e_1(:)
    real(dp), allocatable :: my_f_1(:,:)
    real(dp), allocatable :: my_df_1(:,:)
    real(dp) :: my_virial_1(3,3)

    INIT_ERROR(error)

    if (present(local_e)) allocate(my_local_e_1(size(local_e)))
    if (present(f)) allocate(my_f_1(size(f, 1),size(f, 2)))
    if (present(df)) allocate(my_df_1(size(df, 1),size(df, 2)))

    call calc(this%pot1, at, e, local_e, f, df, virial, args_str, error)
    PASS_ERROR(error)

    if (present(e)) my_e_1 = e
    if (present(local_e)) my_local_e_1 = local_e
    if (present(f)) my_f_1 = f
    if (present(df)) my_df_1 = df
    if (present(virial)) my_virial_1 = virial

    call calc(this%pot2, at, e, local_e, f, df, virial, args_str, error)
    PASS_ERROR(error)

    if (present(e)) e = my_e_1 + e
    if (present(local_e)) local_e = my_local_e_1 + local_e
    if (present(f)) f = my_f_1 + f
    if (present(df)) df = my_df_1 + df
    if (present(virial)) virial = my_virial_1 + virial

    if (present(local_e)) deallocate(my_local_e_1)
    if (present(f)) deallocate(my_f_1)
    if (present(df)) deallocate(my_df_1)
  
  end subroutine Potential_Sum_Calc

  function Potential_Sum_Cutoff(this)
    type(Potential_Sum), intent(in) :: this
    real(dp) :: potential_sum_cutoff

    if(associated(this%pot1) .and. associated(this%pot2)) then
       potential_sum_cutoff = max(cutoff(this%pot1), cutoff(this%pot2))
    else
       potential_sum_cutoff = 0.0_dp
    endif

  end function Potential_Sum_Cutoff

