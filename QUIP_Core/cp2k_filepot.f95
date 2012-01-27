! H0 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! H0 X
! H0 X   libAtoms+QUIP: atomistic simulation library
! H0 X
! H0 X   Portions of this code were written by
! H0 X     Albert Bartok-Partay, Silvia Cereda, Gabor Csanyi, James Kermode,
! H0 X     Ivan Solt, Wojciech Szlachta, Csilla Varnai, Steven Winfield.
! H0 X
! H0 X   Copyright 2006-2010.
! H0 X
! H0 X   These portions of the source code are released under the GNU General
! H0 X   Public License, version 2, http://www.gnu.org/copyleft/gpl.html
! H0 X
! H0 X   If you would like to license the source code under different terms,
! H0 X   please contact Gabor Csanyi, gabor@csanyi.net
! H0 X
! H0 X   Portions of this code were written by Noam Bernstein as part of
! H0 X   his employment for the U.S. Government, and are not subject
! H0 X   to copyright in the USA.
! H0 X
! H0 X
! H0 X   When using this software, please cite the following reference:
! H0 X
! H0 X   http://www.libatoms.org
! H0 X
! H0 X  Additional contributions by
! H0 X    Alessio Comisso, Chiara Gattinoni, and Gianpietro Moras
! H0 X
! H0 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X cp2k_filepot_template program
!X
!% filepot program for CP2K driver using input file template
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

program cp2k_filepot_template

  use libatoms_module
  use cp2k_driver_template_module

  implicit none

    type(Atoms)                   :: my_atoms
    real(dp)                      :: energy
    real(dp), allocatable         :: f0(:,:)
    real(dp), pointer             :: forces_p(:,:)
    type(CInoutput)                :: xyz_io

    character(len=STRING_LENGTH)  :: infile, outfile
    character(len=STRING_LENGTH) :: args_str
    character(len=STRING_LENGTH)  :: arg
    integer :: i, index_insert

    call system_initialise(verbosity=PRINT_SILENT,enable_timing=.true.)
    call verbosity_push(PRINT_NORMAL)
    mainlog%prefix="CP2K_FILEPOT"
    call system_timer('cp2k_filepot_template')

    if (cmd_arg_count() < 2) &
      call system_abort("Usage: cp2k_filepot_template infile outfile [ other_cli_arg1=v1 other_cli_arg2=v2 ...]")

    call get_cmd_arg(1, infile)
    call get_cmd_arg(2, outfile)
    args_str = ""
    if (cmd_arg_count() > 2) then
      do i=3, cmd_arg_count()
	call get_cmd_arg(i, arg)
        !add {} if there is space in the arg
        if (index(trim(arg)," ").ne.0) then
            index_insert = index(trim(arg),"=")
            arg(index_insert+1:len_trim(arg)+2) = "{"//arg(index_insert+1:len_trim(arg))//"}"
            call print('arg: '//trim(arg),PRINT_SILENT)
        endif
	args_str = trim(args_str) // " " // trim(arg)
      end do
    endif

    call read(my_atoms,infile)

    call print('cp2k_filepot_template args_str '//trim(args_str), PRINT_ALWAYS)

    !call CP2K
    allocate(f0(3,my_atoms%N))
    call do_cp2k_calc(at=my_atoms, f=f0, e=energy, args_str=trim(args_str))
    !momentum conservation
    call sum0(f0)

    !write energy and forces
    call set_value(my_atoms%params,'energy',energy)
    call add_property(my_atoms,'force', 0.0_dp, n_cols=3)
    if (.not. assign_pointer(my_atoms, 'force', forces_p)) then
       call system_abort("filepot_read_output needed forces, but couldn't find force in my_atoms ")
    endif
    forces_p = f0

    call initialise(xyz_io,outfile,action=OUTPUT)
    call write(my_atoms,xyz_io,real_format='%20.13f')
    call finalise(xyz_io)

    deallocate(f0)
    call finalise(my_atoms)
    call system_timer('cp2k_filepot_template')
    mainlog%prefix=""
    call verbosity_push(PRINT_SILENT)
    call system_finalise

contains

  !momentum conservation
  !   weighing function: 1 (simply subtract sumF/n)
  !?!   weighing function: m (keeping same acceleration on the atoms)
  subroutine sum0(force)

    real(dp), dimension(:,:), intent(inout) :: force
    integer   :: i
    real(dp)  :: sumF(3)

    do i = 1, size(force,2)
       sumF(1) = sum(force(1,1:size(force,2)))
       sumF(2) = sum(force(2,1:size(force,2)))
       sumF(3) = sum(force(3,1:size(force,2)))
    enddo

    if ((sumF(1).feq.0.0_dp).and.(sumF(2).feq.0.0_dp).and.(sumF(3).feq.0.0_dp)) then
       call print('cp2k_driver: Sum of the forces are zero.')
       return
    endif

    call print('cp2k_driver: Sum of the forces was '//sumF(1:3))
    sumF = sumF / size(force,2)

    do i = 1, size(force,2)
       force(1:3,i) = force(1:3,i) - sumF(1:3)
    enddo

    do i = 1, size(force,2)
       sumF(1) = sum(force(1,1:size(force,2)))
       sumF(2) = sum(force(2,1:size(force,2)))
       sumF(3) = sum(force(3,1:size(force,2)))
    enddo
    call print('cp2k_driver: Sum of the forces after mom.cons.: '//sumF(1:3))

  end subroutine sum0

end program cp2k_filepot_template
