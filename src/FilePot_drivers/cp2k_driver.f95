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
!X cp2k_driver_template program
!X
!% filepot program for CP2K driver using input file template
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#include "error.inc"
program cp2k_driver_template
use libatoms_module
use cp2k_driver_module

implicit none

    integer, parameter :: CP2K_LINE_LENGTH = 1024 !Max line length to be printed into the CP2K input files

    type(Atoms)                   :: my_atoms
    real(dp)                      :: energy
    real(dp), allocatable         :: f0(:,:)
    real(dp), pointer             :: forces_p(:,:)
    type(CInoutput)                :: xyz_io

    character(len=STRING_LENGTH)  :: infile, outfile
    character(len=STRING_LENGTH) :: args_str
    character(len=STRING_LENGTH)  :: arg
    integer :: i, index_insert

    integer :: error = ERROR_NONE

    call system_initialise(verbosity=PRINT_SILENT,enable_timing=.true.)
    call verbosity_push(PRINT_NORMAL)
    call system_timer('cp2k_driver_template')

    if (cmd_arg_count() < 2) &
      call system_abort("Usage: cp2k_driver_template infile outfile [ other_cli_arg1=v1 other_cli_arg2=v2 ...]")

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

    call read(my_atoms,infile, error=error)
    if (error == ERROR_NONE) then
       allocate(f0(3,my_atoms%N))
    else
       ! continue: do_cp2k_calc() should fail, but print first out help message if requested
       CLEAR_ERROR(error)
    endif

    !call CP2K
    call do_cp2k_calc(at=my_atoms, f=f0, e=energy, args_str=trim(args_str), error=error)
    HANDLE_ERROR(error)

    !momentum conservation
!   call sum0(f0)

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
    call system_timer('cp2k_driver_template')
    mainlog%prefix=""
    call verbosity_push(PRINT_SILENT)
    call system_finalise

end program cp2k_driver_template
