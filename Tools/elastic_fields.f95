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

program elastic

  use libAtoms_module
  use elasticity_module, only: elastic_fields
 
  implicit none

  type(Atoms) :: at
  type(CInoutput) :: infile, outfile
  type(Dictionary) :: params
  integer :: i, n, iostat, nargs
  character(len=2048) :: comment, arg1, arg2, missing_params
  real(dp) :: a, C11, C12, C44, cutoff, nneightol

  call system_initialise(PRINT_SILENT)
  call verbosity_push(PRINT_NORMAL)

  call initialise(params)

  call param_register(params, 'a', PARAM_MANDATORY, a, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'C11', PARAM_MANDATORY, C11, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'C12', PARAM_MANDATORY, C12, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'C44', PARAM_MANDATORY, C44, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'cutoff', '5.0', cutoff, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'nneightol', '1.2', nneightol, help_string="No help yet.  This source file was $LastChangedBy$")

  nargs = cmd_arg_count()

  if (nargs < 2) then
     call print('Usage: elastic_fields <infile> <outfile> [params]')
     call print('')
     call print('Parameters and default values are:')
     call param_print(params)
     call verbosity_push(PRINT_SILENT)
     call system_finalise()
     stop
  end if

  call get_cmd_arg(1, arg1)
  call get_cmd_arg(2, arg2)

  call initialise(infile, arg1)
  call initialise(outfile, arg2, action=OUTPUT)


  if (nargs > 2) then
     if (.not. param_read_args(params, (/ (i, i=3,nargs ) /))) &
          call system_abort('Error parsing command line')
  end if

  if (.not. param_check(params, missing_params)) &
       call system_abort('Missing mandatory parameters: '//missing_params)

  call print('Parameters:')
  call param_print(params)
  call print('')


  ! Loop over all frames in input file
  do n=0,infile%n_frame-1
     call read(infile, at, frame=n)

     call print('Frame '//n)

     call set_cutoff(at, cutoff)
     at%nneightol = nneightol

     call elastic_fields(at, a, C11, C12, C44)
     call write(outfile, at)
  end do
  call print('Done '//n//' frames!')

  call verbosity_push(PRINT_SILENT)
  call system_finalise

end program elastic
