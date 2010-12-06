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

program eval
  
  use libAtoms_module
  use potential_module
  use libatoms_misc_utils_module
  use descriptors_module
  
  implicit none
  
  type(Atoms) at
  type(CInOutput) infile
  real(dp) :: rOH(8), rHH(6)
  integer error

  call system_initialise()
  
  
  call initialise(infile, 'dimers_9.xyz')
  ! main loop over frames
  do 
     call read(at, infile, error=error)
     if (error /= 0) then
        if (error == ERROR_IO_EOF) then
           exit
        endif
     endif
     
     call water_dimer(at, (/1,2,3/), (/4,5,6/), 10.0_dp, rOHout = rOH, rHHout = rHH)
     call add_array(at%params, 'rOH', rOH, 8)
     call add_array(at%params, 'rHH', rHH, 6)
     call write(at,'stdout', prefix='AT')

     call finalise(at)
     
  enddo

  call finalise(infile)
  call system_finalise()

end program
