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

!diffusion constant calculating helper program.
!prints the average square displacement as the function of the frames.
!The diffusivity can be calculated by fitting 4*D*t to <r^2> (t).
program diffusion

use libatoms_module

  implicit none

 
    type(Atoms)                           :: structure, reference
    type(CInoutput)                       :: xyzfile
    type(Inoutput)                        :: datafile
    real(dp)                              :: ave_r2, one_r2
    integer                               :: frame_count, frames_processed
    integer                               :: i
    integer                               :: error

    !Input
    type(Dictionary)                      :: params_in
    character(STRING_LENGTH)               :: xyzfilename, datafilename
    integer                               :: from, to
    integer                               :: IO_Rate
    integer                               :: one_atom

    call system_initialise(PRINT_SILENT)
    call verbosity_push(PRINT_NORMAL)

    call initialise(params_in)
    call param_register(params_in, 'xyzfile', param_mandatory, xyzfilename, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params_in, 'datafile', param_mandatory, datafilename, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params_in, 'from', '1', from, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params_in, 'to', '0', to, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params_in, 'one_atom', '0', one_atom, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params_in, 'IO_Rate', '1', IO_Rate, help_string="No help yet.  This source file was $LastChangedBy$")
    if (.not. param_read_args(params_in)) then
       call print('diffusion xyzfile datafile [from=1] [to=0] [one_atom] [IO_Rate]')
       call system_abort('could not parse argument line')
    end if
    call finalise(params_in)

    if (from<1) call system_abort('from must be a positive intger.')
    ! Read the input and output filename and open them
    call initialise(xyzfile,xyzfilename,action=INPUT)
    call initialise(datafile,datafilename,action=OUTPUT)
    call print('#average displacement square of all the atoms',file=datafile)

    call print(' Input file: '//trim(xyzfilename))
    call print('Output file: '//trim(datafilename))
    call print('From step '//from//' to step '//to)
    call print('Also check atom '//one_atom)
    call print('')

    call print('Reading data...')

    call read(structure, xyzfile, error=error)
    !reference = structure

    frame_count = 0
    frames_processed = 0

    do
     
       if (error/=0) exit
  
       frame_count = frame_count + 1
  
       if ((frame_count.gt.to) .and. (to.gt.0)) exit
  
       write(mainlog%unit,'(a,a,i0,$)') achar(13),'Frame ',frame_count
  
       if (frame_count.eq.from) then
          call print('Reference structure is at step '//frame_count)
          reference = structure
          call map_into_cell(reference)
       endif
     
       if (frame_count.ge.from) then
          frames_processed = frames_processed + 1
  
          call map_into_cell(structure)
          one_r2 = 0._dp
          ave_r2 = 0._dp
          do i = 1, structure%N
             ave_r2 = ave_r2 + normsq(realpos(reference,i) - realpos(structure,i))
          enddo
          ave_r2 = ave_r2 / structure%N
          if (one_atom.ne.0) then
             one_r2 = normsq(realpos(reference,one_atom) - realpos(structure,one_atom))
          endif
     
          if (mod(frame_count-1,IO_Rate).eq.0) then
             if (one_atom.eq.0) then
                call print(frame_count//' '//ave_r2,file=datafile)
             else
                call print(frame_count//' '//ave_r2//' '//one_r2,file=datafile)
             endif
          endif
       endif
       call read(structure,xyzfile,error=error)
     
    end do
  
    call print('')
    call print('Read '//frame_count//' frames, processed '//frames_processed//' frames.')
  
    !Free up memory
    call finalise(structure)
    call finalise(reference)
    call finalise(xyzfile)
    call finalise(datafile)
  
    call print('Finished.')

    call verbosity_pop
    call system_finalise

end program diffusion
