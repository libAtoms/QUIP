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

!Program to calculate the EVB E_GAP given
!   -- an XYZ file with the coordinates,
!   -- 2 PSF files with the 2 resonance states
!   -- a form_bond_file containing the pair of atoms that form a bond in PSF1, but not in PSF2
!   -- a break_bond_file containing the pair of atoms that form a bond in PSF2, but not in PSF1
!
program extract_EVB

  use libatoms_module
  use potential_module

  implicit none

  type(Atoms)                              :: at
  type(CInoutput)                          :: coord_file
  type(InOutput)                           :: out_file
  integer                                  :: iframe, step_nr
  real(dp)                                 :: time_passed
  real(dp)                                 :: energy1, energy2
  !real(dp), allocatable                    :: f1(:,:)
  character(len=FIELD_LENGTH), allocatable :: output_lines(:)
  integer                                  :: lines_used, i

  type(Potential)                          :: pot
  type(Potential)                          :: evbpot
  character(len=STRING_LENGTH)             :: args_str
  character(FIELD_LENGTH)                  :: mm_args_str
                                           
  !Input parameters                        
  type(Dictionary)                         :: params_in
  character(len=FIELD_LENGTH)              :: coord_filename, out_filename
  character(len=FIELD_LENGTH)              :: PSF_file1, PSF_file2, &
                                              filepot_program, cp2k_program !, Residue_Library
  integer                                  :: last_frame
  real(dp)                                 :: time_step
  logical                                  :: append_output
  integer                                  :: restart_every


    call system_initialise(verbosity=PRINT_NORMAL,enable_timing=.true.)
    call system_timer('program')

    !INPUT
      call initialise(params_in)
      call param_register(params_in, 'coord_file', PARAM_MANDATORY, coord_filename, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(params_in, 'outfile', PARAM_MANDATORY, out_filename, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(params_in, 'PSF_file1', PARAM_MANDATORY, PSF_file1, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(params_in, 'PSF_file2', PARAM_MANDATORY, PSF_file2, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(params_in, 'filepot_program', PARAM_MANDATORY, filepot_program, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(params_in, 'cp2k_program', PARAM_MANDATORY, cp2k_program, help_string="No help yet.  This source file was $LastChangedBy$")
      !call param_register(params_in, 'Residue_Library', PARAM_MANDATORY, Residue_Library, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(params_in, 'last_frame', '0', last_frame, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(params_in, 'time_step', '0.0', time_step, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(params_in, 'append_output', 'F', append_output, help_string="No help yet.  This source file was $LastChangedBy$")
      call param_register(params_in, 'restart_every', '1', restart_every, help_string="No help yet.  This source file was $LastChangedBy$")

      if (.not. param_read_args(params_in)) then
        !call print_usage
        call print("extract_EVB coord_file outfile PSF_file1 PSF_file2 filepot_program cp2k_program [last_frame] [time_step] [append_output] [restart_every]",PRINT_ALWAYS)
        call system_abort('could not parse argument line')
      end if

      call finalise(params_in)

    !PRINT INPUT PARAMETERS
      call print('Run parameters:')
      call print('  coord_file      '// trim(coord_filename)     )
      call print('  outfile         '// trim(out_filename)       )
      call print('  PSF_file1       '// trim(PSF_file1)          )
      call print('  PSF_file2       '// trim(PSF_file2)          )
      !call print('  Residue_Library '// trim(Residue_Library)    )
      call print('  cp2k_program    '// trim(cp2k_program)       )
      call print('  filepot_program '// trim(filepot_program)    )
      call print('  last_frame      '// last_frame               )
      call print('  time_step       '// time_step                )
      call print('  restart_every   '// restart_every//' steps'  )
      if (restart_every/=1) call print('    Prints after every '//restart_every//' steps, synchronised with the restart frequency.')
      if (append_output) then
         call print("    Will not overwrite out_file (appending probably restart).")
      else
         call print("    Will overwrite out_file!!")
      endif
      call print('---------------------------------------')
      call print('')

  
    !Read coordinates and colvars
    call print('Reading in the coordinates from file '//trim(coord_filename)//'...')
    call initialise(coord_file, coord_filename, action=INPUT)
    call print('Number of frames '//coord_file%n_frame)

    if (last_frame/=0) then
       last_frame = min(last_frame, coord_file%n_frame)
    else
       last_frame = coord_file%n_frame
    endif
    call print('Last frame to be processed '//last_frame)

    !open output file, append if required
    if (append_output) then
       call initialise(out_file, out_filename, action=OUTPUT,append=.true.)
    else
       call initialise(out_file, out_filename, action=OUTPUT,append=.false.)
       call print("#  Step Nr.    Time[fs]   Energy_PSF1[eV]   Energy_PSF2[eV]      EVB(E1-E2)", file=out_file)
    endif

    !initialise EVB potential
    call initialise(pot,'FilePot command='//trim(filepot_program)//' property_list=pos min_cutoff=0.0')
    call initialise(evbpot,args_str='EVB=T form_bond={1 2} break_bond={2 6} offdiagonal_A12=0.0 offdiagonal_mu12=0.0 save_energies=T save_forces=F', &
         pot1=pot)

    !setup topology files
    call system_command('cp '//trim(PSF_file1)//' quip_cp2k_evb1.psf')
    call system_command('cp '//trim(PSF_file2)//' quip_cp2k_evb2.psf')


    step_nr = 0

    !allocate printing buffer
    allocate(output_lines(restart_every))
    lines_used = 0

    do iframe = 0, last_frame-1

       call read(coord_file, at, frame=iframe)

       write(mainlog%unit,'(a,a,i0,$)') achar(13),'Frame ',(iframe+1)

       !get Step_nr and Time
       if (.not.get_value(at%params,"i",step_nr)) then
          step_nr = step_nr + 1
       endif
       if (.not.get_value(at%params,"time",time_passed)) then
          time_passed = step_nr * time_step
       endif

       !Calculate EVB energy

       mm_args_str="Run_Type=MM cp2k_program="//trim(cp2k_program)//" clean_up_files=T save_output_files=F have_silica_potential=F centre_cp2k=F PSF_Print=USE_EXISTING_PSF"

       args_str = "evb_step="//step_nr//" mm_args_str={"//trim(mm_args_str)//"} topology_suffix1=_evb1 topology_suffix2=_evb2"
     !  call set_value(at%params,'Library',trim(Residue_Library))

       call calc(evbpot,at,energy=energy1,args_str=trim(args_str))
       !allocate(f1(1:3,1:at%N))
       !call calc(evbpot,at,energy=energy1,f=f1,args_str=trim(args_str))
       !deallocate(f1)
       !call print_xyz(at,'evb_all.xyz',all_properties=.true.)

       if (.not. get_value(at%params,'EVB1_energy',energy1)) then
          call system_abort('Could not find EVB1_energy')
       endif
       if (.not. get_value(at%params,'EVB2_energy',energy2)) then
          call system_abort('Could not find EVB2_energy')
       endif

       !Printing

       !if not append_output, we want to print time=0, too
       if (.not.append_output .and. iframe==0) then
          write (UNIT=output_lines(1),FMT="(I8,F13.3,3F18.8)") &
                 step_nr, time_passed, energy1, energy2, (energy1-energy2)
          call print(output_lines(1),file=out_file)
          output_lines(1) = ""
          call finalise(at)
          cycle
       endif

       !printing into buffer
       write (UNIT=output_lines(lines_used+1),FMT="(I8,F13.3,3F18.8)") &
              step_nr, time_passed, energy1, energy2, (energy1-energy2)
       lines_used = lines_used + 1

       !printing buffer into file
       if (lines_used == restart_every) then
          do i=1,lines_used
             call print(output_lines(i),file=out_file)
             output_lines(i) = ""
          enddo
          deallocate(output_lines)
          allocate(output_lines(restart_every))
          lines_used = 0
       endif

       call finalise(at)

    enddo

    call finalise(coord_file)
    call finalise(at)

    call print('Finished.')

    call system_timer('program')
    call system_finalise

end program extract_EVB
