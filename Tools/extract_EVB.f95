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

!Program to calculate EVB given
!   -- an XYZ file with the coordinates,
!   -- 2 PSF files with the 2 resonance states
!   -- a form_bond_file containing the pair of atoms that form a bond in PSF1, but not in PSF2
!   -- a break_bond_file containing the pair of atoms that form a bond in PSF2, but not in PSF1
!
program extract_EVB

  use libatoms_module
  use quip_module

  implicit none

  type(Atoms)                   :: at
  type(CInoutput)               :: coord_file
  !real(dp),       pointer       :: avgpos_p(:,:)
  integer                       :: iframe, step_nr
  real(dp)                      :: time_passed
  real(dp)                      :: energy1
  type(Potential)               :: pot
  type(Potential)               :: evbpot
  character(len=STRING_LENGTH)  :: args_str

  !Input parameters
  type(Dictionary)              :: params_in
  character(len=FIELD_LENGTH)   :: coord_filename
  character(len=FIELD_LENGTH)   :: PSF_file1, PSF_file2, filepot_program, cp2k_program !, Residue_Library
  integer                       :: last_frame
  real(dp)                      :: time_step
  character(FIELD_LENGTH)       :: mm_args_str


!    call system_initialise(verbosity=ANAL,enable_timing=.true.)
!    call system_initialise(verbosity=VERBOSE,enable_timing=.true.)
    call system_initialise(verbosity=PRINT_NORMAL,enable_timing=.true.)
    call system_timer('program')

    !INPUT
      call initialise(params_in)
      call param_register(params_in, 'coord_file'         , PARAM_MANDATORY, coord_filename      )
      call param_register(params_in, 'last_frame'         , '0'            , last_frame          )
      !call param_register(params_in, 'time_step'          , '0.0'          , time_step           )
      call param_register(params_in, 'PSF_file1'          , PARAM_MANDATORY, PSF_file1           )
      call param_register(params_in, 'PSF_file2'          , PARAM_MANDATORY, PSF_file2           )
      !call param_register(params_in, 'Residue_Library'    , PARAM_MANDATORY, Residue_Library     )
      call param_register(params_in, 'filepot_program'    , PARAM_MANDATORY, filepot_program     )
      call param_register(params_in, 'cp2k_program'       , PARAM_MANDATORY, cp2k_program        )

      if (.not. param_read_args(params_in, do_check = .true.)) then
        !call print_usage
        call print("extract_EVB coord_file PSF_file1 PSF_file2 filepot_program cp2k_program [last_frame]",PRINT_ALWAYS)
        call system_abort('could not parse argument line')
      end if

      call finalise(params_in)

    !PRINT INPUT PARAMETERS
      call print('Run parameters:')
      call print('  coord_file      '// trim(coord_filename)     )
      call print('  PSF_file1       '// trim(PSF_file1)          )
      call print('  PSF_file2       '// trim(PSF_file2)          )
      !call print('  Residue_Library '// trim(Residue_Library)    )
      call print('  cp2k_program    '// trim(cp2k_program)       )
      call print('  filepot_program '// trim(filepot_program)    )
      call print('  last_frame      '// last_frame               )
      !call print('  time_step       '// time_step                )
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


    !initialise EVB potential
!    call initialise(pot,'FilePot command='//trim(filepot_program)//' property_list=pos:avgpos:mol_id:atom_res_number min_cutoff=0.0')
    call initialise(pot,'FilePot command='//trim(filepot_program)//' property_list=pos min_cutoff=0.0')
    call initialise(evbpot,args_str='EVB=T form_bond={1 2} break_bond={2 6}', &
         pot1=pot)
    !setup topology files
    call system_command('cp '//trim(PSF_file1)//' quip_cp2k_evb1.psf')
    call system_command('cp '//trim(PSF_file2)//' quip_cp2k_evb2.psf')


    step_nr = 0

    do iframe = 0, last_frame-1

       call read(coord_file, at, frame=iframe)

       write(mainlog%unit,'(a,a,i0,$)') achar(13),'Frame ',(iframe+1)

       !if (iframe==0) then !check colvar atoms
       !   if ((.not.input_in_cp2k_units) .and. skip_extra_firsts) cycle !skip very first without force information
       !endif
       !if (iframe==1) then
       !   if ((.not.input_in_cp2k_units) .and. skip_extra_firsts .and. append_output) cycle !skip the overlapping second for the restarts
       !endif

       !get Step_nr and Time
       if (.not.get_value(at%params,"i",step_nr)) then
          step_nr = step_nr + 1
       endif
       if (.not.get_value(at%params,"time",time_passed)) then
          time_passed = step_nr * time_step
       endif

       !calculate EVB energy
          mm_args_str="Run_Type=MM cp2k_program="//trim(cp2k_program)//" clean_up_files=T save_output_files=F have_silica_potential=F centre_cp2k=F PSF_Print=USE_EXISTING_PSF"
!          mm1_args_str="Run_Type=MM cp2k_program="//trim(cp2k_program)//" clean_up_files=T save_output_files=F have_silica_potential=F centre_cp2k=F PSF_Print=DRIVER_PRINT_AND_SAVE topology_suffix=_evb1"

          args_str = "evb_step="//step_nr//" mm_args_str={"//trim(mm_args_str)//"} topology_suffix1=_evb1 topology_suffix2=_evb2"
        !  call calc_connect(at)
        !  call map_into_cell(at)
        !  call calc_dists(at)
        !  call set_value(at%params,'Library',trim(Residue_Library))
        !  call create_residue_labels_arb_pos(at,do_CHARMM=.true.,pos_field_for_connectivity="pos")
          !call add_property(at,"avgpos",0._dp,n_cols=3)
          !if(.not.assign_pointer(at,"avgpos",avgpos_p)) call system_abort("no avgpos")
          !avgpos_p=at%pos(1:3,1:at%N)
          call calc(evbpot,at,e=energy1,args_str=trim(args_str))
   
          call finalise(at)


       call finalise(at)

    enddo

    call finalise(coord_file)
    call finalise(at)

    call print('Finished.')

    call system_timer('program')
    call system_finalise

!contains
!
!
!  subroutine print_usage
!
!    call print("Usage: extract_cv coord_file force_file colvar_file outfile [restart_every=1] [append_output=F] [last_frame=0] [time_step=0.0]")
!    call print("   coord_file     XYZ coordinate file")
!    call print("   velo_file      XYZ file with velocities stored in 'velo' property")
!    call print("   force_file     XYZ file with forces stored in 'frc' property")
!    call print("   atomlist_file  file containing the (colvar) atoms whose properties should be printed")
!    call print("   outfile        output file with positions, velocities and forces on the (colvar) atoms in time.")
!    call print("   restart_every  optionally stores the output in a buffer and prints synchronised with the restart frequency.")
!    call print("                     Useful for restarted runs.")
!    call print("   append_output  optionally append to an existing output file without reprinting the header")
!    call print("                     Useful for restarted runs.")
!    call print("   last_frame     optionally processes only the first last_frame frames.  0 means process all frames.")
!    call print("   time_step      optionally sets the time step between frames.  It is used if 'time' is not stored in the coord_file.")
!
!  end subroutine print_usage

end program extract_EVB
