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

!Program to extract position, velocity and force on certain atoms given
!   -- an XYZ file with the coordinates,
!   -- an XYZ file with the velocities (can be the same file) in a.u.,
!   -- an XYZ file with the forces (can be the same file) in a.u.,
!   -- a file containing the atoms (typically to which the constraint is applied)
!   and the output name specified.
!Optionally restarts can be taken into account by printing synchronised
!  with the restart frequency (e.g. restart_frequency=20).
!Optionally appends the output (of restart files) to the same outfile (e.g. append_output=T).
!
program extract_cv

  use libatoms_module

  implicit none

  type(Atoms)                   :: at, at_force, at_velo
  type(CInoutput)               :: coord_file, velo_file, force_file
  type(InOutput)                :: out_file
  type(Table)                   :: colvar_atoms
  real(dp),       pointer       :: f(:,:), v(:,:)
  character(len=FIELD_LENGTH), allocatable :: output_lines(:)
  integer                       :: lines_used
  integer                       :: iframe, step_nr
  real(dp)                      :: time_passed
  integer                       :: i

  !Input parameters
  type(Dictionary)              :: params_in
  character(len=FIELD_LENGTH)   :: coord_filename
  character(len=FIELD_LENGTH)   :: velo_filename
  character(len=FIELD_LENGTH)   :: force_filename
  character(len=FIELD_LENGTH)   :: out_filename
  character(len=FIELD_LENGTH)   :: colvar_filename
  integer                       :: last_frame
  real(dp)                      :: time_step
  integer                       :: restart_every
  logical                       :: append_output
  logical                       :: input_in_cp2k_units
  logical                       :: skip_extra_firsts


!    call system_initialise(verbosity=PRINT_ANAL,enable_timing=.true.)
!    call system_initialise(verbosity=PRINT_VERBOSE,enable_timing=.true.)
    call system_initialise(verbosity=PRINT_NORMAL,enable_timing=.true.)
    call system_timer('program')

    !INPUT
      call initialise(params_in)
      call param_register(params_in, 'coord_file'         , PARAM_MANDATORY, coord_filename      , help_string="input file containing positions")
      call param_register(params_in, 'velo_file'          , PARAM_MANDATORY, velo_filename       , help_string="input file containing velocities")
      call param_register(params_in, 'force_file'         , PARAM_MANDATORY, force_filename      , help_string="input file containing forces")
      call param_register(params_in, 'atomlist_file'      , PARAM_MANDATORY, colvar_filename     , help_string="file with atomlist whose pos/velo/frc should be extracted")
      call param_register(params_in, 'outfile'            , PARAM_MANDATORY, out_filename        , help_string="output filename")
      call param_register(params_in, 'last_frame'         , '0'            , last_frame          , help_string="last frame to be processed")
      call param_register(params_in, 'time_step'          , '0.0'          , time_step           , help_string="time step in case there is no Time param in the input posfile")
      call param_register(params_in, 'restart_every'      , '1'            , restart_every       , help_string="restart frequency")
      call param_register(params_in, 'append_output'      , 'F'            , append_output       , help_string="whether to append to output file")
      call param_register(params_in, 'input_in_cp2k_units', 'F'            , input_in_cp2k_units , help_string="whether units are in CP2K units or in A, A/fs, eV/A")
      call param_register(params_in, 'skip_extra_firsts'  ,'F'             , skip_extra_firsts   , help_string="whether to skip the extra printed first config (for QUIP outputs)")

      if (.not. param_read_args(params_in)) then
        call print_usage
        call system_abort('could not parse argument line')
      end if

      call finalise(params_in)

    !PRINT INPUT PARAMETERS
      call print('Run parameters:')
      call print('  coord_file      '// trim(coord_filename)     )
      call print('  velo_file       '// trim(velo_filename)      )
      call print('  force_file      '// trim(force_filename)     )
      call print('  outfile         '// trim(out_filename)       )
      call print('  atomlist_file   '// trim(colvar_filename)    )
      call print('  last_frame      '// last_frame               )
      call print('  time_step       '// time_step                )
      call print('  restart_every   '// restart_every//' steps'  )
      if (restart_every/=1) call print('    Prints after every '//restart_every//' steps, synchronised with the restart frequency.')
      call print('  append_output   '//append_output)
      if (append_output) then
         call print("    Will not overwrite out_file (appending probably restart).")
      else
         call print("    Will overwrite out_file!!")
      endif
      call print('  input_in_cp2k_units '//input_in_cp2k_units)
      call print('  skip_extra_firsts   '//skip_extra_firsts//" ONLY FOR cp2k_units=F !")
      call print('---------------------------------------')
      call print('')

  
    !Read coordinates and colvars
    call print('Reading in the coordinates from file '//trim(coord_filename)//'...')
    call initialise(coord_file, coord_filename, action=INPUT)
    call print('Number of frames '//coord_file%n_frame)
    call print('Reading in the velocities from file '//trim(velo_filename)//'...')
    call initialise(velo_file, velo_filename, action=INPUT)
    call print('Number of frames '//velo_file%n_frame)
    call print('Reading in the forces from file '//trim(force_filename)//'...')
    call initialise(force_file, force_filename, action=INPUT)
    call print('Number of frames '//force_file%n_frame)

    if (last_frame/=0) then
       last_frame = min(last_frame, coord_file%n_frame, force_file%n_frame)
    else
       last_frame = min(coord_file%n_frame, velo_file%n_frame, force_file%n_frame)
    endif
    call print('Last frame to be processed '//last_frame)

    call print('Reading in the (colvar) atoms from file '//trim(colvar_filename)//'...')
    call read_integer_table(colvar_atoms,trim(colvar_filename))

    !open output file, append if required
    if (append_output) then
       call initialise(out_file, out_filename, action=OUTPUT,append=.true.)
    else
       call initialise(out_file, out_filename, action=OUTPUT,append=.false.)
       call print("# Step Nr.   Time[fs]  Atom Nr.  Pos_x[A]     Pos_y[A]     Pos_z[A]    Vel_x[A/fs]  Vel_y[A/fs]  Vel_z[A/fs]  Frc_x[eV/A]  Frc_y[eV/A]  Frc_z[eV/A]", file=out_file)
       !call print("#     Step Nr.          Time[fs]        Atom Nr.          Pos.[A]            Vel.[A/fs]        Force[eV/A]",file=out_file)
    endif

    step_nr = 0

    !allocate printing buffer
    allocate(output_lines(restart_every*colvar_atoms%N))
    lines_used = 0

    do iframe = 0, last_frame-1

       call read(coord_file, at, frame=iframe)
       call read(velo_file, at_velo, frame=iframe)
       call read(force_file, at_force, frame=iframe)

       write(mainlog%unit,'(a,a,i0,$)') achar(13),'Frame ',(iframe+1)

       if (iframe==0) then !check colvar atoms
          if (any(colvar_atoms%int(1,1:colvar_atoms%N).gt.at%N)) &
             call system_abort("Constraint atom(s) is >"//at%N)

          if ((.not.input_in_cp2k_units) .and. skip_extra_firsts) cycle !skip very first without force information
       endif
       if (iframe==1) then
          if ((.not.input_in_cp2k_units) .and. skip_extra_firsts .and. append_output) cycle !skip the overlapping second for the restarts
       endif

       !get Step_nr and Time
       if (.not.get_value(at%params,"i",step_nr)) then
          step_nr = step_nr + 1
       endif
       if (.not.get_value(at%params,"time",time_passed)) then
          if (.not.get_value(at%params,"Time",time_passed)) then
             time_passed = step_nr * time_step
          endif
       endif

       !get velocity and force in a.u.
       nullify(v)
       if (.not. assign_pointer(at_velo, 'velo', v)) call system_abort("Could not find velo property.")
       nullify(f)
       if (.not. assign_pointer(at_force, 'frc', f)) then
          if (.not. assign_pointer(at_force, 'force', f)) call system_abort("Could not find frc or force property.")
       endif

       !printing into buffer
       do i=1,colvar_atoms%N
          if (input_in_cp2k_units) then
             write (UNIT=output_lines(lines_used+1),FMT="(I8,F13.3,I8,9F13.8)") &
                    step_nr, time_passed, colvar_atoms%int(1,i), &
                    at%pos(1:3,colvar_atoms%int(1,i)), &
                    (v(1:3,colvar_atoms%int(1,i))*BOHR/AU_FS), &
                    (f(1:3,colvar_atoms%int(1,i))*HARTREE/BOHR)
          else
             write (UNIT=output_lines(lines_used+1),FMT="(I8,F13.3,I8,9F13.8)") &
                    step_nr, time_passed, colvar_atoms%int(1,i), &
                    at%pos(1:3,colvar_atoms%int(1,i)), &
                    (v(1:3,colvar_atoms%int(1,i))), &
                    (f(1:3,colvar_atoms%int(1,i)))
          endif
          lines_used = lines_used + 1
       enddo

       !printing buffer into file
       if (lines_used == restart_every*colvar_atoms%N) then
          !if (mod((iframe+1),restart_every)/=0) call system_abort("lines_used "//lines_used// &
          !   " == restart_every "//restart_every//" * colvar_atoms%N "//colvar_atoms%N// &
          !   ", but not frame "//(iframe+1)//" mod restart_every "//restart_every//" =0")
          do i=1,lines_used
             call print(output_lines(i),file=out_file)
             output_lines(i) = ""
          enddo
          deallocate(output_lines)
          allocate(output_lines(restart_every*colvar_atoms%N))
          lines_used = 0
       endif

       call finalise(at_force)
       call finalise(at_velo)
       call finalise(at)

    enddo

    call finalise(colvar_atoms)
    call finalise(coord_file)
    call finalise(velo_file)
    call finalise(force_file)
    call finalise(out_file)
    call finalise(at_force)
    call finalise(at_velo)
    call finalise(at)
    if (allocated(output_lines)) deallocate(output_lines)

    call print('Finished.')

    call system_timer('program')
    call system_finalise

contains

  subroutine read_integer_table(output_table,input_filename)

    type(Table), intent(out) :: output_table
    character(len=*), intent(in) :: input_filename
    type(InOutput) :: input_file
    integer :: n, num_entries
    integer :: dummy,stat

       call initialise(input_file,trim(input_filename),action=INPUT)
       read (input_file%unit,*) num_entries
       call allocate(output_table,1,0,0,0,num_entries)
       do n=1,num_entries
          dummy = 0
          read (input_file%unit,*,iostat=stat) dummy
          if (stat /= 0) exit
          call append(output_table,dummy)
       enddo
       if (output_table%N.ne.num_entries) call system_abort('read_integer_table: Something wrong with the atomlist file')
       call finalise(input_file)

  end subroutine read_integer_table

  subroutine print_usage

    call print("Usage: extract_cv coord_file force_file colvar_file outfile [restart_every=1] [append_output=F] [last_frame=0] [time_step=0.0]")
    call print("   coord_file     XYZ coordinate file")
    call print("   velo_file      XYZ file with velocities stored in 'velo' property")
    call print("   force_file     XYZ file with forces stored in 'frc' property")
    call print("   atomlist_file  file containing the (colvar) atoms whose properties should be printed")
    call print("   outfile        output file with positions, velocities and forces on the (colvar) atoms in time.")
    call print("   restart_every  optionally stores the output in a buffer and prints synchronised with the restart frequency.")
    call print("                     Useful for restarted runs.")
    call print("   append_output  optionally append to an existing output file without reprinting the header")
    call print("                     Useful for restarted runs.")
    call print("   last_frame     optionally processes only the first last_frame frames.  0 means process all frames.")
    call print("   time_step      optionally sets the time step between frames.  It is used if 'time' is not stored in the coord_file.")

  end subroutine print_usage

end program extract_cv
