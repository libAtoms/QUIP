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

!Program to calculate coordination number as a fn of time, if given is
!   -- an XYZ file with the coordinates,
!   -- the atomic number of the atoms in the coordination sphere
!   -- the atom number of the central atom
!   and the output name specified.
!Optionally restarts can be taken into account by printing synchronised
!  with the restart frequency (e.g. restart_frequency=20).
!Optionally appends the output (of restart files) to the same outfile (e.g. append_output=T).
!
!A Fermi distribution is used to smooth the coordination number
!M. Sprik, Far.Disc. 110, 437-445 (1998)
!
program coordination_number

  use libatoms_module

  implicit none

  type(Atoms)                   :: at
  type(CInoutput)               :: coord_file
  type(InOutput)                :: out_file
  character(len=STRING_LENGTH), allocatable :: output_lines(:)
  integer                       :: lines_used
  integer                       :: iframe, step_nr
  real(dp)                      :: time_passed
  integer                       :: i

  !Input parameters
  type(Dictionary)              :: params_in
  character(len=STRING_LENGTH)   :: coord_filename
  character(len=STRING_LENGTH)   :: out_filename
  integer                       :: last_frame
  real(dp)                      :: time_step
  integer                       :: restart_every
  logical                       :: append_output
  logical                       :: skip_extra_firsts
real(dp) :: dist, coord_num
logical, allocatable :: mask_a(:)
character(STRING_LENGTH) :: mask_str
real(dp) :: kappa, l_cutoff
integer :: iatom, center_atom
logical :: no_smoothing


!    call system_initialise(verbosity=ANAL,enable_timing=.true.)
!    call system_initialise(verbosity=VERBOSE,enable_timing=.true.)
    call system_initialise(verbosity=PRINT_NORMAL,enable_timing=.true.)
    call system_timer('program')

    !INPUT
      call initialise(params_in)
      call param_register(params_in, 'infile'             , PARAM_MANDATORY, coord_filename   , help_string="input XYZ"   )
      call param_register(params_in, 'outfile'            , PARAM_MANDATORY, out_filename     , help_string="output file"   )
      mask_str=''
      call param_register(params_in, 'neighbour_mask'     , ''             , mask_str         , help_string="neighbour atom in the coordination sphere"   )
      call param_register(params_in, 'center_mask'        , '0'            , center_atom      , help_string="center atom, around which the coordination number is calculated"   )
      call param_register(params_in, 'kappa'              , '-1.0'         , kappa            , help_string="smoothing by default: 1/(exp(kappa*(d-cutoff))+1)"   )
      call param_register(params_in, 'cutoff'             , PARAM_MANDATORY, l_cutoff           , help_string="cutoff of the coordination sphere"   )
      call param_register(params_in, 'last_frame'         , '0'            , last_frame       , help_string="last frame to process (0=all)"   )
      call param_register(params_in, 'time_step'          , '0.0'          , time_step        , help_string="time step"   )
      call param_register(params_in, 'restart_every'      , '1'            , restart_every    , help_string="to be consistent with the restart frequency"   )
      call param_register(params_in, 'append_output'      , 'F'            , append_output    , help_string="do not overwrite output file"   )
      call param_register(params_in, 'skip_extra_firsts'  ,'F'             , skip_extra_firsts, help_string="QMMM_md_buf prints the first time step twice"   )
      call param_register(params_in, 'no_smoothing'       ,'F'             , no_smoothing     , help_string="whether to use smoothing, or an abrupt cutoff"   )

      if (.not. param_read_args(params_in)) then
        call print("coordination_number infile outfile cutoff [center_atom] [neighbour_atom] [kappa=-1.0] [last_frame=0] [time_step=0.0] [restart_every=1] [append_output=F] [skip_extra_firsts=F] [no_smoothing=F]", PRINT_ALWAYS)
        call system_abort('could not parse argument line')
      end if

      call finalise(params_in)

      if (kappa<0) call system_abort("kappa<0")
      if (l_cutoff<0) call system_abort("cutoff<0")

    !PRINT INPUT PARAMETERS
      call print('Run parameters:')
      call print('  coord_file      '// trim(coord_filename)     )
      call print('  outfile         '// trim(out_filename)       )
      call print('  neighbour_mask  '// trim(mask_str)           )
      call print('  center_mask     '// center_atom              )
      call print('  last_frame      '// last_frame               )
      call print('  kappa           '// kappa                    )
      call print('  cutoff          '// l_cutoff                   )
      call print('  time_step       '// time_step                )
      call print('  restart_every   '// restart_every//' steps'  )
      if (restart_every/=1) call print('    Prints after every '//restart_every//' steps, synchronised with the restart frequency.')
      call print('  append_output   '//append_output)
      if (append_output) then
         call print("    Will not overwrite out_file (appending probably restart).")
      else
         call print("    Will overwrite out_file!!")
      endif
      call print('  skip_extra_firsts   '//skip_extra_firsts//" ONLY FOR cp2k_units=F !")
      call print('---------------------------------------')
      call print('')

  
    !Read coordinates
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
       call print("# Step Nr.   Time[fs]  n_coord("//trim(mask_str)//")   ((r_cut="//l_cutoff//"))", file=out_file)
    endif

    step_nr = 0

    !allocate printing buffer
    allocate(output_lines(restart_every))
    lines_used = 0

    do iframe = 0, last_frame-1

       call read(coord_file, at, frame=iframe)

       write(mainlog%unit,'(a,a,i0,$)') achar(13),'Frame ',(iframe+1)

       if (iframe==0) then
          if (skip_extra_firsts) cycle !skip very first without force information
       endif
       if (iframe==1) then
          if (skip_extra_firsts .and. append_output) cycle !skip the overlapping second for the restarts
       endif

       !get Step_nr and Time
       if (.not.get_value(at%params,"i",step_nr)) then
          step_nr = step_nr + 1
       endif
       if (.not.get_value(at%params,"time",time_passed)) then
          time_passed = step_nr * time_step
       endif

       if (center_atom<1.or.center_atom>at%N) call system_abort("center_atom="//center_atom//" is <0 or >"//at%N)
       allocate(mask_a(at%N))
       call is_in_mask(mask_a, at, mask_str)

       coord_num=0._dp
       !calc coordination number
       do iatom=1,at%N
          if (iatom==center_atom) cycle
          if (.not.mask_a(iatom)) cycle
          dist = distance_min_image(at,center_atom,iatom)
          if (no_smoothing) then
             if (dist<l_cutoff) coord_num = coord_num +1._dp
          else
             coord_num = coord_num + 1._dp/(exp(kappa*(dist-l_cutoff))+1._dp)
          endif
       enddo
       deallocate(mask_a)

       if (.not.append_output .and. iframe==0) then !print time=0 line
             write (UNIT=output_lines(1),FMT="(I8,F13.3,F13.8)") &
                    step_nr, time_passed, coord_num
             call print(output_lines(1),file=out_file)
             output_lines(1) = ""
             call finalise(at)
             cycle
       endif

       !printing into buffer
       write (UNIT=output_lines(lines_used+1),FMT="(I8,F13.3,F13.8)") &
              step_nr, time_passed, coord_num
       lines_used = lines_used + 1

       !printing buffer into file
       if (lines_used == restart_every) then
          !if (mod((iframe+1),restart_every)/=0) call system_abort("lines_used "//lines_used// &
          !   " == restart_every "//restart_every// &
          !   ", but not frame "//(iframe+1)//" mod restart_every "//restart_every//" =0")
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
    call finalise(out_file)
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

subroutine is_in_mask(mask_a, at, mask_str)
  type(Atoms), intent(in) :: at
  logical, intent(out) :: mask_a(at%N)
  character(len=*), optional, intent(in) :: mask_str

  integer :: i_at, i_Z
  integer :: Zmask
  type(Table) :: atom_indices
  character(len=4) :: species(128)
  integer :: n_species

  if (.not. present(mask_str)) then
    mask_a = .true.
    return
  endif

  if (len_trim(mask_str) == 0) then
    mask_a = .true.
    return
  endif

  mask_a = .false.
  if (mask_str(1:1)=='@') then
    call parse_atom_mask(mask_str,atom_indices)
    do i_at=1, atom_indices%N
      mask_a(atom_indices%int(1,i_at)) = .true.
    end do
  else if (scan(mask_str,'=')/=0) then
    call system_abort("property type mask not supported yet")
  else
    call split_string(mask_str, ' ,', '""', species, n_species)
    do i_Z=1, n_species
      Zmask = Atomic_Number(species(i_Z))
      do i_at=1, at%N
        if (at%Z(i_at) == Zmask) mask_a(i_at) = .true.
      end do
    end do
  end if
end subroutine is_in_mask

end program coordination_number
