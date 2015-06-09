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

!% Atoms_ll_module
!% linked list of Atoms structures, read in from one or more xyz files, with 
!% some automatic time sorting/clean-up capability

#include "error.inc"

module Atoms_ll_module
use error_module
use system_module
use extendable_str_module
use dictionary_module
use atoms_module
use cinoutput_module
implicit none
private

public :: atoms_ll, atoms_ll_entry

type atoms_ll
  type (atoms_ll_entry), pointer :: first => null()
  type (atoms_ll_entry), pointer :: last => null()
end type atoms_ll

type atoms_ll_entry
  type(atoms) :: at
  type(atoms_ll_entry), pointer :: next => null()
  type(atoms_ll_entry), pointer :: prev => null()
  real(dp) :: r_index = -1.0_dp
end type atoms_ll_entry

public :: initialise, finalise, print_xyz, read_xyz, new_entry, remove_last_entry

interface initialise
  module procedure atoms_ll_initialise, atoms_ll_entry_initialise
end interface initialise

interface finalise
  module procedure atoms_ll_finalise, atoms_ll_entry_finalise
end interface finalise

!% print all configs in atoms_ll object
interface print_xyz
  module procedure atoms_ll_print_xyz
end interface

!% read configs into atoms_ll object
interface read_xyz
  module procedure atoms_ll_read_xyz_filename
end interface

!% add new entry to atoms linked list structure
interface new_entry
  module procedure atoms_ll_new_entry
end interface new_entry

!% remove last entry from atoms linked list structure
interface remove_last_entry
  module procedure atoms_ll_remove_last_entry
end interface remove_last_entry

contains

  subroutine atoms_ll_initialise(this)
    type(Atoms_ll), intent(inout) :: this

    call finalise(this)
  end subroutine atoms_ll_initialise

  subroutine atoms_ll_finalise(this)
    type(Atoms_ll), intent(inout) :: this

    type(atoms_ll_entry), pointer :: entry, next

    entry => this%first
    do while (associated(entry))
      next => entry%next
      call finalise(entry)
      entry => next
    end do
    this%first => null()
    this%last => null()
  end subroutine atoms_ll_finalise

  subroutine atoms_ll_entry_initialise(this)
    type(Atoms_ll_entry), intent(inout) :: this

    call finalise(this)
  end subroutine atoms_ll_entry_initialise

  subroutine atoms_ll_entry_finalise(this)
    type(Atoms_ll_entry), intent(inout) :: this

    call finalise(this%at)
    this%next => null()
    this%prev => null()
  end subroutine atoms_ll_entry_finalise
  
  subroutine atoms_ll_new_entry(this, atoms_p, before, after)
    type(atoms_ll), target, intent(inout) :: this    !% this: atoms_ll object to add to
    type(atoms), intent(inout), pointer :: atoms_p   !% atoms_p: pointer to atoms object to insert
    type(atoms_ll_entry), intent(in), target, optional :: before, after !% before, after: if present, put in new entry before/after this one
                                                                        !% before and after are mutually exclusive

    type(atoms_ll_entry), pointer :: my_before, my_after
    type(atoms_ll_entry), pointer :: entry

    if (present(before) .and. present(after)) call system_abort("atoms_ll_new_entry got both before and after")

    if (present(before)) then
      my_before => before
      my_after => before%prev
    else if (present(after)) then
      my_before => after%next
      my_after => after
    else
      my_after => this%last
      my_before => null()
    endif

    allocate(entry)
    if (associated(my_before)) then
      my_before%prev => entry
      entry%next => my_before
    else
      this%last => entry
    endif
    if (associated(my_after)) then
      my_after%next => entry
      entry%prev => my_after
    else
      this%first => entry
    endif

    atoms_p => entry%at
  end subroutine atoms_ll_new_entry

  subroutine atoms_ll_remove_last_entry(this)
    type(atoms_ll), target, intent(inout) :: this

    if (associated(this%last)) then
      call finalise(this%last%at)
      if (associated(this%last%prev)) nullify(this%last%prev%next)
      this%last => this%last%prev
      if (.not. associated (this%last)) nullify(this%first)
    endif

  end subroutine atoms_ll_remove_last_entry


  subroutine atoms_ll_print_xyz(this, xyzfile, comment, properties, real_format, mask)
    type(Atoms_ll),            intent(inout)    :: this     !% Atoms_ll object to print
    type(CInoutput),         intent(inout) :: xyzfile  !% CInOutput object to write to
    character(*), optional, intent(in) :: properties !% List of properties to print 
    character(*), optional, intent(in)    :: comment  !% Comment line (line #2 of xyz file)
    character(len=*), optional, intent(in) :: real_format   !% format of real numbers in output
    logical, optional, intent(in) :: mask(:)                !% mask of which atoms to print


    type(extendable_str) :: my_comment
    type(atoms_ll_entry), pointer :: entry
    integer :: i

    call initialise(my_comment)
    if (present(comment)) call concat(my_comment, comment)

    if (present(mask)) call system_abort('atoms_ll_print_xyz: mask argument not currently supported. if you need it, please reimplement!')

    entry => this%first
    i = 1
    do while (associated(entry)) 
       call set_value(entry%at%params, 'comment', my_comment)
       call set_value(entry%at%params, 'atoms_ll_i', i)
       call write(xyzfile, entry%at, properties, real_format=real_format)
       entry => entry%next
       i = i + 1
    end do
    call finalise(my_comment)

  end subroutine atoms_ll_print_xyz

  !% reads a sequence of configurations with cinoutput, optionally skipping every decimation frames,
  !%  ignoring things outside of min_time--max_time, sorting by Time value and eliminating configs with
  !%  duplicate Time values.
  subroutine atoms_ll_read_xyz_filename(this, filename, file_is_list, decimation, min_time, max_time, sort_Time, no_Time_dups, quiet, no_compute_index, &
                                        properties, all_properties, error)
    type(atoms_ll) :: this !% object to read configs into
    character(len=*), intent(in) :: filename !% file name to read from
    logical, intent(in) :: file_is_list !% is true, file is list of files to read from rather than xyz file
    integer, intent(in), optional :: decimation !% only read one out of every decimation frames
    real(dp), intent(in), optional :: min_time, max_time !% ignore frames before min_time or after max_time
    character(len=*), intent(in), optional :: properties !% properties to include in objects in list
    logical, intent(in), optional :: sort_Time, no_Time_dups, quiet, no_compute_index, all_properties 
    integer, intent(out), optional :: error
       !% sort_time: if true, sort configurations by time
       !% no_Times_dups: if true, eliminate configs that have same time value
       !% quiet: if true, don't output progress bar
       !% no_compute_index: if true, don't compute index for xyz file
       !% do_all_properties: if true, include all properties in save Atoms structures

    integer :: status
    integer :: frame_count, last_file_frame_n
    type(CInOutput) :: cfile
    type(Inoutput) file_list
    character(len=len(filename)) my_filename
    type(Atoms) :: structure_in
    type(Atoms), pointer :: structure
    logical :: skip_frame
    real(dp) :: cur_time, entry_time, entry_PREV_time
    integer :: do_decimation
    real(dp) :: do_min_time, do_max_time
    logical :: do_sort_Time, do_no_Time_dups, do_quiet, do_no_compute_index
    type(Atoms_ll_entry), pointer :: entry
    logical :: is_a_dup, do_all_properties
    character(len=STRING_LENGTH) :: my_properties
    integer :: initial_frame_count
    integer :: l_error

    INIT_ERROR(error)
    my_properties = optional_default("species:pos:Z", properties)
    do_all_properties = optional_default(.false., all_properties)
    do_decimation = optional_default(1, decimation)
    do_min_time = optional_default(-1.0_dp, min_time)
    do_max_time = optional_default(-1.0_dp, max_time)
    do_sort_Time = optional_default(.false., sort_Time)
    do_no_Time_dups = optional_default(.false., no_Time_dups)
    do_quiet = optional_default(.false., quiet)
    do_no_compute_index = optional_default(.false., no_compute_index)

    if (len_trim(my_properties) == 0) then
      call print("WARNING: len_trim(my_properties) == 0, doing all_properties")
      do_all_properties=.true.
    endif

    if (do_no_Time_dups .and. .not. do_sort_Time) then
      RAISE_ERROR("ERROR: atoms_ll_read_xyz no_Times_dups requires sort_Time", error)
    endif

    if (do_decimation /= 1 .and. do_no_Time_dups) then
      RAISE_ERROR("ERROR: atoms_ll_read_xyz decimation="//do_decimation//" /= 1 and no_Time_dups=T conflict", error)
    endif

    if (file_is_list) then
      call initialise(file_list, trim(filename), INPUT)
      my_filename = read_line(file_list)
    else
      my_filename = trim(filename)
    endif

    cur_time = 1.0_dp
    last_file_frame_n = 0
    frame_count = 1
    status = 0
    l_error = ERROR_NONE
    do while (l_error == ERROR_NONE) ! loop over files
      call initialise(cfile,trim(my_filename),action=INPUT, no_compute_index=do_no_compute_index)
      initial_frame_count = frame_count
      do while (l_error == ERROR_NONE) ! loop over frames in this file
	if (.not. do_quiet) write(mainlog%unit,'(4a,i0,a,i0,$)') achar(13), 'Read file ',trim(my_filename), &
	  ' Frame ',frame_count,' which in this file is frame (zero based) ',(frame_count-1-last_file_frame_n)
!	if (.not. do_quiet) write(mainlog%unit,'(3a,i0,a,i0)') 'Read file ',trim(my_filename), &
!	     ' Frame ',frame_count,' which in this file is frame (zero based) ',(frame_count-1-last_file_frame_n)
	call read(cfile, structure_in, frame=frame_count-1-last_file_frame_n, error=l_error)

	if (l_error == ERROR_NONE) then ! we succesfully read a structure
	  skip_frame = .false.
	  is_a_dup = .false.
	  if (do_min_time > 0.0_dp .or. do_max_time > 0.0_dp .or. do_sort_Time .or. no_Time_dups) then ! we need Time value
	    if (get_value(structure_in%params,"Time",cur_time, case_sensitive=.false.)) then
	      if ((do_min_time >= 0.0_dp .and. cur_time < do_min_time) .or. (do_max_time >= 0.0_dp .and. cur_time > do_max_time)) skip_frame = .true.
	    else
              RAISE_ERROR("ERROR: min_time="//do_min_time//" < 0 or max_time="//do_max_time//" < 0 or sort_Time="//do_sort_Time//", but Time field wasn't found in config " // frame_count, error)
	    endif
	  endif
	  if (.not. skip_frame) then ! frame is in appropriate time range
	    if (do_sort_Time) then ! sort by Time value
	      ! look at LAST entry
	      if (associated(this%LAST)) then
		entry => this%LAST
		!NB if (.not. get_value(entry%at%params, "Time", entry_time, case_sensitive=.false.)) call system_abort("atoms_ll_read_xyz sort missing Time for LAST entry")
		entry_time = entry%r_index
		!NB
	      else ! no structure list yet, fake entry_time
		entry_time = -1.0e38_dp
                nullify(entry)
	      endif
	      if (cur_time <= entry_time) then ! new frame is at or BEFORE LAST entry
		do while (associated(entry))
		  !NB if (.not. get_value(entry%at%params, "Time", entry_time, case_sensitive=.false.)) call system_abort("atoms_ll_read_xyz sort missing Time for entry")
		  entry_time = entry%r_index
		  !NB
		  if (do_no_Time_dups .and. (cur_time == entry_time)) then ! entries match in time, i.e. duplicates
		    is_a_dup = .true.
		    exit
		  endif
		  ! check time of PREV entry
		  if (associated(entry%PREV)) then
		    !NB if (.not. get_value(entry%PREV%at%params, "Time", entry_PREV_time, case_sensitive=.false.)) call system_abort("atoms_ll_read_xyz sort missing Time for entry%PREV")
		    entry_PREV_time = entry%PREV%r_index
		    !NB
		    ! if cur_time is <= entry, and cur_time > entry%PREV, we want to insert AFTER entry%prev, so decrement entry and exit now
		    if ((cur_time <= entry_time) .and. (cur_time > entry_PREV_time)) then
		      entry => entry%PREV
		      exit
		    endif
		  endif
		  ! otherwise we keep on going
		  entry => entry%PREV ! entry will be unassociated, and therefore loop will end, if we got to beginning of list
		end do ! associated(entry)
	      end if ! cur frame was BEFORE last entry
	      ! we skipped previous section if cur_time > entry_time for LAST entry, so we want to insert AFTER entry=this%LAST
	    else ! no sorting in time, so always insert AFTER last 
	      entry => this%last
	    endif
	    if (.not. is_a_dup) then
	      if (associated(entry)) then ! we want to insert AFTER entry
		call new_entry(this, structure, AFTER=entry)
		entry%NEXT%r_index = cur_time
	      else ! entry is unassociated, therefore we insert and the BEGINNING of the list
		call new_entry(this, structure, BEFORE=this%FIRST)
		this%FIRST%r_index = cur_time
	      endif
	      ! actually copy the structure
              if (do_all_properties) then
                call atoms_copy_without_connect(structure, structure_in)
              else
                call atoms_copy_without_connect(structure, structure_in, properties=my_properties)
              endif
	    endif
	    if (.not. do_quiet) write (mainlog%unit,'(a,$)') "          "
!	    if (.not. do_quiet) write (mainlog%unit,'(a)') "          "
	  else ! skip_frame was true, we're skipping
	    if (.not. do_quiet) write (mainlog%unit,'(a,$)') " skip     "
!	    if (.not. do_quiet) write (mainlog%unit,'(a)') " skip"
	  endif ! skip_frame
	  frame_count = frame_count + do_decimation
        else
           CLEAR_ERROR(error)
	endif ! l_error == 0 for reading this structure

      end do ! while l_error == 0 for frames in this file
      if (cfile%got_index) then
	last_file_frame_n = last_file_frame_n + cfile%n_frame
      else
	last_file_frame_n = last_file_frame_n + cfile%current_frame
      endif
      !call print("at end of file, frame_count " // frame_count // " last_file_frame_n " // last_file_frame_n)
      call finalise(cfile)

      if (file_is_list) then
	my_filename = read_line(file_list, status)
	if (status /= 0) call finalise(file_list)
      else
	status = 1
      endif

    end do ! while status == 0 for each file

    if (.not. quiet) call print("")
  end subroutine atoms_ll_read_xyz_filename

end module Atoms_ll_module
