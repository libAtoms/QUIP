module read_configs_module
use system_module
use atoms_module
use cinoutput_module
implicit none
private

public :: read_configs

contains

  !% reads a sequence of configurations with cinoutput, optionally skipping every decimation frames,
  !%  ignoring things outside of min_time--max_time, sorting by Time value and eliminating configs with
  !%  duplicate Time values.
  function read_configs(filename, file_is_list, decimation, min_time, max_time, sort_Time, no_Time_dups) result(structure_ll)
    character(len=*), intent(in) :: filename
    logical, intent(in) :: file_is_list
    integer, intent(in) :: decimation
    real(dp), intent(in) :: min_time, max_time
    logical, intent(in), optional :: sort_Time, no_Time_dups
    type(atoms_ll) :: structure_ll

    integer :: status
    integer :: frame_count, last_file_frame_n
    type(CInOutput) :: cfile
    type(Inoutput) file_list
    character(len=len(filename)) my_filename
    type(Atoms) :: structure_in
    type(Atoms), pointer :: structure
    logical :: skip_frame
    real(dp) :: cur_time, entry_time, entry_PREV_time
    logical :: do_sort_Time, do_no_Time_dups
    type(Atoms_ll_entry), pointer :: entry
    logical :: is_a_dup

    do_sort_Time = optional_default(.false., sort_Time)
    do_no_Time_dups = optional_default(.false., no_Time_dups)

    if (do_no_Time_dups .and. .not. do_sort_Time) call system_abort("ERROR: read_configs do_no_Times_dups requires do_sort_Time")

    if (decimation /= 1 .and. do_no_Time_dups) call system_abort("ERROR: read_configs decimation="//decimation//" /= 1 and no_Time_dups=T conflict")

    status = 0
    if (file_is_list) then
      call initialise(file_list, trim(filename), INPUT)
      my_filename = read_line(file_list)
    else
      my_filename = trim(filename)
    endif

    last_file_frame_n = 0
    frame_count = 1
    do while (status == 0) ! loop over files
      call initialise(cfile,trim(my_filename),action=INPUT)
      status = 0
      do while (status == 0) ! loop over frames in this file
	! write(mainlog%unit,'(4a,i0,a,i0,a,$)') achar(13), 'Read file ',trim(my_filename), ' Frame ',frame_count,' which in this file is frame (zero based) ',(frame_count-1-last_file_frame_n),"    "
	write(mainlog%unit,'(3a,i0,a,i0)') 'Read file ',trim(my_filename), ' Frame ',frame_count,' which in this file is frame (zero based) ',(frame_count-1-last_file_frame_n)
	call read(cfile, structure_in, frame=frame_count-1-last_file_frame_n, status=status)

	if (status == 0) then ! we succesfully read a structure
	  skip_frame = .false.
	  is_a_dup = .false.
	  if (min_time > 0.0_dp .or. max_time > 0.0_dp .or. do_sort_Time .or. no_Time_dups) then ! we need Time value
	    if (get_value(structure_in%params,"Time",cur_time)) then
	      if ((min_time >= 0.0_dp .and. cur_time < min_time) .or. (max_time >= 0.0_dp .and. cur_time > max_time)) skip_frame = .true.
	    else
	      call system_abort("ERROR: min_time="//min_time//" < 0 or max_time="//max_time//" < 0 or do_sort_Time="//do_sort_Time//", but Time field wasn't found in config " // frame_count)
	    endif
	  endif
	  if (.not. skip_frame) then ! frame is in appropriate time range
	    if (do_sort_Time) then ! sort by Time value
	      ! look at LAST entry
	      if (associated(structure_ll%LAST)) then
		entry => structure_ll%LAST
		if (.not. get_value(entry%at%params, "Time", entry_time)) call system_abort("read_configs sort missing Time for LAST entry")
	      else ! no structure list yet, fake entry_time
		entry_time = -1.0e38_dp
	      endif
	      if (cur_time <= entry_time) then ! new frame is at or BEFORE LAST entry
		do while (associated(entry))
		  if (.not. get_value(entry%at%params, "Time", entry_time)) call system_abort("read_configs sort missing Time for entry")
		  if (do_no_Time_dups .and. (cur_time == entry_time)) then ! entries match in time, i.e. duplicates
		    is_a_dup = .true.
		    exit
		  endif
		  ! check time of PREV entry
		  if (associated(entry%PREV)) then
		    if (.not. get_value(entry%PREV%at%params, "Time", entry_PREV_time)) call system_abort("read_configs sort missing Time for entry%PREV")
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
	      ! we skipped previous section if cur_time > entry_time for LAST entry, so we want to insert AFTER entry=structure_ll%LAST
	    else ! no sorting in time, so always insert AFTER last 
	      entry => structure_ll%last
	    endif
	    if (.not. is_a_dup) then
	      if (associated(entry)) then ! we want to insert AFTER entry
		call new_entry(structure_ll, structure, AFTER=entry)
	      else ! entry is unassociated, therefore we insert and the BEGINNING of the list
		call new_entry(structure_ll, structure, BEFORE=structure_ll%FIRST)
	      endif
	      ! actually copy the structure
	      call atoms_copy_without_connect(structure, structure_in, properties="species:pos:Z")
	    endif
	  else ! skip_frame was true, we're skipping
	    ! write (mainlog%unit,'(a,$)') " skip     "
	    write (mainlog%unit,'(a)') " skip"
	  endif ! skip_frame
	  frame_count = frame_count + decimation
	endif ! status == 0 for reading this structure

      end do ! while status == 0 for frames in this file
      last_file_frame_n = last_file_frame_n + cfile%n_frame
      ! call print("at end of file, frame_count " // frame_count // " last_file_frame_n " // last_file_frame_n)
      call finalise(cfile)

      if (file_is_list) then
	my_filename = read_line(file_list, status)
	if (status /= 0) call finalise(file_list)
      else
	status = 1
      endif

    end do ! while status == 0 for each file

    call print("")
  end function read_configs

end module read_configs_module
