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

!X
!X ParamReader module
!X
!X James Kermode <jrk33@cam.ac.uk>
!% The ParamReader module provides the facility to read parameters in the
!% form 'key = value' from files or from the command line. In typical
!% usage you would 'Param_Register' a series of parameters with default values,
!% then read some values from a parameter file overriding the defaults,
!% possibly then reading from the command line arguments to override the
!% options set by the parameter file.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

module paramreader_module

  use system_module
  use extendable_str_module
  use dictionary_module
  use table_module

  implicit none
  private

  public :: param_register, PARAM_MANDATORY, param_read_line, param_read_args, param_print_help, param_print, param_check, param_write_string

  integer, parameter, private :: MAX_N_FIELDS = 100       !% Maximum number of fields during parsing

  integer, parameter :: PARAM_NO_VALUE = 0 !% Special parameter type that doesn't get parsed
  integer, parameter :: PARAM_REAL = 1 !% Real (double precision) parameter
  integer, parameter :: PARAM_INTEGER = 2 !% Integer parameter
  integer, parameter :: PARAM_STRING = 3  !% String parameter
  integer, parameter :: PARAM_LOGICAL = 4 !% Logical parameter

  character(len=*), parameter :: PARAM_MANDATORY = '//MANDATORY//' !% If this is passed to 'Param_Register'
                                                                   !% instead of a default value then
                                                                   !% this parameter is not optional.

  type ParamEntry
     !% OMIT

     character(len=STRING_LENGTH):: value
     integer :: N, param_type

     real(dp), pointer :: real => null()
     integer, pointer :: integer => null()
     logical, pointer :: logical => null()
     character(len=STRING_LENGTH), pointer :: string => null()

     real(dp), dimension(:), pointer :: reals => null()
     integer, dimension(:), pointer :: integers => null()

     logical, pointer :: has_value

     character(len=STRING_LENGTH) :: help_string

  end type ParamEntry

  !% Overloaded interface to register a parameter in a Dictionary object.
  !% 'key' is the key name and 'value' the default value. For a mandatory
  !% parameter use a value of 'PARAM_MANDATORY'. The last argument to Register
  !% should be a pointer to the variable that the value of  the parameter should
  !% be copied to after parsing by 'param_read_line', 'param_read_file'
  !% or 'param_read_args'. For a parameter which shouldn't be parsed, do not specify a target.
  interface param_register
     module procedure param_register_single_integer
     module procedure param_register_multiple_integer
     module procedure param_register_single_real
     module procedure param_register_multiple_real
     module procedure param_register_single_string
     module procedure param_register_dontread
     module procedure param_register_single_logical
  end interface

#ifdef POINTER_COMPONENT_MANUAL_COPY
  interface assignment(=)
    module procedure ParamEntry_assign
  end interface
#endif

  contains

    ! Overloaded interface for registering parameter of type single logical
    subroutine param_register_single_logical(dict, key, value, logical_target, help_string, has_value_target)
      type(Dictionary), intent(inout) :: dict
      character(len=*), intent(in) :: key
      character(len=*), intent(in) :: value
      logical, intent(inout), target :: logical_target
      character(len=*), intent(in) :: help_string
      logical, intent(inout), optional, target :: has_value_target

      call param_register_main(dict, key, value, 1, PARAM_LOGICAL, help_string=help_string, &
	 logical_target=logical_target, has_value_target=has_value_target)

    end subroutine param_register_single_logical

    ! Overloaded interface for registering parameter of type single integer
    subroutine param_register_single_integer(dict, key, value, int_target, help_string, has_value_target)
      type(Dictionary), intent(inout) :: dict
      character(len=*), intent(in) :: key
      character(len=*), intent(in) :: value
      integer, intent(inout), target :: int_target
      character(len=*), intent(in) :: help_string
      logical, intent(inout), optional, target :: has_value_target

      call param_register_main(dict, key, value, 1, PARAM_INTEGER, help_string=help_string, &
	 int_target=int_target, has_value_target=has_value_target)

    end subroutine param_register_single_integer


    ! Overloaded interface for registering parameter of type multiple integer
    subroutine param_register_multiple_integer(dict, key, value, int_target_array, help_string, has_value_target)
      type(Dictionary), intent(inout) :: dict
      character(len=*), intent(in) :: key
      character(len=*), intent(in) :: value
      integer, dimension(:), intent(inout), target :: int_target_array
      character(len=*), intent(in) :: help_string
      logical, intent(inout), optional, target :: has_value_target

      call param_register_main(dict, key, value, size(int_target_array), PARAM_INTEGER, help_string=help_string, &
           int_target_array=int_target_array, has_value_target=has_value_target)

    end subroutine param_register_multiple_integer


    ! Overloaded interface for registering parameter of type single real
    subroutine param_register_single_real(dict, key, value, real_target, help_string, has_value_target)
      type(Dictionary), intent(inout) :: dict
      character(len=*), intent(in) :: key
      character(len=*), intent(in) :: value
      real(dp), intent(inout), target :: real_target
      character(len=*), intent(in) :: help_string
      logical, intent(inout), optional, target :: has_value_target

      call param_register_main(dict, key, value, 1, PARAM_REAL, help_string=help_string, &
	 real_target=real_target, has_value_target=has_value_target)

    end subroutine param_register_single_real


    ! Overloaded interface for registering parameter of type multiple real
    subroutine param_register_multiple_real(dict, key, value, real_target_array, help_string, has_value_target)
      type(Dictionary), intent(inout) :: dict
      character(len=*), intent(in) :: key
      character(len=*), intent(in) :: value
      real(dp), dimension(:), intent(inout), target :: real_target_array
      character(len=*), intent(in) :: help_string
      logical, intent(inout), optional, target :: has_value_target

      call param_register_main(dict, key, value, size(real_target_array), PARAM_REAL, help_string=help_string, &
           real_target_array=real_target_array, has_value_target=has_value_target)

    end subroutine param_register_multiple_real


    ! Overloaded interface for registering parameter of type single string
    subroutine param_register_single_string(dict, key, value, char_target, help_string, has_value_target)
      type(Dictionary), intent(inout) :: dict
      character(len=*), intent(in) :: key
      character(len=*), intent(in) :: value
      character(len=*), intent(inout), target :: char_target
      character(len=*), intent(in) :: help_string
      logical, intent(inout), optional, target :: has_value_target

      if (len(char_target) /= STRING_LENGTH) &
	call system_abort('param_register_single_string called for "'//trim(key)// &
	 '" has char_target(len='//len(char_target)//'), must be called with char_target(len=STRING_LENGTH)')

      call param_register_main(dict, key, value, 1, PARAM_STRING, help_string=help_string, &
	 char_target=char_target, has_value_target=has_value_target)

    end subroutine param_register_single_string

    ! Overloaded interface for registering parameters which don't need parsing
    subroutine param_register_dontread(dict, key)
      type(Dictionary), intent(inout) :: dict
      character(len=*), intent(in) :: key

      call param_register_main(dict, key, '', 0, PARAM_NO_VALUE, help_string="NOT PARSED")

    end subroutine param_register_dontread


    ! Helper routine called by all the registration interfaces
    !%OMIT
    subroutine param_register_main(dict, key, value, N, param_type,  &
         help_string, int_target, real_target, char_target,  &
         logical_target, int_target_array, real_target_array, has_value_target)
      type(Dictionary), intent(inout) :: dict

      character(len=*), intent(in) :: key
      character(len=*), intent(in) :: value
      integer, intent(in) :: N, param_type

      character(len=*), intent(in) :: help_string
      integer, intent(inout), optional, target :: int_target
      integer, dimension(:), intent(inout), optional, target :: int_target_array
      real(dp), intent(inout), optional, target :: real_target
      real(dp), dimension(:), intent(inout), optional, target :: real_target_array
      character(len=*), intent(inout), optional, target :: char_target
      logical, intent(inout), optional, target :: logical_target
      logical, intent(inout), optional, target :: has_value_target

      type(ParamEntry) :: entry
      type(DictData) :: data

      if (len_trim(value) > STRING_LENGTH) &
           call system_abort("Param_Register: Value "//trim(value)//" too long")

      entry%value = value
      entry%N = N
      entry%param_type = param_type

      entry%help_string=trim(help_string)

! call print("parser register entry " // trim(key) // " value " // trim(value), PRINT_ALWAYS)

      if (present(has_value_target)) then
	entry%has_value => has_value_target
      else
	nullify(entry%has_value)
      endif

      ! Check target type and size and set pointer to point to it
      select case (param_type)
      case(PARAM_NO_VALUE)
         ! Do nothing

      case (PARAM_LOGICAL)
	if (N == 1) then
	  if (present(logical_target)) then
	    entry%logical => logical_target
	  else
	    call system_abort('Param_Register: no target for single logical parameter')
	  endif
	else
	  call system_abort('Param_Register: no support for logical array parameters')
	endif
      case (PARAM_INTEGER)
         if (N == 1) then
            if (present(int_target)) then
               entry%integer => int_target
            else
               call system_abort('Param_Register: no target for single integer parameter')
            end if
         else
            if (present(int_target_array)) then
               if (size(int_target_array) == N) then
                  entry%integers => int_target_array
               else
                  call system_abort('Param_Register: integer target array wrong size')
               end if
            else
               call system_abort('Param_Register: no target for integer array parameter')
            end if
         end if

      case (PARAM_REAL)
         if (N == 1) then
            if (present(real_target)) then
               entry%real => real_target
            else
               call system_abort('Param_Register: no target for single real parameter')
            end if
         else
            if (present(real_target_array)) then
               if (size(real_target_array) == N) then
                  entry%reals => real_target_array
               else
                  call system_abort('Param_Register: real target array wrong size')
               end if
            else
               call system_abort('Param_Register: no target for real array parameter')
            end if
         end if


      case (PARAM_STRING)
         if (N == 1) then
            if (present(char_target)) then
               entry%string => char_target
            else
               call system_abort('Param_Register: no target for single char parameter')
            end if
         else
            call system_abort('Param_Register: multiple string param type not supported')
         end if

      case default
         write (line, '(a,i0)') 'Param_Register: unknown parameter type ', &
              entry%param_type
         call system_abort(line)

      end select

      ! Parse the default value string
      if (.not. param_parse_value(entry)) then
         call system_abort('Error parsing value '//trim(entry%value))
      end if

      ! Store the entry object in dictionary
      allocate(data%d(size(transfer(entry,data%d))))
      data%d = transfer(entry,data%d)
      call Set_Value(dict, key, data)

    end subroutine param_register_main


    !% Read and parse line and update the key/value pairs stored by
    !% in a Dictionary object. Returns false if an error is encountered parsing
    !% the line.
    function param_read_line(dict, myline, ignore_unknown, check_mandatory, task, did_help) result(status)

      type(Dictionary), intent(inout) :: dict !% Dictionary of registered key/value pairs
      character(len=*), intent(in) :: myline  !% Line to parse
      logical, intent(in), optional :: ignore_unknown !% If true, ignore unknown keys in line
      logical, intent(in), optional :: check_mandatory !% If true, check for missing mandatory parameters
      character(len=*), intent(in), optional :: task
      logical, intent(out), optional :: did_help
      logical :: status

      character(len=STRING_LENGTH) :: field
      integer equal_pos
      character(len=STRING_LENGTH), dimension(MAX_N_FIELDS) :: final_fields
      character(len=STRING_LENGTH) :: key, value
      integer :: i, num_pairs
      type(ParamEntry) :: entry
      type(DictData) :: data
      logical :: my_ignore_unknown
      integer :: entry_i
      logical, allocatable :: my_has_value(:)
      logical :: my_check_mandatory

      my_ignore_unknown=optional_default(.false., ignore_unknown)
      my_check_mandatory=optional_default(.true., check_mandatory)

      if(present(did_help)) did_help = .false.

      call split_string(myline," ,","''"//'""'//'{}', final_fields, num_pairs, matching=.true.)

      allocate(data%d(size(transfer(entry,data%d))))

      allocate(my_has_value(dict%N))
      my_has_value = .false.

      ! zero out has_value pointers for entire dictionary
      do i=1, dict%N
	entry = transfer(dict%entries(i)%d%d,entry)
	if (associated(entry%has_value)) then
	  entry%has_value = .false.
	endif
      end do

      if (present(task)) then
        call print("parser doing "//trim(task), PRINT_VERBOSE)
      else
        call print("parser doing UNKNOWN", PRINT_VERBOSE)
      endif
      call print("parser input line: <"//trim(myline)//">", PRINT_VERBOSE)
      ! Set the new entries
      do i=1,num_pairs
	 field = final_fields(i)
	 equal_pos = index(trim(field),'=')
	 if (equal_pos == 0) then
	  key=field
	  value=''
	 else if (equal_pos == 1) then
	  call system_abort("Malformed field '"//trim(field)//"'")
         else if (equal_pos == len(trim(field))) then
	  key = field(1:equal_pos-1)
	  value = ''
	 else
	  key = field(1:equal_pos-1)
	  value = field(equal_pos+1:len(trim(field)))
	 endif
	 call print("param_read_line key='"//trim(key)//"' value='"//trim(value)//"'", PRINT_NERD)
         if (len_trim(value) > STRING_LENGTH) then
            call print("param_read_line: value "//trim(value)//" too long")
            status = .false.
            return
         end if

	 if (trim(key) == "--help") then
            if(present(did_help)) did_help = .true.
            call param_print_help(dict)
            cycle
         endif

         ! Extract this value
         if (.not. get_value(dict, key, data, i=entry_i)) then
	    if (.not. my_ignore_unknown) then
	      call print("param_read_line: unknown key "//trim(key))
	      status = .false.
	      return
	    endif
	 else
	   entry = transfer(data%d,entry)
	   entry%value = value
           my_has_value(entry_i) = .true.
	   if (associated(entry%has_value)) then
	     entry%has_value = .true.
	   endif
           call print ("parser got: " // trim(paramentry_write_string(key, entry)), PRINT_VERBOSE)
	   if (equal_pos == 0) then
	      status = param_parse_value(entry, key)
	   else
	      status = param_parse_value(entry)
	   endif
	   if (.not. status) then
	      call Print('Error parsing value '//trim(entry%value))
	      return
	   end if

	   ! Put it back in dict
	   data%d = transfer(entry,data%d)
	   call set_value(dict,key,data)
         end if
      end do

      if (current_verbosity() >= PRINT_VERBOSE) then
        do i=1, dict%N
          entry = transfer(dict%entries(i)%d%d,entry)
          if (.not. my_has_value(i)) then
            call print ("parser defaulted: "//trim(paramentry_write_string(string(dict%keys(i)),entry)), PRINT_VERBOSE)
          endif
        end do
      endif

      if (allocated(data%d)) deallocate(data%d)
      status = .true. ! signal success to caller

      if (my_check_mandatory) status = param_check(dict)

    end function param_read_line


    !% Read lines from Inoutput object 'file', in format 'key = value' and set
    !% entries in the Dictionary 'dict'.  Skips lines starting with a '#'. If
    !% 'check_mandatory' is true then finish by checking if all mandatory values have
    !% been specified. Returns false if there was a problem reading the file, or
    !% if the check for mandatory parameters fails.
    function param_read_file(dict, file, check_mandatory,task) result(status)
      type(Dictionary), intent(inout) :: dict !% Dictionary of registered key/value pairs
      type(Inoutput), intent(in) :: file !% File to read from
      logical, intent(in), optional :: check_mandatory !% Should we check if all mandatory parameters have been given
      character(len=*), intent(in), optional :: task
      logical :: status

      integer :: file_status
      character(STRING_LENGTH) :: myline
      logical :: my_check_mandatory

      my_check_mandatory=optional_default(.true., check_mandatory)

      do
         myline = read_line(file, file_status)
         if (file_status /= 0) exit  ! Check for end of file
         ! Discard empty and comment lines, otherwise parse the line
         if (len_trim(myline) > 0 .and. myline(1:1) /= '#') then
            status = param_read_line(dict, myline, task=task)
            if (.not. status) then
               call print('Error reading line: '//myline)
               return
            end if
         end if
      end do

      status = .true.
      ! Should we check if all mandatory params have been specified?
      if (my_check_mandatory) status = param_check(dict)

    end function param_read_file


    !% Read 'key = value' pairs from command line and set entries in this
    !% 'dict'. Array 'args' is a list of the indices of command line
    !% arguments that we should look at, in order, if it's not given we look
    !% at all arguments.  Returns false if fails, or if optional check that
    !% all mandatory values have been specified fails.
    function param_read_args(dict, args, check_mandatory,ignore_unknown,task,command_line, did_help) result(status)
      type(Dictionary), intent(inout) :: dict !% Dictionary of registered key/value pairs
      integer, dimension(:), intent(in), optional :: args !% Argument indices to use
      logical, intent(in), optional :: check_mandatory !% Should we check if all mandatory parameters have been given
      logical, intent(in), optional :: ignore_unknown !% If true, ignore unknown keys in line
      character(len=*), intent(in), optional :: task
      character(len=*), intent(out), optional :: command_line
      logical, intent(out), optional :: did_help
      logical :: status

      integer :: i, nargs
      character(len=STRING_LENGTH) :: this_arg
      character(len=STRING_LENGTH) :: my_command_line
      integer, dimension(:), allocatable :: xargs
      logical :: my_ignore_unknown
      integer :: eq_loc, this_len
      logical :: my_check_mandatory

      my_ignore_unknown=optional_default(.false., ignore_unknown)
      my_check_mandatory=optional_default(.true., check_mandatory)

      nargs = cmd_arg_count()

      ! If args was specified, make a local copy of it
      ! Otherwise make a list of all command line arguments
      if (present(args)) then
         allocate(xargs(size(args)))
         xargs = args
      else

         allocate(xargs(nargs))
         xargs = (/ (i, i=1,size(xargs)) /)
      end if

      ! Concatentate command line options into one string
      my_command_line = ''
      do i=1, size(xargs)
         call get_cmd_arg(xargs(i), this_arg)
         if (index(trim(this_arg),' ') /= 0) then
           if (index(trim(this_arg),'=') /= 0) then
             eq_loc = index(trim(this_arg),'=')
             this_len = len_trim(this_arg)
	     if (scan(this_arg(eq_loc+1:eq_loc+1),"'{"//'"') <= 0) then
	       my_command_line = trim(my_command_line)//' '//this_arg(1:eq_loc)//'"'//this_arg(eq_loc+1:this_len)//'"'
	     else
	       my_command_line = trim(my_command_line)//' '//trim(this_arg)
	     endif
           else
	     if (scan(this_arg(1:1),"{'"//'"') <= 0) then
	       my_command_line = trim(my_command_line)//' "'//trim(this_arg)//'"'
	     else
	       my_command_line = trim(my_command_line)//' '//trim(this_arg)
	     endif
           endif
         else
           my_command_line = trim(my_command_line)//' '//trim(this_arg)
         endif
      end do

      call print("param_read_args got command_line '"//trim(my_command_line)//"'", PRINT_VERBOSE)
      ! Free local copy before there's a chance we return
      deallocate(xargs)

      ! Then parse this string, and return success or failure
      status = param_read_line(dict, my_command_line,ignore_unknown=my_ignore_unknown, task=task, did_help=did_help)
      if (.not. status) return

      status = .true.
      ! Should we check if all mandatory params have been specified?
      if (my_check_mandatory) status = param_check(dict)

      if(present(command_line)) command_line = trim(my_command_line)

    end function param_read_args


    !% Process command line by looking for arguments in 'key=value' form and
    !% parsing their values into 'dict'. Non option arguments
    !% are returned in the array non_opt_args. Optionally check that
    !% all mandatory parameters have been specified.
    function process_arguments(dict, non_opt_args, n_non_opt_args, check_mandatory, task) result(status)
      type(Dictionary) :: dict !% Dictionary of registered key/value pairs
      character(len=*), dimension(:), intent(out) :: non_opt_args !% On exit, contains non option arguments
      integer, intent(out) :: n_non_opt_args !% On exit, number of non option arguments
      logical, intent(in), optional :: check_mandatory !% Should we check if all mandatory parameters have been given
      character(len=*), intent(in), optional :: task
      logical :: status

      character(len=255) :: args(100)
      character(len=1000) :: options
      type(Table) :: opt_inds
      integer :: start, end, i, j, eqpos, n_args
      logical :: my_check_mandatory

      n_args = cmd_arg_count()
      if (n_args > size(args)) call system_abort('Too many command line arguments')

      my_check_mandatory=optional_default(.true., check_mandatory)

      ! Read command line args into array args
      do i=1,n_args
	 call get_cmd_arg(i, args(i))
      end do

      ! Find the key=value pairs, skipping over non-option arguments
      call allocate(opt_inds, 1, 0, 0, 0)
      i = 1
      do
         eqpos = scan(args(i), '=')

         if (eqpos == 0) then
            i = i + 1
            if (i > n_args) exit
            cycle  ! Skip non option arguments
         else if (eqpos == 1) then
            ! Starts with '=', so previous argument is key name
            if (i == 1) call system_abort('Missing option name in command line arguments')
            start = i-1
         else if (eqpos >  1) then
            ! Contains '=' but not at start so key name is in this arg
            start = i
         end if
         end =  i

         ! If last character of arg is '=' then value must be in next argument
         eqpos = len_trim(args(i))
         if (args(i)(eqpos:eqpos) == '=') then
            end = end + 1
         end if

         ! Get whole of quoted string
         if (scan(args(end), '"''') /= 0) then
            do
               end = end + 1
               if (end > n_args) &
                    call system_abort('Mismatched quotes in command line arguments')
               if (scan(args(end), '"''') /= 0) exit
            end do
         end if

         ! Add the indices of the arguments that make up this key/value pair
         call append(opt_inds, (/ (j, j=start,end )/))

         i = end + 1
         if (i > n_args) exit
      end do

      ! Concatenate args and parse line
      options = ''
      do i=1,opt_inds%N
         options = trim(options)//' '//trim(args(opt_inds%int(1,i)))
      end do
      status = param_read_line(dict, options, task=task)
      if (.not. status) return

      status = .true.
      ! Should we check if all mandatory params have been specified?
      if (my_check_mandatory) status = param_check(dict)

      ! Find and return the non-option arguments
      n_non_opt_args = 0
      do i=1,n_args
         if (Find(opt_inds, i) == 0) then
            n_non_opt_args = n_non_opt_args + 1
            if (n_non_opt_args > size(non_opt_args)) call system_abort('Too many non option arguments')

            non_opt_args(n_non_opt_args) = args(i)
         end if
      end do

      call finalise(opt_inds)

    end function process_arguments


    !% Print key/value pairs.
    subroutine param_print(dict, verbosity, out)
      type(Dictionary), intent(in) :: dict !% Dictionary of registered key/value pairs
      integer, intent(in), optional :: verbosity
      type(Inoutput), intent(inout), optional :: out

      type(ParamEntry) :: entry
      integer :: i

      type(DictData) :: data
      allocate(data%d(size(transfer(entry,data%d))))

      do i=1,dict%N
         if (.not. get_value(dict, string(dict%keys(i)), data)) &
              call system_abort('param_print: Key '//string(dict%keys(i))//' missing')
         entry = transfer(data%d,entry)
         ! Print value in quotes if it contains any spaces
         if (index(trim(entry%value), ' ') /= 0) then
            write (line, '(a,a,a,a)')  &
                 string(dict%keys(i)), ' = "', trim(entry%value),'" '
         else
            write (line, '(a,a,a,a)')  &
                 string(dict%keys(i)), ' = ', trim(entry%value), ' '
         end if
         call print(line,verbosity,out)
      end do

      deallocate(data%d)

    end subroutine param_print

    subroutine param_print_help(dict, verbosity, out)
      type(Dictionary), intent(in) :: dict !% Dictionary of registered key/value pairs
      integer, intent(in), optional :: verbosity
      type(Inoutput), intent(inout), optional :: out

      type(ParamEntry) :: entry
      integer :: i

      type(DictData) :: data
      character(len=10) type_str, param_dim_str
      character(len=STRING_LENGTH) out_line

      allocate(data%d(size(transfer(entry,data%d))))

      do i=1,dict%N
         if (.not. get_value(dict, string(dict%keys(i)), data)) &
              call system_abort('param_print: Key '//string(dict%keys(i))//' missing')
         entry = transfer(data%d,entry)

	 out_line=string(dict%keys(i)//" type=")

	 param_dim_str = "unknown"
	 select case(entry%param_type)
	    case(PARAM_NO_VALUE)
	       out_line=trim(out_line)//"NO_VALUE"
	    case(PARAM_REAL)
	       out_line=trim(out_line)//"REAL"
	       if (associated(entry%real)) then
		  out_line=trim(out_line)//" scalar"
	       elseif (associated(entry%reals)) then
		  out_line=trim(out_line)//" dim="//size(entry%reals)
	       else
		  out_line=trim(out_line)//" ERROR-no-value-pointer-associated"
	       endif
	    case(PARAM_INTEGER)
	       out_line=trim(out_line)//"INTEGER"
	       if (associated(entry%integer)) then
		  out_line=trim(out_line)//" scalar"
	       elseif (associated(entry%integers)) then
		  out_line=trim(out_line)//" dim="//size(entry%integers)
	       else
		  out_line=trim(out_line)//" ERROR-no-value-pointer-associated"
	       endif
	    case(PARAM_STRING)
	       out_line=trim(out_line)//"STRING"
	       if (associated(entry%string)) then
		  out_line=trim(out_line)//" scalar"
	       else
		  out_line=trim(out_line)//" ERROR-no-value-pointer-associated"
	       endif
	    case(PARAM_LOGICAL)
	       out_line=trim(out_line)//"LOGICAL"
	       if (associated(entry%logical)) then
		  out_line=trim(out_line)//" scalar"
	       else
		  out_line=trim(out_line)//" ERROR-no-value-pointer-associated"
	       endif
	    case default
	       call system_abort("Unknown param_type "//entry%param_type//" for key " //string(dict%keys(i)))
	 end select

         ! Print value in quotes if it contains any spaces
         if (index(trim(entry%value), ' ') /= 0) then
	    out_line=trim(out_line)//" current_value='"//trim(entry%value)//"'"
         else
	    out_line=trim(out_line)//" current_value="//trim(entry%value)
         end if
	 call print(trim(out_line), verbosity=PRINT_ALWAYS, file=out)

	 call print(linebreak_string(trim(entry%help_string), 80), verbosity=PRINT_ALWAYS, file=out)

	 call print("", verbosity=PRINT_ALWAYS, file=out)
      end do

      deallocate(data%d)

    end subroutine param_print_help

    function paramentry_write_string(key,entry) result(s)
      character(len=*) :: key
      type(ParamEntry), intent(in) :: entry
      character(len=STRING_LENGTH) :: s

      ! Print value in brackets if it contains any spaces
      if (index(trim(entry%value), ' ') /= 0) then
         write (s, '(a,a,a,a)')  &
              trim(key), '={', trim(entry%value),'} '
      else
         if (len(trim(entry%value)) == 0 .and. entry%param_type == PARAM_LOGICAL) then
           write (s, '(a,a,a,a)')  &
                trim(key), '=', 'T', ' '
         else
           write (s, '(a,a,a,a)')  &
                trim(key), '=', trim(entry%value), ' '
         endif
      end if
    end function paramentry_write_string

    !% Return string representation of param dictionary, i.e.
    !% 'key1=value1 key2=value2 quotedkey="quoted value"'
    function param_write_string(dict) result(s)
      type(Dictionary), intent(in) :: dict !% Dictionary of registered key/value pairs
      character(STRING_LENGTH) :: s

      type(ParamEntry) :: entry
      integer :: i

      type(DictData) :: data
      allocate(data%d(size(transfer(entry,data%d))))

      s=""
      do i=1,dict%N
         if (.not. get_value(dict, string(dict%keys(i)), data)) &
              call system_abort('param_print: Key '//string(dict%keys(i))//' missing')
         entry = transfer(data%d,entry)
         s = trim(s)//trim(paramentry_write_string(string(dict%keys(i)),entry))
      end do

      deallocate(data%d)

    end function param_write_string


    ! Parse the value according to its param_type and number of values
    ! and set the targets to the values read. Called by Register and ReadLine
    !%OMIT
    function param_parse_value(entry, key) result(status)
      type(ParamEntry) :: entry
      character(len=*), optional :: key
      logical :: status

      character(len=STRING_LENGTH), dimension(MAX_N_FIELDS) :: fields
      integer :: num_fields, j

      if (entry%param_type == PARAM_NO_VALUE .or. &
           trim(entry%value) == PARAM_MANDATORY) then
         status = .true.
         return
      end if

      call parse_string(entry%value, ' ', fields, num_fields)

      ! jrk33 - commented out these lines as they stop zero length string
      ! values from setting the associated target to an empty string, and
      ! I can't think what good they do!
!!$      if (num_fields == 0 .and. entry%param_type /= PARAM_LOGICAL) then
!!$         status = .true.
!!$         return
!!$      end if

      if ((entry%param_type == PARAM_LOGICAL .and. num_fields /= 0 .and. num_fields /= entry%N) .or. &
          (entry%param_type /= PARAM_STRING .and. entry%param_type /= PARAM_LOGICAL .and. num_fields /= entry%N) ) then
         write (line, '(a,a,a)') 'Param_ParseValue: Number of value fields wrong: ', &
              trim(entry%value)
         call print(line)
         status = .false.
         return
        end if

      select case(entry%param_type)
      case (PARAM_LOGICAL)
         if (entry%N == 1) then
	    if (num_fields == 0) then
	      entry%logical = .true.
	    else
	      entry%logical = String_to_Logical(fields(1))
	    endif
         else
	    call system_abort("param_parse_value no support for logical array yet")
         end if

      case (PARAM_INTEGER)
         if (entry%N == 1) then
            entry%integer = String_to_Int(fields(1))
         else
            do j=1,num_fields
               entry%integers(j) = String_to_Int(fields(j))
            end do
         end if

      case (PARAM_REAL)
         if (entry%N == 1) then
            entry%real = String_to_Real(fields(1))
         else
            do j=1,num_fields
               entry%reals(j) = String_to_Real(fields(j))
            end do
         end if

      case (PARAM_STRING)
	 if (num_fields == 0 .and. present(key)) then
	    entry%string = trim(key)
	 else
	    entry%string = trim(entry%value)
	 endif

      case default
         write (line, '(a,i0)') 'Param_ParseValue: unknown parameter type ', &
              entry%param_type
         call print(line)
         status = .false.
         return
      end select

      status = .true.

    end function param_parse_value


    !% Explicity check if all mandatory values have now been specified
    !% Returns true or false accordingly. Optionally return the list
    !% of missing keys as a string.
    function param_check(dict, missing_keys) result(status)
      type(Dictionary), intent(in) :: dict !% Dictionary of registered key/value pairs
      character(len=*), intent(out), optional :: missing_keys !% On exit, string list of missing mandatory keys,
                                                              !% separated by spaces.
      logical :: status
      integer :: i
      type(ParamEntry) :: entry
      type(DictData) :: data
      allocate(data%d(size(transfer(entry,data%d))))

      if (present(missing_keys)) missing_keys = ''

      status = .true.
      do i=1,dict%N
         if (.not. get_value(dict, string(dict%keys(i)), data)) &
              call system_abort('Param_Check: Key '//string(dict%keys(i))//' missing')
         entry = transfer(data%d,entry)
         if (trim(entry%value) == PARAM_MANDATORY) then
            status = .false.
            if (present(missing_keys)) then
               missing_keys = trim(missing_keys)//' '//string(dict%keys(i))
            else
               call print('missing mandatory parameter '//string(dict%keys(i)))
            endif
         end if
      end do

      deallocate(data%d)

    end function param_check

#ifdef POINTER_COMPONENT_MANUAL_COPY
  subroutine ParamEntry_assign(to, from)
    type(ParamEntry) :: to, from

    to%value = from%value
    to%N = from%N
    to%param_type = from%param_type
    to%real => from%real
    to%integer => from%integer
    to%logical => from%logical
    to%string => from%string
    to%reals => from%reals
    to%integers => from%integers

    to%has_value => from%has_value
  end subroutine ParamEntry_assign
#endif

end module paramreader_module
