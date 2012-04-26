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
!X  System module
!X  
!X  Basic system dependent functionality:
!X  
!X  mpi constants, default output objects, printing
!X  random number generators
!X 
!% The system module contains low-level routines for I/O, timing, random
!% number generation etc. The Inoutput type is used to abstract both
!% formatted and unformatted (i.e. binary) I/O.
!X
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#include "error.inc"

module system_module
  use error_module
!$ use omp_lib
  implicit none

  logical :: system_always_flush = .false.

  logical :: system_use_fortran_random = .false.

#ifdef HAVE_QP
  integer, parameter :: qp = 16 
#else
  integer, parameter :: qp = 8
#endif
  
#ifdef QUAD_PRECISION
  integer, parameter :: dp = 16 ! kind(1.0d0)
#else
  integer, parameter :: dp = 8 ! kind(1.0d0)
#endif

  character :: quip_new_line

  integer, parameter :: INTEGER_SIZE = 4
  integer, parameter :: REAL_SIZE = dp
  integer, parameter :: COMPLEX_SIZE = 2*dp
  logical :: trace_memory = .false.
  integer :: traced_memory = 0

  logical, private :: system_do_timing = .false.
  logical, private :: system_quippy_running = .false.

  type Stack
    integer:: pos
    integer, allocatable :: val(:)
  end type Stack

  type InOutput
     integer:: unit
     character(256)::filename
     character(256)::prefix, postfix
     integer::default_real_precision
     logical::formatted
     logical::append
     logical::active         !% Does it print?
     integer::action
     logical::mpi_all_inoutput_flag = .false.
     logical::mpi_print_id = .false.
     type(Stack) :: verbosity_stack, verbosity_cascade_stack
     logical::initialised = .false.
  end type InOutput

  type allocatable_array_pointers
    integer, allocatable :: i_a(:)
    real(dp), allocatable :: r_a(:)
    complex(dp), allocatable :: c_a(:)
    logical, allocatable :: l_a(:)
  end type allocatable_array_pointers

  public   !standard setting for the module
  integer,private                  :: mpi_n, mpi_myid    ! Number of processes and local process ID
  real(dp),private                 :: start_time         ! Initial time
  integer, parameter, private      :: SYSTEM_STRING_LENGTH = 1024 !max line length read
  integer, parameter, private      :: SYSTEM_STRING_LENGTH_LONG = 102400 !max line length read
  character(10240)                 :: line               ! 'line' is global and is used by other modules
  character(10240),private         :: local_line         ! 'local_line' is private and System should use this instead         
  type(inoutput),target,save      :: mainlog            !% main output, connected to 'stdout' by default
  type(inoutput),target,save      :: errorlog           !% error output, connected to 'stderr' by default
  type(inoutput),target,save      :: mpilog             !% MPI output, written to by each mpi process
  integer,private                  :: idum               ! used in the random generator
  !$omp threadprivate(idum)
  real(dp),parameter               :: NUMERICAL_ZERO = 1.0e-14_dp

  ! system dependent variables 
  integer::RAN_MAX


  ! output labels
  integer,parameter::PRINT_ALWAYS   = -100000
  integer,parameter::PRINT_SILENT  =  -1
  integer,parameter::PRINT_NORMAL  =   0
  integer,parameter::PRINT_VERBOSE =   1
  integer,parameter::PRINT_NERD    =   1000  ! aleph0
  integer,parameter::PRINT_ANAL    =   10000 ! aleph1

  integer,parameter::INPUT=0
  integer,parameter::OUTPUT=1
  integer,parameter::INOUT=2 


  ! random number generator parameters
  integer,parameter::ran_A=16807
  integer,parameter::ran_M=2147483647
  integer,parameter::ran_Q=127773
  integer,parameter::ran_R=2836

  ! System_Timer stack size
  integer, parameter :: TIMER_STACK  = 500

  ! Command argument variables
  integer,              save :: NUM_COMMAND_ARGS  = 0  !% The number of arguments on the command line
  integer, parameter         :: MAX_READABLE_ARGS = 100 !% The maximum number of arguments that will be read
  character(255),       save :: EXEC_NAME              !% The name of the executable
  character(2550), dimension(MAX_READABLE_ARGS), save :: COMMAND_ARG !% The first 'MAX_READABLE_ARGS' command arguments

  private :: inoutput_initialise
  interface initialise
     module procedure inoutput_initialise
  end interface initialise

  private :: inoutput_finalise
  interface finalise
     module procedure inoutput_finalise
  end interface finalise

  private :: inoutput_activate
  interface activate
     module procedure inoutput_activate
  end interface activate

  private :: inoutput_deactivate
  interface deactivate
     module procedure inoutput_deactivate
  end interface deactivate

  interface mpi_all_inoutput
     module procedure inoutput_mpi_all_inoutput
  end interface mpi_all_inoutput

  interface print_mpi_id
     module procedure inoutput_print_mpi_id
  end interface print_mpi_id

  private :: inoutput_print_string
  private :: inoutput_print_integer, inoutput_print_real, inoutput_print_logical

  !% Overloaded interface for printing. With the
  !% 'this' parameter omitted output goes to the default mainlog ('stdout'). The
  !% 'verbosity' parameter controls whether the object is actually printed;
  !% if the verbosity is greater than that currently at the top of the
  !% verbosity stack then output is suppressed. Possible verbosity levels
  !% range from 'ERROR' through 'NORMAL', 'VERBOSE', 'NERD' and 'ANAL'.
  !% Other user-defined types define the Print interface in the same way.
  interface print
     module procedure inoutput_print_string
     module procedure inoutput_print_integer, inoutput_print_real, inoutput_print_logical
     module procedure inoutput_print_char_array
  end interface print

  private :: reada_real_dim1, reada_int_dim1
  interface read_ascii
    module procedure reada_real_dim1, reada_int_dim1
  end interface read_ascii

  private :: inoutput_read_line
  interface read_line
    module procedure inoutput_read_line
  end interface read_line

  private :: inoutput_read_file
  interface read_file
    module procedure inoutput_read_file
  end interface read_file

  private :: inoutput_parse_line
  interface parse_line
     module procedure inoutput_parse_line
  end interface parse_line

  private :: reallocate_int1d, reallocate_int2d, reallocate_real1d, reallocate_real2d
  interface reallocate
     module procedure reallocate_int1d, reallocate_int2d, reallocate_real1d, reallocate_real2d
  end interface reallocate

  interface operator(//)
     module procedure string_cat_logical, string_cat_int, string_cat_real, string_cat_real_array
     module procedure string_cat_complex, string_cat_int_array, string_cat_logical_array
     module procedure string_cat_complex_array, string_cat_string_array
!     module procedure logical_cat_string, logical_cat_logical, logical_cat_int, logical_cat_real
     module procedure int_cat_string!, int_cat_logical, int_cat_int, int_cat_real
     module procedure real_cat_string!, real_cat_logical, real_cat_int, real_cat_real
     module procedure real_array_cat_string
  end interface

  interface system_command
     !%>   call system_command(command)
     !% Interface to a C wrapper to the 'system(3)' system call for
     !% executing external programs. Command can only be up to 1024 characters long.
     !% \begin{description}
     !%   \item['character' --- character(*)]
     !% \end{description}
     subroutine system_command(command,status,error)
       character(*), intent(in) :: command
       integer, optional, intent(out) :: status
       integer, optional, intent(out) :: error
     end subroutine system_command
  end interface

#ifdef NO_FORTRAN_ISNAN
  INTERFACE 	
     elemental function fisnan(r)
       real(8), intent(in)::r
       integer::fisnan
     end function fisnan
  end INTERFACE
#endif

  interface c_mem_info
     subroutine c_mem_info(total_mem,free_mem)
       real(8), intent(out) :: total_mem, free_mem
     endsubroutine c_mem_info
  endinterface c_mem_info

  private :: Stack_Initialise
  interface Initialise
    module procedure Stack_Initialise
  end interface Initialise

  private :: Stack_Finalise
  interface Finalise
    module procedure Stack_Finalise
  end interface Finalise

  private :: Stack_push
  interface push
    module procedure Stack_push
  end interface push

  private :: Stack_pop
  interface pop
    module procedure Stack_pop
  end interface pop

  private :: Stack_value
  interface value
    module procedure Stack_value
  end interface value

  private :: Stack_Print
  interface Print
    module procedure Stack_Print
  end interface Print

  !% takes as arguments a default value and an optional argument, and 
  !% returns the optional argument value if it's present, otherwise
  !% the default value
  private :: optional_default_l, optional_default_i, optional_default_r
  private :: optional_default_c, optional_default_z
  private :: optional_default_ia, optional_default_ra
  interface optional_default
    module procedure optional_default_l, optional_default_i, optional_default_r
    module procedure optional_default_c, optional_default_z
    module procedure optional_default_ia, optional_default_ra
  end interface optional_default

  integer, external :: pointer_to

contains

#ifdef NO_FORTRAN_ISNAN
  elemental function isnan(r)
    real(dp), intent(in)::r
    logical::isnan
    select case(fisnan(r))
       case(0)
          isnan = .false.
       case(1)
          isnan = .true.
       case default
          isnan = .true.
     end select
  end function isnan
#endif

  !% Open a file for reading or writing. The action optional parameter can
  !% be one of 'INPUT' (default), 'OUTPUT' or 'INOUT'.
  !% For unformatted output, the
  !% 'isformatted' optional parameter must
  !% be set to false.
  subroutine inoutput_initialise(this,filename,action,isformatted,append,verbosity,verbosity_cascade,master_only,error)
    type(Inoutput), intent(inout)::this
    character(*),intent(in),optional::filename
    logical,intent(in),optional::isformatted 
    integer,intent(in),optional::action
    logical,intent(in),optional::append
    integer,intent(in),optional :: verbosity, verbosity_cascade
    logical,intent(in),optional :: master_only
    integer,intent(out),optional :: error

    character(32)::formattedstr
    character(32)::position_value
    integer::stat
    logical :: my_master_only

    INIT_ERROR(error)

    ! Default of optional parameters------------------------------- 

    my_master_only = optional_default(.false., master_only)

    if(present(isformatted)) then 
       this%formatted=isformatted
    else 
       this%formatted=.true.
    endif

    if(.NOT.(this%formatted)) then
       formattedstr = 'unformatted'
    else
       formattedstr = 'formatted' 
    end if

    if (present(action)) then
       this%action=action
    else
       this%action=INPUT
    end if

    if (present(append)) then 
       this%append=append
       if(append) then 
          position_value="APPEND"
       else
          position_value="REWIND"
       end if

    else 
       this%append=.false.
       position_value="REWIND"
    end if

    if(this%append .AND. .NOT.this%formatted) then
       RAISE_ERROR(" Append not implemented for unformatted output",error)
    end if

    if (present(filename)) then
       if (filename.eq.'stderr') then
          this%filename=filename
          this%unit=0 !standard error
          this%action=OUTPUT

       else if (filename.eq.'stdout') then
          this%filename=filename
          this%unit = 6 !standard output
          this%action=OUTPUT

       else if (filename.eq.'stdin')  then
          this%filename=filename 
          this%unit = 5 !standard input
          this%action=INPUT

       else
          this%filename=filename
          this%unit=pick_up_unit()! pick up a unit from available(7-99)

          ! actually open the unit
	  if ((.not. my_master_only) .or. mpi_myid == 0) then
	    if (this%action == INPUT) then
	      open(unit=this%unit,file=filename,form=formattedstr,position=position_value,status='OLD',iostat=stat)
	    else
	      open(unit=this%unit,file=filename,form=formattedstr,position=position_value,iostat=stat)
	    endif
	  else
	    stat = 0
	  endif
          if(stat.NE.0)then 
             RAISE_ERROR('IO error opening "'//trim(filename)//'" on unit '//this%unit//', error number: '//stat,error)
          end if

       end if
    else ! no file name passed, so we go to standard output
       this%filename='stdout' 
       this%unit = 6 !standard output 
       this%action = OUTPUT        
    end if  ! end default--------------------------------------------------

    this%prefix = ''
    this%postfix = ''
    this%default_real_precision = 17

    call initialise(this%verbosity_stack)
    if (present(verbosity)) then
      call push(this%verbosity_stack, verbosity)
    else
      call push(this%verbosity_stack, PRINT_NORMAL)
    endif

    call initialise(this%verbosity_cascade_stack)
    if (present(verbosity_cascade)) then
      call push(this%verbosity_cascade_stack, verbosity_cascade)
    else
      call push(this%verbosity_cascade_stack, 0)
    endif

    if ((.not. my_master_only) .or. mpi_myid == 0) then
      call activate(this)  ! now it is active
    endif

    this%initialised = .true.

  end subroutine inoutput_initialise


  !% OMIT
  function pick_up_unit() result(unit)
    integer::unit,i
    logical::iopened
    do i=7,99
       INQUIRE(i,opened=iopened)
       if(.NOT.iopened) then 
          unit=i
          exit
       end if
    end do
  end function pick_up_unit

  !% Deactivate an Inoutput object temporarily.
  subroutine inoutput_deactivate(this)
    type(Inoutput),intent(inout)::this
    this%active=.false.
  end subroutine inoutput_deactivate

  !% Activate an Inoutput object temporarily.
  subroutine inoutput_activate(this)
    type(Inoutput),intent(inout)::this
    this%active=.true.
  end subroutine inoutput_activate


  !% Cleans everything and set members to default 
  subroutine inoutput_finalise(this)
    type(Inoutput), intent(inout)::this
    if (this%unit .ge. 7)  close(this%unit)
    call finalise(this%verbosity_stack)
    call finalise(this%verbosity_cascade_stack)
    call deactivate(this)
    this%initialised = .false.
  end subroutine inoutput_finalise

  !% Close file but don't finalise this Inoutput
  subroutine inoutput_close(this)
    type(Inoutput), intent(inout) :: this

    if (this%unit .ge. 7) close(this%unit)
    call deactivate(this)
  end subroutine inoutput_close

  subroutine inoutput_mpi_all_inoutput(this,value)

    type(inoutput),    intent(inout) :: this
    logical, optional, intent(in)    :: value

    if (present(value)) then
       this%mpi_all_inoutput_flag = value
    else
       this%mpi_all_inoutput_flag = .true.
    end if

  end subroutine inoutput_mpi_all_inoutput

  subroutine inoutput_print_mpi_id(this,value)

    type(inoutput),    intent(inout) :: this
    logical, optional, intent(in)    :: value

    if (present(value)) then
       this%mpi_print_id = value
    else
       this%mpi_print_id = .true.
    end if

  end subroutine inoutput_print_mpi_id

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X
  !X Printing routines for intrinsic types
  !X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine inoutput_print_string(string, verbosity, file, nocr, do_flush)
    character(*),             intent(in) :: string
    integer,  optional,       intent(in) :: verbosity
    type(Inoutput), optional, target, intent(in) :: file
    logical, optional, intent(in) :: nocr, do_flush

    type(Inoutput), pointer :: myoutput
    integer :: myverbosity
    character(len=2) :: nocr_fmt
    logical :: my_do_flush

    my_do_flush = optional_default(system_always_flush, do_flush)

    nocr_fmt = ''
    if (present(nocr)) then
      if (nocr) nocr_fmt = ',$'
    endif

    ! check inoutput object
    myoutput => mainlog
    if(present(file)) myoutput => file

    ! check verbosity request
    myverbosity = PRINT_NORMAL
    if(present(verbosity)) myverbosity = verbosity

    ! if we are not active, do nothing
    if(.not.myoutput%active) return
    ! if request is above threshold, do nothing
    if(myverbosity > value(myoutput%verbosity_stack) ) return

    if(myoutput%action .EQ. INPUT) then
      call system_abort("inoutput_print: you cannot print to an INPUT object")
    end if

    if(.NOT.myoutput%formatted) then
        call system_abort("inoutput_print: this subroutine is not good for unformatted printing")
    end if        

    ! actually write the line, removing trailing blanks
    if (inoutput_do_output(myoutput)) then
      if (len_trim(myoutput%prefix) == 0) then
	if (myoutput%mpi_all_inoutput_flag .and. myoutput%mpi_print_id) then
	  write(myoutput%unit,'(i0,": ",a'//trim(nocr_fmt)//')') mpi_id(), trim(string)//trim(myoutput%postfix)
	else
	  write(myoutput%unit,'(a'//trim(nocr_fmt)//')') trim(string)//trim(myoutput%postfix)
	endif
      else
	if (myoutput%mpi_all_inoutput_flag .and. myoutput%mpi_print_id) then
	  write(myoutput%unit,'(i0,": ",a'//trim(nocr_fmt)//')') mpi_id(), trim(myoutput%prefix)//" "//trim(string) &
	    // " " // trim(myoutput%postfix)
	else
	  write(myoutput%unit,'(a'//trim(nocr_fmt)//')') trim(myoutput%prefix)//" "//trim(string)// " "  // trim(myoutput%postfix)
	endif
      endif
    endif

    if (my_do_flush) call flush(myoutput%unit)
  end subroutine inoutput_print_string

  function inoutput_do_output(this)
    type(inoutput), intent(in) :: this
    logical :: inoutput_do_output

    if (this%mpi_all_inoutput_flag .or. mpi_id() == 0) then
      inoutput_do_output = .true.
    else
      inoutput_do_output = .false.
    end if
  end function inoutput_do_output

  subroutine inoutput_print_char_array(char_a, verbosity, file)
    character(len=*) :: char_a(:)
    integer, optional,        intent(in) :: verbosity
    type(Inoutput), optional, target, intent(in) :: file

    integer i
    character(len=size(char_a)) :: str

    do i=1, size(char_a)
      str(i:i) = char_a(i)
    end do

    call print(str, verbosity, file)
  end subroutine inoutput_print_char_array

  subroutine inoutput_print_logical(log, verbosity, file)
    logical,           intent(in) :: log
    integer, optional, intent(in) :: verbosity
    type(Inoutput), optional, intent(in) :: file

    write(local_line,'(l1)') log
    call print(local_line, verbosity, file)
  end subroutine inoutput_print_logical

  subroutine inoutput_print_integer(int, verbosity, file)
    integer,           intent(in) :: int
    integer, optional, intent(in) :: verbosity
    type(Inoutput), optional, intent(in) :: file

    write(local_line,'(i0)') int
    call print(local_line, verbosity, file)
  end subroutine inoutput_print_integer

  subroutine inoutput_print_real(real, verbosity, file, precision, format)
    real(dp),          intent(in) :: real
    integer, optional, intent(in) :: verbosity
    integer, optional, intent(in) :: precision ! number of decimal places
    character(*), optional, intent(in) :: format
    type(inoutput), optional, intent(in) :: file

    character(7) :: myformat

    if(present(format)) then
       write(local_line, format) real
    else
       if (present(precision)) then
          if (precision > 99) then
             call print_warning('Inoutput_Print_Real: Precision too high. Capping to 99.')
             write(myformat,'(a)')'(f0.99)'
          else
             write(myformat,'(a,i0,a)')'(f0.',precision,')'
          end if
       else
          if(present(file)) then
             write(myformat,'(a,i0,a)')'(f0.',file%default_real_precision,')'
          else
             write(myformat,'(a,i0,a)')'(f0.',mainlog%default_real_precision,')'
          end if
       end if

       write(local_line,myformat) real
    end if

    call print(local_line, verbosity, file)
  end subroutine inoutput_print_real


!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X  Pretty print a title
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  !
  !% Print a centred title, like this:
  !%
  !% '==================================== Title ====================================='
  !
  subroutine print_title(title, verbosity)

    character(*), intent(in) :: title
    integer, intent(in), optional :: verbosity

    character(len(title))    :: my_title
    integer                  :: length, a,b

    my_title = adjustl(title)
    length = len_trim(my_title)

    call print('', verbosity)
    if (length < 76) then
       a = (80 - length) / 2
       b = 80 - length - a - 2
       call print(repeat('=',a)//' '//trim(my_title)//' '//repeat('=',b), verbosity)
    else
       call print(title, verbosity)
    end if
    call print('',verbosity)

  end subroutine print_title

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X  Formatted reading
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  !%Read a line of text from a file (up to a line break, or 1024 characters).
  !%This can then be parsed by the calling routine (using 'parse_line' for example)
  !%
  !%Optionally, a status is returned which is:
  !%
  !%\begin{itemize}
  !% \item $<0$ if the end of the file is reached
  !% \item $=0$ if no problems were encountered
  !% \item $>0$ if there was a read error
  !%\end{itemize}
  !%
  !%The actual number returned is implementation specific
  function inoutput_read_line(this,status)

    type(Inoutput), intent(in)     :: this
    integer, optional, intent(out) :: status
    character(SYSTEM_STRING_LENGTH_LONG) :: inoutput_read_line
    integer                        :: my_status

    if (this%action == OUTPUT) call system_abort('read_line: Cannot read from an output file ('//trim(adjustl(this%filename))//')')

    if (present(status)) then
       read (this%unit,fmt='(a)',iostat=status) inoutput_read_line
    else
       read (this%unit,fmt='(a)',iostat=my_status) inoutput_read_line
       if (my_status < 0) call system_abort('read_line: End of file when reading '//trim(adjustl(this%filename)))
       if (my_status > 0) call system_abort('read_line: Error reading file '//trim(adjustl(this%filename)))
    end if

  end function inoutput_read_line

  subroutine inoutput_read_file(this,line_array, n_lines, status)
    type(Inoutput), intent(in) :: this
    character(len=*), allocatable, intent(inout) :: line_array(:)
    integer, intent(out) :: n_lines
    integer, optional, intent(out) :: status

    integer :: my_status, line_no, line_array_size

    line_no = 0
    my_status = 0

    do while (my_status == 0) 
      line_no = line_no + 1
      if (allocated(line_array)) then
	line_array_size = size(line_array)
      else
	line_array_size = 0
      endif
      if (line_no > line_array_size) call extend_char_array(line_array, 1.5_dp, 10)
      line_array(line_no) = read_line(this, my_status)
    end do

    n_lines = line_no - 1

    if (my_status > 0) then
      if (present(status)) then
	status = my_status
	return
      else
	call system_abort("ERROR inoutput_read_file reading file '"//trim(adjustl(this%filename))//"' status " // my_status)
      endif
    endif

  end subroutine inoutput_read_file

  subroutine extend_char_array(line_array, factor, minlen)
    character(len=*), allocatable, intent(inout) :: line_array(:)
    real(dp), intent(in), optional :: factor
    integer, intent(in), optional :: minlen

    real(dp) :: my_factor
    integer :: my_minlen
    character, allocatable :: t_line_array(:,:)
    integer old_n, i, j

    my_minlen = optional_default(10, minlen)
    my_factor = optional_default(1.5_dp, factor)

    if (.not. allocated(line_array)) then
      allocate(line_array(my_minlen))
      return
    endif

    if (allocated(line_array)) then
      old_n = size(line_array)
    else
      old_n = 0
    end if
    allocate(t_line_array(len(line_array(1)), old_n))
    do i=1, old_n
      do j=1, len(line_array(1))
	t_line_array(j,i) = line_array(i)(j:j)
      end do
    end do
    deallocate(line_array)
    allocate(line_array(int(old_n*my_factor)))
    do i=1, old_n
      do j=1, len(line_array(1))
	line_array(i)(j:j) = t_line_array(j,i)
      end do
    end do
    deallocate(t_line_array)
  end subroutine extend_char_array

  !% Call parse_string on the next line from a file
  subroutine inoutput_parse_line(this,delimiters,fields,num_fields,status)
    type(inoutput),             intent(in)    :: this
    character(*),               intent(in)    :: delimiters
    character(*), dimension(:), intent(inout) :: fields
    integer,                    intent(out)   :: num_fields
    integer, optional,          intent(out)   :: status
    integer                                   :: my_status

    num_fields = 0 ! in case read gives non-zero status
    local_line = read_line(this,my_status)
    if (present(status)) status = my_status
    if (my_status == 0) call split_string_simple(local_line, fields, num_fields, delimiters)

  end subroutine inoutput_parse_line

  !% split a string into fields separated by possible separators
  !% no quoting, matching separators, just a simple split
  subroutine split_string_simple(str, fields, n_fields, separators, error)
    character(len=*), intent(in) :: str !% string to be split
    character(len=*), intent(out) :: fields(:) !% on return, array of fields
    integer, intent(out) :: n_fields !% on return, number of fields
    character(len=*), intent(in) :: separators !% string of possible separators
    integer, intent(out), optional :: error

    integer :: str_len, cur_pos, next_pos, cur_field

    INIT_ERROR(error)

    str_len = len_trim(str)
    cur_pos = 0
    cur_field = 1
    do while (cur_pos <= str_len)
      if (cur_field > size(fields)) then
	RAISE_ERROR("split_string_simple str='"//trim(str)//"' no room for fields size(fields)="//size(fields)//" cur_field "//cur_field, error)
      endif
      next_pos = scan(str(cur_pos+1:str_len),separators)
      if (next_pos > 0) then
	if (next_pos == 1) then ! found another separator, skip it
	   cur_pos = cur_pos + 1
	else ! some stuff between us and separator, must be another field
	   fields(cur_field) = str(cur_pos+1:cur_pos+1+next_pos-2)
	   cur_pos = cur_pos + next_pos
	   cur_field = cur_field + 1
	endif
      else ! end of string, last field
	fields(cur_field) = str(cur_pos+1:str_len)
	cur_field = cur_field + 1
	cur_pos = str_len+1 ! exit loop
      endif
    end do
    n_fields = cur_field - 1

  end subroutine split_string_simple


  !% split a string at separators, making sure not to break up bits that
  !% are in quotes (possibly matching opening and closing quotes), and
  !% also strip one level of quotes off, sort of like a shell would when
  !% tokenizing
  subroutine split_string(this, separators, quotes, fields, num_fields, matching)
    character(len=*), intent(in) :: this
    character(len=*), intent(in) :: separators, quotes
    character(len=*), intent(inout) :: fields(:)
    integer, intent(out) :: num_fields
    logical, intent(in), optional :: matching

    integer :: i, length
    integer :: n_quotes
    character(len=len(quotes)) :: opening_quotes, closing_quotes
    integer :: opening_quote_index, closing_quote_pos
    logical :: do_matching
    character(len=len(fields(1))) :: tmp_field
    character(1) :: c
    integer :: tmp_field_last, t_start, dist
    logical :: in_token

    do_matching = optional_default(.false., matching)

    if (do_matching) then
      if (mod(len(quotes),2) == 0) then
	do i=1, len(quotes)/2
	  opening_quotes(i:i) = quotes(2*(i-1)+1:2*(i-1)+1)
	  closing_quotes(i:i) = quotes(2*(i-1)+2:2*(i-1)+2)
	end do
	n_quotes = len(quotes)/2
      else
	call system_abort("split_string called with matching=.true. but odd number of quotes " // (len(quotes)))
      endif
    else
      n_quotes = len(quotes)
      opening_quotes(1:n_quotes) = quotes(1:n_quotes)
      closing_quotes(1:n_quotes) = quotes(1:n_quotes)
    endif

    length = len_trim(this)
    num_fields = 0
    tmp_field = ""
    tmp_field_last = 0
    in_token = .false.
    i = 1
    do
      if (i > length) then ! last character
	if (in_token) then
	  num_fields = num_fields + 1
	  if (num_fields > size(fields)) call system_abort("split_string on '"//trim(this)//"' ran out of space for fields max " // size(fields))
	  if (tmp_field_last > 0) then
	    if (t_start <= length) then
	      fields(num_fields) = tmp_field(1:tmp_field_last) // this(t_start:length)
	    else
	      fields(num_fields) = tmp_field(1:tmp_field_last)
	    endif
	  else
	    if (t_start <= length) then
	      fields(num_fields) = this(t_start:length)
	    else
	    endif
	  endif
	endif
	exit
      else if (scan(this(i:i),opening_quotes(1:n_quotes)) > 0) then ! found an opening quote
	opening_quote_index = index(opening_quotes,this(i:i))
	closing_quote_pos = find_closing_delimiter(this(i+1:length), closing_quotes(opening_quote_index:opening_quote_index), &
						   opening_quotes(1:n_quotes), closing_quotes(1:n_quotes), do_matching)
	if (closing_quote_pos <= 0) then
	  call print("splitting string '"//trim(this)//"'", PRINT_ALWAYS)
	  call system_abort("split_string on '"//trim(this)//"' couldn't find closing quote matching opening at char " // i)
	endif
	if (in_token) then ! add string from t_start to tmp_field
          if (tmp_field_last > 0) then
            tmp_field = tmp_field(1:tmp_field_last) // this(t_start:i-1)
          else
            tmp_field = this(t_start:i-1)
          endif
	  tmp_field_last = tmp_field_last + (i-1 - t_start + 1)
	endif
	if (tmp_field_last > 0) then ! add contents of quote to tmp_field
	  if (i+closing_quote_pos-1 >= i+1) tmp_field = tmp_field(1:tmp_field_last) // this(i+1:i+closing_quote_pos-1)
	else
	  if (i+closing_quote_pos-1 >= i+1) tmp_field = this(i+1:i+closing_quote_pos-1)
	endif
        ! update tmp_field_last
	tmp_field_last = tmp_field_last + closing_quote_pos - 1
	in_token = .true.
	i = i + closing_quote_pos + 1
        ! reset t_start
	t_start = i
      else if (scan(this(i:i),separators) > 0) then ! found a separator
        if (next_non_separator(this, i+1, length, separators, dist) == '=') then
          if (in_token) then
            ! add string from t_start to tmp_field
            if (tmp_field_last > 0) then
              tmp_field = tmp_field(1:tmp_field_last) // this(t_start:i-1)//'='
            else
              tmp_field = this(t_start:i-1)//'='
            endif
            tmp_field_last = tmp_field_last + (i-1)-t_start+1+1
          else
            tmp_field = '='
            tmp_field_last = 1
          endif
          in_token = .true.
          ! update i and t_start to be after '='
          i = i + dist + 1
          t_start = i 
          if (i <= length) then
            ! look for next non separator
            c = next_non_separator(this, i, length, separators, dist)
            if (dist > 0) then
              i = i + dist-1
              t_start = i
            endif
          endif
        else
          if (in_token) then ! we were in a token before finding this separator
            num_fields = num_fields + 1
            if (num_fields > size(fields)) call system_abort("split_string on '"//trim(this)//"' ran out of space for fields max " // size(fields))
            ! add string from t_start and tmp_field to fields(num_fields)
            if (tmp_field_last > 0) then
              fields(num_fields) = tmp_field(1:tmp_field_last) // this(t_start:i-1)
            else
              fields(num_fields) = this(t_start:i-1)
            endif
            tmp_field = ""
            tmp_field_last = 0
          endif
          in_token = .false.
          i = i + 1
        endif
      else ! plain character
	if (.not. in_token) then
	  t_start = i
	endif
	in_token = .true.
	i = i + 1
      endif
    end do
  end subroutine split_string

  function next_non_separator(this, start, end, separators, dist) result(c)
    character(len=*), intent(in) :: this
    integer, intent(in) :: start, end
    character(len=*) :: separators
    integer, intent(out) :: dist
    character(1) :: c

    integer i

! call print("finding next_non_sep in '"//this(start:end)//"'")
    c=''
    dist = 0
    do i=start, end
      if (scan(this(i:i), separators) == 0) then ! found a non-separator
        c = this(i:i)
        dist = i-start+1
        exit
      endif
    end do
! call print("returning c='"//c//"' and dist "//dist)
  end function next_non_separator

  !% outdated - please use split_string
  !% Parse a string into fields delimited by certain characters. On exit
  !% the 'fields' array will contain one field per entry and 'num_fields'
  !% gives the total number of fields. 'status' will be given the error status
  !% (if present) and so can be used to tell if an end-of-file occurred.
  subroutine parse_string(this, delimiters, fields, num_fields, matching, error)

    character(*),               intent(in)    :: this
    character(*),               intent(in)    :: delimiters
    character(*), dimension(:), intent(inout) :: fields
    integer,                    intent(out)   :: num_fields
    logical, optional,          intent(in)    :: matching
    integer, optional,          intent(out)   :: error

    integer                                   :: field_start, length
    integer :: delim_pos
    integer :: n_delims, i
    character(len=len(delimiters)) :: opening_delims, closing_delims
    character(len=1) :: opening_delim
    integer :: opening_delim_index
    logical :: do_matching


    do_matching = optional_default(.false., matching)

    INIT_ERROR(error)

    field_start = 1
    num_fields = 0
    length = len_trim(this)

    if (do_matching) then
      if (mod(len(delimiters),2) == 0) then
	do i=1, len(delimiters)/2
	  opening_delims(i:i) = delimiters(2*(i-1)+1:2*(i-1)+1)
	  closing_delims(i:i) = delimiters(2*(i-1)+2:2*(i-1)+2)
	end do
	n_delims = len(delimiters)/2
      else
         RAISE_ERROR("parse_string called with matching=.true. but odd number of delimiters " // (len(delimiters)), error)
      endif
    else
      n_delims = len(delimiters)
      opening_delims(1:n_delims) = delimiters(1:n_delims)
      closing_delims(1:n_delims) = delimiters(1:n_delims)
    endif

    do
      delim_pos = scan(this(field_start:length), opening_delims(1:n_delims))
      if (delim_pos == 0) then ! didn't find opening delimiter
	if (len_trim(this(field_start:length)) == 0) then !...and the rest of the string is blank...
	  !... then we've finished
	  exit
	else !otherwise, there's one field left to get
	  if (length >= field_start) then
	    num_fields = num_fields + 1
	    if (num_fields > size(fields)) then
              RAISE_ERROR("parse_string ran out of space for fields", error)
            endif
            fields(num_fields) = this(field_start:length)
	  endif
	  return
	end if
      endif
      ! get here if we found an opening delimiter
      delim_pos = delim_pos + field_start - 1
      if (delim_pos /= field_start) then ! found an opening delimiter after some text
	! save text in a field, and jump over it
	if (delim_pos-1 >= field_start) then
	  num_fields = num_fields + 1
          if (num_fields > size(fields)) then
             RAISE_ERROR("parse_string ran out of space for fields", error)
          endif
	  fields(num_fields) = this(field_start:delim_pos-1)
	end if
	field_start = delim_pos
      endif
      field_start = field_start + 1
      if (do_matching) then
	opening_delim = this(delim_pos:delim_pos)
	opening_delim_index = index(opening_delims(1:n_delims), opening_delim)
	delim_pos = find_closing_delimiter(this(field_start:length), closing_delims(opening_delim_index:opening_delim_index), opening_delims(1:n_delims), closing_delims(1:n_delims), do_matching)
      else
	delim_pos = find_closing_delimiter(this(field_start:length), closing_delims, opening_delims(1:n_delims), closing_delims(1:n_delims), do_matching)
      endif
      if (delim_pos == 0) then ! didn't find closing delimiter
	if (do_matching) then
	  call print("parse_string failed to find closing delimiter to match opening delimiter at position " // (field_start-1), PRINT_ALWAYS)
	  call print("parse_string string='"//this//"'", PRINT_ALWAYS)
          RAISE_ERROR("parse_string failed to find closing delimiter", error)
	else
	  delim_pos = length-field_start+2
	endif
      endif
      delim_pos = delim_pos + field_start - 1
      if (delim_pos-1 >= field_start) then
	num_fields = num_fields + 1
        if (num_fields > size(fields)) then
          RAISE_ERROR("parse_string ran out of space for fields", error)
        endif
	fields(num_fields) = this(field_start:delim_pos-1)
      endif
      field_start = delim_pos+1
      if (field_start > length) return
    end do

  end subroutine parse_string

  recursive function find_closing_delimiter(this, closing_delim, opening_delims, closing_delims, matching) result(pos)
    character(len=*), intent(in) :: this
    character(len=*), intent(in) :: closing_delim
    character(len=*), intent(in) :: opening_delims, closing_delims
    logical :: matching
    integer :: pos

    integer :: length, first_matching_closing_delim, first_opening_delim
    integer :: opening_delim_index
    character(len=1) :: opening_delim
    integer :: substring_end_pos

    pos = 0

    do
      if (matching) then
	first_matching_closing_delim = scan(this, closing_delim)
      else
	first_matching_closing_delim = scan(this, closing_delims)
      endif
      first_opening_delim = scan(this, opening_delims)
      if ((first_opening_delim > 0) .and. (first_opening_delim < first_matching_closing_delim)) then
	length = len(this)
	if (matching) then
	  opening_delim = this(first_opening_delim:first_opening_delim)
	  opening_delim_index = index(opening_delims, opening_delim)
	  substring_end_pos = find_closing_delimiter(this(first_opening_delim+1:length), &
	    closing_delims(opening_delim_index:opening_delim_index), opening_delims, closing_delims, matching)
	else
	  substring_end_pos = find_closing_delimiter(this(first_opening_delim+1:length), &
	    closing_delims, opening_delims, closing_delims, matching)
	endif
	if (substring_end_pos == 0) &
	  call system_abort("find_closing_delimiter failed to find substring closing delimiter '"// &
	    closing_delims(opening_delim_index:opening_delim_index)//"' in string '"//this// &
	    "' for substring starting at "//(first_opening_delim+1))
	substring_end_pos = substring_end_pos + first_opening_delim+1 - 1
	pos = find_closing_delimiter(this(substring_end_pos+1:length), closing_delim, opening_delims, &
	  closing_delims, matching) + substring_end_pos+1 - 1
	return
      else
	pos=first_matching_closing_delim
	return
      endif
    end do

    return

  end function find_closing_delimiter

  !% Parse a string into fields delimited by certain characters. On exit
  !% the 'fields' array will contain one field per entry and 'num_fields'
  !% gives the total number of fields. 'status' will be given the error status
  !% (if present) and so can be used to tell if an end-of-file occurred.
  subroutine parse_string_orig(this, delimiters, fields, num_fields)

    character(*),               intent(in)    :: this
    character(*),               intent(in)    :: delimiters
    character(*), dimension(:), intent(inout) :: fields
    integer,                    intent(out)   :: num_fields

    integer                                   :: field_start,field_end,width,length
    integer                                   :: array_length

    field_start = 1
    num_fields = 0
    length = len_trim(this)
    array_length = size(fields)

    do
       !Try to find a delimiter
       width = scan(this(field_start:length),delimiters)
       !If delimiter not found...
       if (width == 0) then
          !...and the rest of the string is blank...
          if (len_trim(this(field_start:length)) == 0) then
             !... then we've finished
             exit
           else
              !otherwise, there's one field left to get
              field_end = length
           end if
        !On the other hand, if the delimiter is the first character...
        else if (width == 1) then
           !...then move past it and start the do loop again
           field_start = field_start + 1
           cycle
        !Otherwise calculate the end of the field, without the delimiter
        else
           field_end = field_start + width - 2
        end if

        num_fields = num_fields + 1

        if (num_fields > array_length) then
           call print(this)
           call system_abort('inoutput_parse_line: Number of fields is greater than storage array size')
        end if

        fields(num_fields) = adjustl(this(field_start:field_end))
        field_start = field_end + 1

     end do

  end subroutine parse_string_orig

  !% Convert an input string into an integer. If 'err' is present, it is set to true
  !% if an error occurred during the conversion.
  function string_to_int(string,err)
    character(*), intent(in)   :: string
    character(len=len(string)) :: local_string
    logical, optional, intent(out) :: err
    integer                    :: String_To_Int
    character(10)              :: format
    integer                    :: n
    integer stat

    local_string = adjustl(string)
    n = len_trim(local_string)
    write(format,'(a,i0,a)')'(i',n,')'
    read(local_string,format,iostat=stat) string_to_int
    if (present(err)) err = (stat /= 0)

  end function string_to_int

  !% Convert an input string into a logical. If 'err' is present, it is set to true
  !% if an error occurred during the conversion.
  function string_to_logical(string, err)
    character(*), intent(in)   :: string
    logical, optional, intent(out) :: err
    logical                    :: string_to_logical
    integer stat

    read(string,*,iostat=stat) string_to_logical

    if (present(err)) err = (stat /= 0)

  end function string_to_logical


  !% Convert an input string into a real. If 'err' is present, it is set to true
  !% if an error occurred during the conversion.
  function string_to_real(string, err)
    character(*), intent(in)   :: string
    logical, optional, intent(out) :: err
    real(dp)                   :: string_to_real
    integer stat

    if (present(err)) then
       err = .false.
       if (scan(adjustl(string), 'tfTF') == 1) then
          err = .true.
          return
       end if
    end if

    read(string,*,iostat=stat) string_to_real

    if (present(err)) err = (stat /= 0)


  end function string_to_real

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!% Reallocation: used to reduce the need to deallocate and reallocate
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine reallocate_int1d(array, d1, zero, copy)
    integer, allocatable, dimension(:), intent(inout) :: array
    integer,                            intent(in)    :: d1
    logical, optional,                  intent(in)    :: zero, copy

    logical :: do_copy
    integer, allocatable :: tmp(:)

    if (allocated(array)) then
       if (size(array) /= d1) then
	  do_copy = optional_default(.false., copy)
	  if (do_copy) then
	    allocate(tmp(size(array)))
	    tmp = array
	  endif
          deallocate(array)
          allocate(array(d1))
	  if (do_copy) then
	    array = 0
	    array(1:min(size(tmp),size(array))) = tmp(1:min(size(tmp),size(array)))
	    deallocate(tmp)
	  endif
       end if
    else
       allocate(array(d1))
    end if

    if (present(zero)) then
       if (zero) array = 0
    end if

  end subroutine reallocate_int1d

  subroutine reallocate_real1d(array, d1, zero)
    real(dp), allocatable, dimension(:), intent(inout) :: array
    integer,                             intent(in)    :: d1
    logical, optional,                   intent(in)    :: zero

    if (allocated(array)) then
       if (size(array) /= d1) then
          deallocate(array)
          allocate(array(d1))
       end if
    else
       allocate(array(d1))
    end if

    if (present(zero)) then
       if (zero) array = 0.0_dp
    end if

  end subroutine reallocate_real1d

  subroutine reallocate_int2d(array, d1, d2, zero, copy)
    integer, allocatable, dimension(:,:), intent(inout) :: array
    integer,                              intent(in)    :: d1,d2
    logical, optional,                    intent(in)    :: zero, copy

    logical :: do_copy
    integer, allocatable :: tmp(:,:)

    if (allocated(array)) then
       if (.not. all(shape(array) == (/d1,d2/))) then
	  do_copy = optional_default(.false., copy)
	  if (do_copy) then
	    allocate(tmp(size(array,1),size(array,2)))
	    tmp = array
	  endif
          deallocate(array)
          allocate(array(d1,d2))
	  if (do_copy) then
	    array = 0
	    array(1:min(size(tmp,1),size(array,1)),1:min(size(tmp,2),size(array,2))) = tmp(1:min(size(tmp,1),size(array,1)),1:min(size(tmp,2),size(array,2)))
	    deallocate(tmp)
	  endif
       end if
    else
       allocate(array(d1,d2))
    end if

    if (present(zero)) then
       if (zero) array = 0
    end if

  end subroutine reallocate_int2d

  subroutine reallocate_real2d(array, d1, d2, zero, copy)
    real(dp), allocatable, dimension(:,:), intent(inout) :: array
    integer,                               intent(in)    :: d1,d2
    logical, optional,                     intent(in)    :: zero, copy

    logical :: do_copy
    real(dp), allocatable :: tmp(:,:)

    if (allocated(array)) then
       if (.not. all(shape(array) == (/d1,d2/))) then
	  do_copy = optional_default(.false., copy)
	  if (do_copy) then
	    allocate(tmp(size(array,1),size(array,2)))
	    tmp = array
	  endif
          deallocate(array)
          allocate(array(d1,d2))
	  if (do_copy) then
	    array = 0.0_dp
	    array(1:min(size(tmp,1),size(array,1)),1:min(size(tmp,2),size(array,2))) = tmp(1:min(size(tmp,1),size(array,1)),1:min(size(tmp,2),size(array,2)))
	    deallocate(tmp)
	  endif
       end if
    else
       allocate(array(d1,d2))
    end if

    if (present(zero)) then
       if (zero) array = 0.0_dp
    end if

  end subroutine reallocate_real2d

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!% Read scalar and array data from ascii files. These
!% interfaces are not yet heavily overloaded to cater for all intrinsic and most
!% derived types.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine reada_real_dim1(this,da,status)
    type(Inoutput), intent(in)     :: this
    real(dp), intent(out) :: da(:)
    integer, optional, intent(out) :: status

    integer :: my_status

    if (this%action == OUTPUT) call system_abort('read_line: Cannot read from an output file ('//trim(adjustl(this%filename))//')')

    if (present(status)) then
       read (this%unit,fmt=*,iostat=status) da
    else
       read (this%unit,fmt=*,iostat=my_status) da
       if (my_status < 0) call system_abort('read_line: End of file when reading '//trim(adjustl(this%filename)))
       if (my_status > 0) call system_abort('read_line: Error reading file '//trim(adjustl(this%filename)))
    end if
  end subroutine reada_real_dim1

  subroutine reada_int_dim1(this,ia,status)
    type(Inoutput), intent(in)     :: this
    integer, intent(out) :: ia(:)
    integer, optional, intent(out) :: status

    integer :: my_status

    if (this%action == OUTPUT) call system_abort('read_line: Cannot read from an output file ('//trim(adjustl(this%filename))//')')

    if (present(status)) then
       read (this%unit,fmt=*,iostat=status) ia
    else
       read (this%unit,fmt=*,iostat=my_status) ia
       if (my_status < 0) call system_abort('read_line: End of file when reading '//trim(adjustl(this%filename)))
       if (my_status > 0) call system_abort('read_line: Error reading file '//trim(adjustl(this%filename)))
    end if
  end subroutine reada_int_dim1

   !% Rewind to the start of this file. Works for both formatted and unformatted files.
   subroutine rewind(this)
      type(Inoutput), intent(inout) :: this
      rewind this%unit
   end subroutine rewind

   !% Move the file pointer back by 'n' (defaults to 1) records. Works for
   !% formatted and unformatted files.
   subroutine backspace(this,n)
      type(Inoutput), intent(inout) :: this
      integer, optional              :: n
      integer                        :: i
      if (present(n)) then
	  do i=1,n
	     backspace this%unit
	  end do
      else
	  backspace this%unit
      end if
   end subroutine backspace

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!% Concatenation functions.
!% Overloadings for the // operator to make strings from various other types.
!% In each case, we need to work out the exact length of the resultant string
!% in order to avoid printing excess spaces.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  !% Return a string which is the real number 'r' rounded to 'digits' decimal digits
  function round(r,digits)

    real(dp), intent(in) :: r
    integer,  intent(in) :: digits
    ! below we work out the exact length of the resultant string
    !          space for '-' sign or not  + digits in integer part +                              space for . or not             + decimal digits
    character( int(0.5_dp-sign(0.5_dp,r)) + int(log10(max(1.0_dp,abs(r)+0.5_dp*10.0_dp**(-digits)))) + 1 + int(sign(0.5_dp,real(digits,dp)-0.5_dp)+0.5_dp) + max(0,digits)) :: round
    character(8) :: format

    if (digits > 0) then
       write(format,'(a,i0,a)')'(f0.',max(0,digits),')'
       write(round,format) r
    else
       write(round,'(i0)') int(r)
    end if

  end function round

  function string_cat_logical(string, log)
    character(*),      intent(in)  :: string
    logical,           intent(in)  :: log
    character((len(string)+1)) :: string_cat_logical
    write(string_cat_logical,'(a,l1)') string, log
  end function string_cat_logical

  function string_cat_logical_array(string, log)
    character(*),      intent(in)  :: string
    logical,           intent(in)  :: log(:)
    character((len(string)+2*size(log)-1)) :: string_cat_logical_array
    character(len=32) format

    format = '(a,'//size(log)//'(l1,1x),l1)'
    write(string_cat_logical_array,format) string, log
  end function string_cat_logical_array

  elemental function int_format_length(i) result(len)
    integer, intent(in)::i
    integer::len
    len = max(1,(-sign(1, i)+1)/2 + ceiling(log10(abs(real(i,dp))+0.01_dp)))
  end function int_format_length

  function string_cat_int(string, int)
    character(*),      intent(in)  :: string
    integer,           intent(in)  :: int
    ! below we work out the exact length of the resultant string
    character(len(string)+int_format_length(int)) :: string_cat_int

    write(string_cat_int,'(a,i0)') string, int
  end function string_cat_int

  function int_cat_string(int,string)
    character(*),      intent(in)  :: string
    integer,           intent(in)  :: int
    ! below we work out the exact length of the resultant string
    character(len(string)+int_format_length(int)) :: int_cat_string

    write(int_cat_string,'(i0,a)') int,string 
  end function int_cat_string

  function string_cat_int_array(string, values)
    character(*),      intent(in)  :: string
    integer,           intent(in)  :: values(:)
    ! below we work out the exact length of the resultant string
    character(len(string)+size(values)+sum(int_format_length(values)))::string_cat_int_array

      character(32) :: format

      if (size(values) == 1) then
         format = '(a,i0)'
         write(string_cat_int_array,format) string, values
      else if (size(values)>1) then
         format = '(a,' // (size(values)-1) //'(i0,1x),i0)'
         write(string_cat_int_array,format) string, values
      else
         write(string_cat_int_array,'(a)') string
      end if

  end function string_cat_int_array

  pure function real_sci_format_length() result(len)
    integer::len
    !  space sign 0.   fractional part                    E+00
    len = 1 + 1 + 2 + max(0,mainlog%default_real_precision)+4
  end function real_sci_format_length


  function string_cat_real_array(string, values)
    character(*),      intent(in)  :: string
    real(dp),          intent(in)  :: values(:)
    ! we work out the exact length of the resultant string
    character((len(string)+size(values)*real_sci_format_length())) :: string_cat_real_array
    character(32) :: format

    if (size(values)>0) then
       ! replaced concatenation with write... for PGI bug, NB 22/6/2007
       write(format,'("(a,",I0,"e",I0,".",I0,")")') size(values), real_sci_format_length(), &
            mainlog%default_real_precision
       write(string_cat_real_array, format) string, values 
    else
       write(string_cat_real_array, '(a)') string
    end if

  end function string_cat_real_array

  function string_cat_complex_array(string, values)
    character(*),      intent(in)  :: string
    complex(dp),          intent(in)  :: values(:)
    ! we work out the exact length of the resultant string
    character((len(string)+2*size(values)*real_sci_format_length())) :: string_cat_complex_array
    character(32) :: format

    if (size(values)>0) then
       ! replaced concatenation with write... for PGI bug, NB 22/6/2007
       write(format,'("(a,",I0,"e",I0,".",I0,")")') 2*size(values), real_sci_format_length(), &
            mainlog%default_real_precision
       write(string_cat_complex_array, format) string, values 
    else
       write(string_cat_complex_array, '(a)') string
    end if

  end function string_cat_complex_array

  function string_cat_string_array(string, values)
    character(*),      intent(in)  :: string
    character(*),      intent(in)  :: values(:)
    ! we work out the exact length of the resultant string
    character(len(string)+size(values)*len(values(1))) :: string_cat_string_array
    character(32) :: format

    if (size(values)>0) then
       ! replaced concatenation with write... for PGI bug, NB 22/6/2007
       write(format,'("(a",I0,",",I0,"a",I0,")")') len(string), size(values)+1, len(values(1))
       write(string_cat_string_array, format) string, values 
    else
       write(string_cat_string_array, '(a)') string
    end if

  end function string_cat_string_array

  function real_array_cat_string(values, string)
    character(*),      intent(in)  :: string
    real(dp),          intent(in)  :: values(:)
    ! we work out the exact length of the resultant string
    character((len(string)+size(values)*real_sci_format_length())) :: real_array_cat_string
    character(32) :: format

    if (size(values)>0) then
       ! replaced concatenation with write... for PGI bug, NB 22/6/2007
       write(format,'("(",I0,"e",I0,".",I0,",a)")') size(values), real_sci_format_length(), &
            mainlog%default_real_precision
       write(real_array_cat_string, format) values, string
    else
       write(real_array_cat_string, '(a)') string
    end if

  end function real_array_cat_string

  pure function real_format_length(r) result(len)
    real(dp), intent(in)::r
    integer::len

    if(isnan(r)) then
       len = 3
    else       !         sign                           int part         space?          decimal point                                                        fractional part
       len = int(0.5_dp-sign(0.5_dp,r)) + int(log10(max(1.0_dp,abs(r)))) + 1 + & 
           & int(sign(0.5_dp,real(mainlog%default_real_precision,dp)-0.5_dp)+0.5_dp) &
           & + max(0,mainlog%default_real_precision)

#ifdef GFORTRAN_ZERO_HACK
       !gfortran hack - 0.0000... is printed as .00000000
       if (r == 0.0) len = len - 1
#endif

    end if
  end function real_format_length

  pure function complex_format_length(c) result(len)
    complex(dp), intent(in)::c
    integer::len

    len = real_format_length(real(c))+1+real_format_length(imag(c))
  end function complex_format_length

  function real_cat_string(r, string)
    character(*),      intent(in)  :: string
    real(dp),          intent(in)  :: r
    ! we work out the exact length of the resultant string
    character( len(string)+real_format_length(r)) :: real_cat_string
    character(12) :: format

    if (mainlog%default_real_precision > 0) then
       write(format,'(a,i0,a)')'(f0.',max(0,mainlog%default_real_precision),',a)'
       if (isnan(r)) then
          write(real_cat_string,'(a,a)') "NaN", string
       else
          write(real_cat_string,format) r, string
       endif
    else
       write(real_cat_string,'(i0,a)') int(r), string
    end if
  end function real_cat_string

  function string_cat_real(string, r)
    character(*),      intent(in)  :: string
    real(dp),          intent(in)  :: r
    ! we work out the exact length of the resultant string
    character( len(string)+real_format_length(r)) :: string_cat_real
    character(12) :: format

    if (mainlog%default_real_precision > 0) then
       if (isnan(r)) then
	 write(string_cat_real,'(a,a)') string,"NaN"
       else
	 write(format,'(a,i0,a)')'(a,f0.',max(0,mainlog%default_real_precision),')'
	 write(string_cat_real,format) string, r
       endif
    else
       write(string_cat_real,'(a,i0)') string, int(r)
    end if
  end function string_cat_real

  function string_cat_complex(string, c)
    character(*),      intent(in)  :: string
    complex(dp),          intent(in)  :: c
    ! we work out the exact length of the resultant string
    character( len(string)+complex_format_length(c)) :: string_cat_complex
    character(24) :: format

    if (mainlog%default_real_precision > 0) then
       write(format,'(a,i0,a,i0,a)')'(a,f0.',max(0,mainlog%default_real_precision),'," ",f0.', &
 	                                     max(0,mainlog%default_real_precision),')'
       write(string_cat_complex,format) string, c
    else
       write(string_cat_complex,'(i0," ",i0)') string, int(real(c)), int(imag(c))
    end if
  end function string_cat_complex


  !% Return the mpi size and rank for the communicator 'comm'.
  !% this routine aborts of _MPI is not defined
  subroutine get_mpi_size_rank(comm, nproc, rank)
    
    integer, intent(in)  :: comm  !% MPI communicator
    integer, intent(out) :: nproc  !% Total number of processes
    integer, intent(out) :: rank  !% Rank of this process

#ifdef _MPI
    include 'mpif.h'
#endif

#ifdef _MPI
    
    integer::error_code

    call MPI_COMM_SIZE(comm, nproc, error_code)
    if (error_code .ne. MPI_SUCCESS) then
       rank=-1
       nproc=-1
       return
    endif
    call MPI_COMM_RANK(comm, rank, error_code)
    if (error_code .ne. MPI_SUCCESS) then
       rank=-1
       nproc=-1
       return
    endif
#else
    rank = 0
    nproc = 1
#endif
  end subroutine get_mpi_size_rank

  ! System initialiser

  !% Must be called at the start of all programs. Initialises MPI if present,
  !% set the random number seed sets up the default Inoutput objects
  !% logger and errorlog to point to stdout and stderr respectively. Calls
  !% Hello_World to do some of the work and print a friendly welcome. If we're
  !% using MPI, by default we set the same random seed for each process.
  !% This also attempts to read the executable name, the number of command
  !% arguments, and the arguments themselves.
  subroutine system_initialise(verbosity,seed, mpi_all_inoutput, common_seed, enable_timing, quippy_running, mainlog_file)
    integer,intent(in), optional::verbosity           !% mainlog output verbosity
    integer,intent(in), optional::seed                !% Seed for the random number generator.
    logical,intent(in), optional::mpi_all_inoutput    !% Print on all MPI nodes (false by default)
    logical,intent(in), optional::common_seed 
    logical,intent(in), optional::enable_timing           !% Enable system_timer() calls
    logical,intent(in), optional::quippy_running       !% .true. if running under quippy (Python interface)
    character(len=*),intent(in), optional :: mainlog_file
    !% If 'common_seed' is true (default), random seed will be the same for each
    !% MPI process.
    character(30) :: arg
    integer       :: status, i, n

#ifdef _MPI
    integer::PRINT_ALWAYS
    integer :: is_initialised
    include "mpif.h"

    call MPI_initialized(is_initialised, PRINT_ALWAYS)
    call abort_on_mpi_error(PRINT_ALWAYS, "system_initialise, mpi_initialised()")
    if (is_initialised == 0) then
      call MPI_INIT(PRINT_ALWAYS)
      call abort_on_mpi_error(PRINT_ALWAYS, "system_initialise, mpi_init()")
    endif
    call get_mpi_size_rank(MPI_COMM_WORLD, mpi_n, mpi_myid)
    if (mpi_n < 1 .or. mpi_myid < 0) &
      call system_abort("system_initialise Got bad size="//mpi_n// " or rank="//mpi_myid//" from get_mpi_size_rank")
    error_mpi_myid = mpi_myid
#else
    mpi_n  = 1 ! default
    mpi_myid = 0 ! default
#endif

    ! should really be in declaration, but ifort complains
#ifdef NO_F2003_NEW_LINE
    quip_new_line = char(13)
#else
    quip_new_line = new_line(' ')
#endif

    call cpu_time(start_time)
    mainlog%mpi_all_inoutput_flag = optional_default(.false., mpi_all_inoutput)

    ! Initialise the verbosity stack and default logger
    if (present(mainlog_file)) then
      call initialise(mainlog, trim(mainlog_file), verbosity=verbosity, action=OUTPUT)
    else
      call initialise(mainlog, verbosity=verbosity, action=OUTPUT)
    endif
    call print_mpi_id(mainlog)
    call initialise(errorlog,'stderr')
    call print_mpi_id(errorlog)
    error_unit = errorlog%unit

    call hello_world(seed, common_seed)

    RAN_MAX = huge(1)

    ! Query arguments and executable name
    num_command_args = cmd_arg_count()
    call get_cmd_arg(0, arg, status=status)
    if (status == 0) then
       EXEC_NAME = arg
    else
       EXEC_NAME ='<UNKNOWN>'
    end if
    if (NUM_COMMAND_ARGS < MAX_READABLE_ARGS) then
       n = NUM_COMMAND_ARGS
    else
       n = MAX_READABLE_ARGS
    end if
    do i = 1, n
       call get_cmd_arg(i, COMMAND_ARG(i), status = status)
       if (status /= 0) then
          write(line,'(a,i0)')'system_initialise: Problem reading command argument ',i
          call system_abort(line)
       end if
    end do

    system_do_timing = optional_default(.false., enable_timing)
    if (system_do_timing) then
      call print("Calls to system_timer will report times")
    else
      call print("Calls to system_timer will do nothing by default")
    endif
    call print('')

    system_quippy_running = optional_default(.false., quippy_running)

  end subroutine system_initialise

  subroutine get_mainlog_errorlog_ptr(mainlog_ptr, errorlog_ptr)
    type inoutput_ptr
       type(inoutput), pointer :: p
    end type inoutput_ptr
    integer, intent(out), dimension(12) :: mainlog_ptr, errorlog_ptr
    type(inoutput_ptr) :: mainlog_p, errorlog_p

    mainlog_p%p => mainlog
    mainlog_ptr = transfer(mainlog_p, mainlog_ptr)
    errorlog_p%p => errorlog
    errorlog_ptr = transfer(errorlog_p, errorlog_ptr)

  end subroutine get_mainlog_errorlog_ptr

  function cmd_arg_count()
    integer :: cmd_arg_count

#ifndef GETARG_F2003
    integer :: iargc
    external iargc
#endif

#ifndef GETARG_F2003
    cmd_arg_count = iargc()
#else
    cmd_arg_count = command_argument_count()
#endif
  end function cmd_arg_count

  subroutine get_cmd_arg(i,arg, status)
    integer, intent(in) :: i
    character(len=*), intent(out) :: arg
    integer, intent(out), optional :: status

#ifndef GETARG_F2003
    external getarg
#endif

#ifndef GETARG_F2003
    call getarg(i, arg)
    if (present(status)) status = 0
#else
    call get_command_argument(i, arg, status=status)
#endif
  end subroutine get_cmd_arg

  subroutine get_env_var(name, arg, status)
    character(len=*), intent(in) :: name
    character(len=*), intent(out) :: arg
    integer, intent(out), optional :: status

#ifndef GETENV_F2003
    external getenv
#endif

#ifndef GETENV_F2003
    call getenv(trim(name), arg)
    if (present(status)) then
       if( len_trim(arg) > 0 ) then
          status = 0
       else
          status = 1
       endif
    endif
#else
    call get_environment_variable(trim(name), arg, status=status)
#endif
  end subroutine get_env_var

  !% Shut down gracefully, finalising system objects.
  subroutine system_finalise()
    integer :: values(8)
#ifdef _MPI
    integer :: PRINT_ALWAYS
    include "mpif.h"
#endif

    call date_and_time(values=values)
    call print("")
    call print('System::Finalise: '//date_and_time_string(values))
    call print("System::Finalise: Bye-Bye!")
    call finalise(mainlog)
    call finalise(errorlog)
#ifdef _MPI
    call mpi_finalize(PRINT_ALWAYS)
    call abort_on_mpi_error(PRINT_ALWAYS, "system_finalise, mpi_finalise()")
#endif
  end subroutine system_finalise

  !% Print a warning message to default mainlog, but don't quit
  subroutine print_warning(message)
    character(*), intent(in) :: message
    call print('WARNING: '//message)
  end subroutine print_warning

  !% Take the values from 'date_and_time' and make a nice string
  function date_and_time_string(values)
    character(21)       :: date_and_time_string
    integer, intent(in) :: values(8)
    character(2)        :: time(7)
    character(4)        :: year
    integer             :: i

    write(year,'(i0)') values(1)
    do i = 2, 7
       if (i==4) cycle ! We don't use the local adjustment to UTC
       write(time(i),'(i0.2)') values(i)
    end do
    write(date_and_time_string,'(11a)') time(3),'/',time(2),'/',year,'   ',time(5),':',time(6),':',time(7)

  end function date_and_time_string

  subroutine system_set_random_seeds(seed)
#ifdef _OPENMP
    use omp_lib
#endif
    integer, intent(in) :: seed

    integer :: n
    integer, allocatable :: seed_a(:)

#ifdef _OPENMP
    !$omp parallel
    idum=seed*(omp_get_thread_num()+1)
    !$omp end parallel
#else
    idum=seed              !gabor generator
#endif

    call random_seed(size=n)
    allocate(seed_a(n))
    seed_a = seed
    call random_seed(put=seed_a) !fortran 90 generator
    deallocate(seed_a)

  end subroutine system_set_random_seeds

  !% Called by 'system_initialise' to print welcome messages and
  !% seed the random number generator. 
  subroutine hello_world(seed, common_seed)
    integer, optional::seed !% Seed for the random number generator.
    logical, optional :: common_seed
    !% If 'common_seed' is true (default), random seed will be the same for each
    !% MPI process.

    integer:: actual_seed, i, ran_dummy
    integer:: values(20) ! for time inquiry function
    logical :: use_common_seed
#ifdef _MPI
    integer :: PRINT_ALWAYS
    include "mpif.h"    
#endif

    call date_and_time(values=values)

    call print('System::Hello World: '//date_and_time_string(values))
    call print('System::Hello World: SVN version  '//current_version())
    call print('System::Hello World: QUIP_ARCH    '//QUIP_ARCH)
    call print('System::Hello World: compiled on  '//__DATE__//' at '//__TIME__)
#ifdef _MPI
    call print('System::Hello World: MPI parallelisation with '//mpi_n_procs()//' processes')
#endif
!   Open MP stuff
!$OMP parallel
!$OMP master
!$  call print('System::Hello World: OpenMP parallelisation with '//OMP_get_num_threads()//' threads')
!$OMP end master
!$OMP end parallel
    if(present(seed)) then
       actual_seed = seed
    else
       actual_seed=1 + values(8)+values(5)+values(6)+values(7) !hour+minute+seconds+millisecond
       use_common_seed = .true.
       if (present(common_seed)) use_common_seed = common_seed
#ifdef _MPI
       if (.not. use_common_seed) then
          call print('system::Hello World: MPI run with different seeds on each process')
       else
          call print('system::Hello World: MPI run with the same seed on each process')

          ! Broadcast seed from process 0 to all others
          call MPI_Bcast(actual_seed, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, PRINT_ALWAYS)
          call abort_on_mpi_error(PRINT_ALWAYS, 'Hello_World: MPI_Bcast()')
       end if
#endif
    end if

    call print('System::Hello World: Random Seed = '//actual_seed)
    call system_set_random_seeds(actual_seed)

    ! The first seed tends to give very small random numbers. The loop below decorrelates the seed from the initial value so it can be trusted to be uniform.
    do i = 1, 100
       ran_dummy = ran()
    enddo

    call print('System::Hello World: global verbosity = '//value(mainlog%verbosity_stack))
    call print('')
  end subroutine hello_world

  subroutine system_resync_rng
#ifdef _MPI
#ifdef _OPENMP
    use omp_lib
#endif
    integer :: PRINT_ALWAYS
#ifdef _OPENMP
    integer :: idum_buf(omp_get_num_threads())
#endif
    include "mpif.h"    
    call print('Resyncronising random number generator', PRINT_VERBOSE)

#ifdef _OPENMP
    ! Collect all seeds into continuous buffer
    !$omp parallel
    idum_buf(omp_get_thread_num()+1) = idum
    !$omp end parallel
    ! Broadcast seed from process 0 to all others
    call MPI_Bcast(idum_buf, omp_get_num_threads(), MPI_INTEGER, 0, &
         MPI_COMM_WORLD, PRINT_ALWAYS)
    call abort_on_mpi_error(PRINT_ALWAYS, 'system_resync_rng')
    ! Set seeds on individual threads
    !$omp parallel
    idum = idum_buf(omp_get_thread_num()+1)
    !$omp end parallel
#else
    ! Broadcast seed from process 0 to all others
    call MPI_Bcast(idum, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, PRINT_ALWAYS)
    call abort_on_mpi_error(PRINT_ALWAYS, 'system_resync_rng')
#endif
#endif
  end subroutine system_resync_rng

  !
  !% Return the correct ordinal ending (st,nd,rd,th) for the given integer
  !
  elemental function th(n)

    integer, intent(in) :: n
    character(2)        :: th
    integer             :: l,m

    l = mod(n,100)
    m = mod(n,10)

    if (l > 10 .and. l < 20) then
       th = 'th'
    else 
       select case(m)
       case(1)
          th = 'st'
       case(2)
          th = 'nd'
       case(3)
          th = 'rd'
       case default
          th = 'th'
       end select
    end if         

  end function th


  !% Reseed the random number generator. Useful when restarting from check files.
  subroutine system_reseed_rng(new_seed)
    integer, intent(in) :: new_seed

    call print('System::Reseed_RNG: Reseeding random number generator, new seed = '//new_seed,PRINT_VERBOSE)
    call system_set_random_seeds(new_seed)
  end subroutine system_reseed_rng

  !% Return the current random number seed. 
  function system_get_random_seed()
    integer :: system_get_random_seed
    
    system_get_random_seed = idum
  end function system_get_random_seed

  !% Return a random integer
  function ran()
    integer::ran
    integer:: k
    real(dp):: dran

    if (system_use_fortran_random) then
       call random_number(dran)
       ran = dran*RAN_MAX
    else
       if (idum == 0) &
	 call system_abort("function ran(): linear-congruential random number generators fail with seed idum=0")
       k = idum/ran_Q
       idum = ran_A*(idum-k*ran_Q)-ran_R*k
       if(idum < 0) idum = idum + ran_M
       ran =idum
     endif
  end function ran

  !% Return a random real number uniformly distributed in the range [0,1]
  function ran_uniform()
    real(dp)::ran_uniform

    if (system_use_fortran_random) then
       call random_number(ran_uniform)
    else
       ran_uniform = 1.1_dp
       do while (ran_uniform > 1.0_dp) ! generating [0,1]  
	  ran_uniform = real(ran(),dp)/RAN_MAX
       end do
    end if
  end function ran_uniform

  !% Return random real from Normal distribution with mean zero and standard deviation one.
  function ran_normal()
     real(dp) :: ran_normal, r, v1, v2     
     r = 2.0_dp
     do while (r > 1.0_dp)
        v1 = 2.0_dp * ran_uniform() - 1.0_dp
        v2 = 2.0_dp * ran_uniform() - 1.0_dp
        r = v1*v1 + v2*v2
     end do
     ran_normal = v1 * sqrt(-2.0_dp * log(r) / r)
   end function ran_normal

   !% Return a random real distributed exponentially between zero and positive infinity
   !% with mean and variance of unity
   function ran_exp() result(r)    
     real(dp) :: r
     r = -log(ran_uniform())
   end function ran_exp

   !% Return a random string of length l containing the characters A-Z, a-z, and 0-9
   function ran_string(l)
   
     integer, intent(in) :: l
     character(len=l)    :: ran_string    
     integer             :: i, n

     do i = 1, l
        n = mod(ran(),62)
        if (n < 10) then
           ran_string(i:i) = achar(n+48) ! 0-9
        else if (n < 36) then
           ran_string(i:i) = achar(n+55) ! A-Z
        else
           ran_string(i:i) = achar(n+61) ! a-z
        end if           
     end do

   end function ran_string

   subroutine current_times(cpu_t, wall_t, mpi_t)
      real(dp), intent(out), optional :: cpu_t, wall_t, mpi_t
      integer wall_t_count, count_rate, max_count
#ifdef _MPI
     include "mpif.h"
#endif

      if (present(cpu_t)) call cpu_time(cpu_t)
      if (present(wall_t)) then
	 call system_clock(wall_t_count, count_rate, max_count)
	 wall_t = real(wall_t_count,dp) * 1.0_dp/real(count_rate,dp)
      endif
#ifdef _MPI
      if (present(mpi_t)) mpi_t = MPI_Wtime()
#else
      if (present(mpi_t)) mpi_t = 0.0_dp
#endif
   end subroutine current_times

   !% Measure elapsed CPU and wall clock time between pairs of calls with 
   !% matching 'name' parameter. Calls to 'system_timer' must be properly 
   !% nested (i.e. start and stop from different pairs can't overlap), and
   !% maximum depth of calls is set by the 'TIMER_STACK' parameter.
   !%
   !%>   call system_timer(name)  start the clock
   !%>   ...                      do something
   !%>   call system_timer(name)  stop clock and print elapsed time
   !%>
   !%> If optional do_always argument is true, routine will do its thing even 
   !%> if system_do_timing is false.
   !
   subroutine system_timer(name, do_always, time_elapsed, do_print)
#ifdef _OPENMP
     use omp_lib
#endif
     character(len=*), intent(in) :: name !% Unique identifier for this timer
     logical, intent(in), optional :: do_always
     real(dp), intent(out), optional :: time_elapsed
     logical, intent(in), optional :: do_print

     integer, save :: stack_pos = 0
     real(dp), save, dimension(TIMER_STACK) :: wall_t0, cpu_t0
     character(len=255), save, dimension(TIMER_STACK) :: names
     character(len=50) :: out_name

     logical my_do_always, my_do_print

     logical :: found_name
     real(dp) :: cpu_t1, wall_t1
#ifdef _MPI
     real(dp), save ::  mpi_t0(TIMER_STACK)
     real(dp) :: mpi_t1
#endif

     my_do_always = optional_default(.false., do_always)
     my_do_print = optional_default(.true., do_print)

     if (.not. my_do_always .and. .not. system_do_timing) return

#ifdef _OPENMP
     ! stacks are not thread safe, so if we're in an OpenMP parallel region do nothing
     if (omp_get_num_threads() > 1) return
#endif

     found_name = .false.
     if (stack_pos >= 1) found_name = trim(names(stack_pos)) == trim(name)

     if (.not. found_name) then
        ! Start a new timer
        stack_pos = stack_pos + 1
        if (stack_pos > TIMER_STACK) then
           call print_warning('System_Timer: stack overflow, name ' // trim(name))
           return
        end if

        names(stack_pos) = trim(name)

#ifdef _MPI
	call current_times(cpu_t0(stack_pos), wall_t0(stack_pos), mpi_t0(stack_pos))
#else
	call current_times(cpu_t0(stack_pos), wall_t0(stack_pos))
#endif

	if (present(time_elapsed)) time_elapsed = 0.0_dp

     else
        ! Stop the most recently started timer
#ifdef _MPI
	call current_times(cpu_t1, wall_t1, mpi_t1)
#else
	call current_times(cpu_t1, wall_t1)
#endif

        out_name = name
#ifndef _MPI
	if (present(time_elapsed)) time_elapsed = wall_t1-wall_t0(stack_pos)
	if (my_do_print) then
	  call print("TIMER: " // out_name // " done in " // (cpu_t1-cpu_t0(stack_pos)) // &
	     " cpu secs, " // (wall_t1-wall_t0(stack_pos))//" wall clock secs.")
	endif
#else
	if (present(time_elapsed)) time_elapsed = mpi_t1-mpi_t0(stack_pos)
	if (my_do_print) then
	  call print("TIMER: " // out_name // " done in " // (cpu_t1-cpu_t0(stack_pos)) // &
	     " cpu secs, " // (wall_t1-wall_t0(stack_pos))//" wall clock secs, " // &
	     (mpi_t1-mpi_t0(stack_pos)) // " mpi wall secs.")
	endif
#endif

        stack_pos = stack_pos - 1
        if (stack_pos < 0) &
             call system_abort('System_Timer: stack underflow, name ' // trim(name))
     end if

   end subroutine system_timer


   !% Test if the file 'filename' can be accessed.
   function is_file_readable(filename)
     character(len=*), intent(in) :: filename
     logical :: is_file_readable

     integer myunit
     integer stat

     myunit=pick_up_unit()

     open(file=filename,unit=myunit,status="OLD",iostat=stat)

     if(stat == 0)then 
        is_file_readable = .true.
        close(unit=myunit)
     else
        is_file_readable = .false.
     end if

   end function is_file_readable


  subroutine Stack_Initialise(this, value)
    type(stack), intent(inout) :: this
    integer, optional::value

    call Finalise(this)
    allocate(this%val(4))
    if(present(value)) then
       this%val(1) = value
       this%pos = 1
    else
       this%pos = 0
    end if
    
  end subroutine Stack_Initialise

  subroutine Stack_Finalise(this)
    type(stack), intent(inout) :: this

    if (allocated(this%val)) deallocate(this%val)

  end subroutine Stack_Finalise

  subroutine Stack_push(this, val)
    type(Stack), intent(inout) :: this
    integer, intent(in) :: val

    integer, allocatable :: val_t(:)

    if (.not. allocated(this%val)) then
      allocate(this%val(4))
    endif

    if (this%pos+1 > size(this%val)) then
      allocate(val_t(size(this%val)))
      val_t = this%val
      deallocate(this%val)
      allocate(this%val(2*size(val_t)))
      this%val(1:size(val_t)) = val_t
      deallocate(val_t)
    endif

    this%val(this%pos+1) = val
    this%pos = this%pos + 1

  end subroutine Stack_push

  subroutine Stack_pop(this)
    type(Stack), intent(inout) :: this

    if (this%pos > 0) then
      this%pos = this%pos - 1
    else
      call system_abort("Underflow in Stack_pop")
    endif
  end subroutine Stack_pop

  function stack_value(this)
    type(Stack), intent(in) :: this
    integer :: stack_value

    if (this%pos > 0) then
      stack_value = this%val(this%pos)
    else
      call system_abort("Called stack_value on empty stack, pos = " // this%pos)
    endif
  end function stack_value

  subroutine Stack_print(this, verbosity, out)
    type(Stack), intent(in) :: this
    integer, intent(in), optional :: verbosity
    type(inoutput), intent(in), optional :: out

    call Print("Stack:", verbosity, out)
    if (allocated(this%val)) then
      call Print("Stack: size " // size(this%val), verbosity, out)
      call Print("Stack: val " // this%val(1:this%pos), verbosity, out)
    endif
  end subroutine Stack_Print

  !% Map from verbsoity codes to descriptive strings
  function verbosity_to_str(val) result(str)
    integer, intent(in) :: val
    character(10) :: str

    select case(val)
       case(PRINT_ALWAYS)
          str = 'ERROR'
       case(PRINT_SILENT)
          str = 'SILENT'
       case(PRINT_NORMAL)
          str = 'NORMAL'
       case(PRINT_VERBOSE)
          str = 'VERBOSE'
       case(PRINT_NERD)
          str = 'NERD'
       case(PRINT_ANAL)
          str = 'ANAL'
    end select
  end function verbosity_to_str

  !% Map from descriptive verbosity names ('NORMAL', 'VERBOSE' etc.) to numbers
  function verbosity_of_str(str) result(val)
    character(len=*), intent(in) :: str
    integer :: val

    if (trim(str) == 'ERROR') then
       val = PRINT_ALWAYS
    else if (trim(str) == 'SILENT') then
       val = PRINT_SILENT
    else if (trim(str) == 'NORMAL') then
       val = PRINT_NORMAL
    else if (trim(str) == 'VERBOSE') then
       val = PRINT_VERBOSE
    else if (trim(str) == 'NERD') then
       val = PRINT_NERD
    else if (trim(str) == 'ANAL') then
       val = PRINT_ANAL
    else
       call system_abort("verbosity_of_str failed to understand '"//trim(str)//"'")
    end if
  end function verbosity_of_str
    

  !% Push a value onto the verbosity stack
  !% Don't ever lower the verbosity if verbosity minimum is set,
  !%   but always push _something_
  subroutine verbosity_push(val)
    integer, intent(in) :: val 

    if ((value(mainlog%verbosity_cascade_stack) == 0) .or. &
        val > value(mainlog%verbosity_stack)) then
      call push(mainlog%verbosity_stack, val)
    else
      call push(mainlog%verbosity_stack, value(mainlog%verbosity_stack))
    endif
    !call print('verbosity_push now '//current_verbosity(), PRINT_ALWAYS)
  end subroutine verbosity_push

  !% pop the current verbosity value off the stack
  subroutine verbosity_pop()
    call pop(mainlog%verbosity_stack)
    !call print('verbosity_pop now '//current_verbosity(), PRINT_ALWAYS)
  end subroutine verbosity_pop

  !% return the current value of verbosity
  function current_verbosity()
    integer :: current_verbosity
    current_verbosity = value(mainlog%verbosity_stack)
  end function current_verbosity

  !% push the current value + n onto the stack
  subroutine verbosity_push_increment(n)
    integer, intent(in), optional :: n
    integer my_n

    my_n = 1
    if(present(n)) my_n = n
    call verbosity_push(value(mainlog%verbosity_stack)+my_n)
  end subroutine verbosity_push_increment

  !% push the current value - n onto the stack
  subroutine verbosity_push_decrement(n)
    integer, intent(in), optional :: n
    integer my_n

    my_n = 1
    if(present(n)) my_n = n
    call verbosity_push(value(mainlog%verbosity_stack)-my_n)
  end subroutine verbosity_push_decrement

  !% set the minimum verbosity value, by pushing value onto
  !% stack and pushing 1 on to verbosity_cascade_stack
  subroutine verbosity_set_minimum(verbosity)
    integer, intent(in) :: verbosity

    call push(mainlog%verbosity_cascade_stack, 1)
    call verbosity_push(verbosity)
  end subroutine verbosity_set_minimum

  !% unset the minimum verbosity value, by popping value from
  !% stack and popping from verbosity_cascade_stack
  subroutine verbosity_unset_minimum()
    call verbosity_pop()
    call pop(mainlog%verbosity_cascade_stack)
  end subroutine verbosity_unset_minimum

  pure function optional_default_l(def, opt_val)
    logical, intent(in) :: def
    logical, intent(in), optional :: opt_val
    logical :: optional_default_l

    if (present(opt_val)) then
      optional_default_l = opt_val
    else
      optional_default_l = def
    endif

  end function optional_default_l

  pure function optional_default_i(def, opt_val)
    integer, intent(in) :: def
    integer, intent(in), optional :: opt_val
    integer :: optional_default_i

    if (present(opt_val)) then
      optional_default_i = opt_val
    else
      optional_default_i = def
    endif

  end function optional_default_i

  pure function optional_default_ia(def, opt_val)
    integer, intent(in) :: def(:)
    integer, intent(in), optional :: opt_val(size(def))
    integer :: optional_default_ia(size(def))

    if (present(opt_val)) then
      optional_default_ia = opt_val
    else
      optional_default_ia = def
    endif

  end function optional_default_ia


  pure function optional_default_r(def, opt_val)
    real(dp), intent(in) :: def
    real(dp), intent(in), optional :: opt_val
    real(dp) :: optional_default_r

    if (present(opt_val)) then
      optional_default_r = opt_val
    else
      optional_default_r = def
    endif

  end function optional_default_r

  pure function optional_default_ra(def, opt_val)
    real(dp), intent(in) :: def(:)
    real(dp), intent(in), optional :: opt_val(size(def))
    real(dp) :: optional_default_ra(size(def))

    if (present(opt_val)) then
      optional_default_ra = opt_val
    else
      optional_default_ra = def
    endif

  end function optional_default_ra


  pure function optional_default_z(def, opt_val)
    complex(dp), intent(in) :: def
    complex(dp), intent(in), optional :: opt_val
    complex(dp) :: optional_default_z

    if (present(opt_val)) then
      optional_default_z = opt_val
    else
      optional_default_z = def
    endif

  end function optional_default_z

  pure function optional_default_c(def, opt_val)
    character(len=*), intent(in) :: def
    character(len=*), intent(in), optional :: opt_val
    character(SYSTEM_STRING_LENGTH) :: optional_default_c

    if (present(opt_val)) then
      optional_default_c = opt_val
    else
      optional_default_c = def
    endif

  end function optional_default_c

  subroutine enable_timing()
    system_do_timing = .true.
  end subroutine enable_timing

  subroutine disable_timing()
    system_do_timing = .false.
  end subroutine disable_timing

  function get_timing()
    logical :: get_timing
    get_timing = system_do_timing
  end function get_timing

  subroutine set_timing(do_timing)
    logical :: do_timing

    system_do_timing = do_timing
  end subroutine set_timing

  function get_quippy_running()
    logical :: get_quippy_running

    get_quippy_running = system_quippy_running

  end function get_quippy_running

  function increase_stack(stack_size)
    integer, intent(in) :: stack_size
    integer :: increase_stack

    integer, external :: c_increase_stack

    increase_stack = c_increase_stack(stack_size)
  end function increase_stack

  !% Abort with a useful message if an MPI routine returned an error status
  subroutine abort_on_mpi_error(error_code, routine_name)
    integer, intent(in) :: error_code
    character(len=*), intent(in) :: routine_name

#ifdef _MPI
    include 'mpif.h'

    character(MPI_MAX_ERROR_STRING)::error_string
    integer::error_string_length, my_error_code

    if(error_code .ne. MPI_SUCCESS) then
       call  MPI_ERROR_STRING(error_code, error_string, error_string_length, my_error_code)
       if(my_error_code .ne. MPI_SUCCESS) then
	  call system_abort(trim(routine_name) // " returned with error code = " // error_code &
	    // ", which could not be parsed")
       else
	  call system_abort(trim(routine_name) // " had error '"  // trim(error_string) // "'")
       endif
    endif
#else
    call system_abort("abort_on_mpi_error called with routine_name='"//trim(routine_name)//"' " // &
      " error_code " // error_code // " even though MPI is off")
#endif
  end subroutine abort_on_mpi_error


  subroutine parallel_print(lines, comm, verbosity, file)
    character(len=*), intent(in) :: lines(:)
    integer, intent(in) :: comm
    integer, intent(in), optional :: verbosity
    type(inoutput), intent(inout), optional :: file

    integer i
#ifdef _MPI
    include 'mpif.h'

    integer proc
    integer mpi_stat(MPI_STATUS_SIZE)
    integer, allocatable :: lengths(:)
    integer :: lengths_send(2)
    character, allocatable :: in_lines(:)
    integer n_c, line_l, n_lines
    integer err

    if (mpi_id() == 0) then
#endif
      do i=1, size(lines)
	call Print(trim(lines(i)), verbosity, file)
      end do 
#ifdef _MPI
      allocate(lengths(2*mpi_n_procs()))
      lengths_send(1) = size(lines)
      lengths_send(2) = len(adjustr(lines(1)))
      call mpi_gather(lengths_send,2,MPI_INTEGER,lengths,2,MPI_INTEGER,0,comm,err)

      do proc=1, mpi_n_procs()-1
	n_lines = lengths((proc)*2+1)
	line_l = lengths((proc)*2+2)
	allocate(in_lines(n_lines*line_l))
	call mpi_recv(in_lines,n_lines*line_l,MPI_CHARACTER,proc,100+proc,comm,mpi_stat,err)
	do i=1, n_lines*line_l, line_l
	  call Print(in_lines(i:i+line_l-1), verbosity, file)
	end do
	deallocate(in_lines)
      end do

      deallocate(lengths)

    else
      line_l = len(adjustr(lines(1)))
      lengths_send(1) = size(lines)
      lengths_send(2) = line_l
      call mpi_gather(lengths_send,2,MPI_INTEGER,lengths,2,MPI_INTEGER,0,comm,err)
      n_c = line_l*size(lines)
      call mpi_send(lines,n_c,MPI_CHARACTER,0,100+mpi_id(),comm,err)
    endif
#endif

  end subroutine parallel_print

  subroutine ALLOC_TRACE(str,amt)
    character(len=*), intent(in) :: str
    integer, intent(in) :: amt

    if (trace_memory) then
      traced_memory = traced_memory + amt
      call print("TR_ALLOCATE " // str // " " // amt // " " // traced_memory)
    endif
  end subroutine ALLOC_TRACE

  subroutine DEALLOC_TRACE(str,amt)
    character(len=*), intent(in) :: str
    integer, intent(in) :: amt
    if (trace_memory) then
      traced_memory = traced_memory - amt
      call print("TR_ALLOCATE " // str // " " // (-amt) // " " // traced_memory)
    endif
  end subroutine DEALLOC_TRACE

#ifdef _OPENMP
  function system_omp_get_num_threads()
    use omp_lib
    integer :: system_omp_get_num_threads

!$omp parallel
!$omp master
    system_omp_get_num_threads = omp_get_num_threads()
!$omp end master
!$omp end parallel
  end function system_omp_get_num_threads

  subroutine system_omp_set_num_threads(threads)
    use omp_lib
    integer, intent(in) :: threads

    call omp_set_num_threads(threads)
  end subroutine system_omp_set_num_threads
#endif  


!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X These functions provide low-level access to some MPI global variables
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  !% Return this processes' MPI ID
  pure function mpi_id() result (id)
    integer::id
    id = mpi_myid
  end function mpi_id

  !%  Return the total number of MPI processes.
  pure function mpi_n_procs() result (n)
    integer::n
    n = mpi_n
  end function mpi_n_procs

  function reference_true()
    logical :: reference_true
    
    reference_true = .true.

  end function reference_true

  function reference_false()
    logical :: reference_false
    
    reference_false = .false.

  end function reference_false

!% String to character array
function s2a(s) result(a)
  character(len=*), intent(in) :: s
  character(len=1), dimension(len(s)) :: a

  integer i
  
  do i=1,len(s)
     a(i) = s(i:i)
  end do

end function s2a

!% Character array to string
function a2s(a) result(s)
  character(len=1), dimension(:), intent(in) :: a
  character(len=size(a)) :: s

  integer i

  do i=1,size(a)
     s(i:i) = a(i)
  end do
  
end function a2s

!% String to padded character array of length l
function pad(s,l) result(a)
  character(len=*), intent(in) :: s
  integer, intent(in) :: l
  character(len=1), dimension(l) :: a

  integer i

  a = ' '
  do i=1,min(len(s),size(a))
     a(i) = s(i:i)
  end do
end function pad

  function make_run_directory(basename, force_run_dir_i, run_dir_i, error) result(dir)
    character(len=*), optional :: basename
    integer, intent(in), optional :: force_run_dir_i
    integer, intent(out), optional :: run_dir_i
    integer, intent(out), optional :: error
    character(len=SYSTEM_STRING_LENGTH) :: dir

    integer i
    character(len=SYSTEM_STRING_LENGTH) :: use_basename
    integer :: use_force_run_dir_i
   
    logical :: exists
    integer stat

    INIT_ERROR(error)

    use_force_run_dir_i = optional_default(-1, force_run_dir_i)
    use_basename = optional_default("run", basename)
   
    if (use_force_run_dir_i >= 0) then
      i = use_force_run_dir_i
      dir = trim(use_basename)//"_"//i
      call system_command("bash -c '[ -d "//trim(dir)//" ]'", status=stat)
      if (stat /= 0) then
         RAISE_ERROR("make_run_directory got force_run_dir_i="//use_force_run_dir_i//" but "//trim(dir)//" doesn't exist or isn't a directory", error)
      endif
    else
      exists = .true.
      i = 0
      do while (exists)
        i = i + 1
        dir = trim(use_basename)//"_"//i
        call system_command("bash -c '[ -e "//trim(dir)//" ]'", status=stat)
	exists = (stat == 0)
      end do
      call system_command("mkdir "//trim(dir), status=stat)
      if (stat /= 0) then
	 RAISE_ERROR("Failed to mkdir "//trim(dir)//" status " // stat, error)
      endif
    endif

    if (present(run_dir_i)) run_dir_i = i

  end function make_run_directory

  function link_run_directory(sourcename, basename, run_dir_i, error) result(dir)
    character(len=*), intent(in) :: sourcename
    character(len=*), optional :: basename
    integer, intent(out), optional :: run_dir_i
    integer, intent(out), optional :: error
    character(len=SYSTEM_STRING_LENGTH) :: dir

    integer i
    character(len=SYSTEM_STRING_LENGTH) :: use_basename
   
    logical :: exists
    integer stat

    INIT_ERROR(error)

    use_basename = optional_default("run", basename)
   
    exists = .true.
    i = 0
    do while (exists)
      i = i + 1
      dir = trim(use_basename)//"_"//i
      call system_command("bash -c '[ -e "//trim(dir)//" ]'", status=stat)
      exists = (stat == 0)
    end do
    call system_command("ln -s "//trim(sourcename)//" "//trim(dir), status=stat)
    if (stat /= 0) then
       RAISE_ERROR("Failed to link "//trim(dir)//" status " // stat, error)
    endif

    if (present(run_dir_i)) run_dir_i = i

  end function link_run_directory

  function current_version()
    integer :: current_version
    character(len=SYSTEM_STRING_LENGTH) :: string

    integer :: last_num

#ifdef SVN_VERSION    
    string = SVN_VERSION
#else    
    string = "0"
#endif

    do last_num=0, len(trim(string))
      if (last_num == len(trim(string))) exit
      if (scan(string(last_num+1:last_num+1),'0123456789') /= 1) exit
    end do
    if (last_num > 0) then
       current_version = string_to_int(string(1:last_num))
    else
       current_version = 0
    endif

  endfunction current_version

  pure function linebreak_string_length(str, line_len) result(length)
    character(len=*), intent(in) :: str
    integer, intent(in) :: line_len
    integer :: length

    length = len_trim(str)+2*len_trim(str)/line_len+3

  end function linebreak_string_length
 
  function linebreak_string(str, line_len) result(lb_str)
    character(len=*), intent(in) :: str
    integer, intent(in) :: line_len

    character(len=linebreak_string_length(str, line_len)) :: lb_str

    logical :: word_break
    integer :: copy_len, last_space
    character(len=len(lb_str)) :: tmp_str

    lb_str=""
    tmp_str=trim(str)
    do while (len_trim(tmp_str) > 0)
      copy_len = min(len_trim(tmp_str),line_len)
      if (tmp_str(copy_len:copy_len) /= " ") then
	 last_space=scan(tmp_str(1:copy_len), " ", .true.)
	 if ( last_space > 0 .and. (len_trim(tmp_str(1:copy_len)) - last_space) < 4) then
	    copy_len=last_space
	 endif
      endif

      if (len_trim(lb_str) > 0) then ! we already have some text, add newline before concatenating next line
	lb_str = trim(lb_str)//quip_new_line//trim(tmp_str(1:copy_len))
      else ! just concatenate next line
	lb_str = trim(tmp_str(1:copy_len))
      endif
      ! if we broke in mid word, add "-"
      word_break = .true.
      if (tmp_str(copy_len:copy_len) == " ") then ! we broke right after a space, so no wordbreak
	word_break = .false.
      else ! we broke after a character
	if (copy_len < len_trim(tmp_str)) then ! there's another character after this one, check if it's a space
	  if (tmp_str(copy_len+1:copy_len+1) == " ") then
	    word_break = .false.
	  endif
	else ! we broke after the last character
	  word_break = .false.
	endif
      endif
      if (word_break) lb_str = trim(lb_str)//"-"
      tmp_str(1:copy_len) = ""
      tmp_str=adjustl(tmp_str)
    end do

   end function linebreak_string

   subroutine mem_info(total_mem, free_mem)
     real(8), intent(out) :: total_mem, free_mem
     call c_mem_info(total_mem, free_mem)
   end subroutine mem_info

   subroutine wait_for_file_to_exist(filename, max_wait_time, cycle_time, error)
      character(len=*) filename
      real(dp) :: max_wait_time
      real(dp), optional :: cycle_time
      integer, intent(out), optional :: error

      real(dp) :: total_wait_time, use_cycle_time
      integer :: usleep_cycle_time
      logical :: file_exists

      INIT_ERROR(error)

      use_cycle_time = optional_default(0.1_dp, cycle_time)
      usleep_cycle_time = int(use_cycle_time*1000000)

      inquire(file=filename, exist=file_exists)
      total_wait_time = 0.0_dp
      do while (.not. file_exists)
	 call fusleep(usleep_cycle_time)
	 total_wait_time = total_wait_time + use_cycle_time
	 inquire(file=filename, exist=file_exists)
	 if (.not. file_exists .and. total_wait_time > max_wait_time) then
	    RAISE_ERROR("error waiting too long for '"//trim(filename)//"' to exist", error)
	 endif
      end do

   end subroutine wait_for_file_to_exist


end module system_module
