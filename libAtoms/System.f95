!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X     libAtoms: atomistic simulation library
!X     
!X     Copyright 2006-2007.
!X
!X     Authors: Gabor Csanyi, Steven Winfield, James Kermode
!X     Contributors: Noam Bernstein, Alessio Comisso
!X
!X     The source code is released under the GNU General Public License,
!X     version 2, http://www.gnu.org/copyleft/gpl.html
!X
!X     If you would like to license the source code under different terms,
!X     please contact Gabor Csanyi, gabor@csanyi.net
!X
!X     When using this software, please cite the following reference:
!X
!X     http://www.libatoms.org
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
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

! $Id: System.f95,v 1.114 2008-07-14 10:21:15 jrk33 Exp $

! $Log: not supported by cvs2svn $
! Revision 1.113  2008/07/02 13:37:31  jrk33
! Added optional_default_la and optional_default_ra for int and real arrays
!
! Revision 1.112  2008/06/23 20:09:34  jrk33
! Default to opening files read only
!
! Revision 1.111  2008/05/12 21:54:58  jrk33
! Do not abort if system_timer stack mismatched
!
! Revision 1.110  2008/05/07 15:51:07  nb326
! Clean up ifort warnings, mostly unused variables
!
! Revision 1.109  2008/05/01 02:26:53  nb326
! Add time_elapsed optional argument to system_timer()
!
! Revision 1.108  2008/04/25 18:34:39  nb326
! Add master_only argument to inoutput_initialise
!
! Revision 1.107  2008/04/24 12:47:35  jrk33
! Trivial change to formatting
!
! Revision 1.106  2008/03/14 19:55:26  nb326
! Add optional_default_z
!
! Revision 1.105  2008/03/12 09:58:03  jrk33
! Corrected docstring for verbosity_to_str
!
! Revision 1.104  2008/03/12 09:56:49  jrk33
! Added verbosity_to_str and str_to_verbosity to allow verbosity to be specified in a meaninful way in input files
!
! Revision 1.103  2008/02/21 11:27:59  jrk33
! Added some OMIT strings to docs
!
! Revision 1.102  2008/02/12 18:19:26  jrk33
! Shortened length of string returned by string_cat_int_array
!
! Revision 1.101  2008/02/11 18:46:33  nb326
! Increate TIMER_STACK, add memory tracing
!
! Revision 1.100  2008/02/11 16:52:02  jrk33
! Added missing end sub names
!
! Revision 1.99  2008/02/04 13:34:49  saw44
! Added ran_exp function to draw randomly from an exponential distribution
!
! Revision 1.98  2008/01/08 22:34:00  nb326
! Add string_cat_complex_array
!
! Revision 1.97  2007/12/12 22:49:24  nb326
! Don't try to mpi_init if already initialised
!
! Revision 1.96  2007/12/12 22:35:38  nb326
! If inoutput type is INPUT, make sure file exists when you open it
!
! Revision 1.95  2007/11/28 22:51:35  nb326
! Save trivial amount of time in string_cat_real()
!
! Revision 1.94  2007/11/14 14:21:57  nb326
! Nicer system_abort message for inoutput_initialise
!
! Revision 1.93  2007/11/07 17:02:32  nb326
! Add optional_default_c for strings
!
! Revision 1.92  2007/11/07 10:07:31  gc121
! tweaked printing
!
! Revision 1.91  2007/10/22 16:17:10  jrk33
! ends of functions and subroutines
!
! Revision 1.90  2007/10/22 14:40:52  gc121
! put in some OMP stuff
!
! Revision 1.89  2007/10/19 19:25:22  nb326
! Extend timer stack
!
! Revision 1.88  2007/10/15 16:14:56  nb326
! Fix inoutput_print_char_array() and parallel_print()
!
! Revision 1.87  2007/10/12 19:19:17  nb326
! Define i outside ifdef _MPI in parallel_print
!
! Revision 1.86  2007/10/11 18:00:55  nb326
! Remeber to use mpi_module
!
! Revision 1.85  2007/10/11 16:33:21  nb326
! Add abort_on_mpi_error() and parallel_print() (new semantics), replace decode_mpi_error with abort_on_mpi_error
!
! Revision 1.84  2007/10/08 15:20:02  saw44
! Added mpi_print_id member to inoutput type and mpi_all_inoutput / print_mpi_id interfaces to allow each process to write to files and suppress the printing of the process id into those files
!
! Revision 1.83  2007/10/04 15:40:37  nb326
! Make increase_stack a function that returns an error status
!
! Revision 1.82  2007/10/02 21:01:56  nb326
! Reorder *_format_length functions for PGI compilers, and add external pointer_to
!
! Revision 1.81  2007/09/18 18:27:14  saw44
! Fixed bug in length calculation for reals when rounding. Added ran_string function to return a random alphanumeric string of a required length
!
! Revision 1.80  2007/09/17 15:35:22  nb326
! timing activated via parameter to system_initialise or call to enable_timing
!
! Revision 1.79  2007/09/14 13:09:48  jrk33
! Applied checkcode, added missing names to ends of subroutines and functions
!
! Revision 1.78  2007/09/04 14:06:28  nb326
! Finalise verbosity_stack and verbosity_cascade_stack  in inoutput_finalise to prevent mem leaks
!
! Revision 1.77  2007/09/03 16:22:36  gc121
! sync
!
! Revision 1.76  2007/09/03 13:39:40  gc121
! nothing
!
! Revision 1.75  2007/08/30 14:32:45  gc121
! put in a writeb_real_dim2 which takes a filename
!
! Revision 1.74  2007/08/30 14:04:59  nb326
! Add space after outputing prefix
!
! Revision 1.73  2007/08/30 11:50:37  nb326
! For MPI output proc number on system_abort()
!
! Revision 1.72  2007/08/24 11:18:18  nb326
! Docs for optional_argument()
!
! Revision 1.71  2007/08/22 17:21:05  jrk33
! Added system_resync_rng function
!
! Revision 1.70  2007/08/17 09:21:17  nb326
! Create cmd_arg_count() and get_cmd_arg() that hide F2003/older semantics differences
!
! Revision 1.69  2007/08/16 14:01:44  nb326
! Add optional_default_r for reals
!
! Revision 1.68  2007/08/10 16:52:51  gc121
! broke some more lines, added int_format_length, wrapped isnan for systems which dont have it
!
! Revision 1.67  2007/08/10 16:40:18  nb326
! Fix some typos introduced in v 1.65, should be double checked
!
! Revision 1.66  2007/08/10 16:30:21  nb326
! New function optional_default to return the value of an optional argument if specified, otherwise the default
!
! Revision 1.65  2007/08/10 14:17:17  gc121
! started converting functions to use int_format_length
!
! Revision 1.64  2007/08/09 15:12:22  gc121
! fixed bug in inoutput_initialise, check for presence of verbosity_cascade arg was wrong
!
! Revision 1.63  2007/08/08 17:29:21  gc121
! fixed bug in inoutput_print_string, wrong verbosity object was used
!
! Revision 1.62  2007/08/08 16:19:50  nb326
! Move verbosity into type(inoutput), and make verbosity_push type routines operate on mainlog%verbosity
!
! Revision 1.61  2007/08/08 08:50:18  gc121
! put in optional argument to put something on the stack on init. changed system abort not to use our printing
!
! Revision 1.60  2007/08/07 09:56:04  gc121
! made double precision default, with -DQUAD_PRECISION option giving the obvious
!
! Revision 1.59  2007/08/07 09:15:11  gc121
! reversed inadvertent quad precision change to System.f95
!
! Revision 1.58  2007/08/03 18:08:38  gc121
! sync
!
! Revision 1.57  2007/07/31 16:09:18  nb326
! Don't use direct access for unformatted file, so no need for inoutput%current_record or transfer() in writeb_* and readb_*
!
! Revision 1.56  2007/07/20 16:15:05  nb326
! Make error _very_ negative
!
! Revision 1.55  2007/07/19 15:57:00  gc121
! added more string_cat stuff
!
! Revision 1.54  2007/07/18 15:02:05  nb326
! Use new verbosity system (intro.tex v 1.13)
!
! Revision 1.53  2007/07/18 12:47:35  nb326
! Stack type, and new verbosity system routines
!
! Revision 1.52  2007/07/17 10:35:55  gc121
! added current_verbosity() function
!
! Revision 1.51  2007/07/17 08:44:02  gc121
! fixed bug in string_cat_real in case of NaN
!
! Revision 1.50  2007/07/13 15:07:06  jrk33
! Added save attribute to mainlog etc to appease g95
!
! Revision 1.49  2007/07/13 14:22:25  jrk33
! Removed inoutput_print_normal interface for compatibility with not ifort compilers
!
! Revision 1.48  2007/07/11 13:12:19  jrk33
! Removed unused variable
!
! Revision 1.47  2007/07/11 12:36:31  saw44
! fixed bug in string_cat_int/real_array - seg faulted when printing arrays with no elements
!
! Revision 1.46  2007/07/10 13:32:04  jrk33
! Edited verbosity doc comments
!
! Revision 1.45  2007/07/05 15:06:54  gc121
! changed finalise print message
!
! Revision 1.44  2007/07/04 09:56:34  nb326
! Check for NaN in string_cat_real
!
! Revision 1.43  2007/06/28 09:43:52  nb326
! Add string_cat_logical_array
!
! Revision 1.42  2007/06/22 14:23:21  nb326
! proctect F2003 getarg with ifdef, work around some PGI bugs
!
! Revision 1.41  2007/05/31 11:16:48  nb326
! Fix format for string_cat_real_array when all values are 0
!
! Revision 1.40  2007/05/15 13:18:43  jrk33
! Optionally do an integer divide by zero on abort to force a traceback print
!
! Revision 1.39  2007/05/08 12:52:17  nb326
! Properly allocate return value for string_cat_real_array
!
! Revision 1.38  2007/05/08 12:08:31  nb326
! More space for negative signs in string_cat_real_array
!
! Revision 1.37  2007/05/07 09:55:26  nb326
! Add inoutput_print_char_array
!
! Revision 1.36  2007/05/02 15:44:08  nb326
! Make mpi_all_inoutput_flag is now part of type(Inoutput).  Also add poorly tested string_cat_int_array
!
! Revision 1.35  2007/05/02 14:01:05  nb326
! Add string_cat_complex, poorly tested
!
! Revision 1.34  2007/04/27 15:09:04  jrk33
! Added optional verbosity argument to print_title
!
! Revision 1.33  2007/04/27 11:57:45  gc121
! changed a couple of variable names and put in more comments to make the auto-doc nicer
!
! Revision 1.32  2007/04/25 13:46:54  jrk33
! Fixed a couple of bugs in real_cat_string and string_cat_real
!
! Revision 1.31  2007/04/23 13:43:08  saw44
! Added round(r,digits) function to return a correctly sized string containing "r" rounded to "digits" decimal places. "digits" can be zero
!
! Revision 1.30  2007/04/20 17:22:49  jrk33
! Changed // routines to calculate output string length
!
! Revision 1.29  2007/04/20 09:09:50  jrk33
! Updated doc comment to remove ReadB/WriteB names
!
! Revision 1.28  2007/04/18 12:42:58  jrk33
! Fixed some doc strings. Added optional err arguments to string_to_logical and string_to_real to match string_to_int.
!
! Revision 1.27  2007/04/18 01:32:08  gc121
! added more printing routines and changed trimming
!
! Revision 1.26  2007/04/17 17:45:40  jrk33
! Updated some doc strings
!
! Revision 1.25  2007/04/17 11:13:05  jrk33
! Corrected doc comment
!
! Revision 1.24  2007/04/17 09:57:19  gc121
! put copyright statement in each file
!
! Revision 1.23  2007/04/13 13:06:24  saw44
! Added finalise interface for inoutput_finalise
!
! Revision 1.22  2007/04/13 11:58:28  saw44
! Removed trim()s from concatenation functions. Updated some documentation comment lines
!
! Revision 1.21  2007/04/12 15:11:19  gc121
! many convention changes
!
! Revision 1.20  2007/04/03 14:02:01  jrk33
! Updated doc comments
!
! Revision 1.19  2007/03/30 16:44:05  jrk33
! Added verbosity stack. Changed print() arguments
!
! Revision 1.18  2007/03/28 17:47:49  saw44
! Moved #ifdef _MPI around line 2046 to stop unused variable warning
!
! Revision 1.17  2007/03/27 14:20:31  jrk33
! Added complex ReadB and WriteB
!
! Revision 1.16  2007/03/21 16:54:25  nb326
! Optionally call abort() on system_abort() to generate trace for debugger
!
! Revision 1.15  2007/03/12 17:01:28  jrk33
! IN, OUT, INOUT to lowercase. ran_double() -> ran_uniform(). Reformatted documentation
!
! Revision 1.14  2007/03/06 12:05:11  saw44
! Added command line argument reader to System_Initialise
!
! Revision 1.13  2007/03/01 13:51:46  jrk33
! Documentation comments reformatted and edited throughout. Anything starting 
! "!(no space)%" is picked up by the documentation generation script
!
! Revision 1.12  2007/02/16 14:44:17  saw44
! Added Print_Title, to pretty print titles to stdout
!
! Revision 1.11  2007/01/26 15:54:02  saw44
! Fixed bug in th function for numbers over 100
!
! Revision 1.10  2007/01/24 11:24:17  saw44
! Added "th" function to return the correct ordinal ending for a given integer
!
! Revision 1.9  2007/01/09 17:00:28  saw44
! Added System_Reseed_RNG subroutine to reseed the random number generator. Handy if restarting from a checkpoint
!
! Revision 1.8  2007/01/05 18:02:04  nb326
! Print out err code if passed to system_abort
!
! Revision 1.7  2007/01/03 10:38:11  nb326
! MPI is a C++ namespace, replace ifdef MPI with ifdef _MPI
!
! Revision 1.6  2007/01/02 16:52:12  nb326
! 10K line
!
! Revision 1.5  2006/12/19 10:17:24  jrk33
! Call MPI_Abort from System_Abort to kill all MPI processes
!
! Revision 1.4  2006/12/05 16:26:39  nb326
! err in String_To_Int should be logical
!
! Revision 1.3  2006/12/04 18:06:58  nb326
! Return optional error status from String_To_Int, and add String_To_Logical and file_readable
!
! Revision 1.2  2006/12/04 17:16:17  gc121
! got changes from LOTF95 that were checked in late
!
! Revision 1.64  2006/12/04 17:14:02  gc121
! beautified some argument names
!
! Revision 1.63  2006/11/24 12:21:16  saw44
! DP -> dp
!
! Revision 1.62  2006/11/13 16:05:17  saw44
! Wrapped the Inoutput_free routines in a finalise interface
!
! Revision 1.61  2006/11/01 10:09:48  jrk33
! Fixed CVS mangling of System_Timer
!
! Revision 1.60  2006/10/31 09:48:28  jrk33
! Fixed System_Timer blank padding bug
!
! Revision 1.59  2006/10/30 22:21:37  gc121
! kludge to eliminate crap from timer name buffers
!
! Revision 1.58  2006/10/25 12:59:08  jrk33
! Timer name strings now fixed to 30 characters long
!
! Revision 1.57  2006/09/06 13:47:40  jrk33
! Increased command line length limit to 1024 characters
!
! Revision 1.56  2006/07/03 17:06:27  jrk33
! Added MPI timer which calls MPI_WTime() to System_Timer
!
! Revision 1.55  2006/06/26 14:13:31  jrk33
! Typo in System_Timer fixed
!
! Revision 1.54  2006/06/26 13:31:55  jrk33
! Added System_Timer subroutine for measuring elapsed CPU and wall clock time of a section of code
!
! Revision 1.53  2006/06/20 17:23:18  gc121
! added new copyright notice to include James, Gian, Mike and Alessandro
!
! Revision 1.52  2006/06/19 20:03:52  jrk33
! Changed common_seed default to .true. in Hello_World
!
! Revision 1.51  2006/06/19 19:16:06  jrk33
! Added optional parameter common_seed to System_Initialise and Hello_World to use one seed for all MPI nodes
!
! Revision 1.50  2006/06/08 14:03:50  saw44
! Added error checking to Read_Line
!
! Revision 1.49  2006/05/31 11:21:43  jrk33
! Removed debugging line
!
! Revision 1.48  2006/05/30 14:03:02  jrk33
! Changed format string in String_To_Real to simple *. Removed recursive I/O in Hello_World by saving result of date_and_time_string in temp variable
!
! Revision 1.47  2006/05/30 11:11:48  jrk33
! Removed declarations for unused variables
!
! Revision 1.46  2006/05/16 16:47:54  saw44
! Added status returning to Parse_Line
!
! Revision 1.45  2006/05/11 15:46:04  saw44
! Added routines to print single ints and reals to the overloaded Print interface. Removed the unused variable ierr from System_Abort. Added System_Warning, which prints a warning on stderr but does not quit the program.
!
! Revision 1.44  2006/05/02 17:08:42  saw44
! Added extra blank line to System_Finalise and corrected a typo
!
! Revision 1.43  2006/04/21 13:27:05  saw44
! Fixed localisation bug. Thanks Alessio.
!
! Revision 1.42  2006/04/07 14:57:59  saw44
! Fixed bug in reallocate (my bad!): == -> /=
!
! Revision 1.41  2006/04/06 10:56:30  saw44
! moved system_command to System.f95
!
! Revision 1.40  2006/04/06 09:40:33  saw44
! Fixed bug in read_line which would have seg faulted if status was not present. Added more helpful comment and error message.
!
! Revision 1.39  2006/04/05 17:02:20  saw44
! Error checking added to read_line
!
! Revision 1.38  2006/04/04 16:31:41  saw44
! Fixed bug in String_To_Real: Integers (no decimal point) were read incorrectly
!
! Revision 1.37  2006/03/31 16:25:21  saw44
! Added Parse_String and changed parse_line to use parse_string
!
! Revision 1.36  2006/03/09 12:52:55  gc121
! print warning if seeds are different for MPI processes
!
! Revision 1.35  2006/03/01 13:57:04  saw44
! Added mulitple inoutput type freeing (up to 4 at once), cleaned up interface blocks a bit
!
! Revision 1.34  2006/02/28 14:58:47  saw44
! Added reallocate subroutines to perform reallocation only if needed
!
! Revision 1.33  2006/02/27 11:09:23  gc121
! put in MPI hello world message and other MPI stuff. corrected bug in time printing (out-of-bounds array access
!
! Revision 1.32  2006/02/21 14:52:46  saw44
! Added date_and_time_string function, system_finalise now prints the time
!
! Revision 1.31  2006/02/20 14:56:47  saw44
! Fixed incorrect time display in Hello_World
!
! Revision 1.30  2006/02/08 17:51:58  saw44
! Added String_To_Int/Real conversion functions
!
! Revision 1.29  2006/02/08 12:27:48  saw44
! Added Inoutput_Read_Line for formatted reading
!
! Revision 1.28  2006/02/06 16:46:55  saw44
! General Code clean-up: some routine names changed, some error messages changed, some code tweaks
!
! Revision 1.27  2006/01/31 14:41:57  saw44
! finalize -> finalise
!
! Revision 1.26  2006/01/31 14:28:19  saw44
! Updated ReadB and WriteB argument order
!
! Revision 1.25  2006/01/31 13:59:12  gc121
! added bye-bye message
!
! Revision 1.24  2006/01/31 11:52:56  saw44
! Changed formatted output error messages
!
! Revision 1.23  2006/01/30 10:40:51  gc121
! removed logger from print calls, its the default
!
! Revision 1.22  2006/01/26 16:10:44  gc121
! added verbosity to printing, fixed function names
!
! Revision 1.21  2006/01/26 01:56:07  gc121
! beautified printing
!
! Revision 1.20  2006/01/24 12:12:57  saw44
! Added workaround for PathScale compiler bug
!
! Revision 1.19  2006/01/20 09:16:39  saw44
! Updated test routine for new style system_init
!
! Revision 1.18  2006/01/19 14:14:55  saw44
! Errors are now written to stderr and not stdout
!
! Revision 1.17  2006/01/19 13:23:08  saw44
! Made the WriteB/ReadB interface subroutines more nicely named
!
! Revision 1.16  2006/01/18 16:07:59  gc121
! cleanup started by gc121 and saw44
!
! Revision 1.15  2006/01/06 17:14:19  saw44
! Fixed random number generator, changed d0 -> _dp
!
! Revision 1.14  2006/01/05 15:10:51  comisso
! added checks over read iostats
!
! Revision 1.13  2005/12/19 17:00:05  saw44
! Added WriteB/ReadB for 2D arrays, Rewinding, Backspacing
!
! Revision 1.12  2005/12/19 12:05:34  saw44
! Alessios update with direct file access
!
! Revision 1.11  2005/12/08 12:41:00  saw44
! Added normally distributed random generator ran_normal()
!
! Revision 1.10  2005/12/07 15:22:31  saw44
! line(1:80) -> trim(line) in Inoutput_Print
!
! Revision 1.9  2005/12/07 10:33:04  saw44
! Fixed System_Abort. No longer locks up the shell!
!
! Revision 1.8  2005/11/29 17:06:55  saw44
! Fixed extra line after every Print() call
!
! Revision 1.7  2005/11/24 10:40:31  comisso
! Small change to logical comparison apllied to get compiled under ibm xlf95
!
! Revision 1.6  2005/11/22 16:11:18  comisso
! #defines in system have been substituted with parameters
!
! Revision 1.5  2005/11/22 12:13:52  gc121
! removed type length variables
!
! Revision 1.4  2005/11/17 11:59:39  comisso
! test
!
! Revision 1.3  2005/11/15 12:09:35  saw44
! Added double precision parameter dp = kind(1.0d0) and changed real(8) to real(dp)
!
! Revision 1.2  2005/11/11 11:18:59  gc121
! forgot to comment cvs magic variables
!
! Revision 1.1.1.1  2005/11/11 10:22:24  gc121
! starting

module system_module
!$ use omp_lib
  implicit none

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

  integer, parameter :: INTEGER_SIZE = 4
  integer, parameter :: REAL_SIZE = dp
  integer, parameter :: COMPLEX_SIZE = 2*dp
  logical :: trace_memory = .false.
  integer :: traced_memory = 0

  logical, private :: system_do_timing = .false.

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
  end type InOutput

  public   !standard setting for the module
  integer,private                  :: mpi_n, mpi_myid    ! Number of processes and local process ID
  real(dp),private                 :: start_time         ! Initial time
  character(10240)                 :: line               ! 'line' is global and is used by other modules
  character(10240),private         :: local_line         ! 'local_line' is private and System should use this instead         
  type(inoutput),target,save      :: mainlog            !% main output, connected to 'stdout' by default
  type(inoutput),target,save      :: errorlog           !% error output, connected to 'stderr' by default
  type(inoutput),target,save      :: mpilog             !% MPI output, written to by each mpi process
  integer,private                  :: idum               ! used in the random generator
  real(dp),parameter               :: NUMERICAL_ZERO = 1.e-14

  ! system dependent variables 
  integer::RAN_MAX


  ! output labels
  integer,parameter::ERROR   = -100000
  integer,parameter::SILENT  =  -1
  integer,parameter::NORMAL  =   0
  integer,parameter::VERBOSE =   1
  integer,parameter::NERD    =   1000  ! aleph0
  integer,parameter::ANAL    =   10000 ! aleph1

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
  character(255), dimension(MAX_READABLE_ARGS), save :: COMMAND_ARG !% The first 'MAX_READABLE_ARGS' command arguments

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

  private :: writeb_int_dim0, writeb_real_dim0, writeb_log_dim0, writeb_char_dim0, writeb_complex_dim0
  private :: writeb_int_dim1, writeb_real_dim1, writeb_log_dim1, writeb_char_dim1, writeb_complex_dim1
  private :: writeb_int_dim2, writeb_real_dim2, writeb_log_dim2, writeb_char_dim2, writeb_complex_dim2
  private :: writeb_real_dim2_filename
  interface write_binary
    module procedure writeb_int_dim0, writeb_real_dim0, writeb_log_dim0, writeb_char_dim0, writeb_complex_dim0
    module procedure writeb_int_dim1, writeb_real_dim1, writeb_log_dim1, writeb_char_dim1, writeb_complex_dim1
    module procedure writeb_int_dim2, writeb_real_dim2, writeb_log_dim2, writeb_char_dim2, writeb_complex_dim2
    module procedure writeb_real_dim2_filename
  end interface write_binary

  private :: readb_int_dim0, readb_real_dim0, readb_log_dim0, readb_char_dim0, readb_complex_dim0
  private :: readb_int_dim1, readb_real_dim1, readb_log_dim1, readb_char_dim1, readb_complex_dim1
  private :: readb_int_dim2, readb_real_dim2, readb_log_dim2, readb_char_dim2, readb_complex_dim2
  interface read_binary
    module procedure readb_int_dim0, readb_real_dim0, readb_log_dim0, readb_char_dim0, readb_complex_dim0
    module procedure readb_int_dim1, readb_real_dim1, readb_log_dim1, readb_char_dim1, readb_complex_dim1
    module procedure readb_int_dim2, readb_real_dim2, readb_log_dim2, readb_char_dim2, readb_complex_dim2
  end interface read_binary

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
     module procedure string_cat_complex_array
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
     subroutine system_command(command,status)
       character(*), intent(in) :: command
       integer, optional,intent(out) :: status
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
  !% For unformatted output, 
  !% for use with the 'read_binary' and 'write_binary' interfaces, the 
  !% 'isformatted' optional parameter must
  !% be set to false.
  subroutine inoutput_initialise(this,filename,action,isformatted,append,verbosity,verbosity_cascade,master_only)
    type(Inoutput), intent(inout)::this
    character(*),intent(in),optional::filename
    logical,intent(in),optional::isformatted 
    integer,intent(in),optional::action
    logical,intent(in),optional::append
    integer,intent(in),optional :: verbosity, verbosity_cascade
    logical,intent(in),optional :: master_only

    character(32)::formattedstr
    character(32)::position_value
    integer::stat
    logical :: my_master_only

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
       call system_abort(" Append not implemented for unformatted output")
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
             call system_abort('IO error opening "'//trim(filename)//'", error number: '//stat)
          end if

       end if
    else ! no file name passed, so we go to standard output
       this%filename='stdout' 
       this%unit = 6 !standard output 
       this%action = OUTPUT        
    end if  ! end default--------------------------------------------------

    this%prefix = ''
    this%postfix = ''
    this%default_real_precision = 16

    call initialise(this%verbosity_stack)
    if (present(verbosity)) then
      call push(this%verbosity_stack, verbosity)
    else
      call push(this%verbosity_stack, NORMAL)
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

  subroutine inoutput_print_string(string, verbosity, file)
    character(*),             intent(in) :: string
    integer,  optional,       intent(in) :: verbosity
    type(Inoutput), optional, target, intent(in) :: file


    type(Inoutput), pointer :: myoutput
    integer :: myverbosity

    ! check inoutput object
    myoutput => mainlog
    if(present(file)) myoutput => file

    ! check verbosity request
    myverbosity = NORMAL
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
	  write(myoutput%unit,'(i0,": ",a)') mpi_id(), trim(string)//trim(myoutput%postfix)
	else
	  write(myoutput%unit,'(a)') trim(string)//trim(myoutput%postfix)
	endif
      else
	if (myoutput%mpi_all_inoutput_flag .and. myoutput%mpi_print_id) then
	  write(myoutput%unit,'(i0,": ",a)') mpi_id(), trim(myoutput%prefix)//" "//trim(string) &
	    // " " // trim(myoutput%postfix)
	else
	  write(myoutput%unit,'(a)') trim(myoutput%prefix)//" "//trim(string)// " "  // trim(myoutput%postfix)
	endif
      endif
    endif
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
    character(102400)                :: inoutput_read_line
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

    local_line = read_line(this,my_status)
    if(present(status)) status = my_status
    if(my_status == 0) call parse_string(local_line,delimiters,fields,&
      num_fields,status=status)
  end subroutine inoutput_parse_line



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
	  call print("splitting string '"//trim(this)//"'", ERROR)
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

  !% Parse a string into fields delimited by certain characters. On exit
  !% the 'fields' array will contain one field per entry and 'num_fields'
  !% gives the total number of fields. 'status' will be given the error status
  !% (if present) and so can be used to tell if an end-of-file occurred.
  subroutine parse_string(this, delimiters, fields, num_fields, matching, status)

    character(*),               intent(in)    :: this
    character(*),               intent(in)    :: delimiters
    character(*), dimension(:), intent(inout) :: fields
    integer,                    intent(out)   :: num_fields
    logical, optional,          intent(in)    :: matching
    integer, optional,          intent(out)   :: status

    integer                                   :: field_start, length
    integer :: delim_pos
    integer :: n_delims, i
    character(len=len(delimiters)) :: opening_delims, closing_delims
    character(len=1) :: opening_delim
    integer :: opening_delim_index
    logical :: do_matching

    do_matching = optional_default(.false., matching)

    if (present(status)) status = 0

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
        if (present(status)) then 
          call print("ERROR: parse_string called with matching=.true. but odd number of delimiters " // (len(delimiters)), ERROR)
          status = 1
          return
        else
          call system_abort("parse_string called with matching=.true. but odd number of delimiters " // (len(delimiters)))
        endif
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
              if (present(status)) then
                call print("ERROR: parse_string ran out of space for fields", ERROR)
                status = 1
                return
              else
                call system_abort("parse_string ran out of space for fields")
              endif
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
            if (present(status)) then
              call print("ERROR: parse_string ran out of space for fields", ERROR)
              status = 1
              return
            else
              call system_abort("parse_string ran out of space for fields")
            endif
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
	  call print("parse_string failed to find closing delimiter to match opening delimiter at position " // (field_start-1), ERROR)
	  call print("parse_string string='"//this//"'", ERROR)
          if (present(status)) then
            call print("ERROR: parse_string failed to find closing delimiter", ERROR)
            status = 1
            return
          else
            call system_abort("parse_string failed to find closing delimiter")
          endif
	else
	  delim_pos = length-field_start+2
	endif
      endif
      delim_pos = delim_pos + field_start - 1
      if (delim_pos-1 >= field_start) then
	num_fields = num_fields + 1
        if (num_fields > size(fields)) then
          if (present(status)) then
            call print("ERROR: parse_string ran out of space for fields", ERROR)
            status = 1
            return
          else
            call system_abort("parse_string ran out of space for fields")
          endif
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

    read(string,*,iostat=stat) string_to_real

    if (present(err)) err = (stat /= 0)


  end function string_to_real

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!% Reallocation: used to reduce the need to deallocate and reallocate
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine reallocate_int1d(array, d1, zero)
    integer, allocatable, dimension(:), intent(inout) :: array
    integer,                            intent(in)    :: d1
    logical, optional,                  intent(in)    :: zero

    if (allocated(array)) then
       if (size(array) /= d1) then
          deallocate(array)
          allocate(array(d1))
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

  subroutine reallocate_int2d(array, d1, d2, zero)
    integer, allocatable, dimension(:,:), intent(inout) :: array
    integer,                              intent(in)    :: d1,d2
    logical, optional,                    intent(in)    :: zero

    if (allocated(array)) then
       if (.not. all(shape(array) == (/d1,d2/))) then
          deallocate(array)
          allocate(array(d1,d2))
       end if
    else
       allocate(array(d1,d2))
    end if

    if (present(zero)) then
       if (zero) array = 0
    end if

  end subroutine reallocate_int2d

  subroutine reallocate_real2d(array, d1, d2, zero)
    real(dp), allocatable, dimension(:,:), intent(inout) :: array
    integer,                               intent(in)    :: d1,d2
    logical, optional,                     intent(in)    :: zero

    if (allocated(array)) then
       if (.not. all(shape(array) == (/d1,d2/))) then
          deallocate(array)
          allocate(array(d1,d2))
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
!% Write scalar and array data to binary files. These
!% interfaces are heavily overloaded to cater for all intrinsic and most
!% derived types.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


  subroutine writeb_int_dim0(int,file)
    type(inoutput),intent(inout):: file
    integer::int

    if(file%action .EQ. INPUT) then
       call system_abort("WriteB_int_dim0: you cannot write to an input object")
    end if

    if (file%formatted) then 
       call system_abort("WriteB_int_dim0: This function cannot be used on formatted objects")
    end if

    if (.not.inoutput_do_output(file)) return

    write(file%unit) int
  end subroutine writeb_int_dim0

  subroutine writeb_real_dim0(r,file)
    type(inoutput),intent(inout):: file
    real(dp)::r

    if(file%action .EQ. INPUT) then
       call system_abort("WriteB_real_dim0: you cannot write to an input object")
    end if

    if (file%formatted) then 
       call system_abort("WriteB_real_dim0: This function cannot be used on formatted objects")
    end if

    if (.not.inoutput_do_output(file)) return

    write (file%unit) r
  end subroutine writeb_real_dim0


  subroutine writeb_log_dim0(l,file)
    type(inoutput),intent(inout):: file
    logical::l

    if(file%action .EQ. INPUT) then
       call system_abort("WriteB_real_dim0: you cannot write to an input obejct")
    end if

    if (file%formatted) then 
       call system_abort("WriteB_real_dim0: This function cannot be used on formatted objects")
    end if

    if (.not.inoutput_do_output(file)) return

    write(file%unit) l
  end subroutine writeb_log_dim0

   subroutine writeb_char_dim0(c,file)
    type(inoutput),intent(inout):: file
    character(*)::c

    if(file%action .EQ. INPUT) then
       call system_abort("WriteB_char_dim0: you cannot write to an input object")
    end if

    if (file%formatted) then 
       call system_abort("WriteB_char_dim0: This function cannot be used on formatted objects")
    end if

    if (.not.inoutput_do_output(file)) return

    write (file%unit) c

  end subroutine writeb_char_dim0


  subroutine writeb_complex_dim0(r,file)
    type(inoutput),intent(inout):: file
    complex(dp)::r

    if(file%action .EQ. INPUT) then
       call system_abort("WriteB_complex_dim0: you cannot write to an input object")
    end if

    if (file%formatted) then 
       call system_abort("WriteB_complex_dim0: This function cannot be used on formatted objects")
    end if

    if (.not.inoutput_do_output(file)) return

    write (file%unit) r

  end subroutine writeb_complex_dim0


  subroutine writeb_int_dim1(intv,file)
    type(inoutput),intent(inout):: file
    integer,intent(in)::intv(:)

    if(file%action .EQ. INPUT) then
       call system_abort("WriteB_int_dim1: you cannot write to an input object")
    end if

    if (file%formatted) then 
       call system_abort("WriteB_int_dim1: This function cannot be used on formatted objects")
    end if

    if (.not.inoutput_do_output(file)) return

    write (file%unit) intv

  end subroutine writeb_int_dim1

   subroutine writeb_real_dim1(rv,file)
    type(inoutput),intent(inout):: file
    real(dp)::rv(:)

    if(file%action .EQ. INPUT) then
       call system_abort("WriteB_real_dim1: you cannot write to an input object")
    end if

    if (file%formatted) then 
       call system_abort("WriteB_real_dim1: This function cannot be used on formatted objects")
    end if

    if (.not.inoutput_do_output(file)) return

    write (file%unit) rv

  end subroutine writeb_real_dim1


 subroutine writeb_log_dim1(lv,file)
    type(inoutput),intent(inout):: file
    logical::lv(:)

    if(file%action .EQ. INPUT) then
       call system_abort("WriteB_log_dim1: you cannot write to an input object")
    end if

    if (file%formatted) then 
       call system_abort("WriteB_log_dim1: This function cannot be used on formatted objects")
    end if

    if (.not.inoutput_do_output(file)) return

    write (file%unit) lv

  end subroutine writeb_log_dim1


   subroutine writeb_char_dim1(cv,file)
    type(inoutput),intent(inout):: file
    character(*)::cv(:)

    if(file%action .EQ. INPUT) then
       call system_abort("WriteB_char_dim1: you cannot write to an input object")
    end if

    if (file%formatted) then 
       call system_abort("WriteB_char_dim1: This function cannot be used on formatted objects")
    end if

    if (.not.inoutput_do_output(file)) return

    write (file%unit) cv

  end subroutine writeb_char_dim1

   subroutine writeb_complex_dim1(rv,file)
    type(inoutput),intent(inout):: file
    complex(dp)::rv(:)

    if(file%action .EQ. INPUT) then
       call system_abort("WriteB_complex_dim1: you cannot write to an input object")
    end if

    if (file%formatted) then 
       call system_abort("WriteB_complex_dim1: This function cannot be used on formatted objects")
    end if

    if (.not.inoutput_do_output(file)) return

    write (file%unit) rv

  end subroutine writeb_complex_dim1


   subroutine writeb_int_dim2(inta2,file)
      type(Inoutput),          intent(inout) :: file
      integer, dimension(:,:), intent(in)    :: inta2

    if(file%action .EQ. INPUT) then
       call system_abort("WriteB_int_dim2: you cannot write to an input object")
    end if

    if (file%formatted) then 
       call system_abort("WriteB_int_dim2: This function cannot be used on formatted objects")
    end if

    if (.not.inoutput_do_output(file)) return

    write (file%unit) inta2

   end subroutine writeb_int_dim2

   subroutine writeb_real_dim2_filename(ra2, filename)
     real(dp), dimension(:,:), intent(in)    :: ra2
     character(len=*) :: filename
     type(inoutput)::file

     call initialise(file, filename, isformatted=.false., action=OUTPUT)
     call writeb_real_dim2(ra2, file)
     call finalise(file)
   end subroutine writeb_real_dim2_filename

   subroutine writeb_real_dim2(ra2,file)
      type(Inoutput),           intent(inout) :: file
      real(dp), dimension(:,:), intent(in)    :: ra2

    if(file%action .EQ. INPUT) then
       call system_abort("WriteB_real_dim2: you cannot write to an input object")
    end if

    if (file%formatted) then 
       call system_abort("WriteB_real_dim2: This function cannot be used on formatted objects")
    end if

    if (.not.inoutput_do_output(file)) return

    write (file%unit) ra2
   end subroutine writeb_real_dim2

   subroutine writeb_log_dim2(la2,file)
      type(Inoutput),           intent(inout) :: file
      logical,  dimension(:,:), intent(in)    :: la2

    if(file%action .EQ. INPUT) then
       call system_abort("WriteB_log_dim2: you cannot write to an input object")
    end if

    if (file%formatted) then 
       call system_abort("WriteB_log_dim2: This function cannot be used on formatted objects")
    end if

    if (.not.inoutput_do_output(file)) return

    write (file%unit) la2
   end subroutine writeb_log_dim2

   subroutine writeb_char_dim2(ca2,file)
      type(Inoutput),               intent(inout) :: file
      character(*), dimension(:,:), intent(in)    :: ca2

    if(file%action .EQ. INPUT) then
       call system_abort("WriteB_char_dim2: you cannot write to an input object")
    end if

    if (file%formatted) then 
       call system_abort("WriteB_char_dim2: This function cannot be used on formatted objects")
    end if

    if (.not.inoutput_do_output(file)) return

    write (file%unit) ca2
   end subroutine writeb_char_dim2

   subroutine writeb_complex_dim2(ra2,file)
      type(Inoutput),           intent(inout) :: file
      complex(dp), dimension(:,:), intent(in)    :: ra2

    if(file%action .EQ. INPUT) then
       call system_abort("WriteB_complex_dim2: you cannot write to an input object")
    end if

    if (file%formatted) then 
       call system_abort("WriteB_complex_dim2: This function cannot be used on formatted objects")
    end if

    if (.not.inoutput_do_output(file)) return

    write (file%unit) ra2
   end subroutine writeb_complex_dim2



!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!% Read scalar and array data from binary files. These
!% interfaces are heavily overloaded to cater for all intrinsic and most
!% derived types.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

   subroutine readb_int_dim0(int,file)
    type(inoutput),intent(inout):: file
    integer,intent(out)::int
    integer::stat

    if(file%action .EQ. OUTPUT) then
       call system_abort("ReadB_int_dim0: you cannot read from an output object")
    end if

    if (file%formatted) then 
       call system_abort("ReadB_int_dim0: This function cannot be used on formatted objects")
    end if

    read(file%unit,IOSTAT=stat) int

    if (stat .GT. 0) call system_abort("ReadB: integer dim0 problem reading from the file, probably reached end of file")

  end subroutine readb_int_dim0

  subroutine readb_real_dim0(r,file)
    type(inoutput),intent(inout):: file
    real(dp)::r

    integer::stat

    if(file%action .EQ. OUTPUT) then
       call system_abort("ReadB_real_dim0: you cannot read from an output object")
    end if

    if (file%formatted) then 
       call system_abort("ReadB_real_dim0: This function cannot be used on formatted objects")
    end if

    read (file%unit,iostat=stat) r

    if (stat .GT. 0) call system_abort("ReadB: real dim0 problem reading from the file, probably reached end of file")

  end subroutine readb_real_dim0

  subroutine readb_log_dim0(l,file)
    type(inoutput),intent(inout):: file
    logical::l

    integer::stat

    if(file%action .EQ. OUTPUT) then
       call system_abort("ReadB_log_dim0: you cannot read from an output object")
    end if

    if (file%formatted) then 
       call system_abort("ReadB_real_dim0: This function cannot be used on formatted objects")
    end if

    read (file%unit, iostat=stat) l

    if (stat .GT. 0) call system_abort("ReadB: logical dim0 problem reading from the file, probably reached end of file")
  end subroutine readb_log_dim0

  subroutine readb_char_dim0(c,file)
    type(inoutput),intent(inout):: file
    character(*)::c

    integer::stat

    if(file%action .EQ. OUTPUT) then
       call system_abort("ReadB_char_dim0: you cannot read from an output object")
    end if

    if (file%formatted) then 
       call system_abort("readB_char_dim0: This function cannot be used on formatted objects")
    end if

    read (file%unit, iostat=stat) c

    if (stat .GT. 0) call system_abort("ReadB: char dim0 problem reading from the file, probably reached end of file")
  end subroutine readb_char_dim0

  subroutine readb_complex_dim0(r,file)
    type(inoutput),intent(inout):: file
    complex(dp)::r

    integer::stat

    if(file%action .EQ. OUTPUT) then
       call system_abort("ReadB_complex_dim0: you cannot read from an output object")
    end if

    if (file%formatted) then 
       call system_abort("ReadB_complex_dim0: This function cannot be used on formatted objects")
    end if

    read (file%unit, iostat=stat) r

    if (stat .GT. 0) call system_abort("ReadB: complex dim0 problem reading from the file, probably reached end of file")
  end subroutine readb_complex_dim0


  subroutine readb_int_dim1(intv,file)
    type(inoutput),intent(inout):: file
    integer::intv(:)

    integer::stat

    if(file%action .EQ. OUTPUT) then
       call system_abort("ReadB_int_dim1: you cannot read from an output object")
    end if

    if (file%formatted) then 
       call system_abort("ReadB_int_dim1: This function cannot be used on formatted objects")
    end if

    read(file%unit, iostat=stat) intv

    if (stat .GT. 0) call system_abort("ReadB: integer dim1 problem reading from the file, probably reached end of file")
  end subroutine readb_int_dim1


  subroutine readb_real_dim1(rv,file)
    type(inoutput),intent(inout):: file
    real(dp)::rv(:)

    integer::stat

    if(file%action .EQ. OUTPUT) then
       call system_abort("ReadB_real_dim1: you cannot read from an output object")
    end if

    if (file%formatted) then 
       call system_abort("ReadB_real_dim1: This function cannot be used on formatted objects")
    end if

    read (file%unit, iostat=stat) rv
    if (stat .GT. 0) call system_abort("ReadB: real dim1 problem reading from the file, probably reached end of file")
  end subroutine readb_real_dim1

  subroutine readb_log_dim1(lv,file)
    type(inoutput),intent(inout):: file
    logical::lv(:)

    integer::stat

    if(file%action .EQ. OUTPUT) then
       call system_abort("WriteB_log_dim1: you cannot read from an output object")
    end if

    if (file%formatted) then 
       call system_abort("WriteB_log_dim1: This function cannot be used on formatted objects")
    end if

    read (file%unit, iostat=stat) lv
    if (stat .GT. 0) call system_abort("ReadB: logical dim1 problem reading from the file, probably reached end of file")
  end subroutine readb_log_dim1

  subroutine readb_char_dim1(cv,file)
    type(inoutput),intent(inout):: file
    character(*)::cv(:)

    integer::stat

    if(file%action .EQ. OUTPUT) then
       call system_abort("ReadB_char_dim1: you cannot read from an output object")
    end if

    if (file%formatted) then 
       call system_abort("ReadB_char_dim1: This function cannot be used on formatted objects")
    end if

    read (file%unit, iostat=stat) cv
    if (stat .GT. 0) call system_abort("ReadB: character dim1 problem reading from the file, probably reached end of file")
  end subroutine readb_char_dim1

    subroutine readb_complex_dim1(rv,file)
    type(inoutput),intent(inout):: file
    complex(dp)::rv(:)

    integer::stat

    if(file%action .EQ. OUTPUT) then
       call system_abort("ReadB_complex_dim1: you cannot read from an output object")
    end if

    if (file%formatted) then 
       call system_abort("ReadB_complex_dim1: This function cannot be used on formatted objects")
    end if

    read (file%unit, iostat=stat) rv

    if (stat .GT. 0) call system_abort("ReadB: complex dim1 problem reading from the file, probably reached end of file")
  end subroutine readb_complex_dim1



  subroutine readb_int_dim2(inta2,file)
      type(Inoutput),          intent(inout) :: file
      integer, dimension(:,:), intent(inout) :: inta2

    integer::stat

    if(file%action .EQ. OUTPUT) then
       call system_abort("ReadB_int_dim2: you cannot read from an output object")
    end if

    if (file%formatted) then 
       call system_abort("ReadB_int_dim2: This function cannot be used on formatted objects")
    end if

    read (file%unit, iostat=stat) inta2

    if (stat .GT. 0) call system_abort("ReadB: integer dim2 problem reading from the file, probably reached end of file")

   end subroutine readb_int_dim2

   subroutine readb_real_dim2(ra2,file)
      type(Inoutput),           intent(inout) :: file
      real(dp), dimension(:,:), intent(inout) :: ra2

    integer::stat

    if(file%action .EQ. OUTPUT) then
       call system_abort("ReadB_real_dim2: you cannot read from an output object")
    end if

    if (file%formatted) then 
       call system_abort("ReadB_real_dim2: This function cannot be used on formatted objects")
    end if

    read (file%unit, iostat=stat) ra2

    if (stat .GT. 0) call system_abort("ReadB: real dim2 problem reading from the file, probably reached end of file")

   end subroutine readb_real_dim2

   subroutine readb_log_dim2(la2,file)
      type(Inoutput),           intent(inout) :: file
      logical,  dimension(:,:), intent(inout) :: la2

    integer::stat

    if(file%action .EQ. OUTPUT) then
       call system_abort("ReadB_log_dim2: you cannot read from an output object")
    end if

    if (file%formatted) then 
       call system_abort("ReadB_log_dim2: This function cannot be used on formatted objects")
    end if

    read (file%unit, iostat=stat) la2

    if (stat .GT. 0) call system_abort("ReadB: logical dim2 problem reading from the file, probably reached end of file")

   end subroutine readb_log_dim2

   subroutine readb_char_dim2(ca2,file)
      type(Inoutput),               intent(inout) :: file
      character(*), dimension(:,:), intent(inout) :: ca2

    integer::stat

    if(file%action .EQ. OUTPUT) then
       call system_abort("ReadB_char_dim2: you cannot read from an output object")
    end if

    if (file%formatted) then 
       call system_abort("ReadB_char_dim2: This function cannot be used on formatted objects")
    end if

    read(file%unit, iostat=stat) ca2

    if (stat .GT. 0) call system_abort("ReadB: character dim2 problem reading from the file, probably reached end of file")
   end subroutine readb_char_dim2

   subroutine readb_complex_dim2(ra2,file)
      type(Inoutput),           intent(inout) :: file
      complex(dp), dimension(:,:), intent(inout) :: ra2

    integer::stat

    if(file%action .EQ. OUTPUT) then
       call system_abort("ReadB_complex_dim2: you cannot read from an output object")
    end if

    if (file%formatted) then 
       call system_abort("ReadB_complex_dim2: This function cannot be used on formatted objects")
    end if

    read (file%unit, iostat=stat) ra2

    if (stat .GT. 0) call system_abort("ReadB: complex dim2 problem reading from the file, probably reached end of file")
   end subroutine readb_complex_dim2

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

      if (size(values)>0) then
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
  subroutine system_initialise(verbosity,seed, mpi_all_inoutput, common_seed, enable_timing)
    integer,intent(in), optional::verbosity           !% mainlog output verbosity
    integer,intent(in), optional::seed                !% Seed for the random number generator.
    logical,intent(in), optional::mpi_all_inoutput    !% Print on all MPI nodes (false by default)
    logical,intent(in), optional::common_seed 
    logical,intent(in), optional::enable_timing           !% Enable system_timer() calls
    !% If 'common_seed' is true (default), random seed will be the same for each
    !% MPI process.
    character(30) :: arg
    integer       :: status, i, n

#ifdef _MPI
    integer::error
    integer :: is_initialised
    include "mpif.h"

    call MPI_initialized(is_initialised, error)
    call abort_on_mpi_error(error, "system_initialise, mpi_initialised()")
    if (is_initialised == 0) then
      call MPI_INIT(error)
      call abort_on_mpi_error(error, "system_initialise, mpi_init()")
    endif
    call get_mpi_size_rank(MPI_COMM_WORLD, mpi_n, mpi_myid)
    if (mpi_n < 1 .or. mpi_myid < 0) &
      call system_abort("system_initialise Got bad size="//mpi_n// " or rank="//mpi_myid//" from get_mpi_size_rank")
#else
    mpi_n  = 1 ! default
    mpi_myid = 0 ! default
#endif

    call cpu_time(start_time)
    mainlog%mpi_all_inoutput_flag = optional_default(.false., mpi_all_inoutput)

    ! Initialise the verbosity stack and default logger
    call initialise(mainlog, verbosity=verbosity)
    call print_mpi_id(mainlog)
    call initialise(errorlog,'stderr')
    call print_mpi_id(errorlog)

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

  end subroutine system_initialise

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

  !% Shut down gracefully, finalising system objects.
  subroutine system_finalise()
    integer :: values(8)
#ifdef _MPI
    integer :: error
    include "mpif.h"
#endif

    call date_and_time(values=values)
    call print("")
    call print('System::Finalise: '//date_and_time_string(values))
    call print("System::Finalise: Bye-Bye!")
    call finalise(mainlog)
    call finalise(errorlog)
#ifdef _MPI
    call mpi_finalize(error)
    call abort_on_mpi_error(error, "system_finalise, mpi_finalise()")
#endif
  end subroutine system_finalise

#ifndef HAVE_QUIPPY
  !% Quit with an error message. Calls 'MPI_Abort' for MPI programs.
  subroutine system_abort(message)
    character(*),      intent(in) :: message
#ifdef IFORT_TRACEBACK_ON_ABORT
    integer :: j
#endif IFORT_TRACEBACK_ON_ABORT
#ifdef SIGNAL_ON_ABORT
    integer :: status
    integer, parameter :: SIGUSR1 = 30
#endif
#ifdef _MPI
    integer::error
    include "mpif.h"
#endif

#ifdef _MPI
    write(unit=errorlog%unit, fmt='(a,i0," ",a)') 'SYSTEM ABORT: proc=',mpi_id(),trim(message)
#else
    write(unit=errorlog%unit, fmt='(a," ",a)') 'SYSTEM ABORT:', trim(message)
#endif

#ifdef _MPI
    call MPI_Abort(MPI_COMM_WORLD, 1, error)
#endif

#ifdef IFORT_TRACEBACK_ON_ABORT
    ! Cause an integer divide by zero error to persuade
    ! ifort to issue a traceback
    j = 1/0
#endif

#ifdef DUMP_CORE_ON_ABORT
    call fabort()
#else
#ifdef SIGNAL_ON_ABORT
    ! send ourselves a USR1 signal rather than aborting
    call kill(getpid(), SIGUSR1, status)
#else
    stop
#endif
#endif
  end subroutine system_abort
#endif

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
    integer, intent(in) :: seed

    integer :: n
    integer, allocatable :: seed_a(:)

    idum=seed              !gabor generator

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

    integer:: actual_seed
    integer:: values(20) ! for time inquiry function
    logical :: use_common_seed
#ifdef _MPI
    integer :: error
    include "mpif.h"    
#endif

    call date_and_time(values=values)

    call print('System::Hello World: '//date_and_time_string(values))
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
          call MPI_Bcast(actual_seed, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, error)
          call abort_on_mpi_error(error, 'Hello_World: MPI_Bcast()')
       end if
#endif
    end if

    call print('System::Hello World: Random Seed = '//actual_seed)
    call system_set_random_seeds(actual_seed)

    call print('System::Hello World: global verbosity = '//value(mainlog%verbosity_stack))
    call print('')
  end subroutine hello_world

  subroutine system_resync_rng
#ifdef _MPI
    integer :: error
    include "mpif.h"    
    call print('Resyncronising random number generator', VERBOSE)

    ! Broadcast seed from process 0 to all others
    call MPI_Bcast(idum, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, error)
    call abort_on_mpi_error(error, 'system_resync_rng')
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

    call print('System::Reseed_RNG: Reseeding random number generator, new seed = '//new_seed,VERBOSE)
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
    if (idum == 0) &
      call system_abort("function ran(): linear-congruential random number generators fail with seed idum=0")
    k = idum/ran_Q
    idum = ran_A*(idum-k*ran_Q)-ran_R*k
    if(idum < 0) idum = idum + ran_M
    ran =idum
  end function ran

  !% Return a random real number uniformly distributed in the range [0,1]
  function ran_uniform()
    real(dp)::ran_uniform
    ran_uniform = 1.1_dp
    do while (ran_uniform > 1.0_dp) ! generating [0,1]  
       ran_uniform = real(ran(),dp)/RAN_MAX
    end do
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
   subroutine system_timer(name, do_always, time_elapsed)
     character(len=*), intent(in) :: name !% Unique identifier for this timer
     logical, optional :: do_always
     real(dp), intent(out), optional :: time_elapsed

     integer, save :: stack_pos = 0
     integer, save, dimension(TIMER_STACK) :: wall_t0
     real(dp), save, dimension(TIMER_STACK) :: cpu_t0
     character(len=255), save, dimension(TIMER_STACK) :: names
     character(len=50) :: out_name

     logical my_do_always

     logical :: found_name
     integer :: count_rate, count_max, wall_t1
     real(dp) :: cpu_t1
#ifdef _MPI
     include "mpif.h"
     real(dp), save ::  mpi_t0(TIMER_STACK)
     real(dp) :: mpi_t1
#endif

     my_do_always = optional_default(.false., do_always)

     if (.not. my_do_always .and. .not. system_do_timing) return

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

        call cpu_time(cpu_t0(stack_pos))
        call system_clock(wall_t0(stack_pos), count_rate, count_max)
#ifdef _MPI
        mpi_t0(stack_pos) = MPI_Wtime()
#endif       

	if (present(time_elapsed)) time_elapsed = 0.0_dp

     else
        ! Stop the most recently started timer
        call cpu_time(cpu_t1)
        call system_clock(wall_t1, count_rate, count_max)
#ifdef _MPI
        mpi_t1 = MPI_Wtime()
#endif


        out_name = name
#ifndef _MPI
	if (present(time_elapsed)) then
	  time_elapsed = real(wall_t1-wall_t0(stack_pos))
	else
	  write (line, '(a,a,a,f0.3,a,f0.3,a)') 'TIMER: ', out_name, &
	       ' done in ', cpu_t1-cpu_t0(stack_pos), ' cpu secs, ', &
	       real(wall_t1 - wall_t0(stack_pos), dp)*1.0_dp/real(count_rate, dp), &
	       ' wall clock secs.'
	  call Print(line)
	endif
#else
	if (present(time_elapsed)) then
	  time_elapsed = real(mpi_t1-mpi_t0(stack_pos))
	else
	  write (line, '(a,a,a,f0.3,a,f0.3,a,f0.3,a)') 'TIMER: ', out_name, &
	       ' done in ', cpu_t1-cpu_t0(stack_pos), ' cpu secs, ', &
	       real(wall_t1 - wall_t0(stack_pos), dp)*1.0_dp/real(count_rate, dp), &
	       ' wall clock secs, ', real(mpi_t1 - mpi_t0(stack_pos)), ' mpi wall secs.'
	  call Print(line)
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
       case(ERROR)
          str = 'ERROR'
       case(SILENT)
          str = 'SILENT'
       case(NORMAL)
          str = 'NORMAL'
       case(VERBOSE)
          str = 'VERBOSE'
       case(NERD)
          str = 'NERD'
       case(ANAL)
          str = 'ANAL'
    end select
  end function verbosity_to_str

  !% Map from descriptive verbosity names ('NORMAL', 'VERBOSE' etc.) to numbers
  function str_to_verbosity(str) result(val)
    character(len=*), intent(in) :: str
    integer :: val

    if (trim(str) == 'ERROR') then
       val = ERROR
    else if (trim(str) == 'SILENT') then
       val = SILENT
    else if (trim(str) == 'NORMAL') then
       val = NORMAL
    else if (trim(str) == 'VERBOSE') then
       val = VERBOSE
    else if (trim(str) == 'NERD') then
       val = NERD
    else if (trim(str) == 'ANAL') then
       val = ANAL
    end if
  end function str_to_verbosity
    

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
    !call print('verbosity_push now '//current_verbosity(), ERROR)
  end subroutine verbosity_push

  !% pop the current verbosity value off the stack
  subroutine verbosity_pop()
    call pop(mainlog%verbosity_stack)
    !call print('verbosity_pop now '//current_verbosity(), ERROR)
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

  function optional_default_l(def, opt_val)
    logical, intent(in) :: def
    logical, intent(in), optional :: opt_val
    logical :: optional_default_l

    if (present(opt_val)) then
      optional_default_l = opt_val
    else
      optional_default_l = def
    endif

  end function optional_default_l

  function optional_default_i(def, opt_val)
    integer, intent(in) :: def
    integer, intent(in), optional :: opt_val
    integer :: optional_default_i

    if (present(opt_val)) then
      optional_default_i = opt_val
    else
      optional_default_i = def
    endif

  end function optional_default_i

  function optional_default_ia(def, opt_val)
    integer, intent(in) :: def(:)
    integer, intent(in), optional :: opt_val(size(def))
    integer :: optional_default_ia(size(def))

    if (present(opt_val)) then
      optional_default_ia = opt_val
    else
      optional_default_ia = def
    endif

  end function optional_default_ia


  function optional_default_r(def, opt_val)
    real(dp), intent(in) :: def
    real(dp), intent(in), optional :: opt_val
    real(dp) :: optional_default_r

    if (present(opt_val)) then
      optional_default_r = opt_val
    else
      optional_default_r = def
    endif

  end function optional_default_r

  function optional_default_ra(def, opt_val)
    real(dp), intent(in) :: def(:)
    real(dp), intent(in), optional :: opt_val(size(def))
    real(dp) :: optional_default_ra(size(def))

    if (present(opt_val)) then
      optional_default_ra = opt_val
    else
      optional_default_ra = def
    endif

  end function optional_default_ra


  function optional_default_z(def, opt_val)
    complex(dp), intent(in) :: def
    complex(dp), intent(in), optional :: opt_val
    complex(dp) :: optional_default_z

    if (present(opt_val)) then
      optional_default_z = opt_val
    else
      optional_default_z = def
    endif

  end function optional_default_z

  function optional_default_c(def, opt_val)
    character(len=*), intent(in) :: def
    character(len=*), intent(in), optional :: opt_val
    character(1024) :: optional_default_c

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

    system_omp_get_num_threads = omp_get_num_threads()
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
  function mpi_id() result (id)
    integer::id
    id = mpi_myid
  end function mpi_id

  !%  Return the total number of MPI processes.
  function mpi_n_procs() result (n)
    integer::n
    n = mpi_n
  end function mpi_n_procs


end module system_module
