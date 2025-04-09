program do_statistics_decorrelated_err
use libatoms_module, only : dp, system_initialise, system_finalise, system_abort, verbosity_push, verbosity_pop, verbosity_of_str, initialise, finalise, read_file, print, operator(//), param_register, param_read_args
use libatoms_module, only : PRINT_SILENT, PRINT_ALWAYS, INPUT, STRING_LENGTH
use libatoms_module, only : inoutput, dictionary
use statistics_module, only : mean_var_decorrelated_err
implicit none
   integer :: i, N
   real(dp), allocatable :: A(:,:)
   real(dp) :: m, v, decorrelation_t, decorrelated_error
   character(len=100), allocatable :: lines(:)
   type(inoutput) :: io
   type(Dictionary) :: cli
   character(len=STRING_LENGTH) :: verbosity_str
   integer :: n_col

   call system_initialise(verbosity=PRINT_SILENT)
   call initialise(cli)
   call param_register(cli, "verbosity", "NORMAL", verbosity_str, help_string="verbosity level")
   call param_register(cli, "n_col", "1", n_col, help_string="number of columns of data to be analyzed")
   if (.not. param_read_args(cli, ignore_unknown=.false., task="mean_var_decorrelated_err CLI arguments")) then
      call print("Usage: mean_var_decorrelated_err [ verbosity=LEVEL ] [ n_col=n(1) ]", PRINT_ALWAYS)
      call system_abort("Failed to parse command line arguments")
   end if
   call finalise(cli)
   call verbosity_push(verbosity_of_str(trim(verbosity_str)))

   call initialise(io, "stdin", INPUT)

   call read_file(io, lines, N)
   allocate(A(N,n_col))
   do i=1, N
      read (lines(i), *) A(i,:)
   end do
   do i=1, n_col
      call mean_var_decorrelated_err(A(:,i), m, v, decorrelated_error, decorrelation_t)
      call print("col " // i // " mean " // m // " variance " // v // " decorrelation_t " // decorrelation_t // " decorrelated_error " // decorrelated_error)
   end do

   call verbosity_pop()
   call system_finalise()

end program
