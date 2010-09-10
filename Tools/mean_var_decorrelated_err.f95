program do_statistics_decorrelated_err
use libatoms_module, only : dp, system_initialise, system_finalise, verbosity_push, verbosity_pop, initialise, read_file, print, operator(//)
use libatoms_module, only : PRINT_SILENT, PRINT_NORMAL, INPUT
use libatoms_module, only : inoutput
use statistics_module, only : mean_var_decorrelated_err
implicit none
   integer :: i, N
   real(dp), allocatable :: A(:)
   real(dp) :: m, v, decorrelation_t, decorrelated_error
   character(len=100), allocatable :: lines(:)
   type(inoutput) :: io

   call system_initialise(verbosity=PRINT_SILENT)
   call verbosity_push(PRINT_NORMAL)
   call initialise(io, "stdin", INPUT)

   call read_file(io, lines, N)
   allocate(A(N))
   do i=1, N
      read (lines(i), *) A(i)
   end do
   call mean_var_decorrelated_err(A, m, v, decorrelated_error, decorrelation_t)

   call print("mean " // m // " variance " // v // " decorrelation_t " // decorrelation_t // " decorrelated_error " // decorrelated_error)

   call verbosity_pop()
   call system_finalise()

end program
