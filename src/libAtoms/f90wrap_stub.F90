! stub file to include dummy defition of Python error abort function when
! building standalone executable
!
! This stub provides a no-op implementation for standalone programs.
! When building quippy, the real implementation from f90wrap will be used instead.
subroutine f90wrap_abort(message)
  character(*) :: message

  ! do nothing - stub for standalone programs
  ! The quippy_running() guard in error.F90 ensures this is never actually called
  ! when quippy is loaded, so the real f90wrap implementation takes precedence

end subroutine f90wrap_abort
