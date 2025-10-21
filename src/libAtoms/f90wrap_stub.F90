! stub file to include dummy defition of Python error abort function when
! building standalone exectuble
!
! This stub is marked as WEAK so it can be overridden by the real
! implementation from f90wrap when building quippy
subroutine f90wrap_abort(message)
  character(*) :: message
  !GCC$ ATTRIBUTES weak :: f90wrap_abort

  ! do nothing - this stub should be overridden by real implementation

end subroutine f90wrap_abort
