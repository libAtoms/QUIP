module minmod

  use minimization_module, only: my_minim => minim

contains


subroutine minim(nsteps, x, f, df, method, conv_tol, max_steps, &
     linminroutine, hook, hook_print_interval, eps_guess, &
     always_do_test_gradient, optdata)


  implicit none
  real*8, intent(inout),dimension(:) :: x
#ifdef PYF_GEN
  external f
  external df
  external hook
#endif
  character(*) method
  real*8  conv_tol
  integer max_steps
  character(*), optional :: linminroutine
  optional :: hook
  integer, optional :: hook_print_interval
  real*8, optional ::  eps_guess
  logical, optional :: always_do_test_gradient
  character(1), optional, intent(inout), dimension(:) :: optdata
  integer, intent(out) :: nsteps
  INTERFACE 
     subroutine hook(x,dx,E,done,do_print,data)
       use system_module
       real*8::x(:)
       real*8::dx(:)
       real*8::E
       logical :: done
       logical, optional:: do_print
       character(1),optional::data(:)
     end subroutine hook

     function f(x,data)
       use system_module
       real*8::x(:)
       character,optional::data(:)
       real*8::f
     end function f

     function df(x,data)
       use system_module
       real*8::x(:)
       character,optional::data(:)
       real*8::df(size(x))
     end function df
  END INTERFACE

        
#ifdef PYF_GEN
  real*8 :: e
  real*8, dimension(*) :: dx
  logical :: done, do_print, 
  
  e = f(x,optdata)
  dx = df(x,optdata)
  call hook(x,dx,e,done,do_print,optdata)
#endif
  
  nsteps = my_minim(x,f,df,method,conv_tol,max_steps,linminroutine, &
       hook, hook_print_interval, eps_guess, always_do_test_gradient, optdata)
  
end subroutine minim
      
end module minmod
