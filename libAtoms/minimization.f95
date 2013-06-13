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
!X  Minimization module
!X  
!%  This module contains subroutines to perform conjugate gradient and 
!%  damped MD minimisation of an objective function.
!%  The conjugate gradient minimiser can use various different line minimisation routines.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#include "error.inc"
module minimization_module

  use error_module
  use system_module
  use linearalgebra_module
  implicit none
  private
  SAVE

  public :: minim, n_minim, fire_minim, test_gradient, n_test_gradient

  real(dp),parameter:: DXLIM=huge(1.0_dp)     !% Maximum amount we are willing to move any component in a linmin step

  ! parameters from Numerical Recipes */
  real(dp),parameter:: GOLD=1.618034_dp
  integer, parameter:: MAXFIT=5
  real(dp),parameter:: TOL=1e-2_dp
  real(dp),parameter:: CGOLD=0.3819660_dp
  real(dp),parameter:: ZEPS=1e-10_dp
  integer, parameter:: ITER=50

  interface minim
     module procedure minim
  end interface

  interface test_gradient
     module procedure test_gradient
  end interface

  interface n_test_gradient
     module procedure n_test_gradient
  end interface

  ! LBFGS stuff
  external QUIPLB2
  integer::MP,LP
  real(dp)::GTOL,STPMIN,STPMAX
  common /quiplb3/MP,LP,GTOL,STPMIN,STPMAX

CONTAINS


  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  !% Smart line minimiser, adapted from Numerical Recipes.
  !% The objective function is 'func'.
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  function linmin(x0,xdir,y,epsilon,func,data)
    real(dp) :: x0(:)  !% Starting position
    real(dp) :: xdir(:)!% Direction of gradient at 'x0'
    real(dp) :: y(:)   !% Finishing position returned in 'y' after 'linmin'
    real(dp)::epsilon  !% Initial step size
    INTERFACE 
       function func(x,data)
         use system_module
         real(dp)::x(:)
	 character(len=1),optional::data(:)
         real(dp)::func
       end function func
    END INTERFACE
    character(len=1),optional::data(:)
    integer:: linmin

    integer::i,it,sizeflag,N, bracket_it
    real(dp)::Ea,Eb,Ec,a,b,c,r,q,u,ulim,Eu,fallback
    real(dp)::v,w,x,Ex,Ev,Ew,e,d,xm,tol1,tol2,p,etmp

    ! Dynamically allocate to avoid stack overflow madness with ifort
    real(dp), dimension(:), allocatable :: tmpa, tmpb, tmpc, tmpu

    !%RV Number of linmin steps taken, or zero if an error occured

    N=size(x0)

    allocate(tmpa(N), tmpb(N), tmpc(N), tmpu(N))
    tmpa=x0
    y=x0

#ifndef _OPENMP
    call verbosity_push_decrement()
#endif
    Ea=func(tmpa,data)
#ifndef _OPENMP
    call verbosity_pop()
#endif
    call print("  Linmin: Ea = " // Ea // " a = " // 0.0_dp, PRINT_NORMAL)
    Eb= Ea + 1.0_dp !just to start us off

    a=0.0_dp
    b=2.0_dp*epsilon


    ! lets figure out if we can go downhill at all 

    it = 2
    do while( (Eb.GT.Ea) .AND. (.NOT.(Eb.FEQ.Ea)))

       b=b*0.5_dp
       tmpb(:)= x0(:)+b*xdir(:)

#ifndef _OPENMP
       call verbosity_push_decrement()
#endif
       Eb=func(tmpb,data)
#ifndef _OPENMP
       call verbosity_pop()
#endif
       call print("  Linmin: Eb = " // Eb // " b = " // b, PRINT_VERBOSE)
       it = it+1

       if(b.LT.1.0e-20) then

          write(line,*) "  Linmin: Direction points the wrong way\n" ; call print(line, PRINT_NORMAL)

          epsilon=0.0_dp
          linmin=0
          return
       end if
    end do

    ! does it work in fortran?
    !   two:    if(isnan(Eb)) then

    !      write(global_err%unit,*) "linmin: got  a NaN!"
    !      epsilon=0
    !      linmin=0
    !      return
    !   end if two

    if(Eb.FEQ.Ea) then

       epsilon=b
       linmin=1
       call print("  Linmin: Eb.feq.Ea, returning after one step", PRINT_VERBOSE)
       return
    end if

    ! we now have Ea > Eb */

    fallback = b ! b is the best point so far, make that the fallback

    write(line,*)  "  Linmin: xdir is ok.."; call print(line,PRINT_VERBOSE)

    c = b + GOLD*b   !first guess for c */
    tmpc = x0 + c*xdir
#ifndef _OPENMP
    call verbosity_push_decrement()
#endif
    Ec = func(tmpc,data)
#ifndef _OPENMP
    call verbosity_pop()
#endif
    call print("  Linmin: Ec = " // Ec // " c = " // c, PRINT_VERBOSE)
    it = it + 1
    !    ! does it work in fortran?
    !    four:   if(isnan(Ec)) then

    !       write(global_err%unit,*) "linmin: Ec is  a NaN!"
    !       epsilon=0
    !       linmin=0
    !       return
    !    end if four





    !let's bracket the minimum 

    do while(Eb.GT.Ec)

       write(line,*) a,Ea; call print(line,PRINT_VERBOSE)
       write(line,*) b,Eb; call print(line,PRINT_VERBOSE)
       write(line,*) c,Ec; call print(line,PRINT_VERBOSE)




       ! compute u by quadratic fit to a, b, c 


       !inverted ?????????????????????
       r = (b-a)*(Eb-Ec)
       q = (b-c)*(Eb-Ea)
       u = b-((b-c)*q-(b-a)*r)/(2.0_dp*max(abs(q-r), 1.0e-20_dp)*sign(q-r))

       ulim = b+MAXFIT*(c-b)

       write(line,*) "u= ",u ; call print(line,PRINT_VERBOSE)


       if((u-b)*(c-u).GT. 0) then ! b < u < c

          write(line,*)"b < u < c" ; call print(line,PRINT_VERBOSE)
          
          tmpu = x0 + u*xdir
#ifndef _OPENMP
          call verbosity_push_decrement()
#endif
          Eu = func(tmpu,data)
#ifndef _OPENMP
          call verbosity_pop()
#endif
	  call print("linmin got one Eu " // Eu // " " // u, PRINT_NERD)
          it = it + 1

          if(Eu .LT. Ec) then ! Eb > Eu < Ec 
             a = b
             b = u
             Ea = Eb
             Eb = Eu
             exit !break?
          else if(Eu .GT. Eb) then  ! Ea > Eb < Eu 
             c = u
             Ec = Eu
             exit
          end if

          !no minimum found yet

          u = c + GOLD*(c-b)
          tmpu = x0 + u*xdir
#ifndef _OPENMP
          call verbosity_push_decrement()
#endif
          Eu = func(tmpu,data)
#ifndef _OPENMP
          call verbosity_pop()
#endif
	  call print("linmin got second Eu " // Eu // " " // u, PRINT_NERD)
          it = it + 1

       else if((u-c)*(ulim-u) .GT. 0) then ! c < u < ulim

          write(line,*) "  c < u < ulim= ", ulim; call print(line,PRINT_VERBOSE)
          tmpu = x0 + u*xdir
#ifndef _OPENMP
          call verbosity_push_decrement()
#endif
          Eu = func(tmpu,data)
#ifndef _OPENMP
          call verbosity_pop()
#endif
	  call print("linmin got one(2) Eu " // Eu // " " // u, PRINT_NERD)
          it = it + 1

          if(Eu .LT. Ec) then
             b = c
             c = u
             u = c + GOLD*(c-b)
             Eb = Ec
             Ec = Eu
             tmpu = x0 + u*xdir
#ifndef _OPENMP
             call verbosity_push_decrement()
#endif
             Eu = func(tmpu,data)
#ifndef _OPENMP
             call verbosity_pop()
#endif
	     call print("linmin got second(2) Eu " // Eu // " " // u, PRINT_NERD)
             it = it + 1
          end if

       else ! c < ulim < u or u is garbage (we are in a linear regime)
          write(line,*) "  ulim=",ulim," < u or u garbage"; call print(line,PRINT_VERBOSE)
          u = ulim
          tmpu = x0 + u*xdir
#ifndef _OPENMP
          call verbosity_push_decrement()
#endif
          Eu = func(tmpu,data)
#ifndef _OPENMP
          call verbosity_pop()
#endif
          it = it + 1
	  call print("linmin got one(3) Eu " // Eu // " " // u, PRINT_NERD)
       end if

       write(line,*) "  "; call print(line,PRINT_VERBOSE)
       write(line,*) "  "; call print(line,PRINT_VERBOSE)

       ! test to see if we change any component too much



       do i = 1,N

          if(abs(c*xdir(i)) .GT. DXLIM)then
             tmpb = x0+b*xdir
#ifndef _OPENMP
             call verbosity_push_decrement()
#endif
             Eb = func(tmpb,data)
#ifndef _OPENMP
             call verbosity_pop()
#endif
	     call print("linmin got new Eb " // Eb // " " // b, PRINT_NERD)
             it = it + 1
             y = tmpb
             write(line,*) " bracket: step too big", b, Eb
             call print(line, PRINT_VERBOSE)
             write(line,'("I= ",I4," C= ",F16.12," xdir(i)= ",F16.12," DXLIM =",F16.12)')&
                  i, c, xdir(i), DXLIM 
             call print(line, PRINT_VERBOSE)

             epsilon = b

             linmin= it
             return
          end if
       end do


       a = b
       b = c 
       c = u
       Ea = Eb
       Eb = Ec
       Ec = Eu

    end do

    epsilon = b
    fallback = b 

    bracket_it = it
    call print("  Linmin: bracket OK in "//bracket_it//" steps", PRINT_VERBOSE)



    ! ahhh.... now we have a minimum between a and c,  Ea > Eb < Ec 

    write(line,*) " "; call print(line,PRINT_VERBOSE)
    write(line,*) a,Ea; call print(line,PRINT_VERBOSE)
    write(line,*) b,Eb; call print(line,PRINT_VERBOSE)
    write(line,*) c,Ec; call print(line,PRINT_VERBOSE)
    write(line,*) " "; call print(line,PRINT_VERBOSE)



    !********************************************************************
    !   * primitive linmin 
    !   * do a quadratic fit to a, b, c 
    !   *
    !  
    !  r = (b-a)*(Eb-Ec);
    !  q = (b-c)*(Eb-Ea);
    !  u = b-((b-c)*q-(b-a)*r)/(2*Mmax(fabs(q-r), 1e-20)*sign(q-r));
    !  y = x0 + u*xdir;
    !  Eu = (*func)(y);
    !
    !  if(Eb < Eu){ // quadratic fit was bad
    !    if(current_verbosity() > MUMBLE) logger("      Quadratic fit was bad, returning 'b'\n");
    !    u = b;
    !    y = x0 + u*xdir;
    !    Eu = (*func)(y);
    !  }
    !  
    !  if(current_verbosity() > PRINT_SILENT)
    !    logger("      simple quadratic fit: %25.16e%25.16e\n\n", u, Eu);
    !
    !  //return(u);
    !   *
    !   * end primitive linmin
    !**********************************************************************/

    ! now we need a<b as the two endpoints
    b=c
    Eb=Ec

    v = b; w = b; x = b
    Ex = Eb; Ev = Eb; Ew = Eb
    e = 0.0_dp
    d = 0.0_dp
    sizeflag = 0

    call print("linmin got bracket", PRINT_NERD)

    ! main loop for parabolic search
    DO it = 1, ITER

       xm = 0.5_dp*(a+b)
       tol1 = TOL*abs(x)+ZEPS
       tol2 = 2.0_dp*tol1

       ! are we done?
       if((abs(x-xm) < tol2-0.5_dp*(b-a)) .OR.(sizeflag.gt.0)) then
          tmpa = x0 + x*xdir
          y = tmpa
#ifndef _OPENMP
          call verbosity_push_decrement()
#endif
          Ex = func(y,data)
#ifndef _OPENMP
          call verbosity_pop()
#endif
	  call print("linmin got new Ex " // Ex // " " // x, PRINT_NERD)

          if(sizeflag.gt.0) call print('Linmin: DXlim exceeded')

          if( Ex .GT. Ea )then ! emergency measure, linmin screwed up, use fallback 
             call print('  Linmin screwed up! current Ex= '//Ex//' at x= '//x//' is worse than bracket') 
             call print('  Linmin variables: a='//a//' b='//b//' Ev='//Ev//' v='//v//' Ew='//Ew//' Eu='//Eu//' u='//u); 
             x = fallback
             epsilon = fallback
             tmpa = x0 + x*xdir
#ifndef _OPENMP
             call verbosity_push_decrement()
#endif
             Ex = func(tmpa,data)
#ifndef _OPENMP
             call verbosity_pop()
#endif
	     call print("linmin got new Ex " // Ex // " " // x, PRINT_NERD)
             y = tmpa
          else
             epsilon = x
          end if

          call print("  Linmin done "//bracket_it//" bracket and "//it//" steps Ex= "//Ex//" x= "//x)
          linmin=it
          return 
       endif

       ! try parabolic fit on subsequent steps
       if(abs(e) > tol1) then
          r = (x-w)*(Ex-Ev)
          q = (x-v)*(Ex-Ew)
          p = (x-v)*q-(x-w)*r
          q = 2.0_dp*(q-r)

          if(q > 0.0_dp) p = -p

          q = abs(q)
          etmp = e
          e = d

          ! is the parabolic fit acceptable?
          if(abs(p) .GE. abs(0.5_dp*q*etmp) .OR.&
               p .LE. q*(a-x) .OR. p .GE. q*(b-x)) then
             ! no, take a golden section step
             call print('  Linmin: Taking Golden Section step', PRINT_VERBOSE+1)
             if(x .GE. xm) then
                e = a-x
             else
                e = b-x
             end if
             d = CGOLD*e
          else
             call print('  Linmin: Taking parabolic step', PRINT_VERBOSE+1)
             ! yes, take the parabolic step
             d = p/q
             u = x+d
             if(u-a < tol2 .OR. b-u < tol2) d = tol1 * sign(xm-x)
          end if

       else 
          ! on the first pass, do golden section
          if(x .GE. xm) then
             e = a-x
          else
             e = b-x
          end if
          d = CGOLD*e
       end if

       ! construct new step
       if(abs(d) > tol1) then
          u = x+d
       else
          u = x + tol1*sign(d)
       end if

       ! evaluate function
       tmpa = u*xdir
#ifndef _OPENMP
       call verbosity_push_decrement()
#endif
       Eu = func(x0+tmpa,data)
#ifndef _OPENMP
       call verbosity_pop()
#endif
       call print('  Linmin: new point u= '//u//' Eu= '//Eu, PRINT_VERBOSE)
       if(any(abs(tmpa) > DXLIM)) then
          if(sizeflag .EQ. 0) then
             call print('  Linmin: an element of x moved more than '//DXLIM)
          end if
          sizeflag = 1
       endif

       ! analyse new point result
       if(Eu .LE. Ex) then
          call print('  Linmin: new point best so far', PRINT_VERBOSE+1)
          if(u .GE. x) then
             a = x
          else
             b = x
          end if

          v=w;w=x;x=u
          Ev=Ew;Ew=Ex;Ex=Eu
       else 
          call print('  Linmin: new point is no better', PRINT_VERBOSE+1)
          if(u < x) then
             a = u
          else
             b = u
          endif

          if(Eu .LE. Ew .OR. w .EQ. x) then
             v=w;w=u      
             Ev = Ew; Ew = Eu
          else if(Eu .LE. Ev .OR. v .EQ. x .OR. v .EQ. w) then
             v = u
             Ev = Eu
          end if
       end if
       call print('  Linmin: a='//a//' b='//b//' x='//x, PRINT_VERBOSE+1)

    end do

    write(line,*) "  Linmin: too many iterations"; call print(line, PRINT_NORMAL)



    y = x0 + x*xdir
    epsilon = x
    linmin = it

    deallocate(tmpa, tmpb, tmpc, tmpu)
    return
  end function linmin

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  !% Simple (fast, dumb) line minimiser. The 'func' interface is
  !% the same as for 'linmin' above.
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


  function linmin_fast(x0,fx0,xdir,y,epsilon,func,data) result(linmin)
    real(dp)::x0(:)   !% Starting vector
    real(dp)::fx0     !% Value of 'func' at 'x0'
    real(dp)::xdir(:) !% Direction
    real(dp)::y(:)    !%  Result
    real(dp)::epsilon !% Initial step            
    INTERFACE 
       function func(x,data)
         use system_module
         real(dp)::x(:)
	 character(len=1),optional::data(:)
         real(dp)::func
       end function func
    END INTERFACE
    character(len=1),optional::data(:)
    integer::linmin  

    ! Dynamically allocate to avoid stack overflow madness with ifort
    real(dp),allocatable::xb(:),xc(:)
    real(dp)::r,q,new_eps,a,b,c,fxb,fxc,ftol
    integer::n
    logical::reject_quadratic_extrap

    !%RV Number of linmin steps taken, or zero if an error occured

    N=size(x0)
    new_eps=-1.0_dp
    a=0.0_dp
    ftol=1.e-13
    allocate(xb(N))
    allocate(xc(N))

    call print("Welcome to linmin_fast", PRINT_NORMAL)

    reject_quadratic_extrap = .true.
    do while(reject_quadratic_extrap)

       b = a+epsilon 
       c = b+GOLD*epsilon
       xb = x0+b*xdir
       xc = x0+c*xdir

#ifndef _OPENMP
       call verbosity_push_decrement()
#endif
       fxb = func(xb,data)
       fxc = func(xc,data)
#ifndef _OPENMP
       call verbosity_pop()
#endif

       write(line,*) "      abc = ", a, b, c; call print(line, PRINT_NORMAL)
       write(line,*) "      f   = ",fx0, fxb, fxc; call print(line, PRINT_NORMAL)

       if(abs(fx0-fxb) < abs(ftol*fx0))then
          write(line,*) "*** fast_linmin is stuck, returning 0" 
          call print(line,PRINT_SILENT)
          
          linmin=0
          return 
       end if

       !WARNING: find a way to do the following

       ! probably we will need to wrap isnan

       !if(isnan(fxb) .OR. isnan(fxc))then
       !   write(global_err%unit,*)"    Got a NaN!!!"
       !   linmin=0
       !   return 
       !end if
       !  if(.NOT.finite(fxb) .OR. .NOT.finite(fxc))then
       !     write(global_err%unit,*)"    Got a INF!!!"
       !     linmin=0
       !     return 
       !  end if

       r = (b-a)*(fxb-fxc)
       q = (b-c)*(fxb-fx0)
       new_eps = b-((b-c)*q-(b-a)*r)/(2.0*max(abs(q-r), 1.0_dp-20)*sign(1.0_dp,q-r))

       write(line,*) "      neweps = ", new_eps; call print(line, PRINT_NORMAL)

       if (abs(fxb) .GT. 100.0_dp*abs(fx0) .OR. abs(fxc) .GT. 100.0_dp*abs(fx0)) then ! extrapolation gave stupid results
          epsilon = epsilon/10.0_dp
          reject_quadratic_extrap = .true.
       else if(new_eps > 10.0_dp*(b-a))then
          write(line,*) "*** new epsilon > 10.0 * old epsilon, capping at factor of 10.0 increase"
          call print(line, PRINT_NORMAL)
          epsilon = 10.0_dp*(b-a)
          a = c
          fx0 = fxc
          reject_quadratic_extrap = .true.
       else if(new_eps <  0.0_dp) then
          ! new proposed minimum is behind us!
          if(fxb < fx0 .and. fxc < fx0 .and. fxc < fxb) then
             ! increase epsilon
             epsilon = epsilon*2.0_dp
             call print("*** quadratic extrapolation resulted in backward step, increasing epsilon: "//epsilon, PRINT_NORMAL)
          else
             epsilon = epsilon/10.0_dp
             write(line,*)"*** xdir wrong way, reducing epsilon ",epsilon
             call print(line, PRINT_NORMAL)
          endif
          reject_quadratic_extrap = .true.
       else if(new_eps < epsilon/10.0_dp) then
          write(line,*) "*** new epsilon < old epsilon / 10, reducing epsilon by a factor of 2"
          epsilon = epsilon/2.0_dp
          reject_quadratic_extrap = .true.
       else
          reject_quadratic_extrap = .false.
       end if
    end do
    y = x0+new_eps*xdir

    epsilon = new_eps
    linmin=1

    deallocate(xb, xc)

    return 
  end function linmin_fast


  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  !% Line minimizer that uses the derivative and extrapolates
  !% its projection onto the search direction to zero. Again,
  !% the 'func' interface is the same as for 'linmin'.
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  function linmin_deriv(x0, xdir, dx0, y, epsilon, dfunc, data) result(linmin)
    real(dp)::x0(:)  !% Starting vector
    real(dp)::xdir(:)!% Search direction
    real(dp)::dx0(:) !% Initial gradient
    real(dp)::y(:)   !% Result
    real(dp)::epsilon!% Initial step            
    INTERFACE 
       function dfunc(x,data)
         use system_module
         real(dp)::x(:)
	 character(len=1),optional::data(:)
         real(dp)::dfunc(size(x))
       end function dfunc
    END INTERFACE
    character(len=1),optional::data(:)
    integer::linmin  

    ! local vars
    integer::N
    real(dp), allocatable::x1(:), dx1(:)
    real(dp)::dirdx0, dirdx1, gamma, new_eps

    !%RV Number of linmin steps taken, or zero if an error occured

    N=size(x0)
    dirdx1 = 9.9e30_dp
    new_eps = epsilon
    allocate(x1(N), dx1(N))

    linmin = 0
    dirdx0 = xdir .DOT. dx0
    if( dirdx0 > 0.0_dp) then ! should be negative
       call print("WARNING: linmin_deriv: xdir*dx0 > 0 !!!!!", PRINT_ALWAYS)
       return
    endif

    do while(abs(dirdx1) > abs(dirdx0))
       x1 = x0+new_eps*xdir
#ifndef _OPENMP
       call verbosity_push_decrement()
#endif
       dx1 = dfunc(x1,data)
#ifndef _OPENMP
       call verbosity_pop()
#endif
       dirdx1 = xdir .DOT. dx1

       if(abs(dirdx1) > abs(dirdx0)) then ! this eps leads to a point with larger abs gradient
          call print("WARNING: linmin_deriv: |dirdx1| > |dirdx0|, reducing epsilon by factor of 5", PRINT_NORMAL)
          new_eps = new_eps/5.0_dp
          if(new_eps < 1.0e-12_dp) then
             call print("WARNING: linmin_deriv: new_eps < 1e-12 !!!!!", PRINT_NORMAL)
             linmin = 0
             return
          endif
       else
          gamma = dirdx0/(dirdx0-dirdx1)
          if(gamma > 2.0_dp) then
             call print("*** gamma > 2.0, capping at 2.0", PRINT_NORMAL)
             gamma = 2.0_dp
          endif
          new_eps = gamma*epsilon
       endif
       linmin = linmin+1
    end do

    write (line,'(a,e10.2,a,e10.2)') '  gamma = ', gamma, ' new_eps = ', new_eps
    call Print(line, PRINT_NORMAL)

    y = x0+new_eps*xdir
    epsilon = new_eps

    deallocate(x1, dx1)
    return 
  end function linmin_deriv


  !% Iterative version of 'linmin_deriv' that avoid taking large steps
  !% by iterating the extrapolation to zero gradient.
  function linmin_deriv_iter(x0, xdir, dx0, y, epsilon, dfunc,data,do_line_scan) result(linmin)
    real(dp)::x0(:)  !% Starting vector
    real(dp)::xdir(:)!% Search direction
    real(dp)::dx0(:) !% Initial gradient
    real(dp)::y(:)   !% Result
    real(dp)::epsilon!% Initial step            
    logical, optional :: do_line_scan

    INTERFACE 
       function dfunc(x,data)
         use system_module
         real(dp)::x(:)
	 character(len=1),optional::data(:)
         real(dp)::dfunc(size(x))
       end function dfunc
    END INTERFACE
    character(len=1),optional :: data(:)
    integer::linmin  

    ! local vars
    integer::N, extrap_steps,i
    real(dp), allocatable::xn(:), dxn(:)
    real(dp)::dirdx1, dirdx2, dirdx_new,  eps1, eps2, new_eps, eps11, step, old_eps
    integer, parameter :: max_extrap_steps = 50
    logical :: extrap

    !%RV Number of linmin steps taken, or zero if an error occured

    N=size(x0)
    allocate(xn(N), dxn(N))

    if (present(do_line_scan)) then
       if (do_line_scan) then
          call print('line scan:', PRINT_NORMAL)
          new_eps = 1.0e-5_dp
          do i=1,50
             xn = x0 + new_eps*xdir
#ifndef _OPENMP
             call verbosity_push_decrement()
#endif
             dxn = dfunc(xn,data)
#ifndef _OPENMP
             call verbosity_pop()
#endif
             dirdx_new = xdir .DOT. dxn
          
             call print(new_eps//' '//dirdx_new//' <-- LS', PRINT_NORMAL)
       
             new_eps = new_eps*1.15
          enddo
       end if
    end if

    eps1 = 0.0_dp
    eps2 = epsilon
    new_eps = epsilon

    linmin = 0
    dirdx1 = xdir .DOT. dx0
    dirdx2 = dirdx1
    if( dirdx1 > 0.0_dp) then ! should be negative
       call print("WARNING: linmin_deriv_iter: xdir*dx0 > 0 !!!!!", PRINT_ALWAYS)
       return
    endif

    extrap_steps = 0

    do while ( (abs(eps1-eps2) > TOL*abs(eps1)) .and. extrap_steps < max_extrap_steps) 
       do
          xn = x0 + new_eps*xdir
#ifndef _OPENMP
          call verbosity_push_decrement()
#endif
          dxn = dfunc(xn,data)
#ifndef _OPENMP
          call verbosity_pop()
#endif
          dirdx_new = xdir .DOT. dxn

          linmin = linmin + 1

          call print('eps1   = '//eps1//' eps2   = '//eps2//' new_eps   = '//new_eps, PRINT_NORMAL)
          call print('dirdx1 = '//dirdx1//' dirdx2 = '//dirdx2//' dirdx_new = '//dirdx_new, PRINT_NORMAL)

          extrap = .false.
          if (dirdx_new > 0.0_dp) then 
             if(abs(dirdx_new) < 2.0_dp*abs(dirdx1)) then
                ! projected gradient at new point +ve, but not too large. 
                call print("dirdx_new > 0, but not too large", PRINT_NORMAL)
                eps2 = new_eps
                dirdx2 = dirdx_new
                extrap_steps = 0
                ! we're straddling the minimum well, so gamma < 2 and we can interpolate to the next step 

                ! let's try to bring in eps1 
                step = 0.5_dp*(new_eps-eps1)
                dirdx1 = 1.0_dp
                do while (dirdx1 > 0.0_dp)
                   eps11 = eps1+step
                   call print("Trying to bring in eps1: "//eps11, PRINT_NORMAL)
                   xn = x0 + eps11*xdir
#ifndef _OPENMP
                   call verbosity_push_decrement()
#endif
                   dxn = dfunc(xn,data)
#ifndef _OPENMP
                   call verbosity_pop()
#endif
                   dirdx1 = xdir .DOT. dxn
                   step = step/2.0_dp
                end do
                eps1 = eps11

                exit
             else
                ! projected gradient at new point +ve and large
                ! let's decrease the step we take
                call print("*** reducing trial epsilon by factor of 2, making eps2=current new_eps", PRINT_NORMAL)
                eps2 = new_eps
                dirdx2 = dirdx_new
                new_eps = (eps1+new_eps)/2.0_dp
                ! check if we are just fiddling around
                if(new_eps < 1.0e-12_dp) then
                   call print("WARNING: linmin_deriv_iter: total_eps < 1e-12 !!!!!", PRINT_NORMAL)
                   linmin = 0
                   return
                endif
             end if
          else ! projected gradient is -ve
             if(abs(dirdx_new) <= abs(dirdx1)) then
                ! projected gradient smaller than  at x1
                if(dirdx2 > 0.0_dp) then
                   ! we have good bracketing, so interpolate
                   call print("dirdx_new <= 0, and we have good bracketing", PRINT_NORMAL)
                   eps1 = new_eps
                   dirdx1 = dirdx_new
                   extrap_steps = 0

                   ! let's try to bring in eps2
                   step = 0.5_dp*(eps2-new_eps)
                   dirdx2 = -1.0_dp
                   do while(dirdx2 < 0.0_dp)
                      eps11 = eps2-step
                      call print("Trying to bring in eps2: "//eps11, PRINT_NORMAL)
                      xn = x0 + eps11*xdir
#ifndef _OPENMP
                      call verbosity_push_decrement()
#endif
                      dxn = dfunc(xn,data)
#ifndef _OPENMP
                      call verbosity_pop()
#endif
                      dirdx2 = xdir .DOT. dxn
                      step = step/2.0_dp
                   end do
                   eps2 = eps11

                   exit
                else
                   ! we have not bracketed yet, but can take this point and extrapolate to a bigger stepsize
                   old_eps = new_eps
                   if(abs(dirdx_new-dirdx1) .fne. 0.0_dp) new_eps = eps1-dirdx1/(dirdx_new-dirdx1)*(new_eps-eps1)
                   call print("we have not bracketed yet, extrapolating: "//new_eps, PRINT_NORMAL)
                   extrap_steps = extrap_steps + 1
                   if(new_eps > 5.0_dp*old_eps) then
                      ! extrapolation is too large, let's just move closer
                      call print("capping extrapolation at "//2.0_dp*old_eps, PRINT_NORMAL)
                      eps1 = old_eps
                      new_eps = 2.0_dp*old_eps
                   else
                      ! accept the extrapolation
                      eps1 = old_eps
                      dirdx1 = dirdx_new
                   end if
                endif
             else
                if (dirdx2 < 0.0_dp) then
                   ! have not bracketed yet, and projected gradient too big - minimum is behind us! lets move forward
                   call print("dirdx2 < 0.0_dp and projected gradient too big, closest stationary point is behind us!", PRINT_NORMAL)
                   eps1 = new_eps
                   dirdx1 = dirdx_new
                   new_eps = eps1*2.0_dp
                   eps2 = new_eps
                   extrap_steps = extrap_steps+1
                else
                   call print("abs(dirdx_new) > abs(dirdx1) but dirdx2 > 0, should only happen when new_eps is converged. try to bring in eps2", PRINT_NORMAL)
                   eps2 = 0.5_dp*(new_eps+eps2)
                   xn = x0 + eps2*xdir
#ifndef _OPENMP
                   call verbosity_push_decrement()
#endif
                   dxn = dfunc(xn,data)
#ifndef _OPENMP
                   call verbosity_pop()
#endif
                   dirdx2 = xdir .DOT. dxn
                   exit
                endif
             endif
          end if
       end do
       new_eps = eps1 - dirdx1/(dirdx2-dirdx1)*(eps2-eps1)


    end do

    if (extrap_steps == max_extrap_steps) then
       call Print('*** linmin_deriv: max consequtive extrapolation steps exceeded', PRINT_ALWAYS)
       linmin = 0
       return
    end if

    call print('linmin_deriv_iter done in '//linmin//' steps')

    epsilon = new_eps
    y = x0 + epsilon*xdir

    deallocate(xn, dxn)
    return 
  end function linmin_deriv_iter


  !% Simplified Iterative version of 'linmin_deriv' that avoid taking large steps
  !% by iterating the extrapolation to zero gradient. It does not try to reduce
  !% the bracketing interval on both sides at every step
  function linmin_deriv_iter_simple(x0, xdir, dx0, y, epsilon, dfunc,data,do_line_scan) result(linmin)
    real(dp)::x0(:)  !% Starting vector
    real(dp)::xdir(:)!% Search direction
    real(dp)::dx0(:) !% Initial gradient
    real(dp)::y(:)   !% Result
    real(dp)::epsilon!% Initial step            
    logical, optional :: do_line_scan
    INTERFACE 
       function dfunc(x,data)
         use system_module
         real(dp)::x(:)
	 character(len=1),optional::data(:)
         real(dp)::dfunc(size(x))
       end function dfunc
    END INTERFACE
    character(len=1),optional :: data(:)
    integer::linmin  

    ! local vars
    integer::N, extrap_steps, i
    real(dp), allocatable::xn(:), dxn(:)
    real(dp)::dirdx1, dirdx2, dirdx_new,  eps1, eps2, new_eps, old_eps 
    integer, parameter :: max_extrap_steps = 50
    logical :: extrap

    !%RV Number of linmin steps taken, or zero if an error occured



    N=size(x0)
    allocate(xn(N), dxn(N))

    if (present(do_line_scan)) then
       if (do_line_scan) then
          call print('line scan:', PRINT_NORMAL)
          new_eps = 1.0e-5_dp
          do i=1,50
             xn = x0 + new_eps*xdir
#ifndef _OPENMP
             call verbosity_push_decrement()
#endif
             dxn = dfunc(xn,data)
#ifndef _OPENMP
             call verbosity_pop()
#endif
             dirdx_new = xdir .DOT. dxn
             
             call print(new_eps//' '//dirdx_new, PRINT_NORMAL)
          
             new_eps = new_eps*1.15
          enddo
       end if
    end if


    eps1 = 0.0_dp
    eps2 = epsilon
    old_eps = 0.0_dp
    new_eps = epsilon

    linmin = 0
    dirdx1 = xdir .DOT. dx0
    dirdx2 = 0.0_dp
    if( dirdx1 > 0.0_dp) then ! should be negative
       call print("WARNING: linmin_deriv_iter_simple: xdir*dx0 > 0 !!!!!", PRINT_ALWAYS)
       return
    endif



    extrap_steps = 0

    do while ( (abs(old_eps-new_eps) > TOL*abs(new_eps)) .and. extrap_steps < max_extrap_steps) 
       xn = x0 + new_eps*xdir
#ifndef _OPENMP
       call verbosity_push_decrement()
#endif
       dxn = dfunc(xn,data)
#ifndef _OPENMP
       call verbosity_pop()
#endif
       dirdx_new = xdir .DOT. dxn
       
       linmin = linmin + 1
       
       call print('eps1   = '//eps1//' eps2   = '//eps2//' new_eps   = '//new_eps, PRINT_NORMAL)
       call print('dirdx1 = '//dirdx1//' dirdx2 = '//dirdx2//' dirdx_new = '//dirdx_new, PRINT_NORMAL)
       extrap = .false.
       if (dirdx_new > 0.0_dp) then 
          if(abs(dirdx_new) < 10.0_dp*abs(dirdx1)) then
             ! projected gradient at new point +ve, but not too large. 
             call print("dirdx_new > 0, but not too large", PRINT_NORMAL)
             eps2 = new_eps
             dirdx2 = dirdx_new
             extrap_steps = 0
             ! we're straddling the minimum well, so gamma < 2 and we can interpolate to the next step 
          else
             ! projected gradient at new point +ve and large
             ! let's decrease the step we take
             call print("*** reducing trial epsilon by factor of 2, making eps2=current new_eps", PRINT_NORMAL)
             eps2 = new_eps
             dirdx2 = dirdx_new
             old_eps = new_eps
             new_eps = (eps1+new_eps)/2.0_dp
             ! check if we are just fiddling around
             if(new_eps < 1.0e-12_dp) then
                call print("WARNING: linmin_deriv_iter_simple: new_eps < 1e-12 !!!!!", PRINT_NORMAL)
                linmin = 0
                return
             endif
             cycle
          end if
       else ! projected gradient is -ve
          if(abs(dirdx_new) <= abs(dirdx1)) then
             ! projected gradient smaller than  at x1
             if(dirdx2 > 0.0_dp) then
                ! we have good bracketing, so interpolate
                call print("dirdx_new <= 0, and we have good bracketing", PRINT_NORMAL)
                eps1 = new_eps
                dirdx1 = dirdx_new
                extrap_steps = 0
             else
                ! we have not bracketed yet, but can take this point and extrapolate to a bigger stepsize
                old_eps = new_eps
                if((dirdx_new-dirdx1) .fne. 0.0_dp) new_eps = eps1-dirdx1/(dirdx_new-dirdx1)*(new_eps-eps1)
                call print("we have not bracketed yet, extrapolating: "//new_eps, PRINT_NORMAL)
                if(new_eps > 5.0_dp*old_eps) then
                   ! extrapolation is too large, let's just move closer
                   call print("capping extrapolation at "//2.0_dp*old_eps, PRINT_NORMAL)
                   eps1 = old_eps
                   new_eps = 2.0_dp*old_eps
                else
                   ! accept the extrapolation
                   eps1 = old_eps
                   dirdx1 = dirdx_new
                end if
                extrap_steps = extrap_steps + 1
                cycle
             endif
          else
             if (dirdx2 .eq. 0.0_dp) then
                ! have not bracketed yet, and projected gradient too big - minimum is behind us! lets move forward
                call print("dirdx2 < 0.0_dp and projected gradient too big, closest stationary point is behind us!", PRINT_NORMAL)
                eps1 = new_eps
                dirdx1 = dirdx_new
                old_eps = new_eps
                new_eps = eps1*2.0_dp
                eps2 = new_eps
                extrap_steps = extrap_steps+1
                cycle
             else
                call print("dirdx_new < 0, abs(dirdx_new) > abs(dirdx1) but dirdx2 > 0, function not monotonic?", PRINT_NORMAL)
                eps1 = new_eps
                dirdx1 = dirdx_new
             endif
          endif
       end if
       old_eps = new_eps
       new_eps = eps1 - dirdx1/(dirdx2-dirdx1)*(eps2-eps1)
       
    end do

    if (extrap_steps == max_extrap_steps) then
       call Print('*** linmin_deriv_iter_simple: max consequtive extrapolation steps exceeded', PRINT_ALWAYS)
       linmin = 0
       return
    end if

    epsilon = new_eps
    y = x0 + epsilon*xdir

    deallocate(xn, dxn)
    return 
  end function linmin_deriv_iter_simple




  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  !% Damped molecular dynamics minimiser. The objective
  !% function is 'func(x)' and its gradient is 'dfunc(x)'.
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX



  function damped_md_minim(x,mass,func,dfunc,damp,tol,max_change,max_steps,data)

    real(dp)::x(:)  !% Starting vector
    real(dp)::mass(:) !% Effective masses of each degree of freedom
    INTERFACE 
       function func(x,data)
         use system_module
         real(dp)::x(:)
	 character(len=1),optional::data(:)
         real(dp)::func
       end function func
    end INTERFACE
    INTERFACE
       function dfunc(x,data)
         use system_module
         real(dp)::x(:)
	 character(len=1),optional::data(:)
         real(dp)::dfunc(size(x))
       end function dfunc
    END INTERFACE
    real(dp)::damp !% Velocity damping factor
    real(dp)::tol  !% Minimisation is taken to be converged when 'normsq(force) < tol'
    real(dp)::max_change !% Maximum position change per time step
    integer:: max_steps  !% Maximum number of MD steps
    character(len=1),optional::data(:)

    integer :: N,i
    integer :: damped_md_minim
    real(dp),allocatable::velo(:),acc(:),force(:)
    real(dp)::dt,df2, f

    !%RV Returns number of MD steps taken.


    N=size(x)
    allocate(velo(N),acc(N),force(N))

    call print("Welcome to damped md minim()")
    write(line,*)"damping = ", damp            ; call print(line)

    velo=0.0_dp
#ifndef _OPENMP
    call verbosity_push_decrement()
#endif
    acc = -dfunc(x,data)/mass
#ifndef _OPENMP
    call verbosity_pop()
#endif


    dt = sqrt(max_change/maxval(abs(acc)))
    write(line,*)"dt = ", dt ; call print(line)

    do I=0,max_steps
       velo(:)=velo(:) + (0.5*dt)*acc(:)
#ifndef _OPENMP
       call verbosity_push_decrement()
#endif
       force(:)= -dfunc(X,data)
#ifndef _OPENMP
       call verbosity_pop()
#endif
       df2=normsq(force)
       if(df2 .LT. tol) then
          write(line,*) I," force^2 =",df2 ; call print(line)
          write(line,*)"Converged in ",i," steps!!" ; call print(line)
          exit
       else if (mod(i,100) .EQ. 0) then
          write(line,*)i," f = ", func(x,data), "df^2 = ", df2, "max(abs(df)) = ",maxval(abs(force)); call print(line)
       end if
       acc(:)=force(:)/mass(:)
       velo(:)=velo(:) + (0.5*dt)*acc(:)

       velo(:)=velo(:) * (1.0-damp)/(1.0+damp)
       x(:)=x(:)+dt*velo(:)
       x(:)=x(:)+0.5*dt*dt*acc(:)

#ifndef _OPENMP
       call verbosity_push_decrement()
#endif
       f = func(x,data)
#ifndef _OPENMP
       call verbosity_pop()
#endif
       call print("f=" // f, PRINT_VERBOSE)
       call print(x, PRINT_NERD)

    end do
    ! 
    if(i .EQ. max_steps) then
       write(line,*) "Failed to converge in ",i," steps" ; call print(line)
    end if
    damped_md_minim = i
    return 
  end function damped_md_minim

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  !% Gradient descent minimizer, using either the conjugate gradients or
  !% steepest descent methods. The objective function is
  !% 'func(x)' and its gradient is 'dfunc(x)'.
  !% There is an additional 'hook' interface which is called at the
  !% beginning of each gradient descent step. 
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX



  function minim(x_in,func,dfunc,method,convergence_tol,max_steps, linminroutine, hook, hook_print_interval,  &
		 eps_guess, always_do_test_gradient, data, status)
    real(dp),     intent(inout) :: x_in(:) !% Starting position
    INTERFACE 
       function func(x,data)
         use system_module
         real(dp)::x(:)
	 character(len=1),optional::data(:)
         real(dp)::func
       end function func
    end INTERFACE
    INTERFACE
       function dfunc(x,data)
         use system_module
         real(dp)::x(:)
	 character(len=1),optional::data(:)
         real(dp)::dfunc(size(x))
       end function dfunc
    END INTERFACE
    character(*), intent(in)    :: method !% 'cg' for conjugate gradients or 'sd' for steepest descent
    real(dp),     intent(in)    :: convergence_tol !% Minimisation is treated as converged once $|\mathbf{\nabla}f|^2 <$
    !% 'convergence_tol'. 
    integer,      intent(in)    :: max_steps  !% Maximum number of 'cg' or 'sd' steps
    integer::minim
    character(*), intent(in), optional :: linminroutine !% Name of the line minisation routine to use. 
    !% This should be one of  'NR_LINMIN', 'FAST_LINMIN' and
    !% 'LINMIN_DERIV'.
    !% If 'FAST_LINMIN' is used and problems with the line
    !% minisation are detected, 'minim' automatically switches
    !% to the more reliable 'NR_LINMIN', and then switches back
    !% once no more problems have occurred for some time.
    !% the default is NR_LINMIN
    optional :: hook
    INTERFACE 
       subroutine hook(x,dx,E,done,do_print,data)
         use system_module
         real(dp), intent(in) ::x(:)
         real(dp), intent(in) ::dx(:)
         real(dp), intent(in) ::E
         logical, intent(out) :: done
	 logical, optional, intent(in) :: do_print
	 character(len=1),optional, intent(in) ::data(:)
       end subroutine hook
    end INTERFACE
    integer, intent(in), optional :: hook_print_interval
    real(dp), intent(in), optional :: eps_guess
    logical, intent(in), optional :: always_do_test_gradient
    character(len=1), optional, intent(inout) :: data(:)
    integer, optional, intent(out) :: status

    !%RV Returns number of gradient descent steps taken during minimisation

    integer, parameter:: max_bad_cg = 5
    integer, parameter:: max_bad_iter = 5
    integer, parameter:: convergence_window = 3
    integer:: convergence_counter
    integer:: main_counter
    integer:: bad_cg_counter
    integer:: bad_iter_counter
    integer:: resetflag
    integer:: exit_flag
    integer:: lsteps
    integer:: i, extra_report
    real(dp), parameter:: stuck_tol = NUMERICAL_ZERO
    real(dp):: linmin_quality
    real(dp):: eps, alpha
    real(dp):: oldeps
    real(dp), parameter :: default_eps_guess = 0.1_dp ! HACK
    real(dp) :: my_eps_guess
    real(dp):: f, f_new
    real(dp):: gg, dgg, hdirgrad_before, hdirgrad_after
    real(dp):: dcosine, gdirlen, gdirlen_old, normsqgrad_f, normsqgrad_f_old
    real(dp):: obj, obj_new
    logical:: do_sd, do_cg, do_pcg, do_lbfgs, do_sd2
    integer:: fast_linmin_switchback
    logical:: do_fast_linmin
    logical:: do_linmin_deriv
    logical:: dumbool, done
    logical:: do_test_gradient
    integer :: my_hook_print_interval

    ! working arrays
    ! Dynamically allocate to avoid stack overflow madness with ifort
    real(dp),dimension(:), allocatable :: x, y, hdir, gdir, gdir_old, grad_f, grad_f_old, x_old
    ! for lbfgs
    real(dp), allocatable :: lbfgs_work(:), lbfgs_diag(:)
    integer :: lbfgs_flag
    integer, parameter :: lbfgs_M = 40

    if (current_verbosity() >= PRINT_VERBOSE) then
      my_hook_print_interval = optional_default(1, hook_print_interval)
    else if (current_verbosity() >= PRINT_NORMAL) then
      my_hook_print_interval = optional_default(10, hook_print_interval)
    else
      my_hook_print_interval = optional_default(100000, hook_print_interval)
    endif

    call system_timer("minim")

    call system_timer("minim/init")
    allocate(x(size(x_in)))
    x = x_in

    allocate(y(size(x)), hdir(size(x)), gdir(size(x)), gdir_old(size(x)), &
         grad_f(size(x)), grad_f_old(size(x)))

    extra_report = 0
    if (present(status)) status = 0

    call print("Welcome to minim()", PRINT_NORMAL)
    call print("space is "//size(x)//" dimensional", PRINT_NORMAL)
    do_sd = .false.
    do_sd2 = .false.
    do_cg = .false.
    do_pcg = .false.
    do_lbfgs = .false.

    if(trim(method).EQ."sd") then
       do_sd = .TRUE.
       call print("Method: Steepest Descent", PRINT_NORMAL)
    else if(trim(method).EQ."sd2")then
       do_sd2 = .TRUE.
       call print("Method: Two-Point Step Size Gradient Methods, J Barzilai and JM Borwein, IMA J Num Anal (1988) 8, 141-148", PRINT_NORMAL)
       allocate(x_old(size(x)))
       y=x
       x_old = x
    else if(trim(method).EQ."cg")then
       do_cg = .TRUE.
       call print("Method: Conjugate Gradients", PRINT_NORMAL)
    else if(trim(method).EQ."pcg")then
       do_cg = .TRUE.
       do_pcg = .TRUE.
       call print("Method: Preconditioned Conjugate Gradients", PRINT_NORMAL)
    else if(trim(method).EQ."lbfgs") then
       do_lbfgs = .TRUE.
       call print("Method: LBFGS by Jorge Nocedal, please cite D. Liu and J. Nocedal, Mathematical Programming B 45 (1989) 503-528", PRINT_NORMAL)
       allocate(lbfgs_diag(size(x)))
       allocate(lbfgs_work(size(x)*(lbfgs_M*2+1)+2*lbfgs_M))
       call print("Allocating LBFGS work array: "//(size(x)*(lbfgs_M*2+1)+2*lbfgs_M)//" bytes")
       lbfgs_flag = 0
       ! set units
       MP = mainlog%unit
       LP = errorlog%unit
       !
       y=x ! this is reqired for lbfgs on first entry into the main loop, as no linmin is done
    else
       call System_abort("Invalid method in optimize: '"//trim(method)//"'")
    end if


    do_fast_linmin = .FALSE. 
    do_linmin_deriv = .FALSE.
    if(present(linminroutine)) then
       if(do_lbfgs) &
          call print("Minim warning: a linminroutine was specified for use with LBFGS")
       if(do_sd2) &
          call print("Minim warning: a linminroutine was specified for use with two-point steepest descent SD2")
       if(trim(linminroutine) .EQ. "FAST_LINMIN") then
          do_fast_linmin =.TRUE.
          call print("Using FAST_LINMIN linmin", PRINT_NORMAL)
       else if(trim(linminroutine).EQ."NR_LINMIN") then
          call print("Using NR_LINMIN linmin", PRINT_NORMAL)
       else if(trim(linminroutine).EQ."LINMIN_DERIV") then
          do_linmin_deriv = .TRUE.
          call print("Using LINMIN_DERIV linmin", PRINT_NORMAL)
       else
          call System_abort("Invalid linminroutine: "//linminroutine)
       end if
    end if

    if (current_verbosity() .GE. PRINT_NERD .and. .not. do_linmin_deriv) then 
       dumbool=test_gradient(x, func, dfunc,data=data)
    end if

    ! initial function calls
    if (.not. do_linmin_deriv)  then
#ifndef _OPENMP
       call verbosity_push_decrement(2)
#endif
       f = func(x,data)
#ifndef _OPENMP
       call verbosity_pop()
#endif
    else
       f = 0.0_dp
    end if

#ifndef _OPENMP
    call verbosity_push_decrement(2)
#endif
    grad_f = dfunc(x,data)
#ifndef _OPENMP
    call verbosity_pop()
#endif

    if (my_hook_print_interval > 0) then
      if (present(hook)) then 
	call hook(x, grad_f, f, done, .true., data)
	if (done) then
	  call print('hook reports that minim finished, exiting.', PRINT_NORMAL)
	  exit_flag = 1
	end if
      else
	call print("hook is not present", PRINT_VERBOSE)
      end if
    endif


    grad_f_old = grad_f
    gdir = (-1.0_dp)*grad_f
    hdir = gdir
    gdir_old = gdir
    normsqgrad_f = normsq(grad_f)
    normsqgrad_f_old = normsqgrad_f
    gdirlen =  sqrt(normsqgrad_f)
    gdirlen_old = gdirlen
    ! quality monitors
    dcosine = 0.0_dp
    linmin_quality = 0.0_dp
    bad_cg_counter = 0
    bad_iter_counter = 0
    convergence_counter = 0
    resetflag = 0
    exit_flag = 0
    fast_linmin_switchback = -1
    lsteps = 0



    do_test_gradient = optional_default(.false., always_do_test_gradient)
    my_eps_guess = optional_default(default_eps_guess, eps_guess)
    eps = my_eps_guess    
    
    !********************************************************************
    !*
    !*  MAIN CG LOOP
    !*
    !**********************************************************************
    
    
    if(normsqgrad_f .LT. convergence_tol)then  
       call print("Minimization is already converged!")
       call print(trim(method)//" iter = "// 0 //" df^2 = " // normsqgrad_f // " f = " // f &
            &// " "//lsteps//" linmin steps eps = "//eps,PRINT_VERBOSE)
       exit_flag = 1
    end if

    call system_timer("minim/init")
    call system_timer("minim/main_loop")
    
    main_counter=1 ! incremented at the end of the loop
    do while((main_counter .LT. max_steps) .AND. (.NOT.(exit_flag.gt.0)))
      call system_timer("minim/main_loop/"//main_counter)
       
       if ((current_verbosity() >= PRINT_ANAL .or. do_test_gradient) & 
            .and. .not. do_linmin_deriv) then
	  dumbool=test_gradient(x, func, dfunc,data=data)
          if (.not. dumbool) call print("Gradient test failed")
       end if

       !**********************************************************************
       !*
       !*  Print stuff
       !*
       !**********************************************************************/


#ifndef _OPENMP
       if (my_hook_print_interval == 1 .or. mod(main_counter,my_hook_print_interval) == 1) call verbosity_push_increment()
#endif
       call print(trim(method)//" iter = "//main_counter//" df^2 = "//normsqgrad_f//" f = "//f// &
            ' max(abs(df)) = '//maxval(abs(grad_f)),PRINT_VERBOSE)
       if(.not. do_lbfgs) &
            call print(" dcos = "//dcosine//" q = " //linmin_quality,PRINT_VERBOSE)
#ifndef _OPENMP
       if (my_hook_print_interval == 1 .or. mod(main_counter,my_hook_print_interval) == 1) call verbosity_pop()
#endif

       ! call the hook function
       if (present(hook)) then 
          call hook(x, grad_f, f, done, (mod(main_counter-1,my_hook_print_interval) == 0), data)
          if (done) then
             call print('hook reports that minim finished, exiting.', PRINT_NORMAL)
             exit_flag = 1
             cycle
          end if
       else
          call print("hook is not present", PRINT_VERBOSE)
       end if

       !**********************************************************************
       !*
       !* test to see if we've converged, let's quit
       !*
       !**********************************************************************/
       ! 
       if(normsqgrad_f < convergence_tol) then
          convergence_counter = convergence_counter + 1
          !call print("Convergence counter = "//convergence_counter)
       else
          convergence_counter = 0
       end if
       if(convergence_counter == convergence_window) then
          call print("Converged after step " // main_counter)
          call print(trim(method)//" iter = " // main_counter // " df^2= " // normsqgrad_f // " f= " // f)
          exit_flag = 1 ! while loop will quit
          cycle !continue
       end if

       if(.not. do_lbfgs .and. .not. do_sd2) then
          !**********************************************************************
          !*
          !*  do line minimization
          !*
          !**********************************************************************/
          
          oldeps = eps
          ! no output from linmin unless level >= PRINT_VERBOSE
	  call system_timer("minim/main_loop/"//main_counter//"/linmin")
#ifndef _OPENMP
          call verbosity_push_decrement()
#endif
          if(do_fast_linmin) then
             lsteps = linmin_fast(x, f, hdir, y, eps, func,data)
             if (lsteps .EQ.0) then
                call print("Fast linmin failed, calling normal linmin at step " // main_counter)
                lsteps = linmin(x, hdir, y, eps, func,data)
             end if
          else if(do_linmin_deriv) then
             lsteps = linmin_deriv_iter_simple(x, hdir, grad_f, y, eps, dfunc,data)
          else
             lsteps = linmin(x, hdir, y, eps, func,data)
          end if
          if ((oldeps .fne. my_eps_guess) .and. (eps > oldeps*2.0_dp)) then
            eps = oldeps*2.0_dp
          endif
#ifndef _OPENMP
          call verbosity_pop()
#endif
	  call system_timer("minim/main_loop/"//main_counter//"/linmin")

          !**********************************************************************
          !*
          !*  check the result of linmin
          !*
          !**********************************************************************
          
          if(lsteps .EQ. 0) then ! something very bad happenned, gradient is bad?
             call print("*** LINMIN returned 0, RESETTING CG CYCLE and eps  at step " // main_counter)
#ifndef _OPENMP
             call verbosity_push_increment()
#endif
             extra_report = extra_report + 1
             hdir = -1.0 * grad_f
             eps = my_eps_guess

	     if (current_verbosity() >= PRINT_NERD) call line_scan(x, hdir, func, .not. do_linmin_deriv, dfunc, data)
             
             bad_iter_counter = bad_iter_counter + 1
             if(bad_iter_counter .EQ. max_bad_iter) then  
                
                call print("*** BAD linmin counter reached maximum, exiting " // max_bad_iter)
                
                exit_flag = 1
                if (present(status)) status = 1
             end if

             cycle !continue
          end if

       end if ! .not. lbfgs
       
       !**********************************************************************
       !*
       !*  Evaluate function at new position
       !*
       !**********************************************************************
       

       if (.not. do_linmin_deriv) then
#ifndef _OPENMP
          call verbosity_push_decrement(2)
#endif
          f_new = func(y,data)
#ifndef _OPENMP
          call verbosity_pop()
#endif
       else
          f_new = 0.0_dp
       end if
       !       if(!finite(f_new)){
       !	logger(ERROR, "OUCH!!! f_new is not finite!\n");
       !	return -1; // go back screaming
       !       end if

       !********************************************************************
       ! Is everything going normal?
       !*********************************************************************/

       ! let's make sure we are progressing 
       ! obj is the thing we are trying to minimize
       ! are we going down,and enough?

       if(.not. do_lbfgs .and. .not. do_sd2) then
          if (.not. do_linmin_deriv) then
             obj = f
             obj_new = f_new
          else
             ! let's pretend that linmin_deriv never goes uphill!
             obj     =  0.0_dp
             obj_new = -1.0_dp
          end if
          
          if(obj-obj_new > abs(stuck_tol*obj_new)) then
             ! everything is fine, clear some monitoring flags
             bad_iter_counter = 0
             do i=1, extra_report
#ifndef _OPENMP
                call verbosity_pop()
#endif
             end do
             extra_report = 0
             
             if(present(linminroutine)) then
                if(trim(linminroutine) .EQ. "FAST_LINMIN" .and. .not. do_fast_linmin .and. &
		      main_counter >  fast_linmin_switchback) then
                   call print("Switching back to FAST_LINMIN linmin")
                   do_fast_linmin = .TRUE.
                end if
             end if
          !**********************************************************************
          !* Otherwise, diagnose problems
          !**********************************************************************
          else if(obj_new > obj)then !  are we not going down ? then things are very bad
          
             if( abs(obj-obj_new) < abs(stuck_tol*obj_new))then ! are we stuck ?
                call print("*** Minim is stuck, exiting")
                exit_flag = 1
                if (present(status)) status = 0
                cycle !continue 
             end if
             
             
             call print("*** Minim is not going down at step " // main_counter //" ==> eps /= 10")
             eps = oldeps / 10.0_dp

	     if (current_verbosity() >= PRINT_NERD) call line_scan(x, hdir, func, .not. do_linmin_deriv, dfunc, data)
             if(current_verbosity() >= PRINT_NERD .and. .not. do_linmin_deriv) then     
                if(.NOT.test_gradient(x, func, dfunc,data=data)) then
                   call print("*** Gradient test failed!!")
                end if
             end if
             
             
             if(do_fast_linmin) then
                do_fast_linmin = .FALSE.
                fast_linmin_switchback = main_counter+5
                call print("Switching off FAST_LINMIN (back after " //&
                     (fast_linmin_switchback - main_counter) // " steps if all OK")
             end if
             
             call print("Resetting conjugacy")
             resetflag = 1
             
             main_counter=main_counter-1      
             bad_iter_counter=bad_iter_counter+1 ! increment BAD counter
             
          else        ! minim went downhill, but not a lot
          
             call print("*** Minim is stuck at step " // main_counter // ", trying to unstick", PRINT_VERBOSE)
          
#ifndef _OPENMP
             if (current_verbosity() >= PRINT_NORMAL) then
                call verbosity_push_increment()
                extra_report = extra_report + 1
             end if
#endif
      
	     if (current_verbosity() >= PRINT_NERD) call line_scan(x, hdir, func, .not. do_linmin_deriv, dfunc, data)
             !**********************************************************************
             !*
             !* do gradient test if we need to
             !*
             !**********************************************************************/
             ! 
             if(current_verbosity() >= PRINT_NERD .and. .not. do_linmin_deriv) then     
                if(.NOT.test_gradient(x, func, dfunc,data=data)) then
                   call print("*** Gradient test failed!! Exiting linmin!")
                   exit_flag = 1
                   if (present(status)) status = 1
                   cycle !continue
                end if
             end if
          
             bad_iter_counter=bad_iter_counter+1 ! increment BAD counter
             eps = my_eps_guess ! try to unstick
             call print("resetting eps to " // eps,PRINT_VERBOSE)
             resetflag = 1 ! reset CG
          end if

          if(bad_iter_counter == max_bad_iter)then
             call print("*** BAD iteration counter reached maximum " // max_bad_iter // " exiting")
             exit_flag = 1
             if (present(status)) status = 1
             cycle !continue
          end if
       end if ! .not. do_bfgs and .not. do_sd2

       !**********************************************************************
       !*
       !*  accept iteration, get new gradient
       !*
       !**********************************************************************/

       f = f_new
       x = y
      
! Removed resetting 
!       if(mod(main_counter,50) .EQ. 0) then ! reset CG every now and then regardless
!          resetflag = 1
!       end if

       ! measure linmin_quality
       if(.not. do_lbfgs) hdirgrad_before = hdir.DOT.grad_f  

       grad_f_old = grad_f
#ifndef _OPENMP
       call verbosity_push_decrement(2)
#endif
       grad_f = dfunc(x,data)
#ifndef _OPENMP
       call verbosity_pop()
#endif

       if(.not. do_lbfgs) then
          hdirgrad_after = hdir.DOT.grad_f 
          if (hdirgrad_after /= 0.0_dp) then
            linmin_quality = hdirgrad_before/hdirgrad_after
          else
            linmin_quality = HUGE(1.0_dp)
          endif
       end if

       normsqgrad_f_old = normsqgrad_f
       normsqgrad_f = normsq(grad_f)
          

       !**********************************************************************
       !* Choose minimization method
       !**********************************************************************/
       
       if(do_sd) then  !steepest descent
          hdir = -1.0_dp * grad_f 
       elseif(do_sd2) then
          if(main_counter == 1) then
             alpha = 1.0e-6_dp
          else
             alpha = dot_product(x-x_old,grad_f-grad_f_old) / dot_product(grad_f-grad_f_old,grad_f-grad_f_old)
          endif

          x_old = x
          x = x - alpha * grad_f
          y = x
       else if(do_cg)then ! conjugate gradients
          
          if( bad_cg_counter == max_bad_cg .OR.(resetflag > 0)) then ! reset the conj grad cycle
             
             if( bad_cg_counter .EQ. max_bad_cg)  then 
                call print("*** bad_cg_counter == "//  max_bad_cg, PRINT_VERBOSE)
             end if
             call print("*** Resetting conjugacy", PRINT_VERBOSE)
             
             
             hdir = -1.0_dp * grad_f
             bad_cg_counter = 0
             resetflag = 0
          else  ! do CG
             gdirlen = 0.0
             dcosine = 0.0
             
             gg = normsq(gdir)
             
             if(.NOT.do_pcg) then ! no preconditioning
                dgg = max(0.0_dp, (gdir + grad_f).DOT.grad_f) ! Polak-Ribiere formula
                gdir = (-1.0_dp) * grad_f
                ! the original version was this, I had to change because intrinsic does'nt return  allocatables.
                !                dgg = (grad_f + gdir).DOT.grad_f
                !                gdir = (-1.0_dp) * grad_f
             else ! precondition
                call System_abort("linmin: preconditioning not implemented")
                dgg = 0.0 ! STUPID COMPILER
                !  //dgg = (precond^grad_f+gdir)*grad_f;
                !  //gdir = -1.0_dp*(precond^grad_f);
             end if
             if(gg .ne. 0.0_dp) then
                hdir = gdir+hdir*(dgg/gg)
             else
                hdir = gdir
             endif
             ! calculate direction cosine
             dcosine = gdir.DOT.gdir_old
             gdir_old = gdir
             gdirlen = norm(gdir)
             
             if(gdirlen .eq. 0.0_dp .or. gdirlen_old .eq. 0.0_dp) then
                dcosine = 0.0_dp
             else
                dcosine = dcosine/(gdirlen*gdirlen_old)
             endif
             gdirlen_old = gdirlen
             
             if(abs(dcosine) >  0.2) then 
                bad_cg_counter= bad_cg_counter +1
             else 
                bad_cg_counter = 0
             end if
             
          end if
       else if(do_lbfgs)then  ! LBFGS method
          y = x
          call LBFGS(size(x),lbfgs_M,y, f, grad_f, .false., lbfgs_diag, (/-1,0/), 1e-12_dp, 1e-12_dp, lbfgs_work, lbfgs_flag)
!          do while(lbfgs_flag == 2)
!             call LBFGS(size(x),lbfgs_M,y, f, grad_f, .false., lbfgs_diag, (/-1,0/), 1e-12_dp, 1e-12_dp, lbfgs_work, lbfgs_flag)             
!          end do
          if(lbfgs_flag < 0) then ! internal LBFGS error
             call print('LBFGS returned error code '//lbfgs_flag//', exiting')
             exit_flag = 1
             if (present(status)) status = 1
             cycle
          end if
       else
          call System_abort("minim(): c'est ci ne pas une erreur!")
       end if

      call system_timer("minim/main_loop/"//main_counter)
       main_counter=main_counter + 1
    end do
    call system_timer("minim/main_loop")

    if(main_counter >= max_steps) then 
       call print("Iterations exceeded " // max_steps)
    end if

    call print("Goodbye from minim()")
    call print("")
    minim = main_counter-1

    x_in = x

    deallocate(x)
    deallocate(y, hdir, gdir, gdir_old, grad_f, grad_f_old)
    if(allocated(x_old)) deallocate(x_old)

    if(do_lbfgs) then
       deallocate(lbfgs_diag)
       deallocate(lbfgs_work)
    end if

    ! just in case extra pushes weren't popped
#ifndef _OPENMP
    do i=1, extra_report
       call verbosity_pop()
    end do
#endif

    call system_timer("minim")

  end function minim


  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  ! test_gradient
  !
  !% Test a function against its gradient by evaluating the gradient from the
  !% function by finite differnces. We can only test the gradient if energy and force 
  !% functions are pure in that they do not change the input vector 
  !% (e.g. no capping of parameters). The interface for 'func(x)' and 'dfunc(x)'
  !% are the same as for 'minim' above.
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


  function test_gradient(xx,func,dfunc, dir,data)
    real(dp),intent(in)::xx(:) !% Position
    INTERFACE 
       function func(x,data)
         use system_module
         real(dp)::x(:)
	 character(len=1),optional::data(:)
         real(dp)::func
       end function func
    end INTERFACE
    INTERFACE
       function dfunc(x,data)
         use system_module
         real(dp)::x(:)
	 character(len=1),optional::data(:)
         real(dp)::dfunc(size(x))
       end function dfunc
    END INTERFACE
    real(dp), intent(in), optional, target :: dir(:) !% direction along which to test
    character(len=1),optional::data(:)
    !%RV Returns true if the gradient test passes, or false if it fails
    logical :: test_gradient


    integer:: N,I,loopend
    logical::ok,printit,monitor_ratio
    real(dp)::f0, f, tmp, eps, ratio, previous_ratio
    ! Dynamically allocate to avoid stack overflow madness with ifort
    real(dp),dimension(:), allocatable ::x,dx,x_0
    real(dp), allocatable :: my_dir(:)

    N=size(xx)
    allocate(x(N), dx(N), x_0(N))
    x=xx

    if(current_verbosity() > PRINT_VERBOSE) then 
       printit = .TRUE.
       loopend=0
    else
       printit = .FALSE.
       loopend=1
    end if

    if(printit) then
       write(line,*)" ";  call print(line)
       write(line,*)" " ; call print(line)
       write(line,*) "Gradient test";  call print(line)
       write(line,*)" " ; call print(line)
    end if

    if(printit) then
       write(line, *) "Calling func(x)"; call print(line)
    end if
    !f0 = param_penalty(x_0);
#ifndef _OPENMP
    call verbosity_push_decrement()
#endif
    f0 = func(x,data)
#ifndef _OPENMP
    call verbosity_pop()
#endif
    !!//logger("f0: %24.16f\n", f0);

    if(printit) then
       write(line, *) "Calling dfunc(x)"; call print(line)
    end if
#ifndef _OPENMP
    call verbosity_push_decrement()
#endif
    dx = dfunc(x,data)
#ifndef _OPENMP
    call verbosity_pop()
#endif
    !!//dx = param_penalty_deriv(x);

    allocate(my_dir(N))
    if (present(dir)) then
       if (norm(dir) .feq. 0.0_dp) &
	 call system_abort("test_gradient got dir = 0.0, can't use as normalized direction for test")
       my_dir = dir/norm(dir)
    else
       if (norm(dx) == 0.0_dp) &
	 call system_abort("test_gradient got dfunc = 0.0, can't use as normalized direction for test")
       my_dir = (-1.0_dp)*dx/norm(dx)
    endif
    !//my_dir.zero();
    !//my_dir.randomize(0.1);
    !//my_dir.x[4] = 0.1;
    !//logger("dx: ");dx.print(logger_stream);

    tmp = my_dir.DOT.dx
    x_0(:) = x(:)
    !//logger("x0:  "); x_0.print(logger_stream);

    if(printit) then
       call print("GT  eps        (f-f0)/(eps*df)     f")
    end if

    ok = .FALSE.
    ratio = 0.0 

    do i=0,loopend ! do it twice, print second time if not OK
       previous_ratio = 0.0
       monitor_ratio = .FALSE.

       eps=1.0e-1_dp
       do while(eps>1.e-20)


          x = x_0 + eps*my_dir
          !//logger("x:  "); x.print(logger_stream);
          !//f = param_penalty(x);
#ifndef _OPENMP
          call verbosity_push_decrement()
#endif
          f = func(x,data)
#ifndef _OPENMP
          call verbosity_pop()
#endif
          !//logger("f: %24.16f\n", f);
          previous_ratio = ratio
          ratio = (f-f0)/(eps*tmp)

          if(printit) then
             write(line,'("GT ",e8.2,f22.16,e24.16)') eps, ratio, f;
             call print(line)
          end if

          if(abs(ratio-1.0_dp) .LT. 1e-2_dp) monitor_ratio = .TRUE.

          if(.NOT.monitor_ratio) then
             if(abs((f-f0)/f0) .LT. NUMERICAL_ZERO) then ! we ran out of precision, gradient is really bad
		call print("(f-f0)/f0 " // ((f-f0)/f0) // " ZERO " // NUMERICAL_ZERO, PRINT_ANAL)
		call print("ran out of precision, quitting loop", PRINT_ANAL)
                exit
             end if
          end if

          if(monitor_ratio) then
             if( abs(ratio-1.0_dp) > abs(previous_ratio-1.0_dp) )then !  sequence broke
                if(abs((f-f0)/f0*(ratio-1.0_dp)) < 1e-10_dp) then ! lets require 10 digits of precision
                   ok = .TRUE. 
                   !//break;
                end if
             end if
          end if

          eps=eps/10.0
       end do

       if(.NOT.ok) then 
          printit = .TRUE.  ! go back and print it
       else
          exit
       end if

    end do
    if(printit) then
       write(line,*)" "; call print(line, PRINT_NORMAL)
    end if

    if(ok) then
       write(line,*)"Gradient test OK"; call print(line, PRINT_VERBOSE)
    end if
    test_gradient= ok

    deallocate(x, dx, x_0)
    deallocate(my_dir)


  end function test_gradient

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  ! n_test_gradient
  !
  !% Test a function against its gradient by evaluating the gradient from the
  !% function by symmetric finite differnces. We can only test the gradient if 
  !% energy and force functions are pure in that they do not change the input 
  !% vector (e.g. no capping of parameters). The interface for 'func(x)' and 
  !% 'dfunc(x)'are the same as for 'minim' above.
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine n_test_gradient(xx,func,dfunc, dir,data)
    real(dp),intent(in)::xx(:) !% Position
    INTERFACE 
       function func(x,data)
         use system_module
         real(dp)::x(:)
	 character(len=1),optional::data(:)
         real(dp)::func
       end function func
    end INTERFACE
    INTERFACE
       function dfunc(x,data)
         use system_module
         real(dp)::x(:)
	 character(len=1),optional::data(:)
         real(dp)::dfunc(size(x))
       end function dfunc
    END INTERFACE
    real(dp), intent(in), optional, target :: dir(:) !% direction along which to test
    character(len=1),optional::data(:)

    integer N, i
    real(dp) :: eps, tmp, fp, fm
    real(dp), allocatable :: x(:), dx(:)
    real(dp), pointer :: my_dir(:)

    N=size(xx)
    allocate(x(N), dx(N))
    x=xx

    dx = dfunc(x,data)

    if (present(dir)) then
       if (norm(dir) == 0.0_dp) &
	 call system_abort("n_test_gradient got dir = 0.0, can't use as normalized direction for test")
       allocate(my_dir(N))
       my_dir = dir/norm(dir)
    else
       if (norm(dx) == 0.0_dp) &
	 call system_abort("n_test_gradient got dfunc = 0.0, can't use as normalized direction for test")
       allocate(my_dir(N))
       my_dir = (-1.0_dp)*dx/norm(dx)
    endif

    tmp = my_dir.DOT.dx

    do i=2, 10
       eps = 10.0_dp**(-i)
       x = xx + eps*my_dir
       fp = func(x,data)
       x = xx - eps*my_dir
       fm = func(x,data)
       call print ("fp " // fp // " fm " // fm)
       call print("GT_N eps " // eps // " D " // tmp // " FD " // ((fp-fm)/(2.0_dp*eps)) // " diff " // (tmp-(fp-fm)/(2.0_dp*eps)))
    end do

    deallocate(x, dx)
    deallocate(my_dir)

  end subroutine n_test_gradient


  !% FIRE MD minimizer from Bitzek et al., \emph{Phys. Rev. Lett.} {\bfseries 97} 170201.
  !% Beware, this algorithm is patent pending in the US.
  function fire_minim(x, mass, func, dfunc, dt0, tol, max_steps, hook, hook_print_interval, data, dt_max, status)
    real(dp), intent(inout), dimension(:) :: x
    real(dp), intent(in) :: mass
    interface 
       function func(x,data)
         use system_module
         real(dp)::x(:)
	 character(len=1),optional::data(:)
         real(dp)::func
       end function func
       function dfunc(x,data)
         use system_module
         real(dp)::x(:)
	 character(len=1),optional::data(:)
         real(dp)::dfunc(size(x))
       end function dfunc
    end interface
    real(dp), intent(in) :: dt0
    real(dp), intent(in) :: tol
    integer,  intent(in) :: max_steps
    optional :: hook
    interface
       subroutine hook(x,dx,E,done,do_print,data)
         use system_module
         real(dp), intent(in) ::x(:)
         real(dp), intent(in) ::dx(:)
         real(dp), intent(in) ::E
         logical, intent(out) :: done
	 logical, optional, intent(in) :: do_print
	 character(len=1),optional, intent(in) ::data(:)
       end subroutine hook
    end interface
    integer, optional :: hook_print_interval
    character(len=1),optional::data(:)
    real(dp), intent(in), optional :: dt_max
    integer, optional, intent(out) :: status
    integer :: fire_minim

    integer  :: my_hook_print_interval

    real(dp), allocatable, dimension(:) :: velo, acc, force
    real(dp) :: f, df2, alpha_start, alpha, P, dt, my_dt_max
    integer :: i, Pcount
    logical :: done

    if (present(status)) status = 0

    if (current_verbosity() >= PRINT_VERBOSE) then
      my_hook_print_interval = optional_default(1, hook_print_interval)
    else if (current_verbosity() >= PRINT_NORMAL) then
      my_hook_print_interval = optional_default(10, hook_print_interval)
    else if (current_verbosity() >= PRINT_SILENT) then
      my_hook_print_interval = optional_default(100000, hook_print_interval)
    endif

    allocate (velo(size(x)), acc(size(x)), force(size(x)))

    alpha_start = 0.1_dp
    alpha = alpha_start
    dt = dt0
    my_dt_max = optional_default(20.0_dp*dt, dt_max)
    Pcount = 0

    call Print('Welcome to fire_minim', PRINT_NORMAL)

    velo = 0
    acc = -dfunc(x,data)/mass

    do i = 1, max_steps

       velo = velo + (0.5_dp*dt)*acc
       force = dfunc(x,data)  ! keep original sign of gradient for now as hook() expects gradient, not force

       df2 = normsq(force)

       if (df2 < tol) then
          write (line, '(i4,a,e10.2)')  i, ' force^2 = ', df2
          call Print(line, PRINT_NORMAL)

          write (line, '(a,i0,a)') 'Converged in ', i, ' steps.'
          call Print(line, PRINT_NORMAL)
          exit
       else if(mod(i, my_hook_print_interval) == 0) then
          f = func(x,data)
          write (line, '(i4,a,e24.16,a,e24.16,a,f0.3)') i,' f=',f,' df^2=',df2,' dt=', dt
          call Print(line, PRINT_NORMAL)

          if(present(hook)) then
             call hook(x, force, f, done, (mod(i-1,my_hook_print_interval) == 0), data)
             if (done) then
                call print('hook reports that fire_minim is finished, exiting.', PRINT_NORMAL)
                exit
             end if
          end if
       end if

       force = -force ! convert from gradient of energy to force
       acc = force/mass
       velo = velo + (0.5_dp*dt)*acc

       P = velo .dot. force
       velo = (1.0_dp-alpha)*velo + (alpha*norm(velo)/norm(force))*force
       if (P > 0) then
          if(Pcount > 5) then
             dt = min(dt*1.1_dp, my_dt_max)
             alpha = alpha*0.99_dp
          else
             Pcount = Pcount + 1
          end if
       else
          dt = dt*0.5_dp
          velo = 0.0_dp
          alpha = alpha_start
          Pcount = 0
       end if

       x = x + dt*velo
       x = x + (0.5_dp*dt*dt)*acc

       if(current_verbosity() >= PRINT_NERD) then
          write (line, '(a,e24.16)')  'E=', func(x,data)
          call print(line,PRINT_NERD)
          call print(x,PRINT_NERD)
       end if
    end do

    if(i == max_steps) then
       write (line, '(a,i0,a)') 'fire_minim: Failed to converge in ', i, ' steps.'
       call Print(line, PRINT_ALWAYS)
       if (present(status)) status = 1
    end if

    fire_minim = i

    deallocate(velo, acc, force)

  end function fire_minim



!%%%%
!% Noam's line minimizer
!%
!% args
!% x: input vector, flattened into 1-D double array
!% bothfunc: input pointer to function computing value and gradient
!% neg_gradient: input negative of gradient (i.e. force) for input x
!% E: input value at x
!% search_dir: input vector search direction, not necessarily normalized
!% new_x: output x at minimimum
!% new_new_gradient: output negative of gradient at new_x
!% new_E: output value at new_x
!% max_step_size: input max step size (on _normalized_ search dir), output actual step size for this linmin
!% accuracy: input desired accuracy on square of L2 norm of projected gradient
!% N_evals: input initial number of function evals so far, output final number of function evalutions
!% max_N_evals: input max number of function evaluations before giving up
!% hook: input pointer to function to call after each evaluation
!% data: input pointer to other data needed for calc, flattened to 1-D character array with transfer()
!% error: output error state
!%%%%

subroutine n_linmin(x, bothfunc, neg_gradient, E, search_dir, &
		  new_x, new_neg_gradient, new_E, &
		  max_step_size, accuracy, N_evals, max_N_evals, hook, hook_print_interval, &
		  data, error)
    real(dp), intent(inout) :: x(:)
    real(dp), intent(in) :: neg_gradient(:)
    interface 
       subroutine bothfunc(x,E,f,data,error)
         use system_module
         real(dp)::x(:)
         real(dp)::E
         real(dp)::f(:)
	 character(len=1),optional::data(:)
	 integer,optional :: error
       end subroutine bothfunc
    end interface
    real(dp), intent(inout) :: E
    real(dp), intent(inout) :: search_dir(:)
    real(dp), intent(out) :: new_x(:), new_neg_gradient(:)
    real(dp), intent(out) ::  new_E
    real(dp), intent(inout) ::  max_step_size
    real(dp), intent(in) :: accuracy
    integer, intent(inout) :: N_evals
    integer, intent(in) :: max_N_evals
    optional :: hook
    interface 
       subroutine hook(x,dx,E,done,do_print,data)
         use system_module
         real(dp), intent(in) ::x(:)
         real(dp), intent(in) ::dx(:)
         real(dp), intent(in) ::E
         logical, intent(out) :: done
	 logical, optional, intent(in) :: do_print
	 character(len=1),optional, intent(in) ::data(:)
       end subroutine hook
    end interface
    integer, intent(in), optional :: hook_print_interval
    character(len=1),optional::data(:)
    integer, intent(out), optional :: error

    logical :: do_print
    real(dp) search_dir_mag
    real(dp) p0_dot, p1_dot, new_p_dot
    real(dp), allocatable :: p0(:), p1(:), p0_ng(:), p1_ng(:), new_p(:), new_p_ng(:), t_projected(:)
    real(dp) p0_E, p1_E, new_p_E
    real(dp) p0_pos, p1_pos, new_p_pos

    real(dp) tt
    real(dp) t_a, t_b, t_c, Ebar, pbar, soln_1, soln_2

    logical done, use_cubic, got_valid_cubic

    integer l_error
    real(dp) est_step_size

    INIT_ERROR(error)

    do_print = .false.

    allocate(p0(size(x)))
    allocate(p1(size(x)))
    allocate(p0_ng(size(x)))
    allocate(p1_ng(size(x)))
    allocate(new_p(size(x)))
    allocate(new_p_ng(size(x)))
    allocate(t_projected(size(x)))

    call print ("n_linmin starting line minimization", PRINT_NERD)
    call print ("n_linmin initial |x| " // norm(x) // " |neg_gradient| " // &
      norm(neg_gradient) // " |search_dir| " // norm(search_dir), PRINT_NERD)

    search_dir_mag = norm(search_dir)
    ! search_dir = search_dir / search_dir_mag
    search_dir = search_dir / search_dir_mag

    p0_dot = neg_gradient .dot. search_dir

    t_projected = search_dir*p0_dot
    if (normsq(t_projected) .lt. accuracy) then
	call print ("n_linmin initial config is apparently converged " // norm(t_projected) // " " // accuracy, PRINT_NERD)
    endif

    p0_pos = 0.0_dp
    p0 = x
    p0_ng = neg_gradient
    p0_E = E

    p0_dot = p0_ng .dot. search_dir
    call print("initial p0_dot " // p0_dot, PRINT_NERD)

    if (p0_dot .lt. 0.0_dp) then
	p0_ng = -p0_ng
	p0_dot = -p0_dot
    endif

    call print("cg_n " // p0_pos // " " // p0_e // " " // p0_dot // " " // &
	       normsq(p0_ng) // " " // N_evals // " bracket starting", PRINT_VERBOSE)

    est_step_size = 4.0_dp*maxval(abs(p0_ng))**2/p0_dot
    if (max_step_size .gt. 0.0_dp .and. est_step_size .gt. max_step_size) then
	est_step_size = max_step_size
    endif

    p1_pos = est_step_size
    p1 = x + p1_pos*search_dir
    p1_ng = p0_ng
    l_error = 1
    do while (l_error .ne. 0)
	N_evals = N_evals + 1
	l_error=0
	call bothfunc(p1, p1_e, p1_ng, data, error=l_error); p1_ng = -p1_ng
	if (present(hook_print_interval)) do_print = (mod(N_evals, hook_print_interval) == 0)
	if (present(hook)) call hook(p1,p1_ng,p1_E,done,do_print,data)
	if (l_error .ne. 0) then
	    call print("cg_n " // p1_pos // " " // p1_e // " " // 0.0_dp // " " // &
		       0.0_dp // " " // N_evals // " bracket first step ERROR", PRINT_ALWAYS)
	    est_step_size = est_step_size*0.5_dp
	    p1_pos = est_step_size
	    p1 = x + p1_pos*search_dir
	endif
	if (N_evals .gt. max_N_evals) then
	    RAISE_ERROR_WITH_KIND(ERROR_MINIM_NOT_CONVERGED, "n_linmin ran out of iterations", error)
	endif
    end do

    p1_dot = p1_ng .dot. search_dir

    call print("cg_n " // p1_pos // " " // p1_e // " " // p1_dot // " " // &
	       normsq(p1_ng) // " " // N_evals // " bracket first step", PRINT_VERBOSE)

    t_projected = search_dir*p1_dot
!     if (object_norm(t_projected,norm_type) .lt. accuracy) then
! 	new_x = p1
! 	new_neg_gradient = p1_ng
! 	new_E = p1_E
! 	! search_dir = search_dir * search_dir_mag
! 	call scalar_selfmult (search_dir, search_dir_mag)
! 	minimize_along = 0
! call print ("returning from minimize_along, t_projected is_converged")
! 	return
!     endif

    call print ("starting bracketing loop", PRINT_NERD)

    ! bracket solution
    do while (p1_dot .ge. 0.0_dp)
	p0 = p1
	p0_ng = p1_ng
	p0_E = p1_E
	p0_pos = p1_pos
	p0_dot = p1_dot

	p1_pos = p1_pos + est_step_size

	call print ("checking bracketing for " // p1_pos, PRINT_NERD)

	p1 = x + p1_pos*search_dir
	l_error = 1
	do while (l_error .ne. 0)
	    N_evals = N_evals + 1
	    l_error = 0
	    call bothfunc (p1, p1_E, p1_ng, data, error=l_error); p1_ng = -p1_ng
	    if (present(hook_print_interval)) do_print = (mod(N_evals, hook_print_interval) == 0)
	    if (present(hook)) call hook(p1,p1_ng,p1_E,done,do_print,data)
	    if (done) then
	      call print("hook reported done", PRINT_NERD)
	      search_dir = search_dir * search_dir_mag
	      return
	    endif
	    if (l_error .ne. 0) then
	      call print("cg_n " // p0_pos // " " // p0_e // " " // 0.0_dp // " " // &
			 0.0_dp // " " // N_evals // " bracket loop ERROR", PRINT_ALWAYS)
	      call print ("Error in bracket loop " // l_error // " stepping back", PRINT_ALWAYS)
	      p1_pos = p1_pos - est_step_size
	      est_step_size = est_step_size*0.5_dp
	      p1_pos = p1_pos + est_step_size
	      p1 = x + p1_pos*search_dir
	    endif
	    if (N_evals .gt. max_N_evals) then
	      search_dir = search_dir * search_dir_mag
	      RAISE_ERROR_WITH_KIND(ERROR_MINIM_NOT_CONVERGED, "n_linmin ran out of iterations", error)
	    endif
	end do

	p1_dot = p1_ng .dot. search_dir

!	tt = -p0_dot/(p1_dot-p0_dot)
!	if (1.5D0*tt*(p1_pos-p0_pos) .lt. 10.0*est_step_size) then
!	    est_step_size = 1.5_dp*tt*(p1_pos-p0_pos)
!	else
	    est_step_size = est_step_size*2.0_dp
!	end if

      call print("cg_n " // p1_pos // " " // p1_e // " " // p1_dot // " " // &
		 normsq(p1_ng) // " " // N_evals // " bracket loop", PRINT_VERBOSE)
    end do
    
    call print ("bracketed by" // p0_pos // " " // p1_pos, PRINT_NERD)

    done = .false.
    t_projected = 2.0_dp*sqrt(accuracy)
    !new_p_dot = accuracy*2.0_dp
    do while (normsq(t_projected) .ge. accuracy .and. (.not. done))
        call print ("n_linmin starting true minimization loop", PRINT_NERD)

	use_cubic = .false.
	if (use_cubic) then
	    !!!! fit to cubic polynomial
	    Ebar = p1_E-p0_E
	    pbar = p1_pos-p0_pos
	    t_a = (-p0_dot)
	    t_c = (pbar*((-p1_dot)-(-p0_dot)) - 2.0_dp*Ebar + 2*(-p0_dot)*pbar)/pbar**3
	    t_b = (Ebar - (-p0_dot)*pbar - t_c*pbar**3)/pbar**2

	    soln_1 = (-2.0_dp*t_b + sqrt(4.0_dp*t_b**2 - 12.0_dp*t_a*t_c))/(6.0_dp*t_c)
	    soln_2 = (-2.0_dp*t_b - sqrt(4.0_dp*t_b**2 - 12.0_dp*t_a*t_c))/(6.0_dp*t_c)

	    if (soln_1 .ge. 0.0_dp .and. soln_1 .le. pbar) then
		new_p_pos = p0_pos + soln_1
		got_valid_cubic = .true.
	    else if (soln_2 .ge. 0.0_dp .and. soln_2 .le. pbar) then
		new_p_pos = p0_pos + soln_2
		got_valid_cubic = .true.
	    else
		call print ("n_linmin warning: no valid solution for cubic", PRINT_ALWAYS)
		!!!! use only derivative information to find pt. where derivative = 0
		tt = -p0_dot/(p1_dot-p0_dot)
		new_p_pos = p0_pos + tt*(p1_pos-p0_pos)
		done = .false.
	    endif
	else
	    !!!! use only derivative information to find pt. where derivative = 0
	    tt = -p0_dot/(p1_dot-p0_dot)
	    new_p_pos = p0_pos + tt*(p1_pos-p0_pos)
	    ! done = .true.
	    done = .false.
	endif

	new_p = x + new_p_pos*search_dir
	N_evals = N_evals + 1
	call bothfunc (new_p, new_p_E, new_p_ng, data, error); new_p_ng = -new_p_ng
	if (error .ne. 0) then
	    call system_abort("n_linmin: Error in line search " // error)
	endif

	if (N_evals .gt. max_N_evals) done = .true.

!	if (inner_prod(new_p_ng,new_p_ng) .lt. 0.1 .and. got_valid_cubic) done = .true.
!	if (got_valid_cubic) done = .true.

	new_p_dot = new_p_ng .dot. search_dir

	call print("cg_n " // new_p_pos // " " // new_p_E // " " // new_p_dot // " " // &
		   normsq(new_p_ng) // " " // N_evals // " during line search", PRINT_VERBOSE)

	if (new_p_dot .gt. 0) then
	    p0 = new_p
	    p0_pos = new_p_pos
	    p0_dot = new_p_dot
	    p0_ng = new_p_ng
	    p0_E = new_p_E
	else
	    p1 = new_p
	    p1_pos = new_p_pos
	    p1_dot = new_p_dot
	    p1_ng = new_p_ng
	    p1_E = new_p_E
	endif

	t_projected = search_dir*new_p_dot 
    end do

    new_x = new_p
    new_neg_gradient = new_p_ng
    new_E = new_p_E

    max_step_size = new_p_pos
    search_dir = search_dir * search_dir_mag

    call print ("done with line search", PRINT_NERD)
end subroutine n_linmin

!%%%%
!% Noam's minimizer with preconditioning from Cristoph Ortner
!% return value: number of function evaluations
!% args
!% x_i: input vector, flattened into a 1-D double array
!% bothfunc: input pointer to function that returns value and gradient
!%    can apply simple constraints and external (body) forces
!% use_precond: input logical controling preconditioning
!% apply_precond_func: input pointer to function that applies preconditioner to a vector
!% initial_E: output initial value
!% final_E: output final value
!% expected reduction: input expected reduction in value used to estimate initial step
!% max_N_evals: max number of function evaluations before giving up
!% accuracy: desired accuracy on square of L2 norm of gradient
!% hook: pointer to function passed to linmin and called once per CG step
!%       does stuff like print configuration, and can also apply other ending conditions,
!%       although latter capability isn't used
!% hook_print_interval: how often to call hook, depending on verbosity level
!% data: other data both_func will need to actually do calculation, flattened into
!%       character array by transfer() function.  Maybe be replaced with F2003 pointer soon
!% error: output error state
!%%%%
function n_minim(x_i, bothfunc, use_precond, apply_precond_func, initial_E, final_E, &
    expected_reduction, max_N_evals, accuracy, hook, hook_print_interval, data, error) result(N_evals)
    real(dp), intent(inout) :: x_i(:)
    interface 
       subroutine bothfunc(x,E,f,data,error)
         use system_module
         real(dp)::x(:)
         real(dp)::E
         real(dp)::f(:)
	 character(len=1),optional::data(:)
	 integer, optional :: error
       end subroutine bothfunc
    end interface 
    logical :: use_precond
    interface 
       subroutine apply_precond_func(x,g,P_g,data,error)
         use system_module
         real(dp)::x(:),g(:),P_g(:)
	 character(len=1),optional::data(:)
	 integer, optional :: error
       end subroutine apply_precond_func
    end interface 
    real(dp), intent(out) :: initial_E, final_E
    real(dp), intent(inout) :: expected_reduction
    integer, intent(in) :: max_N_evals
    real(dp), intent(in) :: accuracy
    optional :: hook
    integer :: N_evals
    interface 
       subroutine hook(x,dx,E,done,do_print,data)
         use system_module
         real(dp), intent(in) ::x(:)
         real(dp), intent(in) ::dx(:)
         real(dp), intent(in) ::E
         logical, intent(out) :: done
	 logical, optional, intent(in) :: do_print
	 character(len=1),optional, intent(in) ::data(:)
       end subroutine hook
    end interface
    integer, optional :: hook_print_interval
    character(len=1),optional::data(:)
    integer, optional, intent(out) :: error

    real(dp) :: E_i, E_ip1

    logical :: done, hook_done
    integer :: iter
    real(dp) :: max_step_size, initial_step_size

    real(dp), allocatable :: g_i(:)
    real(dp), allocatable :: x_ip1(:), g_ip1(:)
    real(dp), allocatable :: h_i(:)
    real(dp), allocatable :: P_g(:)
   
    real(dp) :: g_i_dot_g_i, g_ip1_dot_g_i, g_ip1_dot_g_ip1
    real(dp) :: gamma_i

    integer :: l_error

    integer :: my_hook_print_interval

    INIT_ERROR(error)

    if (current_verbosity() >= PRINT_VERBOSE) then
      my_hook_print_interval = optional_default(1, hook_print_interval)
    else if (current_verbosity() >= PRINT_NORMAL) then
      my_hook_print_interval = optional_default(10, hook_print_interval)
    else if (current_verbosity() >= PRINT_SILENT) then
      my_hook_print_interval = optional_default(100000, hook_print_interval)
    endif

    allocate(g_i(size(x_i)))
    allocate(x_ip1(size(x_i)))
    allocate(g_ip1(size(x_i)))
    allocate(h_i(size(x_i)))
    N_evals = 1
    call bothfunc(x_i, E_i, g_i, data, error=error); g_i = -g_i
    if (present(hook)) call hook(x_i, g_i, E_i, done, .true., data)
    PASS_ERROR_WITH_INFO("n_minim first evaluation", error)
    initial_E = E_i

    allocate(P_g(size(x_i)))
    ! need to get P into the routine somehow
    if (use_precond) then
      call apply_precond_func(x_i, g_i, P_g, data, error=error)
      PASS_ERROR_WITH_INFO("n_minim initial preconditioning call", error)
    else
      P_g = g_i
    Endif

    call print("#cg_n use_precond="//use_precond)

    call print("#cg_n " // "              x" // &
                           "             val" // &
                           "            -grad.dir" // &
                           "                |grad|^2" // "   n_evals")
    call print("cg_n " // 0.0_dp // " " // E_i // " " // (g_i.dot.g_i) // " " // &
	       normsq(g_i) // " " // N_evals // " INITIAL_VAL")

    if (normsq(g_i) .lt. accuracy) then
	call print ("n_minim initial config is converged " // norm(g_i) // " " // accuracy, PRINT_VERBOSE)
	final_E = initial_E
	return
    endif

    h_i = P_g

    iter = 1
    done = .false.
    do while (N_evals .le. max_N_evals .and. (.not.(done)))

	!! max_step_size = 4.0_dp*expected_reduction / norm(g_i)
	max_step_size = 1.0_dp * expected_reduction / (g_i .dot. (h_i/norm(h_i))) ! dividing by norm(h_i) because n_linmin will normalize h_i
	! if (max_step_size .gt. 1.0_dp) then
	    ! max_step_size = 1.0_dp
	! endif
	call print("max_step_size "//max_step_size, verbosity=PRINT_VERBOSE)

	call print("cg_n " // 0.0_dp // " " // E_i // " " // (g_i.dot.h_i) // " " // &
		  normsq(g_i) // " " // N_evals // " n_minim pre linmin")
	l_error = ERROR_NONE
	call n_linmin(x_i, bothfunc, g_i, E_i, h_i, &
		      x_ip1, g_ip1, E_ip1, &
		      max_step_size, accuracy, N_evals, max_N_evals, hook, hook_print_interval, data, l_error)
	if (l_error == ERROR_MINIM_NOT_CONVERGED) then
	   if (N_evals > max_N_evals) then
	     final_E = E_i
	     RAISE_ERROR_WITH_KIND(l_error, "linmin: n_minim didn't converge", error)
	   endif
	   ! we're just going to continue after an unconverged linmin,
	   CLEAR_ERROR()
	else
	   PASS_ERROR_WITH_INFO("linmin: n_minim error", error)
	endif

	call print("cg_n " // 0.0_dp // " " // E_ip1 // " " // (g_ip1.dot.h_i) // " " // &
		    normsq(g_ip1) // " " // N_evals // " n_minim post linmin")

	if (E_ip1 > E_i) then
	  final_E = E_i
	  call print("WARNING:n_minim: n_limin stepped uphill - forces may not be consistent with energy", verbosity=PRINT_ALWAYS)
	  ! RAISE_ERROR("n_minim: n_limin stepped uphill - forces may not be consistent with energy", error)
	endif

	if (normsq(g_ip1) .lt. accuracy) then
	    call print("n_minim is converged", PRINT_VERBOSE)
	    E_i = E_ip1
	    x_i = x_ip1
	    g_i = g_ip1
	    done = .true.
	endif

	! gamma_i = sum(g_ip1*g_ip1)/sum(g_i*g_i) ! Fletcher-Reeves
	! gamma_i = sum((g_ip1-g_i)*g_ip1)/sum(g_i*g_i) ! Polak-Ribiere
	g_i_dot_g_i = g_i .dot. P_g
	g_ip1_dot_g_i = g_ip1 .dot. P_g
	!! perhaps have some way of telling apply_precond_func to update/not update preconitioner?
	if (use_precond) then
	  call apply_precond_func(x_ip1, g_ip1, P_g, data, error=error)
	  PASS_ERROR_WITH_INFO("n_minim in-loop preconditioning call", error)
	else
	  P_g = g_ip1
	endif
	g_ip1_dot_g_ip1 = g_ip1 .dot. P_g
	gamma_i = (g_ip1_dot_g_ip1 - g_ip1_dot_g_i)/g_i_dot_g_i
	! steepest descent
	! gamma_i = 0.0_dp
	h_i = gamma_i*h_i + P_g

	if (iter .eq. 1) then
	  expected_reduction = abs(E_i - E_ip1)/10.0_dp
	else
	  expected_reduction = abs(E_i - E_ip1)/2.0_dp
	endif

	E_i = E_ip1
	x_i = x_ip1
	g_i = g_ip1
	! P_g is already P_g_ip1 from gamma_i evaluation code

	if (present(hook)) then
	  call hook(x_i,g_i,E_i,hook_done,(mod(iter-1,my_hook_print_interval) == 0), data)
	  if (hook_done) done = .true.
	endif

	call print("cg_n " // 0.0_dp // " " // E_i // " " // (g_i.dot.h_i) // " " // &
		    normsq(g_i) // " " // N_evals // " n_minim with new dir")

	call print("n_minim loop end, N_evals " // N_evals // " max_N_evals " // max_N_evals // &
	  " done " // done, PRINT_VERBOSE)

	iter = iter + 1

    end do

    if (present(hook)) then
      call hook(x_i,g_i,E_i,hook_done,.true.,data)
      if (hook_done) done = .true.
    endif

    final_E = E_i

end function n_minim

subroutine line_scan(x0, xdir, func, use_func, dfunc, data)
  real(dp)::x0(:)  !% Starting vector
  real(dp)::xdir(:)!% Search direction
  INTERFACE 
     function func(x,data)
       use system_module
       real(dp)::x(:)
       character(len=1),optional::data(:)
       real(dp)::func
     end function func
  END INTERFACE
  logical :: use_func
  INTERFACE 
     function dfunc(x,data)
       use system_module
       real(dp)::x(:)
       character(len=1),optional::data(:)
       real(dp)::dfunc(size(x))
     end function dfunc
  END INTERFACE
  character(len=1), optional::data(:)

  integer :: i
  real(dp) :: new_eps
  real(dp) :: fn, dirdx_new
  real(dp), allocatable :: xn(:), dxn(:)

  allocate(xn(size(x0)))
  allocate(dxn(size(x0)))

  fn = 0.0_dp
  call print('line scan:', PRINT_NORMAL)
  new_eps = 1.0e-5_dp
  do i=1,50
     xn = x0 + new_eps*xdir
#ifndef _OPENMP
     call verbosity_push_decrement()
#endif
     if (use_func) fn = func(xn,data)
     dxn = dfunc(xn,data)
#ifndef _OPENMP
     call verbosity_pop()
#endif
     dirdx_new = xdir .DOT. dxn

     call print('LINE_SCAN ' // new_eps//' '//fn// ' '//dirdx_new, PRINT_NORMAL)

     new_eps = new_eps*1.15
  enddo

  deallocate(xn)

end subroutine line_scan

end module minimization_module

