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

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X misc routines to be included in Potential.f95
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  
 subroutine linear_square(x,afunc)
    real(dp) :: x, afunc(:)

    afunc(1) = x*x
    afunc(2) = x
    afunc(3) = 1
  
  end subroutine linear_square

  subroutine potential_bulk_modulus(pot, at, B, V0, minimise_bulk, eps, args_str)
    type(Potential), intent(inout) :: pot
    type(Atoms), intent(in) :: at
    real(dp), intent(out) :: V0, B
    logical, intent(in), optional :: minimise_bulk
    real(dp), intent(in), optional :: eps
    character(len=*), intent(in), optional :: args_str

    integer  :: i,j, steps 
    real(dp) :: E(11), V(11), a(3), chisq, use_eps
    type(Atoms) :: at0, at1
    logical :: do_minimise_bulk

    do_minimise_bulk = optional_default(.true., minimise_bulk)
    use_eps = optional_default(2.0e-3_dp, eps)

    at0 = at
    if (do_minimise_bulk) then
       ! Minimise initial cell
       call verbosity_push_decrement()
       steps = minim(pot, at0, 'cg', 1e-6_dp, 100, &
            'LINMIN_DERIV', do_pos=.true.,do_lat=.true., do_print=.false., args_str=args_str)
       call verbosity_pop()
    end if

    do i = 1, 11 
       at1 = at0
       call set_lattice(at1, at1%lattice*(1.0_dp - real(i-6,dp)*use_eps), scale_positions=.true.)  
       V(i) = cell_volume(at1)
       call set_cutoff (at1, cutoff(pot)+1.0_dp)
       call calc_connect(at1)

       if (do_minimise_bulk) then
          ! Minimise positions with lattice fixed
          call verbosity_push_decrement()
          steps = minim(pot, at1, 'cg', 1e-6_dp, 100, &
               'LINMIN_DERIV', do_pos=.true.,do_lat=.false., do_print=.false., args_str=args_str)
          call verbosity_pop()
       end if

       call calc(pot, at1, energy=E(i), args_str=args_str)
       call print('bulk_modulus V='//V(i)//' E='//E(i), PRINT_NERD)
    end do

    call least_squares(V, E, (/ ( 1.0_dp, j=1,11) /), &
         a, chisq, linear_square)
 
    V0 = -a(2)/(2*a(1))
    B = -1.0*a(2)*GPA

  end subroutine potential_bulk_modulus
 
!% We define a potential energy function $E(r)$, and then a scaled function
!% \[
!% E\'(r) = \beta E(\alpha r)
!% \]
!% The corresponding force in the original coordinate system is
!% \[
!% F\'(r) = -\frac{\partial E\'}{\partial r} = -\beta \alpha \frac{\partial E\'}{\partial r} = \beta \alpha\, F
!% \]
!% The equilibrium lattice constant changes from $a_0$ to $a_0\'$ and the equilbrium cell volume changes from $V_0$ to $V_0\'$ according to
!% \begin{eqnarray*}
!% a_0' & = & \frac{a_0}{\alpha} \\ V_0' & = & \frac{V_0}{\alpha^3}
!% \end{eqnarray*}
!% The scaled bulk modulus is 
!% \begin{eqnarray*}
!% B' & = & V \frac{\partial^2 E'}{\partial V^2} \\ & = & \beta \alpha^3 V \frac{\partial^2E}{\partial V^2} \\ & = & \beta \alpha^3 B
!% \end{eqnarray*}
!% Thus if we want to match a target volume $V_0\'$ and bulk modulus $B\'$ we should use
!% \begin{eqnarray*}
!%   \alpha & = & \left( \frac{V_0}{V_0'} \right)^\frac{1}{3} = \frac{a_0}{a_0'} \\   \beta  & = & \frac{B'}{B \alpha^3}
!% \end{eqnarray*}
!% where $a_0$ and $a_0\'$ are the lattice constants before and after rescaling.
!% 
!% For QM/MM force mixing, where we label the QM potential as region 1 and the classical (MM) potential 
!% as region 2, the aim is to rescale the QM potential to match the MM lattice constant $a_2$ and bulk 
!% modulus $B_2$, so we have
!% \begin{eqnarray*}
!%   \alpha & = & \frac{a_1}{a_2} \\ \beta  & = & \frac{B_2}{B_1 \alpha^3} 
!% \end{eqnarray*}
!% where $a_1$ and $B_1$ are the unmodified QM lattice constant and bulk modulus, respectively.
!% Note that this is equivalent to rescaling the MM \emph{positions} to match the QM lattice constant.

  subroutine do_reference_bulk(reference_bulk, region1_pot, region2_pot, minimise_bulk, do_rescale_r, do_rescale_E, &
			       r_scale_pot1, E_scale_pot1, do_tb_defaults)
    type(Atoms), intent(inout) :: reference_bulk
    type(Potential), intent(inout), target :: region1_pot, region2_pot
    logical, intent(in) :: minimise_bulk, do_rescale_r, do_rescale_E, do_tb_defaults
    real(dp), intent(inout) :: r_scale_pot1, E_scale_pot1

    real(dp) :: e_bulk, B_region2, B_region1, vol

    type(Atoms) :: bulk_region1, bulk_region2
    integer :: it

    call verbosity_push_decrement(PRINT_NERD)
    call print("potential_local_e_mix_initialise called with reference_bulk")
    call write(reference_bulk, "stdout")
    call verbosity_pop()

    ! set up region1 stuff
    bulk_region1 = reference_bulk
    call set_cutoff(bulk_region1, cutoff(region1_pot)+0.5_dp)
    call calc_connect(bulk_region1)
!    call initialise(region1_pot, "Simple", region1_pot)
    if (minimise_bulk) then
      call print("MINIMISING bulk in region1 potential", PRINT_VERBOSE)
      call verbosity_push_decrement()
      it = minim(region1_pot, at=bulk_region1, method='cg', convergence_tol=0.001_dp, max_steps=100, linminroutine='NR_LINMIN', &
	do_print = .false., do_lat = .true., args_str='solver=DIAG SCF_NONE')
      call verbosity_pop()
    endif
    call calc(region1_pot, bulk_region1, energy=e_bulk, args_str='solver=DIAG SCF_NONE')
    call print("region1 Potential energy " // e_bulk)
    call print("region1 Potential lattice")
    call print(bulk_region1%lattice)

#ifdef HAVE_TB
    if (do_tb_defaults .and. associated(region1_pot%simple%tb)) then
      region1_pot%simple%tb%tbsys%tbmodel%default_fermi_E = region1_pot%simple%tb%fermi_E
      region1_pot%simple%tb%tbsys%tbmodel%has_default_fermi_E = .true.
      region1_pot%simple%tb%tbsys%tbmodel%default_band_width = (region1_pot%simple%tb%fermi_E - minval(TB_evals(region1_pot%simple%tb)))*1.2_dp
      region1_pot%simple%tb%tbsys%tbmodel%has_default_band_width = .true.
      region1_pot%simple%tb%tbsys%tbmodel%default_fermi_T = 0.2_dp
      region1_pot%simple%tb%tbsys%tbmodel%has_default_fermi_T = .true.

      region1_pot%simple%tb%tbsys%scf%conv_tol = 1e-8_dp
    endif
#endif

    ! set up region2 stuff
    bulk_region2 = reference_bulk
!    call initialise(region2_pot, "Simple", region2_pot)
    call set_cutoff(bulk_region2, cutoff(region2_pot)+0.5_dp)
    call calc_connect(bulk_region2)
    if (minimise_bulk) then
      call print("MINIMISING bulk in region2 potential", PRINT_VERBOSE)
      call verbosity_push_decrement()
      it = minim(region2_pot, at=bulk_region2, method='cg', convergence_tol=0.001_dp, max_steps=100, linminroutine='NR_LINMIN', &
	do_print = .false., do_lat = .true.)
      call verbosity_pop()
    end if
    call calc(region2_pot, bulk_region2, energy=e_bulk)
    call print("region2 Potential Bulk energy = " // e_bulk)
    call print("region2 Potential lattice")
    call print(bulk_region2%lattice)

#ifdef HAVE_TB
    if (do_tb_defaults .and. associated(region2_pot%simple%tb)) then
      region2_pot%simple%tb%tbsys%tbmodel%default_fermi_E = region2_pot%simple%tb%fermi_E
      region2_pot%simple%tb%tbsys%tbmodel%has_default_fermi_E = .true.
      region2_pot%simple%tb%tbsys%tbmodel%default_band_width = (region2_pot%simple%tb%fermi_E - minval(TB_evals(region2_pot%simple%tb)))*1.2_dp
      region2_pot%simple%tb%tbsys%tbmodel%has_default_band_width = .true.
      region2_pot%simple%tb%tbsys%tbmodel%default_fermi_T = 0.2_dp
      region2_pot%simple%tb%tbsys%tbmodel%has_default_fermi_T = .true.

      region2_pot%simple%tb%tbsys%scf%conv_tol = 1e-8_dp
    endif
#endif
    r_scale_pot1 = 1.0_dp
    if (do_rescale_r) then
      ! compute scale, assuming that lattices are only off by a uniform scale
      r_scale_pot1 = maxval(abs(bulk_region1%lattice(:,1)))/maxval(abs(bulk_region2%lattice(:,1)))
    endif

    E_scale_pot1 = 1.0_dp
    if (do_rescale_E) then
      call bulk_modulus(region2_pot, bulk_region2, B_region2, vol, minimise_bulk)
      call bulk_modulus(region1_pot, bulk_region1, B_region1, vol, minimise_bulk)

      call print('Bulk modulus in region 1 = '//B_region1//' GPa')
      call print('Bulk modulus in region 2 = '//B_region2//' GPa')

      E_scale_pot1 = B_region2/(B_region1*r_scale_pot1**3)
    end if

    call finalise(bulk_region1)
    call finalise(bulk_region2)
  end subroutine do_reference_bulk

  subroutine do_minimise_mm(relax_pot, at, minim_mm_method, minim_mm_tol, minim_mm_max_steps, &
      minim_mm_linminroutine, minim_mm_do_pos, minim_mm_do_lat, minim_mm_do_print, & 
      minim_mm_args_str, minim_mm_eps_guess, minim_inoutput_movie,  &
      minim_cinoutput_movie, &
      constrained_list)
    type(Potential), intent(inout) :: relax_pot
    type(Atoms), intent(inout) :: at
    character(len=*), intent(in) :: minim_mm_method
    real(dp), intent(in) :: minim_mm_tol
    integer, intent(in) :: minim_mm_max_steps
    character(len=*), intent(in) :: minim_mm_linminroutine
    logical, intent(in) :: minim_mm_do_pos, minim_mm_do_lat, minim_mm_do_print
    character(len=*), intent(in) :: minim_mm_args_str
    real(dp), intent(in) :: minim_mm_eps_guess
    type(Inoutput), intent(inout) :: minim_inoutput_movie
    type(CInoutput), intent(inout) ::  minim_cinoutput_movie
    integer, intent(in) :: constrained_list(:)

    integer, pointer :: fixed_pot(:)
    integer :: mm_steps

    if (.not. assign_pointer(at, "fixed_pot", fixed_pot)) then
      call add_property(at, "fixed_pot", 0)
      if (.not. assign_pointer(at, "fixed_pot", fixed_pot)) &
	call system_abort("do_minimise_mm failed to assign pointer for property 'fixed_pot'")
    endif

    if (minval(constrained_list) < 1 .or. maxval(constrained_list) > size(fixed_pot)) then
      call print("do_minimise_mm: size(fixed_pot) " // size(fixed_pot) // &
	" size(constrained_list) " // size(constrained_list), PRINT_ALWAYS)
      call print("do_minimise_mm: constrained_list", PRINT_ALWAYS)
      call print(constrained_list, PRINT_ALWAYS)
      call system_abort("do_minimise_mm with invalid value in constrained list")
    endif

    fixed_pot = 0
    fixed_pot(constrained_list) = 1

    call verbosity_push_decrement()
    mainlog%prefix="NESTED_MINIM"
    mm_steps = minim(relax_pot, at, minim_mm_method, minim_mm_tol, minim_mm_max_steps, &
	linminroutine=minim_mm_linminroutine, do_pos=minim_mm_do_pos, do_lat=minim_mm_do_lat, &
	do_print=minim_mm_do_print, args_str=minim_mm_args_str, &
	eps_guess=minim_mm_eps_guess, print_inoutput=minim_inoutput_movie, &
        print_cinoutput=minim_cinoutput_movie)
    mainlog%prefix=""
    call verbosity_pop()
    call print("Potential_local_e_mix_calc minimise_mm got mm_steps " // mm_steps, PRINT_VERBOSE)

    call calc_connect(at)

    call print("MM minimisation done in "//mm_steps//" steps", PRINT_VERBOSE)

    fixed_pot=0

  end subroutine do_minimise_mm

