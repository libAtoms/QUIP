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
!X misc routines to be included in MetaPotential.f95
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine do_reference_bulk(reference_bulk, region1_pot, region2_pot, minimise_bulk, do_rescale_r, do_rescale_E, &
			       r_scale_pot1, E_scale_pot1, do_tb_defaults)
    type(Atoms), intent(inout) :: reference_bulk
    type(Potential), intent(in), target :: region1_pot, region2_pot
    logical, intent(in) :: minimise_bulk, do_rescale_r, do_rescale_E, do_tb_defaults
    real(dp), intent(inout) :: r_scale_pot1, E_scale_pot1

    real(dp) :: e_bulk !, B_region2, B_region1

    type(Atoms) :: bulk_region1, bulk_region2
    type(MetaPotential) :: region1_metapot, region2_metapot
    integer :: it

    call verbosity_push_decrement(PRINT_NERD)
    call print("metapotential_local_e_mix_initialise called with reference_bulk")
    call print_xyz(reference_bulk, "stdout")
    call verbosity_pop()

    ! set up region1 stuff
    bulk_region1 = reference_bulk
    call calc_connect(bulk_region1)
    call initialise(region1_metapot, "Simple", region1_pot)
    if (minimise_bulk) then
      call print("MINIMISING bulk in region1 potential", PRINT_VERBOSE)
      call verbosity_push_decrement()
      it = minim(region1_metapot, at=bulk_region1, method='cg', convergence_tol=0.001_dp, max_steps=100, linminroutine='NR_LINMIN', &
	do_print = .false., do_lat = .true., args_str='solver=DIAG SCF_NONE')
      call verbosity_pop()
    endif
    call calc(region1_metapot, bulk_region1, e=e_bulk, args_str='solver=DIAG SCF_NONE')
    call print("region1 Potential energy " // e_bulk)
    call print("region1 Potential lattice")
    call print(bulk_region1%lattice)

    if (do_tb_defaults .and. associated(region1_metapot%pot%tb)) then
      region1_metapot%pot%tb%tbsys%tbmodel%default_fermi_E = region1_metapot%pot%tb%fermi_E
      region1_metapot%pot%tb%tbsys%tbmodel%has_default_fermi_E = .true.
      region1_metapot%pot%tb%tbsys%tbmodel%default_band_width = (region1_metapot%pot%tb%fermi_E - minval(TB_evals(region1_metapot%pot%tb)))*1.2_dp
      region1_metapot%pot%tb%tbsys%tbmodel%has_default_band_width = .true.
      region1_metapot%pot%tb%tbsys%tbmodel%default_fermi_T = 0.2_dp
      region1_metapot%pot%tb%tbsys%tbmodel%has_default_fermi_T = .true.

      region1_metapot%pot%tb%tbsys%scf%conv_tol = 1e-8_dp
    endif

    ! set up region2 stuff
    bulk_region2 = reference_bulk
    call initialise(region2_metapot, "Simple", region2_pot)
    call set_cutoff(bulk_region2, cutoff(region2_metapot)+0.5_dp)
    call calc_connect(bulk_region2)
    if (minimise_bulk) then
      call print("MINIMISING bulk in region2 potential", PRINT_VERBOSE)
      call verbosity_push_decrement()
      it = minim(region2_metapot, at=bulk_region2, method='cg', convergence_tol=0.001_dp, max_steps=100, linminroutine='NR_LINMIN', &
	do_print = .false., do_lat = .true.)
      call verbosity_pop()
    end if
    call calc(region2_metapot, bulk_region2, e=e_bulk)
    call print("region2 Potential Bulk energy = " // e_bulk)
    call print("region2 Potential lattice")
    call print(bulk_region2%lattice)

    if (do_tb_defaults .and. associated(region2_metapot%pot%tb)) then
      region2_metapot%pot%tb%tbsys%tbmodel%default_fermi_E = region2_metapot%pot%tb%fermi_E
      region2_metapot%pot%tb%tbsys%tbmodel%has_default_fermi_E = .true.
      region2_metapot%pot%tb%tbsys%tbmodel%default_band_width = (region2_metapot%pot%tb%fermi_E - minval(TB_evals(region2_metapot%pot%tb)))*1.2_dp
      region2_metapot%pot%tb%tbsys%tbmodel%has_default_band_width = .true.
      region2_metapot%pot%tb%tbsys%tbmodel%default_fermi_T = 0.2_dp
      region2_metapot%pot%tb%tbsys%tbmodel%has_default_fermi_T = .true.

      region2_metapot%pot%tb%tbsys%scf%conv_tol = 1e-8_dp
    endif

    if (do_rescale_r) then
      ! compute scale, assuming that lattices are only off by a uniform scale
      r_scale_pot1 = maxval(abs(bulk_region2%lattice(:,1)))/maxval(abs(bulk_region1%lattice(:,1)))
    endif

    if (do_rescale_E) then
      call system_abort("do_reference_bulk: no energy rescaling implemented yet")
      ! B_region2 = bulk_modulus(region2_metapot, bulk_region2)
      ! B_region1 = bulk_modulus(region1_metapot, bulk_region1, r_scale_pot1)
      ! E_scale_pot1 = B_region2/B_region1
    end if

    call finalise(bulk_region1)
    call finalise(bulk_region2)
    call finalise(region1_metapot)
    call finalise(region2_metapot)
  end subroutine do_reference_bulk

  subroutine do_minimise_mm(relax_metapot, at, minim_mm_method, minim_mm_tol, minim_mm_max_steps, &
      minim_mm_linminroutine, minim_mm_do_pos, minim_mm_do_lat, minim_mm_do_print, & 
      minim_mm_args_str, minim_mm_eps_guess, minim_mm_use_n_minim, minim_inoutput_movie,  &
      minim_cinoutput_movie, &
      constrained_list)
    type(MetaPotential), intent(inout) :: relax_metapot
    type(Atoms), intent(inout) :: at
    character(len=*), intent(in) :: minim_mm_method
    real(dp), intent(in) :: minim_mm_tol
    integer, intent(in) :: minim_mm_max_steps
    character(len=*), intent(in) :: minim_mm_linminroutine
    logical, intent(in) :: minim_mm_do_pos, minim_mm_do_lat, minim_mm_do_print
    character(len=*), intent(in) :: minim_mm_args_str
    real(dp), intent(in) :: minim_mm_eps_guess
    logical, intent(in) :: minim_mm_use_n_minim
    type(Inoutput), intent(inout) :: minim_inoutput_movie
    type(CInoutput), intent(inout) ::  minim_cinoutput_movie
    integer, intent(in) :: constrained_list(:)

    integer, pointer :: fixed_metapot(:)
    integer :: mm_steps

    if (.not. assign_pointer(at, "fixed_metapot", fixed_metapot)) then
      call add_property(at, "fixed_metapot", 0)
      if (.not. assign_pointer(at, "fixed_metapot", fixed_metapot)) &
	call system_abort("do_minimise_mm failed to assign pointer for property 'fixed_metapot'")
    endif

    if (minval(constrained_list) < 1 .or. maxval(constrained_list) > size(fixed_metapot)) then
      call print("do_minimise_mm: size(fixed_metapot) " // size(fixed_metapot) // &
	" size(constrained_list) " // size(constrained_list), PRINT_ALWAYS)
      call print("do_minimise_mm: constrained_list", PRINT_ALWAYS)
      call print(constrained_list, PRINT_ALWAYS)
      call system_abort("do_minimise_mm with invalid value in constrained list")
    endif

    fixed_metapot = 0
    fixed_metapot(constrained_list) = 1

    call verbosity_push_decrement()
    mainlog%prefix="NESTED_MINIM"
    mm_steps = minim(relax_metapot, at, minim_mm_method, minim_mm_tol, minim_mm_max_steps, &
	linminroutine=minim_mm_linminroutine, do_pos=minim_mm_do_pos, do_lat=minim_mm_do_lat, &
	do_print=minim_mm_do_print, args_str=minim_mm_args_str, &
	eps_guess=minim_mm_eps_guess, use_n_minim=minim_mm_use_n_minim, print_inoutput=minim_inoutput_movie, &
        print_cinoutput=minim_cinoutput_movie)
    mainlog%prefix=""
    call verbosity_pop()
    call print("MetaPotential_local_e_mix_calc minimise_mm got mm_steps " // mm_steps, PRINT_VERBOSE)

    call calc_connect(at)

    call print("MM minimisation done in "//mm_steps//" steps", PRINT_VERBOSE)

    fixed_metapot=0

  end subroutine do_minimise_mm

