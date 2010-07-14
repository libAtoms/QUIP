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
!X Local Energy Mixing routines to be included in MetaPotential.f95
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine metapotential_local_e_mix_initialise(this, args_str, region1_pot, region2_pot, reference_bulk, mpi_obj)
    type(MetaPotential_Local_E_Mix), intent(inout) :: this
    character(len=*), intent(in) :: args_str
    type(Potential), intent(in), target :: region1_pot, region2_pot
    type(Atoms), optional, intent(inout) :: reference_bulk
    type(MPI_Context), intent(in), optional :: mpi_obj

    type(Dictionary) :: params
    logical :: minimise_bulk
    logical :: do_rescale_r, do_rescale_E, do_tb_defaults

    call initialise(params)
    call param_register(params, "r_scale", "1.0", this%r_scale_pot1)
    call param_register(params, "E_scale", "1.0", this%E_scale_pot1)
    call param_register(params, "terminate", "T", this%terminate)
    call param_register(params, "minimise_bulk", "F", minimise_bulk)

    call param_register(params, "do_rescale_r", "F", do_rescale_r)
    call param_register(params, "do_rescale_E", "F", do_rescale_E)
    call param_register(params, "do_tb_defaults", "T", do_tb_defaults)

    call param_register(params, "minimise_mm", "F", this%minimise_mm)
    call param_register(params, "minim_mm_method", "cg", this%minim_mm_method)
    call param_register(params, 'minim_mm_tol', '1e-6', this%minim_mm_tol)
    call param_register(params, 'minim_mm_eps_guess', '1e-4', this%minim_mm_eps_guess)
    call param_register(params, 'minim_mm_max_steps', '500', this%minim_mm_max_steps)
    call param_register(params, 'minim_mm_linminroutine', 'FAST_LINMIN', this%minim_mm_linminroutine)
    call param_register(params, 'minim_mm_do_pos', 'T', this%minim_mm_do_pos)
    call param_register(params, 'minim_mm_do_lat', 'F', this%minim_mm_do_lat)
    call param_register(params, 'minim_mm_do_print', 'F', this%minim_mm_do_print)
    call param_register(params, 'minim_mm_use_n_minim', 'F', this%minim_mm_use_n_minim)
    call param_register(params, 'minim_mm_args_str', '', this%minim_mm_args_str)

    if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='MetaPotential_Local_E_Mix_initialise args_str') ) &
      call system_abort("MetaPotential_Local_E_Mix_initialise failed to parse args_str='"//trim(args_str)//"'")
    call finalise(params)

    call initialise(this%create_hybrid_weights_params)
    call read_string(this%create_hybrid_weights_params, args_str)
    call remove_value(this%create_hybrid_weights_params, 'r_scale')
    call remove_value(this%create_hybrid_weights_params, 'E_scale')
    call remove_value(this%create_hybrid_weights_params, 'terminate')
    call remove_value(this%create_hybrid_weights_params, 'minimise_bulk')
    call remove_value(this%create_hybrid_weights_params, 'do_rescale_r')
    call remove_value(this%create_hybrid_weights_params, 'do_rescale_E')
    call remove_value(this%create_hybrid_weights_params, 'do_tb_defaults')
    call remove_value(this%create_hybrid_weights_params, 'minimise_mm')
    call remove_value(this%create_hybrid_weights_params, 'minim_mm_method')
    call remove_value(this%create_hybrid_weights_params, 'minim_mm_tol')
    call remove_value(this%create_hybrid_weights_params, 'minim_mm_eps_guess')
    call remove_value(this%create_hybrid_weights_params, 'minim_mm_max_steps')
    call remove_value(this%create_hybrid_weights_params, 'minim_mm_linminroutine')
    call remove_value(this%create_hybrid_weights_params, 'minim_mm_do_pos')
    call remove_value(this%create_hybrid_weights_params, 'minim_mm_do_lat')
    call remove_value(this%create_hybrid_weights_params, 'minim_mm_do_print')
    call remove_value(this%create_hybrid_weights_params, 'minim_mm_use_n_minim')
    call remove_value(this%create_hybrid_weights_params, 'minim_mm_args_str')

    if (this%minimise_mm) then
      call initialise(this%relax_metapot, "Simple", region2_pot)
    endif

    if (do_rescale_r .or. do_rescale_E .or. do_tb_defaults) then
      if (.not. present(reference_bulk)) &
	call system_abort("metapotential_local_e_mix_initialise got do_rescale_r="//do_rescale_r//" do_rescale_E="//do_rescale_E// &
	  " do_tb_defaults="//do_tb_defaults//" but reference_bulk is not present")

      call do_reference_bulk(reference_bulk, region1_pot, region2_pot, minimise_bulk, do_rescale_r, do_rescale_E, &
	this%r_scale_pot1, this%E_scale_pot1, do_tb_defaults)

    endif

    if (this%r_scale_pot1 <= 0.0_dp) this%r_scale_pot1 = 1.0_dp
    if (this%E_scale_pot1 <= 0.0_dp) this%E_scale_pot1 = 1.0_dp

    call print ("Rescaling positions in region1 potential by " // this%r_scale_pot1 // " to match lattice constants")
    call print ("Rescaling energies in region1 potential by " // this%E_scale_pot1 // " to match bulk modulus")

    this%pot_region1 => region1_pot
    this%pot_region2 => region2_pot

  end subroutine

  subroutine metapotential_local_e_mix_finalise(this)
    type(MetaPotential_Local_E_Mix), intent(inout) :: this

    call finalise(this%pot_region1)
    call finalise(this%pot_region2)

    deallocate(this%pot_region1)
    deallocate(this%pot_region2)
  end subroutine metapotential_local_e_mix_finalise


  function metapotential_local_e_mix_cutoff(this)
    type(MetaPotential_Local_E_Mix), intent(in) :: this
    real(dp) :: metapotential_local_e_mix_cutoff

    metapotential_local_e_mix_cutoff = max(cutoff(this%pot_region1), cutoff(this%pot_region2))
  end function metapotential_local_e_mix_cutoff

  subroutine metapotential_local_e_mix_print(this, file)
    type(MetaPotential_Local_E_Mix),          intent(inout) :: this
    type(Inoutput),intent(inout),optional:: file
    
    call print('Local_E_Mix metapotential:', file=file)
    call print('Potential in region 1 (embedded):', file=file)
    call print('=======================', file=file)
    call print(this%pot_region1, file=file)
    call print('')
    call print('Potential in region 2 (surroundings):', file=file)
    call print('=======================', file=file)
    call print(this%pot_region2, file=file)
    call print('do_terminate ' // this%terminate, file=file)
    call print('r_scale_pot1=' // this%r_scale_pot1 // ' E_scale_pot1=' // this%E_scale_pot1, file=file)
    call print('')
  end subroutine metapotential_local_e_mix_print

  subroutine metapotential_local_e_mix_calc(this, at, energy, local_e, force, virial, args_str, err)
    type(MetaPotential_Local_E_Mix), intent(inout) :: this
    type(Atoms), intent(inout) :: at
    real(dp), intent(out), optional :: energy
    real(dp), intent(out), optional :: local_e(:)
    real(dp), intent(out), optional :: force(:,:)
    real(dp), intent(out), optional :: virial(3,3)
    character(len=*), intent(in), optional :: args_str
    integer, intent(out), optional :: err

    logical :: calc_weights
    integer :: core_hops, transition_hops, buffer_hops
    type(Dictionary) :: params, calc_create_hybrid_weights_params
    type(Table) :: region1_table
    integer :: i
    integer, pointer :: hybrid(:), hybrid_mark(:)

    if (.not. associated(this%pot_region1) .or. .not. associated(this%pot_region2)) &
      call system_abort("MetaPotential_Local_E_Mix_calc: this%pot_region1 or this%pot_region2 not initialised")

    call initialise(params)
    call param_register(params, "calc_weights", "F", calc_weights)
    call param_register(params, "core_hops", "0", core_hops)
    if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='MetaPotential_Local_E_Mix_Calc args_str') ) &
      call system_abort("MetaPotential_Local_E_Mix_calc_energy failed to parse args_str='"//trim(args_str)//"'")
    call finalise(params)

    if (calc_weights) then
      call print("MetaPotential_Local_E_Mix_calc got calc_weights core_hops " // core_hops, PRINT_VERBOSE)
      call add_property(at, "weight_region1", 0.0_dp)
      call add_property(at, "hybrid_mark", HYBRID_NO_MARK)
      if (.not. assign_pointer(at, "hybrid", hybrid)) &
        call system_abort("QC_QUIP_calc at doesn't have hybrid property and calc_weights was specified")
      if (.not. assign_pointer(at, "hybrid_mark", hybrid_mark)) &
        call system_abort("QC_QUIP_calc Failed to add hybrid_mark property to at")

      call calc_connect(at)

      call wipe(region1_table)
      hybrid_mark = HYBRID_NO_MARK
      do i=1, size(hybrid)
        if (hybrid(i) /= 0) call append(region1_table, (/ i, 0, 0, 0 /) )
      end do
      call bfs_grow(at, region1_table, core_hops, nneighb_only = .false.)
      hybrid_mark(int_part(region1_table,1)) = HYBRID_ACTIVE_MARK

      call read_string(calc_create_hybrid_weights_params, write_string(this%create_hybrid_weights_params))
      call read_string(calc_create_hybrid_weights_params, trim(args_str), append=.true.)
      call create_hybrid_weights(at, write_string(calc_create_hybrid_weights_params))
    endif

    call calc_local_energy_mix(this, at, energy, local_e, force, virial, args_str)

    if (present(err)) err = 0

  end subroutine metapotential_local_e_mix_calc

  subroutine calc_local_energy_mix(this, at, energy, local_e, force, virial, args_str)
    type(MetaPotential_Local_E_Mix), intent(inout) :: this
    type(Atoms), intent(inout) :: at
    real(dp), intent(out), optional :: energy
    real(dp), intent(out), optional :: local_e(:)
    real(dp), intent(out), optional :: force(:,:)
    real(dp), intent(out), optional :: virial(3,3)
    character(len=*), intent(in), optional :: args_str

    ! local variables
    integer :: i
    integer,  pointer :: hybrid_mark(:), index(:), termindex(:)
    real(dp), pointer :: weight(:), weight_region1(:), rescale(:), c_weight(:)
    real(dp), allocatable :: f_cluster(:,:), local_e_cluster(:)
    real(dp), allocatable :: weight_saved(:), local_e_weight(:), local_e_region2(:), local_e_region1(:), f_region1(:,:)
    real(dp) :: e_cluster
    logical :: dummy
    type(Atoms) :: cluster
    type(Dictionary) :: params
    character(FIELD_LENGTH) :: cc_args_str

    real(dp), pointer :: local_e_pot1(:), local_e_pot2(:)

    allocate(weight_saved(at%N))
    allocate(local_e_weight(at%N))
    allocate(local_e_region2(at%N))
    allocate(local_e_region1(at%N))
    allocate(f_region1(3,at%N))

    call system_timer("calc_local_energy_mix")

    if(.not. (associated(this%pot_region1) .and. associated(this%pot_region2))) then
       call system_abort('local_e_mix_calc_energy_mix: potential for region 1 or region 2 is not initialised')
    end if

    ! check for a compatible weight property, create it if necessary
    if (.not.assign_pointer(at, 'weight', weight)) then
      call add_property(at, 'weight', 1.0_dp)
      if(.not.assign_pointer(at, 'weight', weight)) &
        call system_abort("Impossible failure to create weight property in calc_local_energy_mix")
      weight_saved = 1.0_dp
    else
      weight_saved = weight
    endif
    ! check for a compatible hybrid and weight_region1 property. they must be present
    if (.not.assign_pointer(at, 'hybrid_mark', hybrid_mark)) &
       call system_abort('hybrid_calc_energy_mix: atoms structure has no "hybrid_mark" property')
    if (.not.assign_pointer(at, 'weight_region1', weight_region1)) &
       call system_abort('hybrid_calc_energy_mix: atoms structure has no "weight_region1" property')

    if(.not. all(hybrid_mark == HYBRID_NO_MARK)) then
      if(present(virial)) then
         virial = 0.0_dp
         call system_abort('hybrid_calc_energy_mix: virials are not yet implemented when QM region is active')
      end if
    endif

    ! if hybrid is active, carve cluster and do minimisation if necessary
    if(any(hybrid_mark /= HYBRID_NO_MARK)) then
      call initialise(params)
      call set_value(params, 'terminate', this%terminate)
      call set_value(params, 'cluster_allow_modification', 'F')
      call set_value(params, 'cluster_periodic_x', 'T')
      call set_value(params, 'cluster_periodic_y', 'T')
      call set_value(params, 'cluster_periodic_z', 'T')
      call set_value(params, 'randomise_buffer', 'F')
      cc_args_str = write_string(params)
      call finalise(params)
      cluster =  create_cluster_hybrid_mark(at, cc_args_str)

      if (this%minimise_mm) then
	dummy = assign_pointer(cluster, 'index', index)
	call do_minimise_mm(this%relax_metapot, at, this%minim_mm_method, this%minim_mm_tol, this%minim_mm_max_steps, &
	  this%minim_mm_linminroutine, this%minim_mm_do_pos, this%minim_mm_do_lat, this%minim_mm_do_print, &
	  this%minim_mm_args_str, this%minim_mm_eps_guess, this%minim_mm_use_n_minim, this%minim_inoutput_movie, &
          this%minim_cinoutput_movie, index)
      endif
    endif

    ! do calculation in region 2 first
    weight = weight_saved-weight_region1 ! "negate" the weights

    call system_timer("calc_local_energy_mix/calc_region_2")
    call Calc(this%pot_region2, at, energy, local_e_region2, force, args_str=args_str)
    call system_timer("calc_local_energy_mix/calc_region_2")

    if (assign_pointer(at, "local_e_pot2", local_e_pot2)) local_e_pot2 = local_e_region2

    call system_timer("calc_local_energy_mix/prep_cluster")
    if(present(local_e)) then
      local_e_weight = 1.0_dp-weight_region1
      local_e = local_e_weight * local_e_region2
    endif

    ! if there are no marked atoms, we are done, so return
    if(all(hybrid_mark == HYBRID_NO_MARK)) then
       ! restore weights
       weight = weight_saved
       call print('No REGION 1 marks, not carving REGION 1 cluster')
       return
    end if

    ! Now compute region 1

    ! copy weight_region1 back into cluster for region1 calculation
    if (.not.assign_pointer(cluster, 'weight', c_weight)) &
      call system_abort("calc_local_energy_mix didn't find weight property for cluster")
    if (.not. assign_pointer(cluster, 'index', index)) &
      call system_abort("calc_local_energy_mix didn't find index property for cluster")
    c_weight = weight_region1(index)

    ! rescale cluster
    call set_lattice(cluster, this%r_scale_pot1 * cluster%lattice, scale_positions=.true.)

    call calc_connect(cluster)

    call system_timer("calc_local_energy_mix/prep_cluster")
    call system_timer("calc_local_energy_mix/calc_region_1")
    allocate(f_cluster(3,cluster%N))
    allocate(local_e_cluster(cluster%N))
    if (present(force)) then
      call Calc(this%pot_region1, cluster, e_cluster, local_e_cluster, f_cluster, args_str = trim(args_str) // " solver=DIAG_GF SCF=GLOBAL_U")
      ! scale back forces
      f_cluster = f_cluster*this%E_scale_pot1*this%r_scale_pot1
    else
      call Calc(this%pot_region1, cluster, e_cluster, local_e_cluster, args_str = trim(args_str) // " solver=DIAG_GF SCF=GLOBAL_U")
    endif

    e_cluster = e_cluster*this%E_scale_pot1
    local_e_cluster = local_e_cluster*this%E_scale_pot1
    call system_timer("calc_local_energy_mix/calc_region_1")

    call system_timer("calc_local_energy_mix/combine")
    if(current_verbosity() > PRINT_VERBOSE) &
         call print_cfg(cluster, "cluster.cfg", all_properties = .true.)

    ! redistribute local energies and weights from the cluster to region 1
    if (.not. assign_pointer(cluster, 'index', index)) &
       call system_abort('calc_local_energy_mix: cluster structure has no "index" property')
    if (.not. assign_pointer(cluster, 'rescale', rescale)) &
       call system_abort('calc_local_energy_mix: cluster structure has no "rescale" property')
    if (.not. assign_pointer(cluster, 'termindex', termindex)) &
       call system_abort('calc_local_energy_mix: cluster structure has no "termindex" property')
    ! first add terms arising from the fact that the position of the rescaled termination atoms
    ! depends on the position of the atoms they are terminating
    local_e_region1 = 0.0_dp
    f_region1 = 0.0_dp
    if(present(force)) then
       do i=1,cluster%N
          if(termindex(i) > 0) &
               f_cluster(:,termindex(i)) = f_cluster(:,termindex(i)) + (1.0_dp-rescale(i))*f_cluster(:,i)
       end do
    end if
    ! now do the redistribution
    do i=1,cluster%N
       if(present(force)) f_region1(:,index(i)) = f_region1(:,index(i)) + f_cluster(:,i)*rescale(i)
       local_e_region1(index(i)) = local_e_region1(index(i)) + local_e_cluster(i)
    end do

    if (assign_pointer(at, "local_e_pot1", local_e_pot1)) local_e_pot1 = local_e_region1

    ! finally add region 1 stuff to the output variables
    if(present(energy)) energy = energy + e_cluster
    if(present(local_e)) then
      local_e_weight = weight_region1
      local_e = local_e + local_e_weight * local_e_region1
    endif
    if(present(force)) force = force + f_region1

    weight = weight_saved

    deallocate(f_cluster)
    deallocate(local_e_cluster)
    call finalise(cluster)

    deallocate(weight_saved)
    deallocate(local_e_weight)
    deallocate(local_e_region2)
    deallocate(local_e_region1)
    deallocate(f_region1)

    call system_timer("calc_local_energy_mix/combine")
    call system_timer("calc_local_energy_mix")
  end subroutine calc_local_energy_mix
