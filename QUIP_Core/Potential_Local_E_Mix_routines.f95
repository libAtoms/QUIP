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
!X Local Energy Mixing routines to be included in Potential.f95
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  recursive subroutine potential_local_e_mix_initialise(this, args_str, region1_pot, region2_pot, reference_bulk, mpi_obj, error)
    type(Potential_Local_E_Mix), intent(inout) :: this
    character(len=*), intent(in) :: args_str
    type(Potential), intent(inout), target :: region1_pot, region2_pot
    type(Atoms), optional, intent(inout) :: reference_bulk
    type(MPI_Context), intent(in), optional :: mpi_obj
    integer, intent(out), optional :: error

    type(Dictionary) :: params
    logical :: minimise_bulk
    logical :: do_rescale_r, do_rescale_E, do_tb_defaults

    INIT_ERROR(error)

    call initialise(params)
    call param_register(params, "r_scale", "1.0", this%r_scale_pot1, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, "E_scale", "1.0", this%E_scale_pot1, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, "terminate", "T", this%terminate, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, "minimise_bulk", "F", minimise_bulk, help_string="No help yet.  This source file was $LastChangedBy$")

    call param_register(params, "do_rescale_r", "F", do_rescale_r, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, "do_rescale_E", "F", do_rescale_E, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, "do_tb_defaults", "T", do_tb_defaults, help_string="No help yet.  This source file was $LastChangedBy$")

    call param_register(params, "minimise_mm", "F", this%minimise_mm, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, "minim_mm_method", "cg", this%minim_mm_method, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'minim_mm_tol', '1e-6', this%minim_mm_tol, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'minim_mm_eps_guess', '1e-4', this%minim_mm_eps_guess, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'minim_mm_max_steps', '500', this%minim_mm_max_steps, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'minim_mm_linminroutine', 'FAST_LINMIN', this%minim_mm_linminroutine, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'minim_mm_do_pos', 'T', this%minim_mm_do_pos, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'minim_mm_do_lat', 'F', this%minim_mm_do_lat, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'minim_mm_do_print', 'F', this%minim_mm_do_print, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'minim_mm_use_n_minim', 'F', this%minim_mm_use_n_minim, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'minim_mm_args_str', '', this%minim_mm_args_str, help_string="No help yet.  This source file was $LastChangedBy$")

    if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='Potential_Local_E_Mix_initialise args_str') ) then
      RAISE_ERROR("Potential_Local_E_Mix_initialise failed to parse args_str='"//trim(args_str)//"'", error)
    endif
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
      call initialise(this%relax_pot, "Simple", region2_pot)
    endif

    if (do_rescale_r .or. do_rescale_E .or. do_tb_defaults) then
      if (.not. present(reference_bulk)) then
	RAISE_ERROR("potential_local_e_mix_initialise got do_rescale_r="//do_rescale_r//" do_rescale_E="//do_rescale_E//" do_tb_defaults="//do_tb_defaults//" but reference_bulk is not present", error)
      endif

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

  recursive subroutine potential_local_e_mix_finalise(this)
    type(Potential_Local_E_Mix), intent(inout) :: this

    call finalise(this%pot_region1)
    call finalise(this%pot_region2)

    deallocate(this%pot_region1)
    deallocate(this%pot_region2)
  end subroutine potential_local_e_mix_finalise


  recursive function potential_local_e_mix_cutoff(this)
    type(Potential_Local_E_Mix), intent(in) :: this
    real(dp) :: potential_local_e_mix_cutoff

    potential_local_e_mix_cutoff = max(cutoff(this%pot_region1), cutoff(this%pot_region2))
  end function potential_local_e_mix_cutoff

  recursive subroutine potential_local_e_mix_print(this, file)
    type(Potential_Local_E_Mix),          intent(inout) :: this
    type(Inoutput),intent(inout),optional:: file
    
    call print('Local_E_Mix potential:', file=file)
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
  end subroutine potential_local_e_mix_print

  recursive subroutine potential_local_e_mix_calc(this, at, args_str, error)
    type(Potential_Local_E_Mix), intent(inout) :: this
    type(Atoms), intent(inout) :: at
    character(len=*), intent(in), optional :: args_str
    integer, intent(out), optional :: error

    logical :: calc_weights
    integer :: core_hops, transition_hops, buffer_hops
    type(Dictionary) :: params, calc_create_hybrid_weights_params
    type(Table) :: region1_table
    integer :: i
    integer, pointer :: hybrid(:), hybrid_mark(:)

    INIT_ERROR(error)

    if (.not. associated(this%pot_region1) .or. .not. associated(this%pot_region2)) &
      call system_abort("Potential_Local_E_Mix_calc: this%pot_region1 or this%pot_region2 not initialised")

    call initialise(params)
    call param_register(params, "calc_weights", "F", calc_weights, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, "core_hops", "0", core_hops, help_string="No help yet.  This source file was $LastChangedBy$")
    if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='Potential_Local_E_Mix_Calc args_str') ) &
      call system_abort("Potential_Local_E_Mix_calc_energy failed to parse args_str='"//trim(args_str)//"'")
    call finalise(params)

    if (calc_weights) then
      call print("Potential_Local_E_Mix_calc got calc_weights core_hops " // core_hops, PRINT_VERBOSE)
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

    call calc_local_energy_mix(this, at, args_str, error=error)
    PASS_ERROR(error)

  end subroutine potential_local_e_mix_calc

  recursive subroutine calc_local_energy_mix(this, at, args_str, error)
    type(Potential_Local_E_Mix), intent(inout) :: this
    type(Atoms), intent(inout) :: at
    character(len=*), intent(in), optional :: args_str
    integer, intent(out), optional :: error

    real(dp) :: energy, virial(3,3)
    real(dp), pointer :: at_force_ptr(:,:), at_local_energy_ptr(:)

    ! local variables
    integer :: i
    integer,  pointer :: hybrid_mark(:), index(:), termindex(:)
    real(dp), pointer :: weight(:), weight_region1(:), rescale(:), c_weight(:)
    real(dp), pointer :: f_cluster(:,:), local_e_cluster(:)
    real(dp), allocatable :: weight_saved(:), local_e_weight(:), local_e_region2(:), local_e_region1(:), f_region1(:,:)
    real(dp) :: e_cluster
    logical :: dummy
    type(Atoms) :: cluster
    type(Table) :: cluster_info
    type(Dictionary) :: params
    character(FIELD_LENGTH) :: cc_args_str

    real(dp), pointer :: local_e_pot1(:), local_e_pot2(:)
    character(STRING_LENGTH) :: calc_energy, calc_force, calc_local_energy, calc_virial, calc_local_virial, local_energy_args_str, local_energy_name

    INIT_ERROR(error)

    allocate(weight_saved(at%N))
    allocate(local_e_weight(at%N))
    allocate(local_e_region2(at%N))
    allocate(local_e_region1(at%N))
    allocate(f_region1(3,at%N))

    call system_timer("calc_local_energy_mix")

    call initialise(params)
    call param_register(params, 'energy', '', calc_energy, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'force', '', calc_force, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'local_energy', '', calc_local_energy, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'virial', '', calc_virial, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'local_virial', '', calc_local_virial, help_string="No help yet.  This source file was $LastChangedBy$")
    if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='Potential_Local_E_Mix args_str') ) &
      call system_abort("Potential_Local_E_Mix failed to parse args_str='"//trim(args_str)//"'")
    call finalise(params)

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
      if(len_trim(calc_virial) > 0 .or. len_trim(calc_local_virial) > 0) then
         RAISE_ERROR('hybrid_calc_energy_mix: virial or local_virial not yet implemented when QM region is active', error)
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
      cluster_info =  create_cluster_info_from_mark(at, cc_args_str)
      call carve_cluster(at, cc_args_str, cluster_info, cluster)

      if (this%minimise_mm) then
	dummy = assign_pointer(cluster, 'index', index)
	call do_minimise_mm(this%relax_pot, at, this%minim_mm_method, this%minim_mm_tol, this%minim_mm_max_steps, &
	  this%minim_mm_linminroutine, this%minim_mm_do_pos, this%minim_mm_do_lat, this%minim_mm_do_print, &
	  this%minim_mm_args_str, this%minim_mm_eps_guess, this%minim_mm_use_n_minim, this%minim_inoutput_movie, &
          this%minim_cinoutput_movie, index)
      endif
    endif

    if (len_trim(calc_local_energy) > 0) then
      call assign_property_pointer(at, trim(local_energy_name), at_local_energy_ptr)
      local_energy_name = trim(calc_local_energy)
      local_energy_args_str = ""
    else
      local_energy_name = "LEMIX_local_energy"
      call add_property(at, trim(local_energy_name), 0.0_dp, ptr=at_local_energy_ptr)
      local_energy_args_str = "local_energy="//trim(local_energy_name)
    endif
    if (len_trim(calc_force) > 0) then
       call assign_property_pointer(at, trim(calc_force), at_force_ptr)
    endif

    ! do calculation in region 2 first
    weight = weight_saved-weight_region1 ! "negate" the weights

    call system_timer("calc_local_energy_mix/calc_region_2")
    call Calc(this%pot_region2, at, args_str=args_str//" "//trim(local_energy_args_str))
    call system_timer("calc_local_energy_mix/calc_region_2")
    if (len_trim(calc_energy) > 0) call get_param_value(at, trim(calc_energy), energy)
    local_e_region2 = at_local_energy_ptr

    if (assign_pointer(at, "local_e_pot2", local_e_pot2)) local_e_pot2 = local_e_region2

    call system_timer("calc_local_energy_mix/prep_cluster")
    if(len_trim(calc_local_energy) > 0) then
      local_e_weight = 1.0_dp-weight_region1
      at_local_energy_ptr = local_e_weight * local_e_region2
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
    call Calc(this%pot_region1, cluster, args_str = args_str//" solver=DIAG_GF SCF=GLOBAL_U "//trim(local_energy_args_str))
    if (len_trim(calc_force) > 0) then
       call assign_property_pointer(cluster, trim(calc_force), f_cluster, error=error)
       PASS_ERROR_WITH_INFO("Potential_Local_E_Mix failed to assign pointer for region 1 force property '"//trim(calc_force)//"'", error)
      ! scale back forces
      f_cluster = f_cluster*this%E_scale_pot1*this%r_scale_pot1
    endif
    if (len_trim(calc_energy) > 0) call get_param_value(at, trim(calc_energy), e_cluster)
    call assign_property_pointer(at, trim(local_energy_name), local_e_cluster, error=error)
    PASS_ERROR_WITH_INFO("Potential_Local_E_Mix failed to assign pointer for region 1 local_energy property '"//trim(calc_local_energy)//"'", error)

    if (len_trim(calc_energy) > 0) e_cluster = e_cluster*this%E_scale_pot1
    if (len_trim(calc_local_energy) > 0) local_e_cluster = local_e_cluster*this%E_scale_pot1
    call system_timer("calc_local_energy_mix/calc_region_1")

    call system_timer("calc_local_energy_mix/combine")
    if(current_verbosity() > PRINT_VERBOSE) call write(cluster, "cluster.cfg")

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
    if(len_trim(calc_force) > 0) then
       do i=1,cluster%N
          if(termindex(i) > 0) &
               f_cluster(:,termindex(i)) = f_cluster(:,termindex(i)) + (1.0_dp-rescale(i))*f_cluster(:,i)
       end do
    end if
    ! now do the redistribution
    do i=1,cluster%N
       if(len_trim(calc_force) > 0) f_region1(:,index(i)) = f_region1(:,index(i)) + f_cluster(:,i)*rescale(i)
       if(len_trim(calc_local_energy) > 0) local_e_region1(index(i)) = local_e_region1(index(i)) + local_e_cluster(i)
    end do

    if (assign_pointer(at, "local_e_pot1", local_e_pot1)) local_e_pot1 = local_e_region1

    ! finally add region 1 stuff to the output variables
    if(len_trim(calc_energy) > 0) then
      energy = energy + e_cluster
      call set_param_value(at, trim(calc_energy), energy)
    endif
    if(len_trim(calc_local_energy) > 0) then
      local_e_weight = weight_region1
      at_local_energy_ptr = at_local_energy_ptr + local_e_weight * local_e_region1
    endif
    if(len_trim(calc_force) > 0) at_force_ptr = at_force_ptr + f_region1

    weight = weight_saved

    call finalise(cluster)

    deallocate(weight_saved)
    deallocate(local_e_weight)
    deallocate(local_e_region2)
    deallocate(local_e_region1)
    deallocate(f_region1)

    if (len_trim(calc_local_energy) == 0) call remove_property(at, trim(local_energy_name))

    call system_timer("calc_local_energy_mix/combine")
    call system_timer("calc_local_energy_mix")
  end subroutine calc_local_energy_mix
