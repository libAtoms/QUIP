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
!X ONIOM routines to be included in Potential.f95
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  recursive subroutine potential_ONIOM_initialise(this, args_str, region1_pot, region2_pot, reference_bulk, mpi_obj, error)
    type(Potential_ONIOM), intent(inout) :: this
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
    call param_register(params, "hybrid_mark_postfix", "", this%hybrid_mark_postfix, help_string="No help yet.  This source file was $LastChangedBy$")

    call param_register(params, "r_scale", "1.0", this%r_scale_pot1, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, "E_scale", "1.0", this%E_scale_pot1, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, "terminate", "T", this%terminate, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, "minimise_bulk", "F", minimise_bulk, help_string="No help yet.  This source file was $LastChangedBy$")

    call param_register(params, "do_rescale_r", "F", do_rescale_r, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, "do_rescale_E", "F", do_rescale_E, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, "do_tb_defaults", "F", do_tb_defaults, help_string="No help yet.  This source file was $LastChangedBy$")

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

    if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='Potential_ONIOM_initialise args_str') ) then
      RAISE_ERROR("Potential_ONIOM_initialise failed to parse args_str='"//trim(args_str)//"'",error)
    endif
    call finalise(params)

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

  recursive subroutine potential_ONIOM_finalise(this)
    type(Potential_ONIOM), intent(inout) :: this

    call finalise(this%pot_region1)
    call finalise(this%pot_region2)

    deallocate(this%pot_region1)
    deallocate(this%pot_region2)
  end subroutine potential_ONIOM_finalise


  recursive function potential_ONIOM_cutoff(this)
    type(Potential_ONIOM), intent(in) :: this
    real(dp) :: potential_ONIOM_cutoff

    potential_ONIOM_cutoff = max(cutoff(this%pot_region1), cutoff(this%pot_region2))
  end function potential_ONIOM_cutoff

  recursive subroutine potential_ONIOM_print(this, file)
    type(Potential_ONIOM),          intent(inout) :: this
    type(Inoutput),intent(inout),optional:: file
    
    call print('ONIOM potential:', file=file)
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
  end subroutine potential_ONIOM_print


  recursive subroutine potential_ONIOM_calc(this, at, args_str, error)
    type(Potential_ONIOM), intent(inout) :: this
    type(Atoms), intent(inout) :: at
    character(len=*), intent(in), optional :: args_str
    integer, intent(out), optional :: error

    logical :: calc_weights
    integer :: core_hops
    type(Dictionary) :: params
    type(Table) :: region1_table
    integer :: i
    integer, pointer :: hybrid(:), hybrid_mark(:)
    character(len=STRING_LENGTH) :: hybrid_mark_postfix

    INIT_ERROR(error)

    if (.not. associated(this%pot_region1) .or. .not. associated(this%pot_region2)) then
      RAISE_ERROR("Potential_ONIOM_calc: this%pot_region1 or this%pot_region2 not initialised", error)
    endif

    call initialise(params)
    call param_register(params, "hybrid_mark_postfix", trim(this%hybrid_mark_postfix), hybrid_mark_postfix, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, "calc_weights", "F", calc_weights, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, "core_hops", "0", core_hops, help_string="No help yet.  This source file was $LastChangedBy$")
    if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='Potential_ONIOM_Calc args_str') ) then
      RAISE_ERROR("Potential_ONIOM_calc_energy failed to parse args_str='"//trim(args_str)//"'",error)
    endif
    call finalise(params)
    if (calc_weights) then
      call print("Potential_ONIOM_calc got calc_weights core_hops " // core_hops, PRINT_VERBOSE)
      call add_property(at, "weight_region1"//trim(hybrid_mark_postfix), 0.0_dp)
      call add_property(at, "hybrid_mark"//trim(hybrid_mark_postfix), HYBRID_NO_MARK)
      if (.not. assign_pointer(at, "hybrid"//trim(hybrid_mark_postfix), hybrid)) then
        RAISE_ERROR("QC_QUIP_calc at doesn't have hybrid"//trim(hybrid_mark_postfix)//" property and calc_weights was specified", error)
      endif
      if (.not. assign_pointer(at, "hybrid_mark"//trim(hybrid_mark_postfix), hybrid_mark)) then
        RAISE_ERROR("QC_QUIP_calc Failed to add hybrid_mark"//trim(hybrid_mark_postfix)//" property to at", error)
      endif

      call calc_connect(at)

      call wipe(region1_table)
      hybrid_mark = HYBRID_NO_MARK
      do i=1, size(hybrid)
        if (hybrid(i) /= 0) call append(region1_table, (/ i, 0, 0, 0 /) )
      end do
      call bfs_grow(at, region1_table, core_hops, nneighb_only = .false.)
      hybrid_mark(int_part(region1_table,1)) = HYBRID_ACTIVE_MARK

      call create_hybrid_weights(at, args_str="transition_hops=0 buffer_hops=0")
    endif

    call calc_oniom(this, at, args_str, error)
    PASS_ERROR(error)

  end subroutine potential_ONIOM_calc


  recursive subroutine calc_oniom(this, at, args_str, error)
    type(Potential_ONIOM), intent(inout) :: this
    type(Atoms), intent(inout) :: at
    character(len=*), intent(in), optional :: args_str
    integer, intent(out), optional :: error

    real(dp) :: energy, virial(3,3)
    real(dp), pointer :: at_force_ptr(:,:), at_local_energy_ptr(:), at_local_virial_ptr(:,:)

    type(Atoms) :: cluster
    type(Table) :: cluster_info
    type(Dictionary) :: params
    character(FIELD_LENGTH) :: cc_args_str

    real(dp) :: cluster_energy_1, cluster_energy_2
    real(dp), allocatable :: cluster_local_e_1(:), cluster_local_e_2(:)
    real(dp), pointer :: cluster_force_ptr(:,:), cluster_local_energy_ptr(:), cluster_local_virial_ptr(:,:)
    real(dp), allocatable :: cluster_force_1(:,:), cluster_force_2(:,:), cluster_local_virial_1(:,:), cluster_local_virial_2(:,:)
    real(dp) :: cluster_virial_1(3,3), cluster_virial_2(3,3)

    integer, pointer :: hybrid_mark(:), index(:), termindex(:)
    real(dp), pointer :: rescale(:)
    logical :: dummy
    integer i

    character(STRING_LENGTH) :: calc_energy, calc_force, calc_virial, calc_local_energy, calc_local_virial, hybrid_mark_postfix

    INIT_ERROR(error)

    call system_timer("calc_oniom")

    call initialise(params)
    call param_register(params, "hybrid_mark_postfix", trim(this%hybrid_mark_postfix), hybrid_mark_postfix, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, "energy", "", calc_energy, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, "force", "", calc_force, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, "virial", "", calc_virial, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, "local_energy", "", calc_local_energy, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, "local_virial", "", calc_local_virial, help_string="No help yet.  This source file was $LastChangedBy$")
    if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='Calc_ONIOM args_str')) then
      RAISE_ERROR("Calc_ONIOM failed to parse args_str='"//trim(args_str)//"'", error)
    endif
    call finalise(params)

    if(.not. (associated(this%pot_region1) .and. associated(this%pot_region2))) then
       RAISE_ERROR('calc_oniom: potential for region 1 or region 2 is not initialised', error)
    end if

    ! Check for a compatible hybrid_mark property. It must be present.
    if (.not.assign_pointer(at, 'hybrid_mark'//trim(hybrid_mark_postfix), hybrid_mark)) then
       RAISE_ERROR('calc_oniom: atoms structure has no "hybrid_mark'//trim(hybrid_mark_postfix)//'" property', error)
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
      cluster_info = create_cluster_info_from_mark(at, cc_args_str)
      call carve_cluster(at, cc_args_str, cluster_info, cluster)


      if (this%minimise_mm) then
	dummy = assign_pointer(cluster, 'index', index)
	call do_minimise_mm(this%relax_pot, at, this%minim_mm_method, this%minim_mm_tol, this%minim_mm_max_steps, &
	  this%minim_mm_linminroutine, this%minim_mm_do_pos, this%minim_mm_do_lat, this%minim_mm_do_print, &
	  this%minim_mm_args_str, this%minim_mm_eps_guess, this%minim_mm_use_n_minim, this%minim_inoutput_movie, &
          this%minim_cinoutput_movie, index)
      endif
    endif

    if (len_trim(calc_force) > 0) call assign_property_pointer(at, trim(calc_force), at_force_ptr)
    if (len_trim(calc_local_energy) > 0) call assign_property_pointer(at, trim(calc_local_energy), at_local_energy_ptr)
    if (len_trim(calc_local_virial) > 0) call assign_property_pointer(at, trim(calc_local_virial), at_local_virial_ptr)

    ! do calculation in whole system first using pot_region2
    call system_timer("calc_oniom/calc_whole_sys")
    call Calc(this%pot_region2, at, args_str=args_str)
    call system_timer("calc_oniom/calc_whole_sys")
    if (len_trim(calc_energy) > 0) call get_param_value(at, trim(calc_energy), energy)
    if (len_trim(calc_virial) > 0) call get_param_value(at, trim(calc_virial), virial)

    call system_timer("calc_oniom/prep_cluster")

    ! if there are no marked atoms, we are done, so return
    if(all(hybrid_mark == HYBRID_NO_MARK)) then
       call print('WARNING: ONIOM with No REGION 1 marks, not carving REGION 1 cluster')
       return
    end if

    ! rescale cluster
    call set_lattice(cluster, this%r_scale_pot1 * cluster%lattice, scale_positions=.true.)

    call calc_connect(cluster)

    call system_timer("calc_oniom/prep_cluster")

    call system_timer("calc_oniom/calc_cluster_1")
    allocate(cluster_force_1(3,cluster%N))
    allocate(cluster_local_e_1(cluster%N))
    allocate(cluster_local_virial_1(9,cluster%N))
    call Calc(this%pot_region1, cluster, args_str=args_str)
    if (len_trim(calc_force) > 0) then
      call assign_property_pointer(cluster, trim(calc_force), cluster_force_ptr, error=error)
      PASS_ERROR_WITH_INFO("Calc_ONIOM failed to get value for force property '"//trim(calc_force)//"'", error)
      cluster_force_1 = cluster_force_ptr*this%E_scale_pot1*this%r_scale_pot1
    endif
    if (len_trim(calc_virial) > 0) then
      call get_param_value(cluster, trim(calc_virial), cluster_virial_1)
      cluster_virial_1 = cluster_virial_1*this%E_scale_pot1
    endif
    if (len_trim(calc_energy) > 0) then
      call get_param_value(cluster, trim(calc_energy), cluster_energy_1)
      cluster_energy_1 = cluster_energy_1*this%E_scale_pot1
    endif
    if (len_trim(calc_local_energy) > 0) then
      call assign_property_pointer(cluster, trim(calc_local_energy), cluster_local_energy_ptr, error=error)
      RAISE_ERROR("Calc_ONIOM failed to get value for local_energy property '"//trim(calc_local_energy)//"'", ERROR)
      cluster_local_e_1 = cluster_local_energy_ptr*this%E_scale_pot1
    endif
    if (len_trim(calc_local_virial) > 0) then
      call assign_property_pointer(cluster, trim(calc_local_virial), cluster_local_virial_ptr, error=error)
      PASS_ERROR_WITH_INFO("Calc_ONIOM failed to get value for local_virial property '"//trim(calc_local_virial)//"'", error)
      cluster_local_virial_1 = cluster_local_virial_ptr*this%E_scale_pot1
    endif
    call system_timer("calc_oniom/calc_cluster_1")

    ! unrescale cluster
    call set_lattice(cluster, cluster%lattice / this%r_scale_pot1, scale_positions=.true.)

    call calc_connect(cluster)

    call system_timer("calc_oniom/calc_cluster_2")
    allocate(cluster_force_2(3,cluster%N))
    allocate(cluster_local_e_2(cluster%N))
    allocate(cluster_local_virial_2(9,cluster%N))
    call Calc(this%pot_region2, cluster, args_str=args_str)
    if (len_trim(calc_force) > 0) cluster_force_2 = cluster_force_ptr
    if (len_trim(calc_virial) > 0) call get_param_value(cluster, trim(calc_virial), cluster_virial_2)
    if (len_trim(calc_energy) > 0) call get_param_value(cluster, trim(calc_energy), cluster_energy_2)
    if (len_trim(calc_local_energy) > 0) cluster_local_e_2 = cluster_local_energy_ptr
    if (len_trim(calc_local_virial) > 0) cluster_local_virial_2 = cluster_local_virial_ptr
    call system_timer("calc_oniom/calc_cluster_2")

    call system_timer("calc_oniom/combine")
    if(current_verbosity() > PRINT_VERBOSE) &
         call write(cluster, "cluster.cfg")

    ! redistribute local energies and weights from the cluster to region 1
    dummy = assign_pointer(cluster, 'index', index)
    dummy = assign_pointer(cluster, 'rescale', rescale)
    dummy = assign_pointer(cluster, 'termindex', termindex)
    ! first add terms arising from the fact that the position of the rescaled termination atoms
    ! depends on the position of the atoms they are terminating
    if(len_trim(calc_force) > 0) then
       do i=1,cluster%N
          if(termindex(i) > 0) then
             cluster_force_1(:,termindex(i)) = cluster_force_1(:,termindex(i)) + (1.0_dp-rescale(i))*cluster_force_1(:,i)
             cluster_force_2(:,termindex(i)) = cluster_force_2(:,termindex(i)) + (1.0_dp-rescale(i))*cluster_force_2(:,i)
          endif
       end do
    end if
    ! now do the redistribution
    if(len_trim(calc_force) > 0) then
      do i=1,cluster%N
        at_force_ptr(:,index(i)) = at_force_ptr(:,index(i)) + cluster_force_1(:,i)*rescale(i) - cluster_force_2(:,i)*rescale(i)
      end do
    end if
    if(len_trim(calc_virial) > 0) then
      virial = virial + cluster_virial_1 - cluster_virial_2
      call set_param_value(at, trim(calc_virial), virial)
    endif
    if(len_trim(calc_local_energy) > 0) then
      do i=1,cluster%N
        at_local_energy_ptr(index(i)) = at_local_energy_ptr(index(i)) + cluster_local_e_1(i) - cluster_local_e_2(i)
      end do
    end if
    if(len_trim(calc_energy) > 0) then
      energy = energy + cluster_energy_1 - cluster_energy_2
      call set_param_value(at, trim(calc_energy), energy)
    endif
    if(len_trim(calc_local_virial) > 0) then
      do i=1,cluster%N
        at_local_virial_ptr(:,index(i)) = at_local_virial_ptr(:,index(i)) + cluster_local_virial_1(:,i) - cluster_local_virial_2(:,i)
      end do
    end if
    call system_timer("calc_oniom/combine")

    call finalise(cluster)
    deallocate(cluster_force_1)
    deallocate(cluster_local_e_1)
    deallocate(cluster_local_virial_1)
    deallocate(cluster_force_2)
    deallocate(cluster_local_e_2)
    deallocate(cluster_local_virial_2)

    call system_timer("calc_oniom")
  end subroutine calc_oniom
