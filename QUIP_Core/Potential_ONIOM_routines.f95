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

  subroutine potential_ONIOM_initialise(this, args_str, region1_pot, region2_pot, reference_bulk, mpi_obj)
    type(Potential_ONIOM), intent(inout) :: this
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
    call param_register(params, "do_tb_defaults", "F", do_tb_defaults)

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

    if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='Potential_ONIOM_initialise args_str') ) &
      call system_abort("Potential_ONIOM_initialise failed to parse args_str='"//trim(args_str)//"'")
    call finalise(params)

    if (this%minimise_mm) then
      call initialise(this%relax_pot, "Simple", region2_pot)
    endif

    if (do_rescale_r .or. do_rescale_E .or. do_tb_defaults) then
      if (.not. present(reference_bulk)) &
	call system_abort("potential_local_e_mix_initialise got do_rescale_r="//do_rescale_r//" do_rescale_E="//do_rescale_E// &
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

  subroutine potential_ONIOM_finalise(this)
    type(Potential_ONIOM), intent(inout) :: this

    call finalise(this%pot_region1)
    call finalise(this%pot_region2)

    deallocate(this%pot_region1)
    deallocate(this%pot_region2)
  end subroutine potential_ONIOM_finalise


  function potential_ONIOM_cutoff(this)
    type(Potential_ONIOM), intent(in) :: this
    real(dp) :: potential_ONIOM_cutoff

    potential_ONIOM_cutoff = max(cutoff(this%pot_region1), cutoff(this%pot_region2))
  end function potential_ONIOM_cutoff

  subroutine potential_ONIOM_print(this, file)
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


  subroutine potential_ONIOM_calc(this, at, energy, local_e, force, virial, args_str, err)
    type(Potential_ONIOM), intent(inout) :: this
    type(Atoms), intent(inout) :: at
    real(dp), intent(out), optional :: energy
    real(dp), intent(out), optional :: local_e(:)
    real(dp), intent(out), optional :: force(:,:)
    real(dp), intent(out), optional :: virial(3,3)
    character(len=*), intent(in), optional :: args_str
    integer, intent(out), optional :: err

    logical :: calc_weights
    integer :: core_hops
    type(Dictionary) :: params
    type(Table) :: region1_table
    integer :: i
    integer, pointer :: hybrid(:), hybrid_mark(:)

    if (.not. associated(this%pot_region1) .or. .not. associated(this%pot_region2)) &
      call system_abort("Potential_ONIOM_calc: this%pot_region1 or this%pot_region2 not initialised")

    call initialise(params)
    call param_register(params, "calc_weights", "F", calc_weights)
    call param_register(params, "core_hops", "0", core_hops)
    if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='Potential_ONIOM_Calc args_str') ) &
      call system_abort("Potential_ONIOM_calc_energy failed to parse args_str='"//trim(args_str)//"'")
    call finalise(params)
    if (calc_weights) then
      call print("Potential_ONIOM_calc got calc_weights core_hops " // core_hops, PRINT_VERBOSE)
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

      call create_hybrid_weights(at, trans_width=0, buffer_width=0)
    endif

    call calc_oniom(this, at, energy, local_e, force, virial, args_str)

    if (present(err)) err = 0

  end subroutine potential_ONIOM_calc


  subroutine calc_oniom(this, at, energy, local_e, force, virial, args_str)
    type(Potential_ONIOM), intent(inout) :: this
    type(Atoms), intent(inout) :: at
    real(dp), intent(out), optional :: energy
    real(dp), intent(out), optional :: local_e(:)
    real(dp), intent(out), optional :: force(:,:)
    real(dp), intent(out), optional :: virial(3,3)
    character(len=*), intent(in), optional :: args_str

    type(Atoms) :: cluster
    type(Dictionary) :: params
    character(FIELD_LENGTH) :: cc_args_str

    real(dp) :: cluster_energy_1, cluster_energy_2
    real(dp), allocatable :: cluster_local_e_1(:), cluster_local_e_2(:)
    real(dp), allocatable :: cluster_force_1(:,:), cluster_force_2(:,:)
    real(dp) :: cluster_virial_1(3,3), cluster_virial_2(3,3)

    integer, pointer :: hybrid_mark(:), index(:), termindex(:)
    real(dp), pointer :: rescale(:)
    logical :: dummy
    integer i

    call system_timer("calc_oniom")

    if(.not. (associated(this%pot_region1) .and. associated(this%pot_region2))) then
       call system_abort('calc_oniom: potential for region 1 or region 2 is not initialised')
    end if

    ! Check for a compatible hybrid_mark property. It must be present.
    if (.not.assign_pointer(at, 'hybrid_mark', hybrid_mark)) &
       call system_abort('calc_oniom: atoms structure has no "hybrid_mark" property')

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
      cluster =  create_cluster_hybrid_mark(at, cc_args_str)

      if (this%minimise_mm) then
	dummy = assign_pointer(cluster, 'index', index)
	call do_minimise_mm(this%relax_pot, at, this%minim_mm_method, this%minim_mm_tol, this%minim_mm_max_steps, &
	  this%minim_mm_linminroutine, this%minim_mm_do_pos, this%minim_mm_do_lat, this%minim_mm_do_print, &
	  this%minim_mm_args_str, this%minim_mm_eps_guess, this%minim_mm_use_n_minim, this%minim_inoutput_movie, &
          this%minim_cinoutput_movie, index)
      endif
    endif

    ! do calculation in whole system first using pot_region2
    call system_timer("calc_oniom/calc_whole_sys")
    call Calc(this%pot_region2, at, energy, local_e, force, virial, args_str=args_str)
    call system_timer("calc_oniom/calc_whole_sys")

    call system_timer("calc_oniom/prep_cluster")

    ! if there are no marked atoms, we are done, so return
    if(all(hybrid_mark == HYBRID_NO_MARK)) then
       call print('No REGION 1 marks, not carving REGION 1 cluster')
       return
    end if

    ! rescale cluster
    call set_lattice(cluster, this%r_scale_pot1 * cluster%lattice, scale_positions=.true.)

    call calc_connect(cluster)

    call system_timer("calc_oniom/prep_cluster")

    call system_timer("calc_oniom/calc_cluster_1")
    allocate(cluster_force_1(3,cluster%N))
    allocate(cluster_local_e_1(cluster%N))
    if (present(force)) then
      if (present(virial)) then
        call Calc(this%pot_region1, cluster, cluster_energy_1, cluster_local_e_1, cluster_force_1, cluster_virial_1, args_str=args_str)
        ! scale back forces
        cluster_force_1 = cluster_force_1*this%E_scale_pot1*this%r_scale_pot1
        cluster_virial_1 = cluster_virial_1*this%E_scale_pot1
      else
        call Calc(this%pot_region1, cluster, cluster_energy_1, cluster_local_e_1, cluster_force_1, args_str=args_str)
        ! scale back forces
        cluster_force_1 = cluster_force_1*this%E_scale_pot1*this%r_scale_pot1
      endif
    else
      if (present(virial)) then
        call Calc(this%pot_region1, cluster, cluster_energy_1, cluster_local_e_1, virial=cluster_virial_1, args_str=args_str)
        ! scale back virial
        cluster_virial_1 = cluster_virial_1*this%E_scale_pot1
      else
        call Calc(this%pot_region1, cluster, cluster_energy_1, cluster_local_e_1, args_str=args_str)
      endif
    endif
    cluster_energy_1 = cluster_energy_1*this%E_scale_pot1
    cluster_local_e_1 = cluster_local_e_1*this%E_scale_pot1
    call system_timer("calc_oniom/calc_cluster_1")

    ! unrescale cluster
    call set_lattice(cluster, cluster%lattice / this%r_scale_pot1, scale_positions=.true.)

    call calc_connect(cluster)

    call system_timer("calc_oniom/calc_cluster_2")
    allocate(cluster_force_2(3,cluster%N))
    allocate(cluster_local_e_2(cluster%N))
    if (present(force)) then
      if (present(virial)) then
        call Calc(this%pot_region2, cluster, cluster_energy_2, cluster_local_e_2, cluster_force_2, cluster_virial_2, args_str=args_str)
      else
        call Calc(this%pot_region2, cluster, cluster_energy_2, cluster_local_e_2, cluster_force_2, args_str=args_str)
      endif
    else
      if (present(virial)) then
        call Calc(this%pot_region2, cluster, cluster_energy_2, cluster_local_e_2, virial=cluster_virial_2, args_str=args_str)
      else
        call Calc(this%pot_region2, cluster, cluster_energy_2, cluster_local_e_2, args_str=args_str)
      endif
    endif
    call system_timer("calc_oniom/calc_cluster_2")

    call system_timer("calc_oniom/combine")
    if(current_verbosity() > PRINT_VERBOSE) &
         call print_cfg(cluster, "cluster.cfg", all_properties = .true.)

    ! redistribute local energies and weights from the cluster to region 1
    dummy = assign_pointer(cluster, 'index', index)
    dummy = assign_pointer(cluster, 'rescale', rescale)
    dummy = assign_pointer(cluster, 'termindex', termindex)
    ! first add terms arising from the fact that the position of the rescaled termination atoms
    ! depends on the position of the atoms they are terminating
    if(present(force)) then
       do i=1,cluster%N
          if(termindex(i) > 0) then
             cluster_force_1(:,termindex(i)) = cluster_force_1(:,termindex(i)) + (1.0_dp-rescale(i))*cluster_force_1(:,i)
             cluster_force_2(:,termindex(i)) = cluster_force_2(:,termindex(i)) + (1.0_dp-rescale(i))*cluster_force_2(:,i)
          endif
       end do
    end if
    ! now do the redistribution
    if(present(force)) then
      do i=1,cluster%N
        force(:,index(i)) = force(:,index(i)) + cluster_force_1(:,i)*rescale(i) - cluster_force_2(:,i)*rescale(i)
      end do
    end if
    if(present(virial)) then
      virial = virial + cluster_virial_1 - cluster_virial_2
    endif
    if (present(local_e)) then
      do i=1,cluster%N
        local_e(index(i)) = local_e(index(i)) + cluster_local_e_1(i) - cluster_local_e_2(i)
      end do
    end if
    if(present(energy)) then
      energy = energy + cluster_energy_1 - cluster_energy_2
    endif
    call system_timer("calc_oniom/combine")

    call finalise(cluster)
    deallocate(cluster_force_1)
    deallocate(cluster_local_e_1)
    deallocate(cluster_force_2)
    deallocate(cluster_local_e_2)

    call system_timer("calc_oniom")
  end subroutine calc_oniom
