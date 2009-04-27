
  !*************************************************************************
  !*
  !*  MetaPotential_ForceMixing routines
  !*
  !*************************************************************************

  subroutine MetaPotential_ForceMixing_initialise(this, args_str, mmpot, qmpot, reference_bulk, mpi)
    type(MetaPotential_ForceMixing), intent(inout) :: this
    character(len=*), intent(in) :: args_str
    type(Potential), intent(in), target :: mmpot, qmpot
    type(Atoms), optional, intent(inout) :: reference_bulk
    type(MPI_Context), intent(in), optional :: mpi

    type(Dictionary) :: params
    logical :: minimise_bulk, do_tb_defaults
    real(dp) :: dummy_E

    call finalise(this)

    this%init_args_str = args_str

    call initialise(params)
    call param_register(params, 'minimise_mm', 'F', this%minimise_mm)
    call param_register(params, 'calc_weights', 'T', this%calc_weights)
    call param_register(params, 'method', 'conserve_momentum', this%method)
    call param_register(params, 'mm_reweight', '1.0', this%mm_reweight)
    call param_register(params, 'conserve_momentum_weight_method', 'uniform', this%conserve_momentum_weight_method)
    call param_register(params, 'mm_args_str', '', this%mm_args_str)
    call param_register(params, 'qm_args_str', '', this%qm_args_str)
    call param_register(params, 'buffer_hops', '0', this%buffer_hops)
    call param_register(params, 'hysteretic_buffer', 'F', this%hysteretic_buffer)
    call param_register(params, 'hysteretic_buffer_inner_radius', '5.0', this%hysteretic_buffer_inner_radius)
    call param_register(params, 'hysteretic_buffer_outer_radius', '7.0', this%hysteretic_buffer_outer_radius)
    call param_register(params, 'fit_hops', '3', this%fit_hops)
    call param_register(params, 'randomise_buffer', 'F', this%randomise_buffer)
    call param_register(params, 'transition_hops', '0', this%transition_hops)
    call param_register(params, 'weight_interpolation', 'hop_ramp', this%weight_interpolation)
    call param_register(params, 'nneighb_only', 'T', this%nneighb_only)
    call param_register(params, 'save_forces', 'T', this%save_forces)
    call param_register(params, 'lotf_spring_hops', '2', this%lotf_spring_hops)
    call param_register(params, 'lotf_interp_order', 'linear', this%lotf_interp_order)
    call param_register(params, 'lotf_interp_space', 'F', this%lotf_interp_space)


    ! Parameters for do_reference_bulk calculation
    call param_register(params, "minimise_bulk", "F", minimise_bulk)
    call param_register(params, "do_rescale_r", "F", this%do_rescale_r)
    call param_register(params, "do_tb_defaults", "F", do_tb_defaults)

    ! Parameters for the MM minimisation
    call param_register(params, 'minim_mm_method', 'cg', this%minim_mm_method)
    call param_register(params, 'minim_mm_tol', '1e-6', this%minim_mm_tol)
    call param_register(params, 'minim_mm_eps_guess', '1e-4', this%minim_mm_eps_guess)
    call param_register(params, 'minim_mm_max_steps', '1000', this%minim_mm_max_steps)
    call param_register(params, 'minim_mm_linminroutine', 'FAST_LINMIN', this%minim_mm_linminroutine)
    call param_register(params, 'minim_mm_do_pos', 'T', this%minim_mm_do_pos)
    call param_register(params, 'minim_mm_do_lat', 'F', this%minim_mm_do_lat)
    call param_register(params, 'minim_mm_do_print', 'F', this%minim_mm_do_print)
    call param_register(params, 'minim_mm_use_n_minim', 'F', this%minim_mm_use_n_minim)
    call param_register(params, 'minim_mm_args_str', '', this%minim_mm_args_str)

    if (.not. param_read_line(params, args_str, ignore_unknown=.true.)) then
       call system_abort('MetaPotential_ForceMixing_initialise failed to parse args_str="'//trim(args_str)//'"')
    end if

    call finalise(params)

    this%r_scale_pot1 = 1.0_dp
    if (this%do_rescale_r) then
       if (.not. present(reference_bulk)) &
            call system_abort("metapotential_forcemixing_initialise got do_rescale_r=T do_tb_defaults="//&
            do_tb_defaults//" but reference_bulk is not present")

       call do_reference_bulk(reference_bulk, mmpot, qmpot, minimise_bulk, .true., .false., &
            this%r_scale_pot1, dummy_E, do_tb_defaults)
      
       if (this%r_scale_pot1 <= 0.0_dp) this%r_scale_pot1 = 1.0_dp
       call print ("Rescaling positions in region1 potential by " // this%r_scale_pot1 // " to match lattice constants")
    end if

    this%mmpot => mmpot
    this%qmpot => qmpot

    if (this%minimise_mm) then
      call initialise(this%relax_metapot, "Simple", mmpot)
    endif

    if (present(mpi)) this%mpi = mpi

  end subroutine MetaPotential_ForceMixing_initialise


  subroutine MetaPotential_ForceMixing_finalise(this)
    type(MetaPotential_ForceMixing), intent(inout) :: this
    
    nullify(this%mmpot)
    nullify(this%qmpot)
    call finalise(this%embedlist)
    call finalise(this%fitlist)

  end subroutine MetaPotential_ForceMixing_finalise


  subroutine MetaPotential_ForceMixing_print(this, file)
    type(MetaPotential_ForceMixing), intent(inout) :: this
    type(Inoutput), intent(inout), optional :: file

    if (current_verbosity() < NORMAL) return

    call Print('MetaPotential_ForceMixing:',file=file)
    call Print(' minimise_mm='//this%minimise_mm,file=file)
    call Print(' calc_weights='//this%calc_weights,file=file)
    call Print(' method='//trim(this%method),file=file)
    call Print(' mm_reweight='//this%mm_reweight,file=file)
    call Print(' conserve_momentum_weight_method='//this%conserve_momentum_weight_method)
    call Print(' mm_args_str='//trim(this%mm_args_str), file=file)
    call Print(' qm_args_str='//trim(this%qm_args_str), file=file)
    call Print(' buffer_hops='//this%buffer_hops, file=file)
    call Print(' fit_hops='//this%fit_hops, file=file)
    call Print(' randomise_buffer='//this%randomise_buffer, file=file)
    call Print(' transition_hops='//this%transition_hops, file=file)
    call Print(' weight_interpolation='//trim(this%weight_interpolation), file=file)
    call Print(' nneighb_only='//this%nneighb_only, file=file)
    call Print(' save_forces='//this%save_forces, file=file)
    call Print(' do_rescale_r='//this%do_rescale_r, file=file)
    call Print(' lotf_spring_hops='//this%lotf_spring_hops, file=file)
    call Print(' lotf_interp_order='//this%lotf_interp_order, file=file)
    call Print(' lotf_interp_space='//this%lotf_interp_space, file=file)
    call Print('',file=file)
    call Print(' minim_mm_method='//this%minim_mm_method,file=file)
    call Print(' minim_mm_tol='//this%minim_mm_tol,file=file)
    call Print(' minim_mm_eps_guess='//this%minim_mm_eps_guess,file=file)
    call Print(' minim_mm_max_steps='//this%minim_mm_max_steps,file=file)
    call Print(' minim_mm_linminroutine='//this%minim_mm_linminroutine,file=file)
    call Print(' minim_mm_do_pos='//this%minim_mm_do_pos,file=file)
    call Print(' minim_mm_do_lat='//this%minim_mm_do_lat,file=file)
    call Print(' minim_mm_do_print='//this%minim_mm_do_print,file=file)
    call Print(' minim_mm_use_n_minim='//this%minim_mm_use_n_minim,file=file)
    call Print(' minim_mm_args_str='//trim(this%minim_mm_args_str),file=file)
    call Print('',file=file)
    if (associated(this%mmpot)) then
       call Print('MM Potential',file=file)
       call Print(this%mmpot,file=file)
       call Print('',file=file)
    end if
    if (associated(this%qmpot)) then
       call Print('QM Potential',file=file)
       call Print(this%qmpot,file=file)
    end if
    call Print('',file=file)

  end subroutine MetaPotential_ForceMixing_print


  subroutine MetaPotential_ForceMixing_calc(this, at, e, local_e, f, virial, args_str, err)
    type(MetaPotential_ForceMixing), intent(inout) :: this
    type(Atoms), intent(inout) :: at
    real(dp), intent(out), optional :: e
    real(dp), intent(out), optional :: local_e(:)
    real(dp), intent(out), optional :: f(:,:)
    real(dp), intent(out), optional :: virial(3,3)
    integer, intent(out), optional :: err
    character(*), intent(in), optional :: args_str

    real(dp), allocatable, dimension(:,:) :: df, df_fit
    real(dp), allocatable, dimension(:,:) :: f_mm, f_qm
    real(dp), pointer, dimension(:,:) :: force_ptr
    real(dp), pointer, dimension(:) :: weight_region1, conserve_momentum_weight
    integer, pointer, dimension(:) :: hybrid, hybrid_mark
    integer :: i
    type(Dictionary) :: params
    logical :: dummy
    logical  :: minimise_mm, calc_weights, nneighb_only, save_forces, lotf_do_init, &
         lotf_do_map, lotf_do_fit, lotf_do_interp, lotf_do_qm, lotf_interp_space, do_rescale_r, &
         randomise_buffer, hysteretic_buffer
    character(FIELD_LENGTH) :: method, mm_args_str, qm_args_str, conserve_momentum_weight_method, &
         AP_method, weight_interpolation, lotf_interp_order
    real(dp) :: mm_reweight, dV_dt, f_tot(3), w_tot, weight, lotf_interp
    real(dp) :: hysteretic_buffer_inner_radius, hysteretic_buffer_outer_radius

    integer :: weight_method, buffer_hops, transition_hops, fit_hops, lotf_spring_hops
    integer,      parameter   :: UNIFORM_WEIGHT=1, MASS_WEIGHT=2, MASS2_WEIGHT=3, USER_WEIGHT=4
    integer, allocatable, dimension(:) :: embed, fit

    if (present(err)) err = 0

    ! Override parameters with those given in args_str
    call initialise(params)
    call param_register(params, "minimise_mm", ''//this%minimise_mm, minimise_mm)
    call param_register(params, "calc_weights", ''//this%calc_weights, calc_weights)
    call param_register(params, "method", this%method, method)
    call param_register(params, "mm_reweight", ''//this%mm_reweight, mm_reweight)
    call param_register(params, 'conserve_momentum_weight_method', ''//this%conserve_momentum_weight_method, &
         conserve_momentum_weight_method)
    call param_register(params, 'mm_args_str', this%mm_args_str, mm_args_str)
    call param_register(params, 'qm_args_str', this%qm_args_str, qm_args_str)
    call param_register(params, 'buffer_hops', ''//this%buffer_hops, buffer_hops)
    call param_register(params, 'hysteretic_buffer', ''//this%hysteretic_buffer, hysteretic_buffer)
    call param_register(params, 'hysteretic_buffer_inner_radius', ''//this%hysteretic_buffer_inner_radius, hysteretic_buffer_inner_radius)
    call param_register(params, 'hysteretic_buffer_outer_radius', ''//this%hysteretic_buffer_outer_radius, hysteretic_buffer_outer_radius)
    call param_register(params, 'transition_hops', ''//this%transition_hops, transition_hops)
    call param_register(params, 'fit_hops', ''//this%fit_hops, fit_hops)
    call param_register(params, 'randomise_buffer', ''//this%randomise_buffer, randomise_buffer)
    call param_register(params, 'weight_interpolation', this%weight_interpolation, weight_interpolation)
    call param_register(params, 'nneighb_only', ''//this%nneighb_only, nneighb_only)
    call param_register(params, 'lotf_spring_hops', ''//this%lotf_spring_hops, lotf_spring_hops)
    call param_register(params, 'lotf_interp_order', this%lotf_interp_order, lotf_interp_order)
    call param_register(params, 'lotf_interp_space', ''//this%lotf_interp_space, lotf_interp_space)
    call param_register(params, 'save_forces', ''//this%save_forces, save_forces)
    call param_register(params, 'do_rescale_r', ''//this%do_rescale_r, do_rescale_r)

    call param_register(params, 'lotf_do_init', 'T', lotf_do_init)
    call param_register(params, 'lotf_do_map', 'F', lotf_do_map)
    call param_register(params, 'lotf_do_qm', 'T', lotf_do_qm)
    call param_register(params, 'lotf_do_fit', 'T', lotf_do_fit)
    call param_register(params, 'lotf_do_interp', 'F', lotf_do_interp)
    call param_register(params, 'lotf_interp', '0.0', lotf_interp)

    
    ! Apply the shortcuts
    if (trim(method) == 'force_mixing_abrupt') then
       buffer_hops = 0
       transition_hops = 0
       weight_interpolation = 'hop_ramp'
    else if (trim(method) == 'force_mixing_smooth') then
       weight_interpolation = 'hop_ramp'
    else if (trim(method) == 'force_mixing_super_smooth') then
       weight_interpolation = 'distance_ramp'
    end if

    if (.not. param_read_line(params, args_str, ignore_unknown=.true.) ) &
      call system_abort("MetaPotential_ForceMixing_calc_energy failed to parse args_str='"//trim(args_str)//"'")
    call finalise(params)

    if (present(e) .or. present(local_e) .or. present(virial) .or. .not. present(f)) &
         call system_abort('MetaPotential_ForceMixing_calc: supports only forces, not energy, virial or local_e')

    allocate(f_mm(3,at%N),f_qm(3,at%N))
    f_mm = 0.0_dp
    f_qm = 0.0_dp

    if (calc_weights) then 

       if (.not. has_property(at, 'hybrid_mark')) &
            call add_property(at, 'hybrid_mark', HYBRID_NO_MARK)

       if (.not. assign_pointer(at, "hybrid", hybrid)) &
            call system_abort("MetaPotential_ForceMixing_calc: at doesn't have hybrid property and calc_weights was specified")

       if (.not. assign_pointer(at, 'hybrid_mark', hybrid_mark)) &
            call system_abort('MetaPotential_ForceMixing_Calc: hybrid_mark property missing')

       hybrid_mark = HYBRID_NO_MARK
       where(hybrid /= 0) hybrid_mark = HYBRID_ACTIVE_MARK

       call Print('MetaPotential_ForceMixing_calc: got '//count(hybrid /= 0)//' active atoms.')

       if (count(hybrid_mark == HYBRID_ACTIVE_MARK) == 0) &
            call system_abort('MetatPotential_ForceMixing_Calc: zero active atoms and calc_weights was specified')

       call create_hybrid_weights(at, trans_width=transition_hops, buffer_width=buffer_hops, &
            weight_interpolation=weight_interpolation, nneighb_only=nneighb_only, min_images_only=.true., &
            mark_buffer_outer_layer=randomise_buffer, hysteretic_buffer=hysteretic_buffer, &
            hysteretic_buffer_inner_radius=hysteretic_buffer_inner_radius, &
            hysteretic_buffer_outer_radius=hysteretic_buffer_outer_radius)

       ! reassign pointers
       
       if (.not. assign_pointer(at, "hybrid", hybrid)) &
            call system_abort("MetaPotential_ForceMixing_calc: at doesn't have hybrid property and calc_weights was specified")

       if (.not. assign_pointer(at, 'hybrid_mark', hybrid_mark)) &
            call system_abort('MetaPotential_ForceMixing_Calc: hybrid_mark property missing')

    end if

    if (.not. assign_pointer(at, 'hybrid_mark', hybrid_mark)) &
         call system_abort('MetaPotential_ForceMixing_Calc: hybrid_mark property missing')

    ! Do the MM minimisation, freezing the atoms marked as active
    if (minimise_mm) then
       call do_minimise_mm(this%relax_metapot, at, this%minim_mm_method, this%minim_mm_tol, this%minim_mm_max_steps, &
            this%minim_mm_linminroutine, this%minim_mm_do_pos, this%minim_mm_do_lat, this%minim_mm_do_print, &
            this%minim_mm_args_str, this%minim_mm_eps_guess, this%minim_mm_use_n_minim, this%minim_inoutput_movie, &
            this%minim_cinoutput_movie, find(hybrid_mark == HYBRID_ACTIVE_MARK))
    endif

    ! Do the classical calculation
    call calc(this%mmpot, at, f=f_mm, err=err, args_str=mm_args_str)

    if (.not. any(hybrid_mark /= HYBRID_NO_MARK)) then
       f = f_mm
       return
    end if

    call print('MetaPotential_ForceMixing_Calc: reweighting classical forces by '//mm_reweight//' in active region', VERBOSE)
    do i=1,at%N
       if (hybrid_mark(i) /= HYBRID_ACTIVE_MARK) cycle
       f_mm(:,i) = mm_reweight*f_mm(:,i)
    end do

    ! Make embed and fit lists if we need them
    if (method(1:4) == 'lotf' .or. trim(method) == 'conserve_momentum') then

       if ((method(1:4) == 'lotf' .and. lotf_do_init) .or. trim(method) == 'conserve_momentum') then
          call create_embed_and_fit_lists(at, fit_hops, this%embedlist, this%fitlist, &
               nneighb_only=nneighb_only, min_images_only=.true.)
       end if

       ! Make some convenient arrays of embed and fit list
       allocate(embed(this%embedlist%N))
       embed = int_part(this%embedlist,1)
       allocate(fit(this%fitlist%N))
       fit = int_part(this%fitlist,1)
    end if

    ! Do the expensive QM calculation. For the LOTF methods we only do this if lotf_do_qm is true
    if (method(1:4) /= 'lotf' .or. (method(1:4) == 'lotf' .and. lotf_do_qm)) then

       call initialise(params)
       call read_string(params, qm_args_str)

       !! TODO - possibly want to set more default options in the qm_args_str here
       if (.not. has_key(params, 'buffer_hops')) call set_value(params, 'buffer_hops', buffer_hops)

       ! We may have to rescale the cluster
       if (do_rescale_r) then
          call set_value(params, 'do_rescale_r', .true.)
          call set_value(params, 'r_scale', this%r_scale_pot1)
       end if

       qm_args_str = write_string(params)
       call finalise(params)

       ! Do the QM. If qm_args_str contatins 'small_cluster' or 'little_clusters' options
       ! then potential calc() will do cluster carving. 
       ! We pass along our mpi context object to allow little clusters to be split over
       ! nodes.
       call calc(this%qmpot, at, f=f_qm, err=err, args_str=qm_args_str, mpi_obj=this%mpi)
    end if
    
    if (method(1:4) == 'lotf') then
       
#ifdef HAVE_LOTF
       
       if (trim(method) == 'lotf_adj_pot_svd' .or. trim(method) == 'lotf_adj_pot_minim') then

          if (trim(method) == 'lotf_adj_pot_svd') then
             AP_method = 'SVD'
          else
             AP_method = 'minim'
          end if
          
          ! Optimise the adjustable potential
          allocate(df(3,this%fitlist%N))
          df = 0.0_dp
          df(:,1:this%embedlist%N) = f_qm(:,embed) - f_mm(:,embed)
    
          if (lotf_do_init) then
             call print('Initialising adjustable potential with map='//lotf_do_map, VERBOSE)
             call adjustable_potential_init(at, this%fitlist, directionN=this%embedlist%N, &
                  method=AP_method, nnonly=nneighb_only, &
                  spring_hops=lotf_spring_hops, map=lotf_do_map)
          end if
            
          if (lotf_do_qm .and. lotf_do_fit) &
               call adjustable_potential_optimise(at, df, method=AP_method)

          if (lotf_do_interp) then
             call adjustable_potential_force(at, df, interp=lotf_interp, &
                  interp_space=lotf_interp_space, interp_order=lotf_interp_order,power=dV_dt)
          else
             call adjustable_potential_force(at, df, power=dV_dt)
          end if
       
          f = f_mm                   ! Start with classical forces
          f(:,fit) = f(:,fit) + df   ! Add correction in fit region
       
          deallocate(df)

       else if (trim(this%method) == 'lotf_adj_pot_sw') then

          call system_abort('lotf_adj_pot_sw no longer supported.')

          ! old style adj pot: SW with variable parameters

!!$          allocate(df(3,this%fitlist%N))
!!$       
!!$          ! old style method: we fit to QM force in embed region and MM force in outer fit region,
!!$          ! not to force difference
!!$          df = f_mm(:,fit)
!!$          df(:,1:this%embedlist%N) = f_qm(:,embed)
!!$
!!$          if (lotf_do_init) &
!!$               call adjustable_potential_sw_init(at, this%fitlist)
!!$
!!$          if (lotf_do_qm .and. lotf_do_fit) &
!!$               call adjustable_potential_sw_optimise(at, df)
!!$
!!$          if (lotf_do_interp) then
!!$             call adjustable_potential_force(at, df, interp=lotf_interp)
!!$          else
!!$             call adjustable_potential_force(at, df)
!!$          end if
!!$       
!!$          f = f_mm           ! start with classical forces
!!$          f(:,fit) = df      ! replace with fitted force
!!$
!!$          deallocate(df)

       end if

#else
       call system_abort('MetaPotential_ForceMixing_Calc: support for method '//trim(method)//' not compiled in')
#endif

    else if (trim(method) == 'conserve_momentum') then
       
       call verbosity_push(NORMAL)

       call print('Conserving momentum using fit list with '//this%fitlist%N//' atoms')

       select case(conserve_momentum_weight_method)
       case ('uniform')
          weight_method = UNIFORM_WEIGHT
          call print('conserve_momentum: uniform weighting', VERBOSE)
       case ('mass')
          weight_method = MASS_WEIGHT
          call print('conserve_momentum: using mass weighting', VERBOSE)
       case ('mass^2')
          weight_method = MASS2_WEIGHT
          call print('conserve_momentum: using mass squared weighting', VERBOSE)
       case ('user')
          weight_method = USER_WEIGHT
          call print('conserve_momentum: using user defined weighting', VERBOSE)
       case default
          call system_abort('MetaPotential_ForceMixing_Calc: unknown conserve_momentum_weight method: '//&
               trim(conserve_momentum_weight_method))
       end select

       if (weight_method == USER_WEIGHT) then
          if (.not. assign_pointer(at, 'conserve_momentum_weight', conserve_momentum_weight)) &
               call system_abort('MetaPotential_ForceMixing_Calc: missing property conserve_momentum_weight')
       end if

       allocate(df(3,this%fitlist%N),df_fit(3,this%fitlist%N))

       df = 0.0_dp
       df(:,1:this%embedlist%N) = f_qm(:,embed) - f_mm(:,embed)

       f_tot = sum(df,dim=2)
       call print('conserve_momentum: norm(sum of target forces) = '//round(norm(f_tot),15))
    
       w_tot = 0.0_dp
       do i = 1, this%fitlist%N
          select case(weight_method)
          case(UNIFORM_WEIGHT)
             weight = 1.0_dp
          case(MASS_WEIGHT)
             weight = at%mass(fit(i))
          case(MASS2_WEIGHT)
             weight = at%mass(fit(i))*at%mass(fit(i))
          case(USER_WEIGHT)
             weight = conserve_momentum_weight(fit(i))
          end select
          df_fit(:,i) = -weight * f_tot
          w_tot = w_tot + weight
       end do
       df_fit = (df_fit / w_tot) + df

       call print('conserve_momentum: norm(sum of    fit forces) = '//round(norm(sum(df_fit,dim=2)),15))

       ! Final forces are classical forces plus corrected QM forces
       f = f_mm
       f(:,fit) = f(:,fit) + df_fit

       call verbosity_pop()

       deallocate(df, df_fit)
       
    else if (method(1:12) == 'force_mixing') then

       if (.not. assign_pointer(at, 'weight_region1', weight_region1)) &
            call system_abort('MetaPotential_ForceMixing_Calc: missing weight_region1 property - try setting calc_weights=T in args_str')

       ! Straight forward force mixing using weight_region1 created by create_hybrid_weights() 
       do i=1,at%N
          f(:,i) = weight_region1(i)*f_qm(:,i) + (1.0_dp - weight_region1(i))*f_mm(:,i)
       end do

    else
       call system_abort('MetaPotential_ForceMixing_calc: unknown method '//trim(method))
    end if
       
    ! Save QM and MM forces and total force as properties of Atoms object
    if (save_forces) then
       if (.not. has_property(at, 'qm_force')) &
            call add_property(at, 'qm_force', 0.0_dp, n_cols=3)
       dummy = assign_pointer(at, 'qm_force', force_ptr)
       force_ptr = f_qm

       if (.not. has_property(at, 'mm_force')) &
            call add_property(at, 'mm_force', 0.0_dp, n_cols=3)
       dummy = assign_pointer(at, 'mm_force', force_ptr)
       force_ptr = f_mm

       if (.not. has_property(at, 'force')) &
            call add_property(at, 'force', 0.0_dp, n_cols=3)
       dummy = assign_pointer(at, 'force', force_ptr)
       force_ptr = f
    end if

    deallocate(f_mm,f_qm)

    if (allocated(embed)) deallocate(embed)
    if (allocated(fit))   deallocate(fit)
    
  end subroutine MetaPotential_ForceMixing_calc


  function MetaPotential_ForceMixing_cutoff(this)
    type(MetaPotential_ForceMixing), intent(in) :: this
    real(dp) :: metapotential_forcemixing_cutoff

    ! Return larger of QM and MM cutoffs
    metapotential_forcemixing_cutoff = max(cutoff(this%mmpot), cutoff(this%qmpot))

  end function MetaPotential_ForceMixing_cutoff


