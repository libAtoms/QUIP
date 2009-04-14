program makecrack

!% The 'makecrack' program is used to prepare a crack system
!% in the thin strip geometry periodic in the crack-front
!% direction and with vacuum surrounding the in-plane edges.
!% Since the cell thickness in the periodic direction is fixed,
!% plane strain conditions apply.
!% 
!% \begin{figure}cr
!%   \centering
!% 
!%   \includegraphics[scale=.4]{loading.eps}
!% 
!%   \caption[Displacement loading scheme] {\label{fig:loading}
!%   Displacement loading scheme. (a) An initial load is applied to an
!%   unstrained slab. The atoms in Region I are vertically shifted by
!%   $\delta$ since there is zero strain behind the tip of a relaxed
!%   crack. In Region III the strain is constant and equal to the far
!%   field value of $\delta/h$. The seed crack length is $c$. Across
!%   Region II the strain increases linearly. The lower panels show the
!%   atomistic configuration of a $1200\times400\times3.68$ \AA{}$^3$
!%   slab containing 90720 atoms with a seed crack of length 600 \AA{}
!%   (b) before and (c) after classical relaxation, with the insets
!%   showing Region II in more detail.}
!% \end{figure}
!% 
!% The external load on the model crack is applied in a way that is
!% equivalent to the `fixed grips' displacement boundary conditions.
!% To exactly mimic experiment, the loading would be applied simply by
!% displacing the top and bottom rows of atoms, then fixing these atoms
!% and allowing the system to relax.
!% However the time between loadings in experiments is of the order of
!% several minutes which is completely inaccessible to simulation.
!% We therefore need to perform the loading more smoothly to reduce the
!% time required for equilibration.
!% We also need to insert a seed crack into the system, so that the crack
!% tip is far from the edges of the simulated system and the thin strip
!% approximation applies.
!% 
!% Both of these requirements can be met by considering the equilibrium
!% stress distribution in a slab containing a stationary crack, then
!% choosing simple initial conditions that approximate this
!% configuration: far behind the tip the strain is zero, and far ahead it
!% approaches the applied loading.
!% 
!% Fig.~\ref{fig:loading} illustrates these conditions, and shows the
!% changes that take place in a classical geometry optimisation as a
!% result of this stress distribution.
!% The optimisation is quick since the initial conditions are a
!% reasonable approximation to the relaxed strain field.
!% 
!% Prior to straining the slab, care has to be taken to align the
!% vertical coordinate origin with the centre of an vertical bond
!% to ensure that the crack opens cleanly.
!% 
!% The amount of equilibration required after each loading can be reduced
!% by performing classical relaxations at strains $\epsilon$ and
!% $\epsilon+\alpha$ then computing the displacement field as
!% \begin{equation}
!%   \mathbf{u}_i = \mathbf{r}^{(\epsilon+\alpha)}_i - \mathbf{r}^{(\epsilon)}_i
!% \end{equation}
!% where $\mathbf{r}^{(\epsilon)}_i$ denotes the relaxed position of atom
!% $i$ under an external strain of $\epsilon$.
!% The load can then be incremented using the equation
!% \begin{equation} \label{eq:loading2}
!%   \mathbf{r}'_i = \mathbf{r}_i + \frac{\beta}{\alpha} \mathbf{u}_i
!% \end{equation}
!% where $\beta$ is the desired load increment.
!% The $\{\mathbf{u}_i\}$ only need to be evaluated once for a given
!% crack configuration.
!%
!% 
!%   The parameters used to control loading (all in the 'crack' stanza of the XML file)
!%   are as follows:
!%   \begin{itemize}
!%    \item {\bf G} --- the energy release rate $G$ in units of J/m$^2$ to use for the
!%    initial loading of the crack slab (converted to a strain corresponding to $\epsilon$
!%    in the equations above).
!%    \item {\bf loading} --- type of loading to apply; can be 'uniform', 'ramp', 'kfield'
!%    or 'interp_kield_uniform'.
!%    \item {\bf load_interp_length} --- if loading type is 'interp_kfield_uniform', the
!%    length in Angstrom over which to interpolate from $K$-field to uniform
!%    loading.
!%    \item {\bf ramp_length} --- if loading type is 'ramp', then length of the loading ramp
!%    \item {\bf ramp_end_G} --- value of $G$ at the end of the loading ramp
!%    \item {\bf initial_loading_strain} --- the initial strain increment used to
!%    generate the load field ($\alpha$ in the equations above)
!%    \item {\bf relax_loading_field} --- whether or not to relax the initial loading field
!%    \item {\bf G_increment} --- the increment of loading that is applied to the
!%    crack slab at the end of each load cycle by the 'crack' program.
!%  \end{itemize}
!% 
!%  When 'makecrack' is run it first strains the uniform crack slab using the value of the 'G'
!%  parameter. This configuration is then classically relaxed if 'relax_loading_field' is true.
!%  The strain is then incremented by an amount 'initial_loading_strain', and the configuration
!%  is re-relaxed if 'relax_loading_field' is true. The 'load' field is now generated as the difference
!%  between these two configurations, and this is written to the 'crack.xyz' (and 'crack.nc') output
!%  files. 
!% 
!%  When the 'crack' program wants to increment the strain, it takes this saved 'load' field
!%  and scales it by an appropriate prefactor so that after the loading
!%  the the new $G_1$ will be given by $G_1 = G_0 + \Delta G$, where $\Delta G$ is equal to 
!%  the 'G_increment' parameter.

  use libAtoms_module
  use QUIP_module

  use cracktools_module
  use crackparams_module
  use elasticity_module
  
  implicit none


  ! Objects
  type(Atoms), target :: bulk, crack_slab, crack_layer
  type(CrackParams) :: params
  type(Potential) :: classicalpot
  type(Metapotential) :: simple
  type(MPI_context) :: mpi_glob
  type(Inoutput) :: xmlfile, infile
  type(CInoutput) :: movie, netcdf
  type(Inoutput) :: checkfile

  ! Pointers into Atoms data structure
  real(dp), pointer, dimension(:,:) :: load, k_disp, u_disp
  integer, pointer, dimension(:) :: move_mask, nn, changed_nn, edge_mask
  
  real(dp), allocatable :: pos1(:,:), pos2(:,:), f(:,:)
  real(dp), dimension(6,6) :: c, c0
  real (dp), dimension(3,3) :: axes, lattice

  real (dp) :: maxy, miny, maxx, minabsy, shift, ydiff, mindiff, &
       l_crack_pos, r_crack_pos, width, height, strain, a, E, v, v2, uij(3), &
       orig_height, G1, K1, r, energy

  integer :: i, j, k, n_fixed, atom1, atom2, n, Z(2), steps, nargs

  character(len=STRING_LENGTH) :: stem, xmlfilename, xyzfilename

!  call initialise(mpi_glob)

  if (mpi_glob%active) then
     ! Same random seed for each process, otherwise velocites will differ
     call system_initialise (common_seed = .true., enable_timing=.true., mpi_all_inoutput=.false.)
     call print('MPI run with '//mpi_glob%n_procs//' processes')
  else
     call system_initialise( enable_timing=.true.)
     call print('Serial run')
  end if

  call initialise(params)

  nargs = cmd_arg_count()

  ! Print usage information if no arguments given
  if (nargs /= 1) then
     call Print('Usage: makecrack <stem>')
     call print('Where <stem>.xml is parameter XML file and <stem>.xyz is the crack slab output XYZ file')
     call print('')
     call print('Available parameters and their default values are:')
     call Print('')
     call Print(params)

     call system_finalise()
     stop
  end if

  call get_cmd_arg(1, stem)

  xmlfilename = trim(stem)//'.xml'
  xyzfilename = trim(stem)//'.xyz'

  call print_title('Initialisation')
  call print('Reading parameters from file '//trim(xmlfilename))
  call initialise(xmlfile,xmlfilename,INPUT)
  call read_xml(params,xmlfile)
  call verbosity_push(params%io_verbosity)   ! Set base verbosity
  call print(params)

  call print ("Initialising classical potential with args " // trim(params%classical_args) &
       // " from file " // trim(xmlfilename))
  call rewind(xmlfile)
  call initialise(classicalpot, params%classical_args, xmlfile, mpi_obj=mpi_glob)
  call finalise(xmlfile)
  call Print(classicalpot)

  call print('Initialising metapotential')
  call initialise(simple, 'Simple', classicalpot, mpi_obj=mpi_glob)
  call print(simple)

  ! Crack structure specific code
  if (trim(params%crack_structure) == 'graphene') then

     call print_title('Graphene Crack')
     
!!$     call graphene_elastic(simple, a, v, E)

     a = 1.42_dp
     v = 0.144_dp
     E = 344.0_dp

     call Print('graphene sheet, lattice constant a='//a)

     bulk = graphene_cubic(a)
    
     if (trim(params%crack_slab_filename).ne.'') then
        call Initialise(infile, trim(params%crack_slab_filename))
        call print('Reading atoms from input file')
        call read_xyz(crack_slab, infile)
     else
        call graphene_slab(crack_layer, a, params%crack_graphene_theta, &
             params%crack_width, params%crack_height)
        crack_layer%Z = 6

        if (abs(params%crack_graphene_theta - 0.0_dp) < 1e-3_dp) then
           call print('armchair sheet')
           shift = 0.61567821_dp
        else if (abs(params%crack_graphene_theta - pi/6.0_dp) < 1e-3_dp) then
           call print('zigzag sheet')
           shift = 0.53319064_dp
        end if

        lattice = crack_layer%lattice
        lattice(1,1) = lattice(1,1) + params%crack_vacuum_size
        lattice(2,2) = lattice(2,2) + params%crack_vacuum_size
        call atoms_set_lattice(crack_layer, lattice)

        do i=1,crack_layer%N
           crack_layer%pos(2,i) = crack_layer%pos(2,i) + shift
        end do
     endif

     ! Actual width and height differ a little from requested
     width = maxval(crack_layer%pos(1,:))-minval(crack_layer%pos(1,:))
     height = maxval(crack_layer%pos(2,:))-minval(crack_layer%pos(2,:))
     call Print('Actual slab dimensions '//width//' A x '//height//' A')

     ! Cut notch
     if (params%crack_graphene_notch_width  > 0.0_dp .and. &
         params%crack_graphene_notch_height > 0.0_dp) then
        call Print('Cutting notch with width '//params%crack_graphene_notch_width// &
             ' A, height '//params%crack_graphene_notch_height//' A.')
     
        i = 1
        do 
           if ((crack_layer%pos(2,i) < &
                -(0.5_dp*params%crack_graphene_notch_height/ &
                params%crack_graphene_notch_width*(crack_layer%pos(1,i)+width/2.0_dp)) + &
                params%crack_graphene_notch_height/2.0_dp) .and. &
                (crack_layer%pos(2,i) > &
                (0.5_dp*params%crack_graphene_notch_height/ &
                params%crack_graphene_notch_width*(crack_layer%pos(1,i)+width/2.0_dp)) - &
                params%crack_graphene_notch_height/2.0_dp)) then
              call remove_atoms(crack_layer, i)

              i = i - 1 ! retest
           end if
           if (i == crack_layer%N) exit
           i = i + 1
        end do
     end if

     crack_layer%lattice(3,3) = 10.0_dp
     crack_slab = crack_layer

     call atoms_set_cutoff(crack_slab, cutoff(classicalpot)+params%md_crust)
     call calc_connect(crack_slab)

     call Print('Graphene sheet contains '//crack_slab%N//' atoms.')

  else if (trim(params%crack_structure) == 'diamond'.or.trim(params%crack_structure) == 'bcc'.or.trim(params%crack_structure) == 'fcc') then

     call crack_parse_atomic_numbers(params%crack_element, Z)
     if(trim(params%crack_structure) == 'diamond') then 
       call print_title('Diamond Structure Crack')
       call diamond(bulk, params%crack_lattice_guess, Z)
     elseif(trim(params%crack_structure) == 'bcc') then
       call print_title('BCC Structure Crack')
       call bcc(bulk, params%crack_lattice_guess, Z(1))
       bulk%cutoff = cutoff(simple)
       bulk%use_uniform_cutoff = .true. 
     elseif(trim(params%crack_structure) == 'fcc') then
       call print_title('FCC Structure Crack')
       call fcc(bulk, params%crack_lattice_guess, Z(1))
       bulk%cutoff = cutoff(simple)
       bulk%use_uniform_cutoff = .true. 
     endif 
     call calc_elastic_constants(simple, bulk, c=c, c0=c0, relax_initial=.true., return_relaxed=.true.)

     call print('Relaxed elastic constants (GPa):')
     call print(c*GPA)
     call print('')
     call print('Unrelaxed elastic constants (GPa):')
     call print(c0*GPA)
     call print('')

     call print('Relaxed lattice')
     call print(bulk%lattice)
     a = bulk%lattice(1,1)

     call Print(trim(params%crack_element)//' crack: atomic number Z='//Z//&
          ', lattice constant a = '//a)
     call Print('Crack name '//params%crack_name)

     ! Parse crack name and make crack slab
     call crack_parse_name(params%crack_name, axes)

     ! Get elastic constants relevant for a pull in y direction
     E = Youngs_Modulus(C, axes(:,2))*GPA
     v = Poisson_Ratio(C, axes(:,2), axes(:,1))
     v2 = Poisson_Ratio(C, axes(:,2), axes(:,3))

     call Print('')
     call print('Youngs modulus E_y = '//E) 
     call print('Poisson ratio v_yx = '//v) 
     call print('Poisson ratio v_yz = '//v2) 
     call Print('')

     if (trim(params%crack_slab_filename).ne.'') then
        call Initialise(infile, trim(params%crack_slab_filename))
        call print('Reading atoms from input file')
        call read_xyz(crack_slab, infile)
     else
        call slab(crack_layer, axes,  a, params%crack_width, params%crack_height, 1, Z, &
             trim(params%crack_structure))
        call atoms_set_cutoff(crack_layer, cutoff(classicalpot)+params%md_crust)
        call supercell(crack_slab, crack_layer, 1, 1, params%crack_num_layers)
     endif

     call calc_connect(crack_slab)

     call Print('Slab contains '//crack_slab%N//' atoms.')

     ! Actual width and height differ a little from requested
     width = maxval(crack_slab%pos(1,:))-minval(crack_slab%pos(1,:))
     height = maxval(crack_slab%pos(2,:))-minval(crack_slab%pos(2,:))
     call Print('Actual slab dimensions '//width//' A x '//height//' A')

  else
     ! Add code here for other structures...
     
     call system_abort("Don't (yet!) know how to make cracks with structure "//trim(params%crack_structure))
     
  end if ! select on crack_structure 

  ! Save bulk cube (used for qm_rescale_r parameter in crack code)
  call print_xyz(bulk, trim(stem)//'_bulk.xyz')

  call set_value(crack_slab%params, 'OrigWidth', width)
  call set_value(crack_slab%params, 'OrigHeight', height)

  call set_value(crack_slab%params, 'YoungsModulus', E)
  call set_value(crack_slab%params, 'PoissonRatio', v)

  ! Open x and y surfaces, remain periodic in z direction (normal to plane)
  if (trim(params%crack_structure) /= 'graphene') then
     lattice = crack_slab%lattice
     lattice(1,1) = lattice(1,1) + params%crack_vacuum_size
     lattice(2,2) = lattice(2,2) + params%crack_vacuum_size
     call atoms_set_lattice(crack_slab, lattice)
  end if

  ! Add various properties to crack_slab
  call add_property(crack_slab, 'hybrid', 0)
  call add_property(crack_slab, 'hybrid_mark', HYBRID_NO_MARK)
  call add_property(crack_slab, 'changed_nn', 0)
  call add_property(crack_slab, 'move_mask', 0)
  call add_property(crack_slab, 'nn', 0)
  call add_property(crack_slab, 'old_nn', 0)
  call add_property(crack_slab, 'md_old_changed_nn', 0)
  call add_property(crack_slab, 'edge_mask', 0)
  call add_property(crack_slab, 'load', 0.0_dp, n_cols=3)
  call add_property(crack_slab, 'k_disp', 0.0_dp, n_cols=3)
  call add_property(crack_slab, 'uniform_disp', 0.0_dp, n_cols=3)

  call fix_pointers()

  call print_title('Fixing Atoms')

  ! Fix top and bottom edges - anything within crack_edge_fix_tol of ymax or ymin is fixed
  maxy = 0.0; miny = 0.0; maxx = 0.0
  do i=1, crack_slab%N
     if (crack_slab%pos(2,i) > maxy) maxy = crack_slab%pos(2,i)
     if (crack_slab%pos(2,i) < miny) miny = crack_slab%pos(2,i)
  end do

  move_mask = 1
  n_fixed = 0
  do i=1, crack_slab%N
     if ((abs(crack_slab%pos(2,i)-maxy) < params%crack_edge_fix_tol) .or. &
         (abs(crack_slab%pos(2,i)-miny) < params%crack_edge_fix_tol)) then
        move_mask(i) = 0
        n_fixed = n_fixed + 1
     end if
  end do
  call Print(crack_slab%N//' atoms. '//n_fixed//' fixed atoms.')
  
  call print_title('Aligning Seed Crack at y=0')

  ! Determine position of seed crack
  l_crack_pos = -width ! single ended crack
  r_crack_pos = -width/2.0_dp + params%crack_seed_length

  if (trim(params%crack_structure) == 'graphene') then
     r_crack_pos = -width
  end if

  ! Find an atom close to y=0
  minabsy = 1000.0_dp
  atom1 = -1; atom2 = -1
  mindiff = 1000.0_dp
  do i=1, crack_slab%N
     if (abs(crack_slab%pos(2,i)) < minabsy) then
        minabsy = abs(crack_slab%pos(2,i))
        atom1 = i
     end if
  end do

  ! Apply shift to centre the seed crack in the right place
  if (trim(params%crack_name) == '(111)[11b0]') then
     ! Find atom1's closest neighbour vertically above or below it (x and z equal, not y)
     do n = 1, atoms_n_neighbours(crack_slab, atom1)
        j = atoms_neighbour(crack_slab, atom1, n, diff=uij) ! nth neighbour of atom1
        if (abs(uij(1)) < 1e-4_dp .and. & 
            abs(uij(2)) > 1e-4_dp .and. &
            abs(uij(3)) < 1e-4_dp) then

	    ydiff = abs(crack_slab%pos(2,atom1)-crack_slab%pos(2,j))
	    if (ydiff < mindiff) then
	      mindiff = ydiff
	      atom2 = j
	    end if
         end if
      end do
      
      if (atom1 == -1 .or. atom2 == -1) &
           call system_abort('Failed to find a pair of atoms vertically aligned!')

      ! Align y=0 to centre line of atom1-atom2 bond
      shift = (crack_slab%pos(2,atom1) + crack_slab%pos(2,atom2))/2.0_dp

      call Print('Centering on (atom '//atom1//')--(atom '//atom2//') bond')
      call print('  Atom 1 pos = '//crack_slab%pos(:,atom1))
      call print('  Atom 2 pos = '//crack_slab%pos(:,atom2))
      call Print('Shifting atoms vertically by '//shift)
      do i=1,crack_slab%N
         crack_slab%pos(2,i) = crack_slab%pos(2,i) + shift
      end do

  else if(trim(params%crack_name) == '(110)[11b0]') then
     ! Align y=0 to atom1
     shift = -crack_slab%pos(2,atom1)
     
     call Print('Centering on atom '//atom1)
     call print('  Atom 1 pos = '//crack_slab%pos(:,atom1))
     call Print('Shifting atoms vertically by '//shift)
     do i=1,crack_slab%N
        crack_slab%pos(2,i) = crack_slab%pos(2,i) + shift
     end do

  else if (trim(params%crack_name) == '(110)[001b]') then
     ! Do nothing - correctly aligned already
  else if (trim(params%crack_name) == '(100)(010)') then
     ! Do nothing - correctly aligned already
  else if (trim(params%crack_structure) == 'graphene') then
     ! Do nothing - correctly aligned already
  else
     ! Get shift from params
     do i=1,crack_slab%N
        crack_slab%pos(2,i) = crack_slab%pos(2,i) + params%crack_y_shift
     end do
  end if


  if (params%io_netcdf) then
     call initialise(movie, 'makecrack_movie.nc', action=OUTPUT)
  else
     call initialise(movie, 'makecrack_movie.xyz', action=OUTPUT)
  end if

  if (mpi_glob%active .and. params%crack_relax_loading_field) then
     allocate(f(3,crack_slab%N))
     call setup_parallel(classicalpot, crack_slab, e=energy, f=f, args_str=params%classical_args_str)
     deallocate(f)
  end if

  ! Apply initial load
  if (trim(params%crack_loading) == 'uniform') then
     if (params%crack_G > 0.0_dp) then

        call print_title('Applying Uniform Load')
        
        call crack_uniform_load(crack_slab, l_crack_pos, r_crack_pos, &
             params%crack_strain_zone_width, params%crack_G, apply_load=.true.)

        if(params%crack_rescale_x_z) then
       !  Rescale in x direction by v and in z direction by v2
          if (.not. get_value(crack_slab%params,'OrigHeight',orig_height)) orig_height = 0.0_dp
          strain = crack_g_to_strain(params%crack_G, E, v, orig_height)
          crack_slab%pos(1,:) = crack_slab%pos(1,:)*(1.0_dp-v*strain)
          crack_slab%pos(3,:) = crack_slab%pos(3,:)*(1.0_dp-v2*strain)
          crack_slab%lattice(3,3) = crack_slab%lattice(3,3)*(1.0_dp-v2*strain)
          call atoms_set_lattice(crack_slab, crack_slab%lattice)
        elseif(params%crack_rescale_x) then 
       !  Rescale in x direction by v 
          if (.not. get_value(crack_slab%params,'OrigHeight',orig_height)) orig_height = 0.0_dp
          strain = crack_g_to_strain(params%crack_G, E, v, orig_height)
          crack_slab%pos(1,:) = crack_slab%pos(1,:)*(1.0_dp-v*strain)
        endif

     end if

     call Print_title('Initialising QM region for unrelaxed system')
     call setup_crack_marks(crack_slab, params)
     call crack_print(crack_slab, movie, params, mpi_glob)
        
     ! Save positions
     pos1 = crack_slab%pos

     call Print_title('Generating Load Field')

     call Print('Initial strain = '//params%crack_initial_loading_strain)

     allocate(pos1(3,crack_slab%N))

     if (params%crack_relax_loading_field) then
        ! Geometry optimise
        steps = minim(simple, crack_slab, params%minim_mm_method, params%minim_mm_tol, params%minim_mm_max_steps, &
             params%minim_mm_linminroutine, do_print=.true., do_pos=.true.,do_lat=.false., &
             use_fire=trim(params%minim_mm_method)=='fire', &
             print_cinoutput=movie)

        call Print_title('Initialising QM region for relaxed system')
        call setup_crack_marks(crack_slab, params)
        call crack_print(crack_slab, movie, params, mpi_glob)

     end if

     ! strain it a bit
     do i=1,crack_slab%N
        crack_slab%pos(2,i) = crack_slab%pos(2,i)*(1.0_dp+params%crack_initial_loading_strain)
     end do
        
     if (params%crack_relax_loading_field) then

        ! now re-relax
        steps = minim(simple, crack_slab, params%minim_mm_method, params%minim_mm_tol, params%minim_mm_max_steps, &
             params%minim_mm_linminroutine, do_print=.true., do_pos=.true.,do_lat=.false., &
             use_fire=trim(params%minim_mm_method)=='fire', &
             print_cinoutput=movie)
     end if

     ! work out displacement field
     do i=1,crack_slab%N
        load(:,i) = crack_slab%pos(:,i) - pos1(:,i)
     end do
     
     ! Restore pos1 - relaxed positions at initial load
     crack_slab%pos = pos1

     call print('Displacement field generated. Max disp: '//maxval(load))
     call print('                              RMS disp: '//sqrt(norm2(reshape(load,(/3*crack_slab%N/)))/(3.0_dp*crack_slab%N)))

  else if (trim(params%crack_loading) == 'ramp') then

     call print_title('Applying Loading Ramp')

     call crack_apply_strain_ramp(crack_slab, params%crack_G, params%crack_ramp_end_G, r_crack_pos, &
          r_crack_pos+params%crack_strain_zone_width, &
          r_crack_pos+params%crack_strain_zone_width+params%crack_ramp_length)

     call Print_title('Initialising QM region for unrelaxed system')
     call setup_crack_marks(crack_slab, params)
     call crack_print(crack_slab, movie, params, mpi_glob)

  else if (trim(params%crack_loading) == 'kfield') then
     
     call print_title('Applying Irwin K-field Loading')

     if (.not. get_value(crack_slab%params,'OrigHeight',orig_height)) orig_height = 0.0_dp

     call print('Initial stress intesity factor K_0 = '//crack_g_to_k(params%crack_G, E, v)/1e6_dp//' MPa.sqrt(m)')
     
     call set_value(crack_slab%params, 'CrackPos', r_crack_pos + 0.85_dp*params%crack_strain_zone_width)
     call crack_k_field(crack_slab, crack_g_to_k(params%crack_G, E, v))

     call fix_pointers()

     allocate(pos1(3,crack_slab%N),pos2(3,crack_slab%N))

     ! Save bulk positions
     pos1 = crack_slab%pos

     do i=1,crack_slab%N
        crack_slab%pos(:,i) = crack_slab%pos(:,i) + k_disp(:,i)
     end do

     call Print_title('Initialising QM region for unrelaxed system')
     call setup_crack_marks(crack_slab, params)
     call crack_print(crack_slab, movie, params, mpi_glob)

     ! relax initial loading
     if (params%crack_relax_loading_field) then
        steps = minim(simple, crack_slab, params%minim_mm_method, params%minim_mm_tol, params%minim_mm_max_steps, &
             params%minim_mm_linminroutine, do_print=.true., do_pos=.true.,do_lat=.false., &
             use_fire=trim(params%minim_mm_method)=='fire', &
             print_cinoutput=movie)
        call Print_title('Initialising QM region for relaxed system')
        call setup_crack_marks(crack_slab, params)
        call crack_print(crack_slab, movie, params, mpi_glob)

     end if

     pos2 = crack_slab%pos

     ! apply load increment
     K1 = crack_g_to_k(crack_strain_to_g( &
          crack_g_to_strain(params%crack_G, E, v, orig_height) + &
          params%crack_initial_loading_strain, E, v, orig_height),E,v)

     call print('Stress Intensity Factor K_1 = '//(K1/1e6_dp)//' MPa.sqrt(m)')
     call crack_k_field(crack_slab, K1)
     call fix_pointers()

     do i=1,crack_slab%N
        crack_slab%pos(:,i) = pos1(:,i) + k_disp(:,i)
     end do

     ! relax again at new load
     if (params%crack_relax_loading_field) then
        steps = minim(simple, crack_slab, params%minim_mm_method, params%minim_mm_tol, params%minim_mm_max_steps, &
             params%minim_mm_linminroutine, do_print=.true., do_pos=.true.,do_lat=.false., &
             use_fire=trim(params%minim_mm_method)=='fire', &
             print_cinoutput=movie)
     end if
       
     do i=1,crack_slab%N
        load(:,i) = crack_slab%pos(:,i) - pos2(:,i)
     end do

     call print('Displacement field generated. Max disp: '//maxval(load))
     call print('                              RMS disp: '//sqrt(norm2(reshape(load,(/3*crack_slab%N/)))/(3.0_dp*crack_slab%N)))

     crack_slab%pos = pos2

     deallocate(pos1,pos2)

  else if (trim(params%crack_loading) == 'interp_kfield_uniform') then
     
     ! interpolate linearly between K field (near tip) and uniform loading (near edge)

     call print_title('Applying Combined K-field and Uniform Loading')

     call print('Interpolation length '//params%crack_load_interp_length//' A')

     if (.not. get_value(crack_slab%params,'OrigHeight',orig_height)) orig_height = 0.0_dp

     call print('Initial energy release rate    G_0 = '//params%crack_G//' J/m^2')
     call print('Initial stress intesity factor K_0 = '//crack_g_to_k(params%crack_G, E, v)/1e6_dp//' MPa.sqrt(m)')
     
     call set_value(crack_slab%params, 'CrackPos', r_crack_pos + 0.85_dp*params%crack_strain_zone_width)
     call crack_k_field(crack_slab, crack_g_to_k(params%crack_G, E, v))

     call crack_uniform_load(crack_slab, l_crack_pos, r_crack_pos, &
          params%crack_strain_zone_width, params%crack_G, apply_load=.false.)

     call fix_pointers()

     allocate(pos1(3,crack_slab%N),pos2(3,crack_slab%N))

     ! Save bulk positions
     pos1 = crack_slab%pos

     do i=1,crack_slab%N
        r = sqrt((crack_slab%pos(1,i) - (r_crack_pos + 0.85*params%crack_strain_zone_width))**2.0_dp + &
             crack_slab%pos(2,i)**2.0_dp)
        if (r > params%crack_load_interp_length) then
           crack_slab%pos(:,i) = pos1(:,i) + u_disp(:,i)
        else
           do k=1,3
              crack_slab%pos(k,i) = pos1(k,i) +  &
                   linear_interpolate(0.0_dp, k_disp(k,i), params%crack_load_interp_length, u_disp(k,i), r)
           end do
        end if
     end do

     call Print_title('Initialising QM region for unrelaxed system')
     call setup_crack_marks(crack_slab, params)
     call crack_print(crack_slab, movie, params, mpi_glob)

     ! relax initial loading
     if (params%crack_relax_loading_field) then
        steps = minim(simple, crack_slab, params%minim_mm_method, params%minim_mm_tol, params%minim_mm_max_steps, &
             params%minim_mm_linminroutine, do_print=.true., do_pos=.true.,do_lat=.false., &
             use_fire=trim(params%minim_mm_method)=='fire', &
             print_cinoutput=movie)
        call Print_title('Initialising QM region for relaxed system')
        call setup_crack_marks(crack_slab, params)
        call crack_print(crack_slab, movie, params, mpi_glob)
     end if

     pos2 = crack_slab%pos

     ! apply load increment to u_disp
     G1 = crack_strain_to_g( &
          crack_g_to_strain(params%crack_G, E, v, orig_height) + &
          params%crack_initial_loading_strain, E, v, orig_height)

     call crack_uniform_load(crack_slab, l_crack_pos, r_crack_pos, &
          params%crack_strain_zone_width, G1, apply_load=.false.)

     ! apply load increment to K_disp
     K1 = crack_g_to_k(G1,E,v)

     call print('Energy release rate     G_1 = '//G1//' J/m^2')
     call print('Stress Intensity Factor K_1 = '//(K1/1e6_dp)//' MPa.sqrt(m)')
     call crack_k_field(crack_slab, K1, do_sig=.false., do_disp=.true.)
     call fix_pointers()

     do i=1,crack_slab%N
        r = sqrt((crack_slab%pos(1,i) - (r_crack_pos + 0.85*params%crack_strain_zone_width))**2.0_dp + &
             crack_slab%pos(2,i)**2.0_dp)
        if (r > params%crack_load_interp_length) then
           crack_slab%pos(:,i) = pos1(:,i) + u_disp(:,i)
        else
           do k=1,3
              crack_slab%pos(k,i) = pos1(k,i) +  &
                   linear_interpolate(0.0_dp, k_disp(k,i), params%crack_load_interp_length, u_disp(k,i), r)
           end do
        end if
     end do

     ! relax again at new load
     if (params%crack_relax_loading_field) then
        steps = minim(simple, crack_slab, params%minim_mm_method, params%minim_mm_tol, params%minim_mm_max_steps, &
             params%minim_mm_linminroutine, do_print=.true., do_pos=.true.,do_lat=.false., &
             use_fire=trim(params%minim_mm_method)=='fire', &
             print_cinoutput=movie)
     end if
       
     do i=1,crack_slab%N
        load(:,i) = crack_slab%pos(:,i) - pos2(:,i)
     end do

     call print('Displacement field generated. Max disp: '//maxval(load))
     call print('                              RMS disp: '//sqrt(norm2(reshape(load,(/3*crack_slab%N/)))/(3.0_dp*crack_slab%N)))

     crack_slab%pos = pos2

     deallocate(pos1,pos2)

  else
     call system_abort('Unknown loading type '//trim(params%crack_loading))
  end if

  call Print_title('Initialising QM region')

  call setup_crack_marks(crack_slab, params)
  call crack_print(crack_slab, xyzfilename, params, mpi_glob)

  if (params%io_netcdf) then
     call initialise(netcdf, trim(stem)//'.nc', action=OUTPUT)
     call crack_print(crack_slab, netcdf, params, mpi_glob)
     call finalise(netcdf)
  else
     if (.not. mpi_glob%active .or. (mpi_glob%active .and. mpi_glob%my_proc == 0)) then
        call initialise(checkfile, trim(stem)//'.check', OUTPUT, isformatted=.false.)
        call write_binary(crack_slab, checkfile)
        call finalise(checkfile)
     end if
  end if

  call verbosity_pop() ! Return to default verbosity

  call finalise(bulk)
  call finalise(crack_slab)
  call finalise(crack_layer)
  call finalise(classicalpot)
  call finalise(simple)
  call system_finalise()

contains

  subroutine setup_crack_marks(crack_slab, params)
    type(Atoms) :: crack_slab
    type(CrackParams) :: params

    integer :: i
    real(dp) :: crack_pos

    ! Setup edge_mask to allow easy exclusion of edge atoms
    do i=1,crack_slab%N
       if (crack_is_edge_atom(crack_slab, i, params%selection_edge_tol)) &
	    edge_mask(i) = 1
    end do

    ! Calculate connectivity and numbers of nearest neighbours
    call crack_update_connect(crack_slab, params)

    ! Find rightmost undercoordinated atoms in bulk - this is initial crack tip position
    crack_pos = crack_find_crack_pos(crack_slab, params)

    ! Artificially set changed_nn to 1 for atoms near to crack tip
    do i = 1, crack_slab%N
       if (distance_min_image(crack_slab, i, (/crack_pos,0.0_dp,0.0_dp/)) < params%crack_seed_embed_tol) &
	    changed_nn(i) = 1
    end do

    call Print('Seeded embed region with '//count(changed_nn /= 0)//' atoms.')

    call crack_update_selection(crack_slab, params)

  end subroutine


  subroutine fix_pointers()

    if (.not. assign_pointer(crack_slab, 'nn', nn)) &
         call system_abort('nn pointer assignment failed')

    if (.not. assign_pointer(crack_slab, 'changed_nn', changed_nn)) &
         call system_abort('changed_nn pointer assignment failed')
  
    if (.not. assign_pointer(crack_slab, 'load', load)) &
         call system_abort('load pointer assignment failed')

    if (.not. assign_pointer(crack_slab, 'move_mask', move_mask)) &
         call system_abort('move_mask pointer assignment failed')

    if (.not. assign_pointer(crack_slab, 'edge_mask', edge_mask)) &
         call system_abort('edge_mask pointer assignment failed')

    if (.not. assign_pointer(crack_slab, 'uniform_disp', u_disp)) &
         call system_abort('u_disp pointer assignment failed')

    if (.not. assign_pointer(crack_slab, 'k_disp', k_disp)) &
         call system_abort('k_disp pointer assignment failed')

  end subroutine fix_pointers


end program makecrack
