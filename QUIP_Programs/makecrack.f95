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
  type(Inoutput) :: xmlfile
  type(CInoutput) :: netcdf
  type(Inoutput) :: checkfile

  ! Pointers into Atoms data structure
  real(dp), pointer, dimension(:,:) :: load
  integer, pointer, dimension(:) :: move_mask, nn, changed_nn, edge_mask, md_old_changed_nn, &
       old_nn, hybrid, hybrid_mark

  
  real(dp), allocatable :: f(:,:), u_disp(:,:), k_disp(:,:)
  real(dp), dimension(6,6) :: c
  real (dp), dimension(3,3) :: lattice

  real (dp) :: maxy, miny, maxx, minabsy, shift, ydiff, mindiff, G, &
       width, height, a, E, v, v2, uij(3), energy

  integer :: i, j,  n_fixed, atom1, atom2, n, Z(2), nargs

  character(len=STRING_LENGTH) :: stem, xmlfilename, xyzfilename

  call initialise(mpi_glob)

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
  call crack_make_slab(params, classicalpot, simple, crack_slab, width, height, E, v, v2, bulk)

  ! Save bulk cube (used for qm_rescale_r parameter in crack code)
  call print_xyz(bulk, trim(stem)//'_bulk.xyz')

  call set_value(crack_slab%params, 'OrigWidth', width)
  call set_value(crack_slab%params, 'OrigHeight', height)

  call set_value(crack_slab%params, 'YoungsModulus', E)
  call set_value(crack_slab%params, 'PoissonRatio_yx', v)
  call set_value(crack_slab%params, 'PoissonRatio_yz', v2)

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

  call crack_fix_pointers(crack_slab, nn, changed_nn, load, move_mask, edge_mask, md_old_changed_nn, &
       old_nn, hybrid, hybrid_mark)

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

  call crack_make_seed(crack_slab, params, u_disp, k_disp)
  !call crack_setup_marks(crack_slab, params) 

  if (params%crack_apply_initial_load) then
     if (mpi_glob%active .and. params%crack_relax_loading_field) then
        allocate(f(3,crack_slab%N))
        call setup_parallel(classicalpot, crack_slab, e=energy, f=f, args_str=params%classical_args_str)
        deallocate(f)
     end if
     
     ! Apply initial load
     call crack_calc_load_field(crack_slab, params, classicalpot, simple, params%crack_loading, overwrite_pos=.true., &
          mpi=mpi_glob)
  end if

  call Print_title('Initialising QM region')

  call crack_setup_marks(crack_slab, params)
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

end program makecrack
