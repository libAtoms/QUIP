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

program crack


  !% The 'crack' program can either load a crack configuration that has
  !% previosuly been created with 'makecrack' or continue a previous
  !% simulation, either from an XYZ file or from a binary checkpoint file.
  !% 
  !% There are several distinct modes of operation for the crack code,
  !% corresponding to different simulation tasks, as listed in the table
  !% below. Each of these is described below.
  !%
  !% \begin{center}
  !% \begin{tabular}{cc}
  !% \hline
  !% \hline
  !% 'simulation_task' & Description \\
  !% \hline
  !% 'md' & Molecular dynamics \\
  !% 'minim' & Hybrid structural relaxation \\
  !% 'force_integraion' & Hybrid force integration \\
  !% 'quasi_static' & Quasi-static loading \\
  !% \hline
  !% \hline
  !% \end{tabular}
  !% \end{center}
  !%
  !% \subsection{Molecular Dynamics}
  !%
  !% We start with a relaxed configuration at a load $G$ below the critical
  !% Griffith load $G_c$.
  !% Depending on the range of lattice trapping, it should be possible to
  !% find a load that satisfies $G_{-} < G < G_c$ so the relaxed crack is
  !% lattice trapped and will not move until the load is increased above
  !% $G_+$.
  !% 
  !% Molecular dynamics is then carried out at a temperature of 300~K, with
  !% a weak Langevin thermostat to correct for the energy drift caused by
  !% the time dependence of the LOTF Hamiltonian.
  !% The dynamics can be accelerated using the predictor-corrector scheme.
  !% This is specified using the 'md_time_step' and 'md_extrapolate_steps'
  !% parameters; if the latter is set to one then no extrapolation is carried
  !% out.
  !% As an example, extrapolation for 10 steps of $\Delta t=1$~fs is possible for
  !% silicon crack systems.
  !%
  !% \begin{figure}[b]
  !%   \centering
  !%   \includegraphics[width=12cm]{flow-chart}
  !%   \caption[Molecular dynamics methodology flow chart]
  !%   {\label{fig:flow-chart} Flow chart illustrating molecular dynamics
  !%   methodology used for the fracture simulations. See text for a
  !%   description of each state and the conditions that have to met for
  !%   transitions to take place.}
  !% \end{figure}
  !% 
  !% Unfortunately it is not possible to run dynamics for long enough to
  !% fully explore the environment at each load and to cross barriers which
  !% the real system would have sufficient thermal energy to pass over at
  !% 300~K.
  !% %
  !% Instead we allow the dynamics to proceed for some fixed amount of time
  !% 'md_wait_time', and then periodically check for rebonding near the
  !% crack tip using the time-averaged coordinates with an averaging time of 
  !% 'md_avg_time'.
  !% 
  !% Fig.~\ref{fig:flow-chart} schematically illustrates the details of the
  !% molecular dynamics methodology used for the fracture simulations.
  !% If no rebonding occurs for some time then we increase the load.
  !% After each loading increment there is a thermalisation period in which
  !% a stronger thermostat is used to dissipate the energy produced by the
  !% rescaling. 
  !% The thermalisation continues until the fluctuations in temperature
  !% are small, defined by the inequality
  !% \begin{equation}
  !% \frac{T - \left<T\right>}{T} < \frac{1}{\sqrt{N}}
  !% \end{equation}
  !% where $T$ and $\left<T\right>$ are the instantaneous and average
  !% temperatures and $N$ is the total number of atoms in the simulation.
  !% Once this condition is satisfied the thermostat is turned down and we
  !% return to almost microcanonical molecular dynamics.
  !% After the crack has started to move, the rebonding is automatically
  !% detected and the load is not increased further.
  !%
  !%\subsection{Molecular Dynamics -- smooth loading}
  !%
  !% There is now an alternative way to smoothly increase the load during
  !% an MD simulation. To use this method set the 'md_smooth_loading_rate'
  !% parameter to a non-zero value. This parameter causes the load to be 
  !% increased smoothly by an amount 'md_smooth_loading_rate*md_time_step'
  !% after each Verlet integration step. Once the crack starts to move
  !% (defined by the tip position changing by more than 'md_smooth_loading_tip_move_tol'
  !% from the original position), then loading will stop. If the crack arrests
  !% (defined by a movement of less than 'md_smooth_loading_tip_move_tol' in
  !% a time 'md_smooth_loading_arrest_time' then loading will recommence.
  !% 
  !% At the moment loading always uses the same 'load' field, but it is
  !% planned to allow the loading field to be recomputed from time to time
  !% if the crack moves and subsequently arrests.                                          
  !% 
  !% \subsection{QM Selection Algorithm}
  !% 
  !% A major advantage of LOTF is that it allows us to keep the QM region
  !% small.
  !% This requires a robust selection algorithm to follow the crack tip as
  !% it moves and identify the atoms that need to be treated with quantum
  !% mechanical accuracy.
  !% This is a difficultproblem since the timescales of thermal vibration 
  !% and crack motion are not well separated.
  !% The hysteretic selection algorithm described in Section 4.6 of my thesis
  !% \footnote{Available at 'http://www.srcf.ucam.org/\~jrk33/Publications'}.
  !% provides an effective solution to this problem.
  !% 
  !% We identify atoms as active when they change their bonding topology,
  !% and then construct embedding ellipses around each active atom.
  !% The set of active atoms is seeded with a few atoms near to the crack
  !% tip at the start of the simulation.
  !% Ellipses are used rather than spheres to allow the embedding region to be biased
  !% forwards so that the QM region always extends ahead of the crack tip.
  !% Fig.~\ref{fig:qm-selection-crack} illustrates how the algorithm works
  !% in a simple case with only two active atoms --- in reality there could
  !% be several hundred.
  !% Ellipses with different radii are used to define inner and outer
  !% selection regions, and then the hysteretic algorithm ensures that
  !% atoms near the edges of the QM region do not oscillate in and out of
  !% the active region.
  !% 
  !% \begin{figure}
  !%   \centering
  !%   \includegraphics[width=120mm]{qm-selection-crack}
  !%   \caption[Hysteretic QM selection algorithm for crack tip]{
  !%     \label{fig:qm-selection-crack} Hysteretic QM selection algorithm
  !%     applied to crack tip region. The red and blue atoms are considered
  !%     `active', and are used to define inner (left panel) and outer
  !%     (right panel) selection regions. The atom indicated with the black
  !%     arrow remains selected despite oscillating in and out of the inner
  !%     region providing that it stays inside the outer region.  }
  !%     
  !% \end{figure}
  !% 
  !% As the crack moves on, we can stop treating atoms behind the crack tip
  !% quantum mechanically.
  !% We cap the size of the QM region at 'selection_max_qm_atoms'
  !% based on our computational capability --- this can be several hundred
  !% atoms for a tight binding simulation, or of the order of a hundred for
  !% an \emph{ab initio} simulation.
  !% By keeping track of the order in which atoms became active, we can
  !% remove them from the QM region in a consistent fashion.
  !% An additional condition prevents atoms further than a threshold
  !% distance 'selection_cutoff_plane' away from the centre of mass of the 
  !% current QM region from becoming active.
  !% 
  !% \subsection{Hybrid Structural Relaxation}
  !%
  !% Classical geometry optimisation can be performed rapidly using the
  !% conjugate gradients algorithm, and this provides a good first
  !% approximation to the hybrid relaxed geometry, which in turn
  !% approximates the atomic configuration that would be found by a fully
  !% quantum mechanical optimisation of the entire system.
  !% Relaxation using the hybrid forces is slightly more involved.
  !% Whenever forces are needed, a QM calculation is performed and the
  !% adjustable potential parameters are optimised to reproduce the QM
  !% forces.
  !% The forces used for the geometry optimisation are the sum of the
  !% classical and adjustable potential forces: as for MD, this ensures
  !% that there is no mechanical incompatibility at the boundary between
  !% classical and quantum regions.
  !% 
  !% The standard formulation of the conjugate gradients algorithm requires
  !% both an objective function and its derivative.
  !% For geometry optimisation, the former is the total energy.
  !% In the LOTF hybrid scheme, there is no well defined total energy, so
  !% this approach is not possible.
  !% It is possible to modify the line minimisation step of the
  !% conjugate gradient algorithm to work with derivative information only.
  !% This is done by extrapolating the projection of the derivative along
  !% the search direction to zero. This is invoked by setting 'minim_linmin_routine' 
  !% to 'LINMIN_DERIV'.
  !% To avoid errors associated with large extrapolation lengths, a maximum
  !% step size is specified and then the procedure is iterated until the
  !% projected derivative is zero.
  !% 
  !% \subsection{Quasi-static loading}
  !% 
  !% As well as MD, the 'crack' code can perform quasi-static
  !% simulations where the system is fully relaxed after each increment of
  !% load.
  !% This approach can be used to estimate the fracture toughness $G\_+$,
  !% the load at which the lattice trapping barrier first falls to zero.
  !% For this purpose, we consider fracture to have occured when the
  !% crack tip moves more than 'quasi_static_tip_move_tol' from its
  !% original position.
  !% 
  !% \subsection{Hybrid Force Integration}
  !% 
  !% Within the LOTF scheme, there is not a meaningful total energy for the
  !% hybrid classical/quantum system; the method matches forces between the
  !% QM and classical regions, rather than energies, so the solution is to
  !% use these accurate forces to evaluate the energy difference between
  !% two configurations by force integration:
  !% \begin{equation} \label{eq:force-integration}
  !%   \Delta E = \int_{\gamma} \mathbf{F}\left(\mathbf{R}\right)\cdot\mathrm{d}\mathbf{R}
  !% \end{equation}
  !% where $\mathbf{R}$ denotes all atomic positions and the integration
  !% contour $\gamma$ can be any path between the two configurations of
  !% interest.
  !% The start and end configurations should be obtained by hybrid minimisation, 
  !% so that $\mathbf{F} = 0$ at both integration limits.
  !% The ending configuation is taken from the file specified in the
  !% 'force_integration_end_file' parameter.
  !% 
  !% The simplest contour connecting the two minima is used: linear interpolation 
  !% between the relaxed unreconstructed state with atomic coordinates $\mathbf{R}_U$ and 
  !% the reconstructed state with coordinates $\mathbf{R}_R$.
  !% The QM region is fixed during the integration process.
  !% The forces $\mathbf{F}(\mathbf{R})$ are calculated in the standard
  !% LOTF fashion: classical and quantum forces are evaluated for the two
  !% regions, then the adjustable potential is optimised to reproduce the
  !% quantum forces.
  !% The forces used for the integration are the sum of the classical and
  !% corrective forces to ensure that there is no mechanical mismatch at
  !% the boundary.
  !% The integration path is discretised into $N+1$ samples according to
  !% \begin{eqnarray}
  !%   \Delta \mathbf{R} & = & \frac{1}{N}\left(\mathbf{R}_R - \mathbf{R}_U\right) \\
  !%   \mathbf{F}_i & = & \mathbf{F}(\mathbf{R}_U + i\,\Delta \mathbf{R}), \;\; 0 \le i \le N
  !% \end{eqnarray}
  !% and then Eq.~\ref{eq:force-integration} can be evaluated using
  !% Simpson's Rule:
  !% \begin{eqnarray}
  !%   \Delta E & \approx & \frac{\Delta \mathbf{R}}{3} \cdot \left[
  !%     \mathbf{F}_0 + 
  !%     2 \sum_{j=1}^{N/2-1} \mathbf{F}_{2j} +
  !%     4 \sum_{j=1}^{N/2} \mathbf{F}_{2j-1} +
  !%     \mathbf{F}_N \right]
  !% \end{eqnarray}
  !% The step-size $\Delta \mathbf{R}$ required for accurate integration of
  !% the energy difference can be calibrated using force integration with the
  !% classical potential alone, where the energy difference can be
  !% calculated exactly --- typically $N=20$ gives good results.
  !% The method has been validated by confirming that perturbing the integration
  !% path does not affect the value of $\Delta E$ obtained.
  !% $N$ is specified using the 'force_integration_n_steps' parameter.


  use libAtoms_module
  use Potential_module

  ! Crack includes
  use CrackTools_module
  use CrackParams_module

  implicit none

  ! Constants
  integer, parameter :: STATE_THERMALISE = 1
  integer, parameter :: STATE_MD = 2
  integer, parameter :: STATE_MD_LOADING = 3
  integer, parameter :: STATE_MD_CRACKING = 4
  integer, parameter :: STATE_DAMPED_MD = 5
  integer, parameter :: STATE_MICROCANONICAL = 6
  integer, parameter :: STATE_MD_CONSTANT = 7
  integer, parameter :: STATE_MD_NVE = 8
  character(len=14), dimension(8), parameter :: STATE_NAMES = &
       (/"THERMALISE    ", "MD            ", "MD_LOADING    ", "MD_CRACKING   ", "DAMPED_MD     ", "MICROCANONICAL", "MD_CONSTANT   ", "MD_NVE        "/)

  ! Objects
  type(InOutput) :: xmlfile
  type(CInOutput) :: movie, crackin, movie_backup
  type(DynamicalSystem), target :: ds, ds_save
  type(Atoms) :: crack_slab, fd_start, fd_end, bulk
  type(CrackParams) :: params
  type(Potential) :: classicalpot, qmpot
  type(Potential) :: hybrid_pot, forcemix_pot
  type(Dictionary) :: pot_params
  type(mpi_context) :: mpi_glob
  type(Table) :: crack_tips, old_crack_tips

  ! Pointers into Atoms data table
  real(dp), pointer, dimension(:,:) :: load, force
  integer, pointer, dimension(:) :: move_mask, nn, changed_nn, edge_mask, load_mask, md_old_changed_nn, &
       old_nn, hybrid, hybrid_mark

  ! Big arrays
  real(dp), allocatable, dimension(:,:) :: f_fm, dr
  real(dp), pointer :: dr_prop(:,:)

  ! Scalars
  integer :: movie_n, nargs, i, state, steps, iunit, k, md_stanza_idx, random_seed
  logical :: mismatch, movie_exist, periodic_clusters(3), dummy, texist
  real(dp) :: fd_e0, f_dr, integral, energy, last_state_change_time, last_print_time, &
       last_checkpoint_time, last_calc_connect_time, &
       last_md_interval_time, time, temp, crack_pos(2), orig_crack_pos, &
       G, last_update_selection_time, last_stanza_change_time, last_update_crack_tip_time
  character(FIELD_LENGTH) :: stem, movie_name, xmlfilename, suffix, checkfile_name
  character(FIELD_LENGTH) :: state_string

  real(dp) :: time_elapsed, total_time_elapsed


  !** Initialisation Code **

  nargs = cmd_arg_count()

  call initialise(mpi_glob)

  if (mpi_glob%active) then
     ! Same random seed for each process, otherwise velocites will differ
     call system_initialise (common_seed = .true., enable_timing=.true., mpi_all_inoutput=.false.)
     call print('MPI run with '//mpi_glob%n_procs//' processes')
  else
     call system_initialise( enable_timing=.true.)
     call print('Serial run')
  end if

  call system_timer('initialisation')
  call initialise(params)

  if (.not. mpi_glob%active) then
     ! Print usage information if no arguments given
     if (nargs /= 1) then
        call print('Usage: crack <stem>')
        call print('')
        call print('Where stem.xyz and stem.xml are the crack slab XYZ file')
        call print('and parameter file respectively')
        call print('')
        call print('Available parameters and their default values are:')
        call print('')
        call print(params)

        call system_finalise
        stop
     end if
  end if

  ! 1st argument contains stem for input xyz or binary checkfile and parameter files
  if (.not. mpi_glob%active) then
     call get_cmd_arg(1, stem)
  else
     stem = 'crack' ! Some MPIs can't handle command line arguments
  end if

  xmlfilename = trim(stem)//'.xml'

  call print_title('Initialisation')
  call print('Reading parameters from file '//trim(xmlfilename))
  call initialise(xmlfile,xmlfilename,INPUT)
  call read_xml(params,xmlfile)
  call verbosity_push(params%io_verbosity)   ! Set base verbosity
  call print(params)

  if (params%io_mpi_print_all) then
     mainlog%mpi_all_inoutput_flag = .true.
     mainlog%mpi_print_id = .true.
  end if

  call Print('Reading bulk cell from file '//trim(stem)//'_bulk.xyz')
  call read(bulk, trim(stem)//'_bulk.xyz')

  call print ("Initialising classical potential with args " // trim(params%classical_args) &
       // " from file " // trim(xmlfilename))
  call rewind(xmlfile)
  call initialise(classicalpot, params%classical_args, xmlfile, mpi_obj=mpi_glob, bulk_scale=bulk)
  call Print(classicalpot)

  if (.not. params%simulation_classical) then
     call print ("Initialising QM potential with args " // trim(params%qm_args) &
          // " from file " // trim(xmlfilename))
     call rewind(xmlfile)
     call initialise(qmpot, trim(params%qm_args)//' little_clusters='//params%qm_little_clusters, xmlfile, mpi_obj=mpi_glob, bulk_scale=bulk)
     call finalise(xmlfile)
     call Print(qmpot)
  end if

  if (params%qm_force_periodic) then
     periodic_clusters = (/ .false., .false., .true. /)
  else
     periodic_clusters = (/ .false., .false., .not. params%qm_little_clusters/)
  endif

  if (.not. params%simulation_classical) then

     call initialise(pot_params)
     call set_value(pot_params, 'method', trim(params%fit_method))
     call set_value(pot_params, 'buffer_hops', params%qm_buffer_hops)
     call set_value(pot_params, 'transition_hops', params%qm_transition_hops)
     call set_value(pot_params, 'fit_hops', params%fit_hops)
     call set_value(pot_params, 'minimise_mm', params%minim_minimise_mm)
     call set_value(pot_params, 'randomise_buffer', params%qm_randomise_buffer)
     call set_value(pot_params, 'mm_reweight', params%classical_force_reweight)
     call set_value(pot_params, 'minim_mm_method', trim(params%minim_mm_method))
     call set_value(pot_params, 'minim_mm_tol', params%minim_mm_tol)
     call set_value(pot_params, 'minim_mm_eps_guess', params%minim_mm_eps_guess)
     call set_value(pot_params, 'minim_mm_max_steps', params%minim_mm_max_steps)
     call set_value(pot_params, 'minim_mm_linminroutine', trim(params%minim_mm_linminroutine))
     call set_value(pot_params, 'minim_mm_args_str', trim(params%minim_mm_args_str) )
     call set_value(pot_params, 'minim_mm_use_n_minim', params%minim_mm_use_n_minim)
     call set_value(pot_params, 'lotf_spring_hops', params%fit_spring_hops)
     call set_value(pot_params, 'do_rescale_r', params%qm_rescale_r)
     call set_value(pot_params, 'minimise_bulk', params%qm_rescale_r)
     call set_value(pot_params, 'hysteretic_buffer', params%qm_hysteretic_buffer)
     call set_value(pot_params, 'hysteretic_buffer_inner_radius', params%qm_hysteretic_buffer_inner_radius)
     call set_value(pot_params, 'hysteretic_buffer_outer_radius', params%qm_hysteretic_buffer_outer_radius)

     call set_value(pot_params, 'hysteretic_connect', params%qm_hysteretic_connect)
     call set_value(pot_params, 'nneighb_only', (.not. params%qm_hysteretic_connect))
     call set_value(pot_params, 'hysteretic_connect_cluster_radius', params%qm_hysteretic_connect_cluster_radius)
     call set_value(pot_params, 'hysteretic_connect_inner_factor', params%qm_hysteretic_connect_inner_factor)
     call set_value(pot_params, 'hysteretic_connect_outer_factor', params%qm_hysteretic_connect_outer_factor)


     call set_value(pot_params, 'mm_args_str', params%classical_args_str)

     call set_value(pot_params, 'qm_args_str', &
          ' little_clusters='//(params%qm_clusters .and. params%qm_little_clusters)// &
          ' single_cluster='//(params%qm_clusters .and. .not. params%qm_little_clusters)// &
          ' terminate='//params%qm_terminate// &
          ' even_electrons='//params%qm_even_electrons// &
          ' cluster_vacuum='//params%qm_vacuum_size// &
          ' cluster_periodic_x=F cluster_periodic_y=F cluster_periodic_z='//periodic_clusters(3)// &
          ' cluster_calc_connect='//(cutoff(qmpot) /= 0.0_dp)// &
          ' buffer_hops='//params%qm_buffer_hops//&
          ' transition_hops='//params%qm_transition_hops//&
          ' randomise_buffer='//params%qm_randomise_buffer//&
          ' hysteretic_connect='//params%qm_hysteretic_connect//&
          ' nneighb_only='//(.not. params%qm_hysteretic_connect)//&
          ' cluster_nneighb_only='//(.not. params%qm_hysteretic_connect)//&
	  ' '//trim(params%qm_args_str))

     call initialise(hybrid_pot, 'ForceMixing '//write_string(pot_params), &
          pot1=classicalpot, pot2=qmpot, bulk_scale=bulk, mpi_obj=mpi_glob)

     call print_title('Hybrid Potential')
     call print(hybrid_pot)

     call set_value(pot_params, 'method', 'force_mixing')
     if (params%qm_rescale_r) then
        call initialise(forcemix_pot, 'ForceMixing '//write_string(pot_params), &
             pot1=classicalpot, pot2=qmpot, bulk_scale=bulk, mpi_obj=mpi_glob)
        call finalise(bulk)
     else
        call initialise(forcemix_pot, 'ForceMixing '//write_string(pot_params), &
             pot1=classicalpot, pot2=qmpot, mpi_obj=mpi_glob)
     end if
     call finalise(pot_params)

  end if

  ! Look for input file. Check for the following, in order
  ! <stem>_check.nc
  ! <stem>_check.xyz
  ! <stem>.nc
  ! <stem>.xyz

  texist=.false.
  if (params%io_netcdf) then
     inquire(file=trim(stem)//'_check.nc', exist=texist)
     if (texist) then
        call print('Reading atoms from input file '//trim(stem)//'_check.nc')
        call initialise(crackin, trim(stem)//'_check.nc', action=INPUT, mpi=mpi_glob)
        call read(crackin, crack_slab)
        call finalise(crackin)
     endif
  end if

  if (.not. texist) then
     inquire(file=trim(stem)//'_check.xyz', exist=texist)
     if (texist) then
        call print('Reading atoms from input file '//trim(stem)//'_check.xyz')
        call initialise(crackin, trim(stem)//'_check.xyz', action=INPUT, mpi=mpi_glob)
        call read(crackin, crack_slab)
        call finalise(crackin)
     endif
  end if

  if (params%io_netcdf .and. .not. texist) then
     inquire(file=trim(stem)//'.nc', exist=texist)
     if (texist) then
        call print('Reading atoms from input file '//trim(stem)//'.nc')
        call initialise(crackin, trim(stem)//'.nc', action=INPUT, mpi=mpi_glob)
        call read(crackin, crack_slab)
        call finalise(crackin)
     endif
  end if

  if (.not. texist) then
     inquire(file=trim(stem)//'.xyz', exist=texist)
     if (texist) then
        call print('Reading atoms from input file '//trim(stem)//'.xyz')
        call initialise(crackin, trim(stem)//'.xyz', action=INPUT, mpi=mpi_glob)
        call read(crackin, crack_slab)
        call finalise(crackin)
     endif
  end if

  if (.not. texist) call system_abort('No input file found - checked for <stem>_check.nc, <stem>_check.xyz, <stem>.nc and <stem>.xyz with stem="'//trim(stem)//'"')

  call initialise(ds, crack_slab)
  if (any(params%md(1:params%num_md_stanza)%extrapolate_steps /= 1)) call initialise(ds_save, crack_slab)
  call finalise(crack_slab)

  call print('Initialised dynamical system with '//ds%N//' atoms')


  call crack_fix_pointers(ds%atoms, nn, changed_nn, load, move_mask, edge_mask, load_mask, md_old_changed_nn, &
       old_nn, hybrid, hybrid_mark, force)

  ds%atoms%damp_mask = 1
  ds%atoms%thermostat_region = 1
!!$  where (ds%atoms%move_mask == 0)
!!$     ds%atoms%thermostat_region = 0
!!$     ds%atoms%damp_mask = 0
!!$  end where

  ! Set number of degrees of freedom correctly
  ds%Ndof = 3*count(ds%atoms%move_mask == 1)

  ! Reseed random number generator if necessary
  if (get_value(ds%atoms%params, 'random_seed', random_seed)) then
     call system_reseed_rng(random_seed)
  else if (params%simulation_seed /= 0) then
     call system_reseed_rng(params%simulation_seed)
  end if

#ifndef HAVE_NETCDF
  if (params%io_netcdf) &
       call system_abort('io_netcdf = .true. but NetCDF support not compiled in')
#endif 

  if (params%io_netcdf) then
     suffix = '.nc'
  else
     suffix = '.xyz'
  end if

  if (.not. mpi_glob%active .or. (mpi_glob%active .and.mpi_glob%my_proc == 0)) then
     ! Avoid overwriting movie file from previous runs by suffixing a number
     movie_n = 1
     movie_exist = .true.
     do while (movie_exist)
        write (movie_name, '(a,i0,a)') trim(stem)//'_movie_', movie_n, trim(suffix)
        inquire (file=movie_name, exist=movie_exist)
        movie_n = movie_n + 1
     end do

     call print('Setting up movie output file '//movie_name)
     call initialise(movie, movie_name, action=OUTPUT)
     if (params%io_backup) then 
        write (movie_name, '(a,i0,a)') trim(stem)//'_movie_backup_', movie_n, trim(suffix)
        call initialise(movie_backup, movie_name, action=OUTPUT)
     end if
  endif

  call Print('Setting neighbour cutoff to '//(cutoff(classicalpot)+params%md(params%md_stanza)%crust)//' A.')
  call set_cutoff(ds%atoms, cutoff(classicalpot)+params%md(params%md_stanza)%crust)
  call print('Neighbour crust is '//params%md(params%md_stanza)%crust// ' A.')

  call calc_connect(ds%atoms, store_is_min_image=.true.)

  if (params%qm_calc_force_error) allocate(f_fm(3,ds%atoms%N))

  ! Allocate various flags 

  call Print('Setting nneightol to '//params%md(params%md_stanza)%nneigh_tol)
  ds%atoms%nneightol = params%md(params%md_stanza)%nneigh_tol

  call crack_update_connect(ds%atoms, params)

  ! Initialise QM region
  if (.not. params%simulation_classical) then
     if (trim(params%selection_method) /= 'static') then
        call print('Initialising dynamic QM region')

        ! See if we read changed_nn from file
        if (trim(params%selection_method) == 'coordination' .and. all(changed_nn == 0)) &
             call system_abort('No seed atoms found - rerun makecrack')
        call print('count(changed_nn /= 0) = '//count(changed_nn /= 0))
        call crack_update_selection(ds%atoms, params)

     else
        call print('Static QM region')
        ! Load QM mask from file, then fix it for entire simulation

        call print('Loaded hybrid mask from XYZ file')
!!$        call crack_update_selection(ds%atoms, params, embedlist=embedlist, fitlist=fitlist, &
!!$             update_embed=.false., num_directionality=directionN)
     end if
  end if

  ! Print a frame before we start
  if (.not. mpi_glob%active .or. (mpi_glob%active .and.mpi_glob%my_proc == 0)) then
     call crack_print(ds%atoms, movie, params)
     if (params%io_backup) &
          call crack_print(ds%atoms, movie_backup, params)
  end if

  if (.not. params%simulation_classical) then
     if (count(hybrid == 1) == 0) call system_abort('Zero QM atoms selected')
  end if

  call setup_parallel(classicalpot, ds%atoms, args_str=trim(params%classical_args_str)//" energy force")

  call crack_fix_pointers(ds%atoms, nn, changed_nn, load, move_mask, edge_mask, load_mask, md_old_changed_nn, &
       old_nn, hybrid, hybrid_mark, force)

  if (all(abs(load) < 1.0e-7_dp)) then
     call print_title('Applying Initial Load')
     call crack_calc_load_field(ds%atoms, params, classicalpot, params%crack_loading, & 
          .true., mpi_glob)

     call crack_fix_pointers(ds%atoms, nn, changed_nn, load, move_mask, edge_mask, load_mask, md_old_changed_nn, &
          old_nn, hybrid, hybrid_mark, force)
  end if

  if (params%simulation_force_initial_load_step) then
     if (.not. has_property(ds%atoms, 'load')) &
          call system_abort('simulation_force_initial_load_step is true but crack slab has no load field - set crack_apply_initial_load = T to regenerate load')
     call print_title('Force_load_step is true,  applying load')
     call crack_apply_load_increment(ds%atoms, params%crack_G_increment)
  end if

  call system_timer('initialisation', time_elapsed=time_elapsed)
  total_time_elapsed = time_elapsed

  !** End of initialisation **

  !****************************************************************
  !*                                                              *
  !*  MOLECULAR DYNAMICS                                          *
  !*                                                              *
  !*                                                              *
  !****************************************************************    
  if (trim(params%simulation_task) == 'md' .or. trim(params%simulation_task) == 'damped_md') then

     md_stanza: do md_stanza_idx = 1, params%num_md_stanza
        params%md_stanza = md_stanza_idx

        call system_timer('md_initialisation')

        call print_title('Molecular Dynamics (md_stanza='//params%md_stanza//')')

        if (.not. get_value(ds%atoms%params, 'Temp', temp)) temp = 2.0_dp*params%md(params%md_stanza)%sim_temp

        ! If velocities are zero, randomise them at correct temperature
        if (maxval(abs(ds%atoms%velo)) < 1e-5_dp) then
           call rescale_velo(ds, temp)
           call zero_momentum(ds)
           ds%cur_temp = temperature(ds, instantaneous=.true.)
           ds%avg_temp = temperature(ds)
        end if

        if (.not. get_value(ds%atoms%params, 'Time', time))  time = 0.0_dp
        ds%t = time

        if (.not. get_value(ds%atoms%params, 'LastStateChangeTime', last_state_change_time)) &
             last_state_change_time = ds%t

        if (.not. get_value(ds%atoms%params, 'LastMDIntervalTime', last_md_interval_time)) &
             last_md_interval_time = ds%t

        if (.not. get_value(ds%atoms%params, 'LastPrintTime', last_print_time)) &
             last_print_time = ds%t

        if (.not. get_value(ds%atoms%params, 'LastCheckpointTime', last_checkpoint_time)) &
             last_checkpoint_time = ds%t

        if (.not. get_value(ds%atoms%params, 'LastCalcConnectTime', last_calc_connect_time)) &
             last_calc_connect_time = ds%t

        last_stanza_change_time = ds%t
        last_update_selection_time = ds%t
        last_update_crack_tip_time = ds%t

        ! Special cases for first time
        if (all(md_old_changed_nn == 0)) md_old_changed_nn = changed_nn
        if (all(old_nn == 0)) old_nn = nn

        ds%avg_time = params%md(params%md_stanza)%avg_time

        if (trim(params%simulation_task) == 'md') then
           if (.not. get_value(ds%atoms%params, "State", state_string)) then
              if (params%md(params%md_stanza)%smooth_loading_rate .fne. 0.0_dp) then 
                 state_string = 'MD_LOADING'
              else
                 state_string = 'THERMALISE'
              end if
           end if
        else 
           state_string = 'DAMPED_MD'
        end if

        ! Allow initial state to be overridden in XML file
        if (trim(params%simulation_initial_state) /= '') state_string = params%simulation_initial_state

        if (trim(state_string) == "MD" .and. (params%md(params%md_stanza)%smooth_loading_rate .fne. 0.0_dp)) state_string = 'MD_LOADING'

        ! Remove any existing thermostats prior to re-adding them
        call finalise(ds%thermostat)
        
        ! Special thermostat for damping
        allocate(ds%thermostat(0:0)); call initialise(ds%thermostat(0),NONE,0.0_dp)

        if (state_string(1:10) == 'THERMALISE') then
           state = STATE_THERMALISE
           call disable_damping(ds)
           call ds_add_thermostat(ds, LANGEVIN, params%md(params%md_stanza)%sim_temp, tau=params%md(params%md_stanza)%thermalise_tau)

        else if (state_string(1:10) == 'MD') then
           state = STATE_MD
           call disable_damping(ds)
           call ds_add_thermostat(ds, LANGEVIN, params%md(params%md_stanza)%sim_temp, tau=params%md(params%md_stanza)%tau)

        else if (state_string(1:10) == 'MD_LOADING') then
           state = STATE_MD_LOADING
           call disable_damping(ds)
           call ds_add_thermostat(ds, LANGEVIN, params%md(params%md_stanza)%sim_temp, tau=params%md(params%md_stanza)%tau)
           call crack_find_tip(ds%atoms, params, old_crack_tips)
           crack_tips = old_crack_tips

        else if (state_string(1:11) == 'MD_CRACKING') then
           state = STATE_MD_CRACKING
           call disable_damping(ds)
           call ds_add_thermostat(ds, LANGEVIN, params%md(params%md_stanza)%sim_temp, tau=params%md(params%md_stanza)%tau)
           call crack_find_tip(ds%atoms, params, old_crack_tips)
           crack_tips = old_crack_tips

        else if (state_string(1:9) == 'DAMPED_MD') then
           state = STATE_DAMPED_MD
           call enable_damping(ds, params%md(params%md_stanza)%damping_time)

        else if (state_string(1:14) == 'MICROCANONICAL') then
           state = STATE_MICROCANONICAL
           call disable_damping(ds)

        else if (state_string(1:11) == 'MD_CONSTANT') then
           state = STATE_MD_CONSTANT
           call disable_damping(ds)
           call ds_add_thermostat(ds, LANGEVIN, params%md(params%md_stanza)%sim_temp, tau=params%md(params%md_stanza)%tau)

        else if (state_string(1:10) == 'MD_NVE') then
           state = STATE_MD_NVE
           call disable_damping(ds)

        else  
           call system_abort("Don't know how to resume in molecular dynamics state "//trim(state_string))
        end if

        call print('Thermostats')
        call print(ds%thermostat)

        call print('Starting in state '//STATE_NAMES(state))

        ! Bootstrap the adjustable potential if we're doing predictor/corrector dynamics
        if (params%md(params%md_stanza)%extrapolate_steps /= 1 .and. .not. params%simulation_classical) then
           call calc(hybrid_pot, ds%atoms, args_str="force=force")
        end if

        call system_timer('md_initialisation', time_elapsed=time_elapsed)
        total_time_elapsed = total_time_elapsed + time_elapsed

        !****************************************************************
        !*  Main MD Loop                                                *
        !*                                                              *
        !****************************************************************    
        do
           call system_timer('step')
           call crack_fix_pointers(ds%atoms, nn, changed_nn, load, move_mask, edge_mask, load_mask, md_old_changed_nn, &
                old_nn, hybrid, hybrid_mark, force)


           select case(state)
           case(STATE_THERMALISE)
              if (ds%t - last_state_change_time >= params%md(params%md_stanza)%thermalise_wait_time .and. &
                   abs(ds%avg_temp - temperature(ds))/temperature(ds) < params%md(params%md_stanza)%thermalise_wait_factor/sqrt(real(ds%atoms%N,dp))) then
                 ! Change to state MD
                 call print('STATE changing THERMALISE -> MD')
                 state = STATE_MD
                 last_state_change_time = ds%t
                 call disable_damping(ds)
                 call initialise(ds%thermostat(1), LANGEVIN, params%md(params%md_stanza)%sim_temp, &
                      gamma=1.0_dp/params%md(params%md_stanza)%tau)
                 md_old_changed_nn = changed_nn
              end if

           case(STATE_MD)
              if ((ds%t - last_state_change_time >= params%md(params%md_stanza)%wait_time) .and. &
                   (ds%t - last_md_interval_time  >= params%md(params%md_stanza)%interval_time)) then

                 mismatch = .false.
                 do i = 1,ds%atoms%N
                    if ((changed_nn(i) == 0 .and. md_old_changed_nn(i) /= 0) .or. &
                         (changed_nn(i) /= 0 .and. md_old_changed_nn(i) == 0)) then
                       mismatch = .true.
                       exit
                    end if
                 end do

                 if (.not. mismatch) then 
                    ! changed_nn hasn't changed for a while so we can increase strain
                    ! Rescale and change to state THERMALISE
                    call print('STATE changing MD -> THERMALISE')
                    state = STATE_THERMALISE
                    last_state_change_time = ds%t

                    call disable_damping(ds)
                    call initialise(ds%thermostat(1), LANGEVIN, params%md(params%md_stanza)%sim_temp, &
                         gamma=1.0_dp/params%md(params%md_stanza)%thermalise_tau)
                 end if

                 ! Apply loading field
                 if (has_property(ds%atoms, 'load')) then
                    call print_title('Applying load')
                    call crack_apply_load_increment(ds%atoms, params%crack_G_increment)
                 else
                    call print('No load field found - not increasing load.')
                 end if

                 md_old_changed_nn = changed_nn
                 last_md_interval_time = ds%t
              end if

           case(STATE_DAMPED_MD)

           case(STATE_MICROCANONICAL)

           case(STATE_MD_LOADING)
              ! If tip has moved by more than smooth_loading_tip_move_tol then
              ! turn off loading. 
              call crack_find_tip(ds%atoms, params, crack_tips)

              if (crack_tips%N /= old_crack_tips%N) &
                   call system_abort('State MD_LOADING: number of crack tips changed from '//old_crack_tips%N//' to '//crack_tips%N)

              if (crack_tips%real(1,crack_tips%N) - old_crack_tips%real(1,crack_tips%N) > params%md(params%md_stanza)%smooth_loading_tip_move_tol) then
                 call print_title('Crack Moving')
                 call print('STATE changing MD_LOADING -> MD_CRACKING')
                 state = STATE_MD_CRACKING
                 last_state_change_time = ds%t
                 old_crack_tips = crack_tips
              else
                 call print('STATE: crack is not moving (crack_pos='//crack_tips%real(1,crack_tips%N)//')')
              end if

           case(STATE_MD_CRACKING)
              ! Monitor tip and if it doesn't move by more than smooth_loading_tip_move_tol in
              ! time smooth_loading_arrest_time then switch back to loading
              if (ds%t - last_state_change_time >= params%md(params%md_stanza)%smooth_loading_arrest_time) then

                 call crack_find_tip(ds%atoms, params, crack_tips)

                 if (crack_tips%N == 0) then
                    call print_title('Cracked Through')
                    exit
                 end if

                 if (crack_tips%real(1,crack_tips%N) - old_crack_tips%real(1,crack_tips%N) < params%md(params%md_stanza)%smooth_loading_tip_move_tol) then
                    call print_title('Crack Arrested')
                    call crack_calc_load_field(ds%atoms, params, classicalpot, params%crack_loading, &
                         .false., mpi_glob)
                    call print('STATE changing MD_CRACKING -> MD_LOADING')
                    state = STATE_MD_LOADING
                 else
                    call print('STATE: crack is moving, crack_tips=')
                    call print(crack_tips)
                 end if

                 last_state_change_time = ds%t
                 old_crack_tips = crack_tips
              end if

           case(STATE_MD_CONSTANT)
              
           case(STATE_MD_NVE)

           case default
              call system_abort('Unknown molecular dynamics state!')
           end select

           ! Are we doing predictor/corrector dynamics?
           if (params%md(params%md_stanza)%extrapolate_steps /= 1) then

              !****************************************************************
              !*  Quantum Selection                                           *
              !*                                                              *
              !****************************************************************    
              call system_timer('selection')
              call print_title('Quantum Selection')
              if (trim(params%selection_method) /= 'static') call crack_update_selection(ds%atoms, params)
              call system_timer('selection')

              !****************************************************************
              !*  Extrapolation                                               *
              !*                                                              *
              !****************************************************************    
              if (.not. params%simulation_classical) then
                 call print_title('Extrapolation')
                 call system_timer('extrapolation')
                 call ds_save_state(ds_save, ds)
              else
                 call system_timer('md_time')
              end if

              do i = 1, params%md(params%md_stanza)%extrapolate_steps

                 if (params%simulation_classical) then
                    call calc(classicalpot, ds%atoms, energy=energy, args_str=trim(params%classical_args_str)//' energy=energy force=force')
                 else
                    if (i== 1) then
                       call calc(hybrid_pot, ds%atoms, args_str="force=force lotf_do_qm=F lotf_do_init=T lotf_do_map=T")
                    else
                       call calc(hybrid_pot, ds%atoms, args_str="force=force lotf_do_qm=F lotf_do_init=F")
                    end if
                    if (params%qm_calc_force_error) call calc(forcemix_pot, ds%atoms, force=f_fm)

                    if (params%hack_qm_zero_z_force) then
                       ! Zero z forces in embed region
                       force(3,find(hybrid == 1)) = 0.0_dp 
                       if (params%qm_calc_force_error) f_fm(3, find(hybrid == 1)) = 0.0_dp
                    end if
                 end if

                 ! advance the dynamics
                 call system_timer('advance_verlet')
                 call advance_verlet(ds, params%md(params%md_stanza)%time_step, force, do_calc_dists=(state /= STATE_MD_LOADING))
                 call system_timer('advance_verlet')
                 if (params%simulation_classical) then
                    call ds_print_status(ds, 'E', epot=energy)
                 else
                    call ds_print_status(ds, 'E')
                 end if
                 if (params%qm_calc_force_error) call print('E err '//ds%t//' '//rms_diff(force, f_fm)//' '//maxval(abs(f_fm-force)))

                 if (state == STATE_MD_LOADING) then
                    ! increment the load
                    call system_timer('apply_load_increment')
                    if (has_property(ds%atoms, 'load')) then
                       call system_timer('load_increment')
                       call crack_apply_load_increment(ds%atoms, params%md(params%md_stanza)%smooth_loading_rate*params%md(params%md_stanza)%time_step)
                       call system_timer('load_increment')
                       if (.not. get_value(ds%atoms%params, 'G', G)) call system_abort('No G in ds%atoms%params')
                    else
                       call print('No load field found - not increasing load.')
                    end if
                    call system_timer('apply_load_increment')
                 end if

              end do

              if (.not. params%simulation_classical) then
                 call system_timer('extrapolation')
              else
                 call system_timer('md_time')
              endif

              if (.not. params%simulation_classical) then

                 !****************************************************************
                 !*  QM Force Computation                                        *
                 !*  and optimisation of Adjustable Potential                    *
                 !*                                                              *
                 !****************************************************************    

                 call print_title('Computation of forces')
                 call system_timer('force computation')
                 call calc(hybrid_pot, ds%atoms, args_str="force=force lotf_do_qm=T lotf_do_init=F lotf_do_fit=T")
                 call system_timer('force computation')


                 !****************************************************************
                 !*  Interpolation                                               *
                 !*                                                              *
                 !****************************************************************    
                 call print_title('Interpolation')
                 call system_timer('interpolation')

                 ! revert to the saved positions etc.
                 call ds_restore_state(ds, ds_save)
                 call crack_fix_pointers(ds%atoms, nn, changed_nn, load, move_mask, edge_mask, load_mask, md_old_changed_nn, &
                      old_nn, hybrid, hybrid_mark, force)

                 do i = 1, params%md(params%md_stanza)%extrapolate_steps

                    call calc(hybrid_pot, ds%atoms, args_str="force=force lotf_do_qm=F lotf_do_init=F lotf_do_interp=T lotf_interp="&
                         //(real(i-1,dp)/real(params%md(params%md_stanza)%extrapolate_steps,dp)))

                    if (params%qm_calc_force_error) call calc(forcemix_pot, ds%atoms, force=f_fm)

                    if (params%hack_qm_zero_z_force) then
                       ! Zero z forces in embed region
                       force(3,find(hybrid == 1)) = 0.0_dp 
                       if (params%qm_calc_force_error) f_fm(3, find(hybrid == 1)) = 0.0_dp
                    end if

                    ! advance the dynamics
                    call advance_verlet(ds, params%md(params%md_stanza)%time_step, force, do_calc_dists=(state /= STATE_MD_LOADING))
                    call ds_print_status(ds, 'I')
                    if (params%qm_calc_force_error) call print('I err '//ds%t//' '//rms_diff(force, f_fm)//' '//maxval(abs(f_fm-force)))

                    if (trim(params%simulation_task) == 'damped_md') &
                         call print('Damped MD: normsq(force) = '//normsq(reshape(force,(/3*ds%N/)))//&
                         ' max(abs(force)) = '//maxval(abs(force)))

                    if (state == STATE_MD_LOADING) then
                       ! increment the load
                       if (has_property(ds%atoms, 'load')) then
                          call crack_apply_load_increment(ds%atoms, params%md(params%md_stanza)%smooth_loading_rate*params%md(params%md_stanza)%time_step)
                          if (.not. get_value(ds%atoms%params, 'G', G)) call system_abort('No G in ds%atoms%params')
                       else
                          call print('No load field found - not increasing load.')
                       end if
                    end if

                 end do
                 call system_timer('interpolation')

              end if ! .not. params%simulation_classical

           else ! params%md(params%md_stanza)%extrapolate_steps /= 1

              !****************************************************************
              !*  Non-Predictor/Corrector Dynamics                            *
              !*                                                              *
              !****************************************************************    

              if (ds%t - last_update_selection_time >= params%selection_update_interval) then
                 last_update_selection_time = ds%t
                 call print_title('Quantum Selection')
                 call system_timer('selection')
                 if (trim(params%selection_method) /= 'static') call crack_update_selection(ds%atoms, params)
                 call system_timer('selection')
              end if

              if (ds%t - last_update_crack_tip_time >= params%md(params%md_stanza)%crack_find_tip_interval) then
                 call system_timer('crack_find_tip')
                 last_update_crack_tip_time = ds%t
                 mainlog%prefix = 'CRACK_TIP'
                 call crack_find_tip(ds%atoms, params, crack_tips)
                 call print(crack_tips)
                 mainlog%prefix = ''
                 call system_timer('crack_find_tip')
              end if

              call print_title('Force Computation')
              call system_timer('force computation/optimisation')
              if (params%simulation_classical) then
                 call calc(classicalpot, ds%atoms, energy=energy, args_str=trim(params%classical_args_str)//' energy=energy force=force')
              else
                 call calc(hybrid_pot, ds%atoms, args_str="force=force")
              end if
              call system_timer('force computation/optimisation')

              if (params%hack_qm_zero_z_force) then
                 ! Zero z forces in embed region
                 force(3,find(hybrid == 1)) = 0.0_dp 
              end if

              call print_title('Advance Verlet')
              call system_timer('advance_verlet')
              call advance_verlet(ds, params%md(params%md_stanza)%time_step, force, do_calc_dists=(state /= STATE_MD_LOADING))
              call system_timer('advance_verlet')
              if (params%simulation_classical) then
                 call ds_print_status(ds, 'D', epot=energy)
              else
                 call ds_print_status(ds, 'D')
              end if

              if (trim(params%simulation_task) == 'damped_md') &
                   call print('Damped MD: normsq(force) = '//normsq(reshape(force,(/3*ds%N/)))//&
                   ' max(abs(force)) = '//maxval(abs(force)))

              call system_timer('load_increment')
              if (state == STATE_MD_LOADING) then
                 ! increment the load
                 if (has_property(ds%atoms, 'load')) then
                    call crack_apply_load_increment(ds%atoms, params%md(params%md_stanza)%smooth_loading_rate*params%md(params%md_stanza)%time_step)
                    if (.not. get_value(ds%atoms%params, 'G', G)) call system_abort('No G in ds%atoms%params')
                 else
                    call print('No load field found - not increasing load.')
                 end if
              end if
              call system_timer('load_increment')

           end if ! params%extrapolate_steps /= 1


           ! Do I/O only on master node
           if (.not. mpi_glob%active .or. (mpi_glob%active .and.mpi_glob%my_proc == 0)) then

              call set_value(ds%atoms%params, 'Time', ds%t)
              call set_value(ds%atoms%params, 'Temp', temperature(ds))
              call set_value(ds%atoms%params, 'LastStateChangeTime', last_state_change_time)
              call set_value(ds%atoms%params, 'LastMDIntervalTime', last_md_interval_time)
              call set_value(ds%atoms%params, 'LastCalcConnectTime', last_calc_connect_time)
              call set_value(ds%atoms%params, 'State', STATE_NAMES(state))

              ! Print movie
              if (ds%t - last_print_time >=  params%io_print_interval) then

                 last_print_time = ds%t
                 call set_value(ds%atoms%params, 'LastPrintTime', last_print_time)

                 if (params%io_backup) then
                    k=k+1           
                    if (mod(k,2).eq.0) then
                       call crack_print(ds%atoms, movie, params)
                       call print('writing .nc file '//trim(stem)//'.nc')
                    else
                       call crack_print(ds%atoms, movie_backup, params)
                       call print('writing .nc file '//trim(stem)//'_backup.nc')
                    endif
                 else
                    call crack_print(ds%atoms, movie, params)
                 end if
              end if

              ! Write checkpoint file
              if (ds%t - last_checkpoint_time >=  params%io_checkpoint_interval) then
                 last_checkpoint_time = ds%t
                 call set_value(ds%atoms%params, 'LastCheckpointTime', last_checkpoint_time)
                 call set_value(ds%atoms%params, 'random_seed', system_get_random_seed())

                 checkfile_name = trim(params%io_checkpoint_path)//trim(stem)//'_check'//suffix
                 inquire (file=checkfile_name,exist=texist)
                 if (texist) then
                    call system_command('mv '//trim(checkfile_name)//' '//trim(checkfile_name)//'.backup')
                 end if
                 call write(ds%atoms, checkfile_name)
              endif
           end if

           ! Recalculate connectivity and nearest neighbour tables
           if (ds%t - last_calc_connect_time >= params%md(params%md_stanza)%calc_connect_interval) then
              last_calc_connect_time = ds%t
              call crack_update_connect(ds%atoms, params)
           end if

           ! check if time exceeds time for this MD stanza
           if (params%md(params%md_stanza)%stanza_time >= 0.0_dp) then
              if (ds%t - last_stanza_change_time >= params%md(params%md_stanza)%stanza_time) cycle md_stanza
           end if

           ! Exit cleanly if file 'stop_run' exists
           inquire (file='stop_run',exist=texist)
           if (texist) then
              iunit = pick_up_unit()
              open (iunit,file='stop_run',status='old')
              close (iunit,status='delete')
              exit md_stanza
           endif

           ! exit cleanly if we exceeded the max run time
           call system_timer('step', time_elapsed=time_elapsed)
           total_time_elapsed = total_time_elapsed + time_elapsed
           if (params%md(params%md_stanza)%max_runtime >= 0.0_dp) then
              call print("elapsed time="//total_time_elapsed//" max_runtime="//params%md(params%md_stanza)%max_runtime)
              if (total_time_elapsed >= params%md(params%md_stanza)%max_runtime) then
                 call print("Exceeded max_runtime, exiting cleanly", PRINT_ALWAYS)
                 exit md_stanza
              endif
           endif

           call flush(mainlog%unit)

        end do

     end do md_stanza


     !****************************************************************
     !*                                                              *
     !*  FORCE INTEGRATION                                           *
     !*                                                              *
     !*                                                              *
     !****************************************************************    
  else if (trim(params%simulation_task) == 'force_integration') then

!!$     params%io_print_properties = trim(params%io_print_properties) // ":dr:forces"

     call print_title('Force Integration')

     fd_start = ds%atoms
     call read(fd_end, params%force_integration_end_file)

     allocate (dr(3,ds%atoms%N))
     dr = (fd_end%pos - fd_start%pos)/real(params%force_integration_n_steps,dp)

     call add_property(ds%atoms, 'dr', 0.0_dp, 3)
     if (.not. assign_pointer(ds%atoms, 'dr', dr_prop)) &
          call system_abort("failed to add dr property to ds%atoms in force_integration task")
     dr_prop = (fd_end%pos - fd_start%pos)

     integral = 0.0_dp

     call print('Force.dr integration')
     write (line, '(a15,a15,a15,a15)') 'Step', 'Energy', 'F.dr', 'Integral(F.dr)'
     call print(line)

     do i=0,params%force_integration_n_steps

        ds%atoms%pos = fd_start%pos + dr*real(i,dp)
        call calc_connect(ds%atoms, store_is_min_image=.true.)

        if (params%simulation_classical) then
           call calc(classicalpot, ds%atoms, energy=energy, args_str=trim(params%classical_args_str)//' force=force')
           if (i == 0) fd_e0 = energy
        else
           call calc(hybrid_pot, ds%atoms, args_str='force=force')
        end if

        f_dr = force .dot. dr

        ! Simpson's rule
        if (i == 0 .or. i == params%force_integration_n_steps) then
           integral = integral + f_dr/3.0_dp
        else if (mod(i,2) == 0) then
           integral = integral + 2.0_dp/3.0_dp*f_dr
        else
           integral = integral + 4.0_dp/3.0_dp*f_dr
        end if

        if (params%simulation_classical) then
           write (line, '(i15,e15.4,f15.4,f15.4)') i, energy, f_dr, integral
        else
           write (line, '(i15,a15,f15.4,f15.4)') i, '----', f_dr, integral
        end if
        call print(line)
        call crack_print(ds%atoms, movie, params)
     end do

     write (line, '(a,f15.8,a)') 'Energy difference = ', fd_e0 - energy, ' eV'
     call print(line)
     write (line, '(a,f15.8,a)') 'Integrated F.dr   = ', integral, ' eV'
     call print(line)

     deallocate(dr)
     call finalise(fd_start, fd_end)


     !****************************************************************
     !*                                                              *
     !*  GEOMETRY OPTIMISATION                                       *
     !*                                                              *
     !*                                                              *
     !****************************************************************    
  else if (trim(params%simulation_task) == 'minim') then

     call print_title('Geometry Optimisation')

     call print('Starting geometry optimisation...')

     if (params%simulation_classical) then
        steps = minim(classicalpot, ds%atoms, method=params%minim_method, convergence_tol=params%minim_tol, &
             max_steps=params%minim_max_steps, linminroutine=params%minim_linminroutine, &
             do_pos=.true., do_lat=.false., do_print=.true., use_fire=trim(params%minim_method)=='fire', &
             print_cinoutput=movie, args_str=params%classical_args_str, eps_guess=params%minim_eps_guess, hook_print_interval=params%minim_print_output)
     else
        steps = minim(hybrid_pot, ds%atoms, method=params%minim_method, convergence_tol=params%minim_tol, &
             max_steps=params%minim_max_steps, linminroutine=params%minim_linminroutine, &
             do_pos=.true., do_lat=.false., do_print=.true., use_fire=trim(params%minim_method)=='fire', &
             print_cinoutput=movie, &
             eps_guess=params%minim_eps_guess, hook_print_interval=params%minim_print_output)
     end if

     if (.not. mpi_glob%active .or. (mpi_glob%active .and.mpi_glob%my_proc == 0)) then
        call crack_update_connect(ds%atoms, params)
        if (params%simulation_classical) then
           call crack_find_tip(ds%atoms, params, crack_tips)
        else
           call crack_update_selection(ds%atoms, params)
        end if
        call crack_print(ds%atoms, movie, params)
     end if

     !****************************************************************
     !*                                                              *
     !*  QUASI-STATIC LOADING                                        *
     !*                                                              *
     !*                                                              *
     !****************************************************************    
  else if (trim(params%simulation_task) == 'quasi_static') then

     call print_title('Quasi Static Loading')

     if (.not. has_property(ds%atoms, 'load')) then
        call print('No load field found. Regenerating load.')
        call crack_calc_load_field(ds%atoms, params, classicalpot, params%crack_loading, &
             .true., mpi_glob)
     end if

     call crack_update_connect(ds%atoms, params)

     dummy = get_value(ds%atoms%params, 'CrackPosx', orig_crack_pos)
     crack_pos(1) = orig_crack_pos

     do while (abs(crack_pos(1) - orig_crack_pos) < params%quasi_static_tip_move_tol)

        if (params%simulation_classical) then
           steps = minim(classicalpot, ds%atoms, method=params%minim_method, convergence_tol=params%minim_tol, &
                max_steps=params%minim_max_steps, linminroutine=params%minim_linminroutine, &
                do_pos=.true., do_lat=.false., do_print=.true., &
                print_cinoutput=movie, &
                args_str=params%classical_args_str, eps_guess=params%minim_eps_guess,use_fire=trim(params%minim_method)=='fire', hook_print_interval=params%minim_print_output)
        else
           steps = minim(hybrid_pot, ds%atoms, method=params%minim_method, convergence_tol=params%minim_tol, &
                max_steps=params%minim_max_steps, linminroutine=params%minim_linminroutine, &
                do_pos=.true., do_lat=.false., do_print=.true., &
                print_cinoutput=movie, eps_guess=params%minim_eps_guess, &
                use_fire=trim(params%minim_method)=='fire', hook_print_interval=params%minim_print_output)
        end if

        call crack_update_connect(ds%atoms, params)
        if (params%simulation_classical) then
           call crack_find_tip(ds%atoms, params, crack_tips)
        else
           call crack_update_selection(ds%atoms, params)
           dummy = get_value(ds%atoms%params, 'CrackPosx', crack_pos(1))
        end if

        ! Apply loading field
        call print_title('Applying load')
        call crack_apply_load_increment(ds%atoms, params%crack_G_increment)

        call crack_print(ds%atoms, movie, params)
     end do

  else

     call system_abort('Unknown task: '//trim(params%simulation_task))

  end if ! switch on simulation_task


  !** Finalisation Code **

  call finalise(movie)
  call finalise(ds)

  call finalise(hybrid_pot)
  call finalise(forcemix_pot)
  call finalise(classicalpot)
  call finalise(qmpot)

  if (allocated(f_fm)) deallocate(f_fm)

  call system_finalise()

end program crack
