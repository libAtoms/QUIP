  !*************************************************************************
  !*
  !*  MetaPotential_FM header
  !*
  !*************************************************************************

  public :: MetaPotential_FM
  type MetaPotential_FM

     type(Potential), pointer :: mmpot => null() 
     type(Potential), pointer :: qmpot => null() 
     type(MPI_context) :: mpi

     character(1024) :: init_args_str

     logical  :: minimise_mm      !% Should classical degrees of freedom be minimised in each calc?
     logical  :: calc_weights     !% Should weights be recalculated on each call to 'calc()'
     character(FIELD_LENGTH) :: method  !% What fit method to use. Options are:
     !% \begin{itemize}
     !%  \item 'lotf_adj_pot_svd' --- LOTF using SVD to optimised the Adj Pot
     !%  \item 'lotf_adj_pot_minim' --- LOTF using conjugate gradients to optimise the Adj Pot
     !%  \item 'conserve_momentum' --- divide the total force on QM region over the fit atoms to conserve momentum 
     !%  \item 'force_mixing' --- force mixing with details depending on values of
     !%      'buffer_hops', 'transtion_hops' and 'weight_interpolation'
     !%  \item 'force_mixing_abrupt' --- simply use QM forces on QM atoms and MM forces on MM atoms 
     !%      (shorthand for 'method=force_mixing buffer_hops=0 transition_hops=0')
     !%  \item 'force_mixing_smooth' --- use QM forces in QM region, MM forces in MM region and 
     !%    linearly interpolate in buffer region  (shorthand for 'method=force_mixing weight_interpolation=hop_ramp')
     !%  \item 'force_mixing_super_smooth' --- as above, but weight forces on each atom by distance from 
     !%    centre of mass of core region (shorthand for 'method=force_mixing weight_interpolation=distance_ramp')
     !% \end{itemize}
     !% Default method is 'conserve_momentum'.
     real(dp) :: mm_reweight      !% Factor by which to reweight classical forces in embed zone
     character(FIELD_LENGTH) :: conserve_momentum_weight_method   !% Weight method to use with 'method=conserve_momentum'. Should be one of
                                                                  !% 'uniform' (default), 'mass', 'mass^2' or 'user',
                                                                  !% with the last referring to a 'conserve_momentum_weight'
                                                                  !% property in the Atoms object.
     character(FIELD_LENGTH) :: mm_args_str !% Args string to be passed to 'calc' method of 'mmpot'
     character(FIELD_LENGTH) :: qm_args_str !% Args string to be passed to 'calc' method of 'qmpot'
     integer :: buffer_hops !% Number of bond hops used for buffer region. Applies to 'force_mixing' methods only.
!     logical :: hysteretic_buffer !% Turn on hysteretic selection algorithm for buffer region
     real(dp) :: buffer_inner_radius,  buffer_outer_radius !% Turn on hysteretic selection algorithm for buffer region
     logical :: hysteretic_connect !% Turn on hysteretic selection algorithm for connection used for buffer region
     real(dp) :: hysteretic_connect_cluster_radius, hysteretic_connect_inner_factor,  hysteretic_connect_outer_factor
     logical :: construct_buffer_use_only_heavy_atoms !% Whether to consider only the heavy atoms radius when constructing the hysteretic buffer

     integer :: transition_hops !% Number of bond hops used for transition region. Applies to 'force_mixing' methods only.
     logical :: use_buffer_for_fitting !% Whether to generate the fit region or just use the buffer as the fit region. Only for method=conserve_momentum
     integer :: fit_hops !% Number of bond hops used for fit region. Applies to 'conserve_momentum' and 'lotf_*' methods only.
     logical :: add_cut_H_in_fitlist !% Whether to extend the fit region where a cut hydrogen is cut after the fitlist selection.
                                     !% This will ensure to only include whole water molecules in the fitlist.
     logical :: randomise_buffer !% If true, then positions of outer layer of buffer atoms will be randomised slightly. Default false.
     character(FIELD_LENGTH) :: weight_interpolation !% Weigh interpolation method for transition region, for 'force_mixing' methods.
       !% Should be one of 'hop_ramp' or 'distance_ramp'; default is 'hop_ramp'.
     logical :: nneighb_only !% Should bond hopping be done to nearest neighbours only when building transition, 
                             !% buffer and fit regions. Default true. Applies to all methods.
     logical :: save_forces !% If true, save MM, QM and total forces as properties in the Atoms object (default true)

     integer :: lotf_spring_hops  !% Maximum lengths of springs for LOTF 'adj_pot_svd' and 'adj_pot_minim' methods (default is 2).
     character(FIELD_LENGTH) :: lotf_interp_order !% Interpolation order: should be one of 'linear', 'quadratic', or 'cubic'. Default is 'linear'.
     logical :: lotf_interp_space !% Do spatial rather than temporal interpolation of adj pot parameters. Default is false.
     logical :: lotf_nneighb_only !% If true (which is the default), uses nearest neigbour hopping to determine fit atoms
     logical :: do_rescale_r !% Should positions in QM region be rescaled
     real(dp) :: r_scale_pot1 !% Rescale positions in QM region by this factor
     
     character(FIELD_LENGTH) :: minim_mm_method 
     real(dp)      :: minim_mm_tol, minim_mm_eps_guess
     integer       :: minim_mm_max_steps
     character(FIELD_LENGTH) :: minim_mm_linminroutine 
     logical       :: minim_mm_do_pos, minim_mm_do_lat
     logical       :: minim_mm_do_print, minim_mm_use_n_minim

     character(FIELD_LENGTH) :: minim_mm_args_str

     type(MetaPotential) :: relax_metapot
     type(Inoutput), pointer :: minim_inoutput_movie
     type(CInoutput), pointer :: minim_cinoutput_movie

     type(Table) :: embedlist, fitlist

  end type MetaPotential_FM

  interface Initialise
     module procedure MetaPotential_FM_initialise
  end interface

  interface Finalise
     module procedure MetaPotential_FM_Finalise
  end interface

  interface Print
     module procedure MetaPotential_FM_Print
  end interface

  interface Cutoff
     module procedure MetaPotential_FM_Cutoff
  end interface

  interface Calc
     module procedure MetaPotential_FM_Calc
  end interface

