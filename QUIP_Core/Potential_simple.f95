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

!X
!X Potential_simple module 
!X
!% The Potential_simple module handles the first-level selection of the
!% desidered force field (tight-binding \texttt{TB_type}, empirical potential \texttt{IP_type} or 
!% hybrid description \texttt{FilePot_type}).
!% It contains the interfaces \texttt{Initialise}, \texttt{Finalise}, \texttt{cutoff}, 
!% \texttt{Calc}, \texttt{Print}, which have the only role to 
!% re-addressing the calls to the corresponding
!% modules. 
!% A Potential object simply contains a pointer to the desired force field type
!% (only one of the three can be selected). 
!% It is initialised with
!%>    call Initialise(pot,arg_str,io_obj,[mpi_obj])
!% where, \texttt{arg_str} is a string defining the potential type and
!% possible some additional options.
!% For interatomic potentials \texttt{arg_str} will start with 'IP'; there are several different IPs defined:
!% see documentation for the 'IP_module' module. For tight binding, 'arg_str' should start with 'TB';
!% see docuementation for 'TB_module' for more details.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#include "error.inc"
module Potential_simple_module

  use libAtoms_module
  use QUIP_Common_module
  use MPI_context_module
  use IP_module
#ifdef HAVE_TB
  use TB_module
#endif
  use FilePot_module
  use CallbackPot_module

#define MAX_CUT_BONDS 4

  implicit none
  private

  public :: potential_simple
  type Potential_Simple
     type(MPI_context) :: mpi

     character(len=124) :: type_args_str

     ! only one of these can be allocated at any given time
#ifdef HAVE_TB
     type(TB_type), pointer     :: tb => null()
#endif
     type(IP_type), pointer     :: ip => null()
     type(FilePot_type), pointer     :: filepot => null()
     type(CallbackPot_type), pointer     :: callbackpot => null()
     logical :: is_wrapper
  end type Potential_Simple

!% Initialise a Potential object (selecting the force field) and, if necessary, the  input file for potential parameters.
  public :: Initialise
  public :: Potential_Simple_Filename_Initialise
  interface Initialise
     module procedure Potential_Simple_Initialise_inoutput, Potential_Simple_Initialise_str
  end interface Initialise

!% Finalise the Potential_Simple object
  public :: Finalise
  interface Finalise
     module procedure Potential_Simple_Finalise
  end interface Finalise

!% Print potential details 
  public :: Print
  interface Print
     module procedure Potential_Simple_Print
  end interface Print

!% Set potential cutoff
  public :: cutoff
  interface cutoff
     module procedure Potential_Simple_cutoff
  end interface cutoff

!% Potential_Simple calculator for energy, forces and virial 
  public :: Calc
  interface Calc
     module procedure Potential_Simple_Calc
  end interface Calc

!% Set up what you need for parallel calculation
  public :: setup_parallel
  interface setup_parallel
     module procedure Potential_Simple_setup_parallel
  end interface setup_parallel

  public :: set_callback
  interface set_callback
     module procedure Potential_Simple_set_callback
  end interface

contains

  subroutine Potential_Simple_Filename_Initialise(this, args_str, filename, mpi_obj, no_parallel, error)
    type(Potential_simple), intent(inout) :: this
    character(len=*), intent(in) :: args_str
    character(len=*), intent(in) :: filename
    type(MPI_context), intent(in), optional :: mpi_obj
    logical, intent(in), optional :: no_parallel
    integer, intent(out), optional :: error

    type(inoutput) :: io

    INIT_ERROR(error)

    ! WARNING: setting master_only=.true. may lead to failure if not all processes are calling the initialise
    call Initialise(io, filename, INPUT, master_only = .true.)
    call Initialise(this, args_str, io, mpi_obj, no_parallel, error=error)
    PASS_ERROR(error)
    call Finalise(io)

  end subroutine Potential_Simple_Filename_Initialise

  subroutine Potential_simple_Initialise_inoutput(this, args_str, io_obj, mpi_obj, no_parallel, error)
    type(Potential_simple), intent(inout) :: this
    character(len=*), intent(in) :: args_str
    type(Inoutput), intent(inout) :: io_obj
    type(MPI_context), intent(in), optional :: mpi_obj
    logical, intent(in), optional :: no_parallel
    integer, intent(out), optional :: error

    type(extendable_str) :: es
    logical my_no_parallel

    INIT_ERROR(error)

    my_no_parallel = optional_default(.false., no_parallel)
    call Initialise(es)
    if (present(mpi_obj)) then
      call read(es, io_obj%unit, convert_to_string=.true., mpi_comm = mpi_obj%communicator)
    else
      call read(es, io_obj%unit, convert_to_string=.true.)
    endif
    if (my_no_parallel) then
      call Initialise(this, args_str, string(es), error=error)
    else
      call Initialise(this, args_str, string(es), mpi_obj, error=error)
    endif
    PASS_ERROR(error)

    call Finalise(es)

  end subroutine Potential_Simple_Initialise_inoutput

  subroutine Potential_Simple_Initialise_str(this, args_str, param_str, mpi_obj, error)
    type(Potential_Simple), intent(inout) :: this
    character(len=*), intent(in) :: args_str
    character(len=*), intent(in), optional :: param_str
    type(MPI_context), intent(in), optional :: mpi_obj
    integer, intent(out), optional :: error

    logical is_TB, is_IP, is_FilePot, is_wrapper, is_callbackpot
    type(Dictionary) :: params

    INIT_ERROR(error)

    call Finalise(this)

    call initialise(params)
    call param_register(params, 'TB', 'false', is_TB)
    call param_register(params, 'IP', 'false', is_IP)
    call param_register(params, 'FilePot', 'false', is_FilePot)
    call param_register(params, 'wrapper', 'false', is_wrapper)
    call param_register(params, 'CallbackPot', 'false', is_CallbackPot)
    if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='Potential_Simple_Initialise_str args_str')) then
      call system_abort("Potential_Simple_Initialise_str failed to parse args_str='"//trim(args_str)//"'")
    endif
    call finalise(params)

    if (count( (/is_TB, is_IP, is_FilePot, is_wrapper, is_callbackpot /) ) /= 1) then
      call system_abort("Potential_Simple_Initialise_str found too few or too many Potential_Simple types args_str='"//trim(args_str)//"'")
    endif

    if (is_IP) then
       allocate(this%ip)
       if(present(param_str)) then
          call Initialise(this%ip, args_str, param_str, mpi_obj, error=error)
	  PASS_ERROR(error)
       else
          RAISE_ERROR('Potential_Simple_initialise: no param_str present during IP init', error)
       endif
#ifdef HAVE_TB
    else if (is_TB) then
       allocate(this%tb)
       if(present(param_str)) then
          call Initialise(this%tb, args_str, param_str, mpi_obj = mpi_obj, error=error)
	  PASS_ERROR(error)
       else
          RAISE_ERROR('Potential_Simple_initialise: no param_str present during TB init', error)
       endif
#else
    else if (is_TB) then
       RAISE_ERROR('Potential_Simple_initialise: TB support not compiled in', error)
#endif
    else if (is_FilePot) then
       allocate(this%filepot)
       call Initialise(this%filepot, args_str, mpi_obj, error=error)
       PASS_ERROR(error)
    else if (is_CallbackPot) then
       allocate(this%CallbackPot)
       call Initialise(this%callbackpot, args_str, mpi=mpi_obj, error=error)
       PASS_ERROR(error)
    else if(is_wrapper) then
       this%is_wrapper = .true.
    endif

    if (present(mpi_obj)) this%mpi = mpi_obj

  end subroutine Potential_Simple_Initialise_str

  subroutine Potential_Simple_Finalise(this, error)
    type(Potential_Simple), intent(inout) :: this
    integer, intent(out), optional :: error

    INIT_ERROR(error)

    if(associated(this%ip)) then
       call Finalise(this%ip)
       deallocate(this%ip)
#ifdef HAVE_TB
    else if(associated(this%tb)) then
       call Finalise(this%tb)
       deallocate(this%tb)
#endif
    elseif(associated(this%filepot)) then
       call Finalise(this%filepot)
       deallocate(this%filepot)
    elseif(associated(this%callbackpot)) then
       call Finalise(this%callbackpot)
       deallocate(this%callbackpot)
    end if
    this%is_wrapper = .false.
  end subroutine Potential_Simple_Finalise

  function Potential_Simple_cutoff(this)
    type(Potential_Simple), intent(in) :: this
    real(dp) :: Potential_Simple_cutoff               

    if(associated(this%ip)) then
       Potential_Simple_cutoff = cutoff(this%ip)
#ifdef HAVE_TB
    else if(associated(this%tb)) then
       Potential_Simple_cutoff = cutoff(this%tb)
#endif
    elseif(associated(this%filepot)) then
       Potential_Simple_cutoff = cutoff(this%filepot)
    elseif(associated(this%callbackpot)) then
       Potential_Simple_cutoff = cutoff(this%callbackpot)
    else
       Potential_Simple_cutoff = 0.0_dp
    end if
  end function Potential_Simple_cutoff

  recursive subroutine Potential_Simple_Calc(this, at, e, local_e, f, df, virial, args_str, err, error)
    type(Potential_Simple), intent(inout) :: this
    type(Atoms), intent(inout) :: at     !% The atoms structure to compute energy and forces
    real(dp), intent(out), optional, target :: e                   !% Total energy
    real(dp), intent(out), optional, target :: local_e(:)          !% Energy per atom    
    real(dp), intent(out), optional, target :: f(:,:)              !% Forces, dimensioned \texttt{(3,at%N)}
    real(dp), intent(out), optional, target :: df(:,:)             !% Finite difference forces, \texttt{(3,at%N)}
    real(dp), intent(out), optional, target :: virial(3,3)         !% Virial
    character(len=*), intent(in), optional :: args_str
    integer, intent(out), optional :: err
    integer, intent(out), optional :: error

    integer:: i,k,n, zero_loc(1)
    real(dp):: e_plus, e_minus, pos_save, r_scale
    real(dp), parameter::delta = 1.0e-4_dp
    type(Dictionary) :: params
    logical :: single_cluster, little_clusters, dummy, do_rescale_r
    character(len=10240) :: my_args_str, new_args_str
    integer, pointer, dimension(:) :: hybrid_mark, cluster_index, termindex, modified_hybrid_mark
    real(dp), pointer, dimension(:) :: weight_region1
    integer, allocatable, dimension(:) :: hybrid_mark_saved
    real(dp), allocatable, dimension(:) :: weight_region1_saved
    real(dp), allocatable, dimension(:,:) :: f_cluster
    logical :: do_carve_cluster
    type(Table) :: cluster_info, cut_bonds
    integer, pointer :: cut_bonds_p(:,:)
    integer :: i_inner, i_outer, n_non_term
    type(Atoms) :: cluster
    real(dp), pointer :: force_ptr(:,:), df_ptr(:,:), local_e_ptr(:)
    real(dp), pointer :: e_ptr, virial_ptr(:,:)
    real(dp), target :: my_e, my_virial(3,3)
    logical :: calc_force, calc_energy, calc_local_e, calc_df, calc_virial, do_calc_force, do_calc_energy, do_calc_local_e, do_calc_df, do_calc_virial
    integer, pointer :: cluster_mark_p(:)
    integer, pointer :: cluster_mark_p_postfix(:)
    integer, pointer :: old_cluster_mark_p(:)
    character(len=FIELD_LENGTH) :: cluster_mark_postfix 
    logical save_mpi

    INIT_ERROR(error)

    if (at%N <= 0) then
      RAISE_ERROR("Potential_Simple_Calc called with at%N <= 0", error)
    endif

    if (present(err)) err = 0

    if (present(args_str)) then
      call print('Potential_Simple_calc got args_str "'//trim(args_str)//'"', PRINT_VERBOSE)
      my_args_str = args_str
    else
      call print('Potential_Simple_calc got not args_str', PRINT_VERBOSE)
      my_args_str = ""
    endif

    call initialise(params)
    call param_register(params, 'single_cluster', 'F', single_cluster)
    cluster_mark_postfix=""
    call param_register(params, 'cluster_mark_postfix', '', cluster_mark_postfix)
    call param_register(params, 'carve_cluster', 'T', do_carve_cluster)
    call param_register(params, 'little_clusters', 'F', little_clusters)
    call param_register(params, 'do_rescale_r', 'F', do_rescale_r)
    call param_register(params, 'r_scale', '1.0', r_scale)

    call param_register(params, 'calc_force', 'F', calc_force)
    call param_register(params, 'calc_energy', 'F', calc_energy)
    call param_register(params, 'calc_local_e', 'F', calc_local_e)
    call param_register(params, 'calc_df', 'F', calc_df)
    call param_register(params, 'calc_virial', 'F', calc_virial)

    if (.not. param_read_line(params, my_args_str, ignore_unknown=.true.,task='Potential_Simple_Calc_str args_str') ) then
      RAISE_ERROR("Potential_Simple_calc failed to parse args_str='"//trim(my_args_str)//"'", error)
    endif
    call finalise(params)

    if (single_cluster .and. little_clusters) then
         RAISE_ERROR('Potential_Simple_calc: single_cluster and little_clusters options are mutually exclusive', error)
    endif

    if (little_clusters) then

       if (present(e) .or. present(local_e) .or. present(df) .or. present(virial) .or. .not. present(f)) then
            RAISE_ERROR('Potential_Simple_calc: little_clusters option only supports calcualtion of forces, not energies, local energies or virials', error)
       endif

       ! must remove "little_clusters" from args_str so that recursion terminates
       call initialise(params)
       call read_string(params, my_args_str)
       call remove_value(params, 'little_clusters')
       new_args_str = write_string(params)
       call finalise(params)
 
       if (.not. assign_pointer(at, 'hybrid_mark', hybrid_mark)) then
            RAISE_ERROR('Potential_Simple_calc: cannot assign pointer to hybrid_mark property ', error)
       endif

       if (.not. any(hybrid_mark == HYBRID_ACTIVE_MARK)) then
          f = 0.0_dp
          return
       end if

       ! Save hybrid_mark and weight_region1, because we're going to overwrite them
       allocate(hybrid_mark_saved(at%N))
       hybrid_mark_saved = hybrid_mark

       if (has_property(at, 'weight_region1')) then
          dummy = assign_pointer(at, 'weight_region1', weight_region1)
          allocate(weight_region1_saved(at%N))
          weight_region1_saved = weight_region1
       end if
       
       ! call ourselves once for each active or transition atom 
       f = 0.0_dp
       n = 0
       do i=1,at%N
          if (hybrid_mark_saved(i) /= HYBRID_ACTIVE_MARK .and. hybrid_mark_saved(i) /= HYBRID_TRANS_MARK) cycle

          if (this%mpi%active) then
             n = n + 1
             call print('Potential_Simple_calc: cluster '//n//' around atom '//i//'  assigned to proc '//mod(n-1,this%mpi%n_procs)//' of '//(this%mpi%n_procs), PRINT_VERBOSE)
             if (mod(n-1, this%mpi%n_procs) .ne. this%mpi%my_proc) cycle
          end if
          call print('Potential_Simple_calc: constructing little_cluster around atom '//i, PRINT_VERBOSE)
          hybrid_mark = HYBRID_NO_MARK
          hybrid_mark(i) = HYBRID_ACTIVE_MARK
          call create_hybrid_weights(at, new_args_str)
          cluster_info = create_cluster_info_from_hybrid_mark(at, new_args_str)
	  cluster = carve_cluster(at, new_args_str, cluster_info)
	  call finalise(cluster_info)
          allocate(f_cluster(3,cluster%N))

          ! Reassign pointers - create_cluster_info_from_hybrid_mark() might have broken them
          if (has_property(at, 'hybrid_mark')) &
               dummy = assign_pointer(at, 'hybrid_mark', hybrid_mark)
          if (has_property(at, 'weight_region1')) &
               dummy = assign_pointer(at, 'weight_region1', weight_region1)

	  if (current_verbosity() >= PRINT_NERD) then
	    call write(cluster, 'stdout', prefix='LITTLE_CLUSTER')
	  endif
          call print('ARGS0 | '//new_args_str,PRINT_VERBOSE)

          ! Disable MPI for duration of calc() call, since we're doing parallelisatin at level of clusters
          save_mpi = this%mpi%active
          this%mpi%active = .false.
          call calc(this, cluster, f=f_cluster, args_str=new_args_str)
          this%mpi%active = save_mpi

          if (do_rescale_r)  f_cluster = f_cluster*r_scale
          f(:,i) = f_cluster(:,1)
          deallocate(f_cluster)
	   call finalise(cluster)
       end do

       if (this%mpi%active) then
          call sum_in_place(this%mpi, f)
       end if
       hybrid_mark = hybrid_mark_saved
       deallocate(hybrid_mark_saved)

       if (allocated(weight_region1_saved)) then
          weight_region1 = weight_region1_saved
          deallocate(weight_region1_saved)
       end if

    else if (single_cluster) then

       if (present(e) .or. present(local_e) .or. present(df) .or. present(virial) .or. .not. present(f)) then
            RAISE_ERROR('Potential_Simple_calc: single_cluster option only supports calcualtion of forces, not energies, local energies or virials', error)
       endif

       if (.not. assign_pointer(at, 'hybrid_mark', hybrid_mark)) then
            RAISE_ERROR('Potential_Simple_calc: cannot assign pointer to hybrid_mark property ', error)
       endif

       if (.not. any(hybrid_mark == HYBRID_ACTIVE_MARK)) then
          f = 0.0_dp
          return
       end if

       ! must remove "single_cluster" from args_str so that recursion terminates
       call initialise(params)
       call read_string(params, my_args_str)
       call remove_value(params, 'single_cluster')
       new_args_str = write_string(params)
       call finalise(params)

       ! call ourselves on a cluster formed from marked atoms
       call print('Potential_Simple_calc: constructing single_cluster', PRINT_VERBOSE)

       if (do_carve_cluster) then
	 call print('Potential_Simple_calc: carving cluster', PRINT_VERBOSE)
	 cluster_info = create_cluster_info_from_hybrid_mark(at, new_args_str, error=error)
	 PASS_ERROR_WITH_INFO("potential_calc: creating cluster info from hybrid_mark", error)

         ! Check there are no repeated indices among the non-termination atoms in the cluster
         n_non_term = count(cluster_info%int(6,1:cluster_info%n) == 0)
         if (multiple_images(int_subtable(cluster_info,(/ (i,i=1,n_non_term) /),(/1/)))) then
              RAISE_ERROR('Potential_Simple_calc: single_cluster=T not yet implemented when cluster contains repeated periodic images', error)
	 endif

	 cluster = carve_cluster(at, new_args_str, cluster_info, error)
	 PASS_ERROR_WITH_INFO("potential_calc: carving cluster", error)
	 call finalise(cluster_info)
	 if (current_verbosity() >= PRINT_NERD) then
           ! prefix should be CLUSTER
	   call write(cluster, 'stdout')
	 endif
	 if (.not. assign_pointer(cluster, 'index', cluster_index)) then
	      RAISE_ERROR('Potential_Simple_calc: cluster is missing index property', error)
	 endif
	 if (.not. assign_pointer(cluster, 'termindex', termindex)) then
	      RAISE_ERROR('Potential_Simple_calc: cluster is missing termindex property', error)
	 endif
	 allocate(f_cluster(3,cluster%N))
         call print('ARGS1 | '//new_args_str,PRINT_VERBOSE)

	 call calc(this, cluster, f=f_cluster, args_str=new_args_str, error=error)
	 PASS_ERROR_WITH_INFO('potential_calc after calc in carve_cluster', error)
	 if (do_rescale_r)  f_cluster = f_cluster*r_scale

         ! Reassign pointers - create_cluster_info_from_hybrid_mark() might have broken them
         if (has_property(at, 'hybrid_mark')) &
              dummy = assign_pointer(at, 'hybrid_mark', hybrid_mark)

	 ! copy forces for all active and transition atoms
	 f = 0.0_dp
	 do i=1,cluster%N
	    if (termindex(i) /= 0) cycle ! skip termination atoms
	    if (hybrid_mark(cluster_index(i)) == HYBRID_ACTIVE_MARK .or. &
		hybrid_mark(cluster_index(i)) == HYBRID_TRANS_MARK) &
		  f(:,cluster_index(i)) = f_cluster(:,i)
	 end do
	 deallocate(f_cluster)   
	 call finalise(cluster)
       else ! not do_carve_cluster
	 call print('Potential_Simple_calc: not carving cluster', PRINT_VERBOSE)
	 cluster_info = create_cluster_info_from_hybrid_mark(at, trim(new_args_str) // " cluster_same_lattice", cut_bonds, error)
	 PASS_ERROR_WITH_INFO('potential_calc creating cluster info from hybrid mark with carve_cluster=F', error)

         !save cluster in cluster_mark property and optionally cluster_mark_postfix property
         call add_property(at,'cluster_mark',HYBRID_NO_MARK)
         if (trim(cluster_mark_postfix)/="") then
           call print('Add cluster_mark'//trim(cluster_mark_postfix),PRINT_ANAL)
           call add_property(at,'cluster_mark'//trim(cluster_mark_postfix),HYBRID_NO_MARK)
         else
           call print('NOT Add cluster_mark'//trim(cluster_mark_postfix),PRINT_ANAL)
         endif
         !save the previous cluster_mark[_postfix] into old_cluster_mark[_postfix]
         call add_property(at,'old_cluster_mark'//trim(cluster_mark_postfix),HYBRID_NO_MARK)
	 call print('Add old_cluster_mark'//trim(cluster_mark_postfix),PRINT_ANAL)
	 if (.not. assign_pointer(at, 'cluster_mark', cluster_mark_p)) then
	   RAISE_ERROR("Potential_Simple_calc failed to assing pointer for cluster_mark"//trim(cluster_mark_postfix)//" pointer", error)
	 endif
         if (trim(cluster_mark_postfix)/="") then
	   if (.not. assign_pointer(at, 'cluster_mark'//trim(cluster_mark_postfix), cluster_mark_p_postfix)) then
	     RAISE_ERROR("Potential_Simple_calc failed to assing pointer for cluster_mark pointer", error)
	   endif
           call print('Assign cluster_mark'//trim(cluster_mark_postfix),PRINT_ANAL)
         else
           call print('NOT Assign cluster_mark'//trim(cluster_mark_postfix),PRINT_ANAL)
         endif
	 if (.not. assign_pointer(at, 'old_cluster_mark'//trim(cluster_mark_postfix), old_cluster_mark_p)) then
	   RAISE_ERROR("Potential_Simple_calc failed to assing pointer for old_cluster_mark pointer", error)
	 endif

         !save old from cluster_mark_postfix
         if (trim(cluster_mark_postfix)/="") then
           old_cluster_mark_p = cluster_mark_p_postfix
         else
           old_cluster_mark_p = cluster_mark_p
         endif

         !zero cluster_mark_[postfix]
         cluster_mark_p = HYBRID_NO_MARK
         if (trim(cluster_mark_postfix)/="") cluster_mark_p_postfix = HYBRID_NO_MARK

         !save modified_hybrid_mark into cluster_mark[_postfix]
	 if (.not. assign_pointer(at, 'modified_hybrid_mark', modified_hybrid_mark)) then
	   cluster_mark_p = hybrid_mark
	 else
	   cluster_mark_p = modified_hybrid_mark
	 endif
         if (trim(cluster_mark_postfix)/="") then
           call print('set value cluster_mark'//trim(cluster_mark_postfix))
           cluster_mark_p_postfix = cluster_mark_p
         else
           call print('NOT set value cluster_mark'//trim(cluster_mark_postfix))
         endif

         !save cut bonds in cut_bonds property
	 call add_property(at, 'cut_bonds', 0, n_cols=MAX_CUT_BONDS)
	 if (.not. assign_pointer(at, 'cut_bonds', cut_bonds_p)) then
	   RAISE_ERROR("Potential_Simple_calc failed to assing pointer for cut_bonds pointer", error)
	 endif
         !zero it
         cut_bonds_p = 0
	 do i=1, cut_bonds%N
	   i_inner = cut_bonds%int(1,i)
	   i_outer = cut_bonds%int(2,i)
	   zero_loc = minloc(cut_bonds_p(:,i_inner))
	   if (cut_bonds_p(zero_loc(1),i_inner) == 0) then ! free space for a cut bond
	     cut_bonds_p(zero_loc(1),i_inner) = i_outer
	   else
	     call print("cut_bonds table:", PRINT_VERBOSE)
	     call print(cut_bonds, PRINT_VERBOSE)
	     call print("ERROR: Potential_Simple_calc ran out of space to store cut_bonds information", PRINT_ALWAYS)
	     call print("ERROR: inner atom " // i_inner // " already has cut_bonds to " // cut_bonds_p(:,i_inner) // &
	      " no space to add cut bond to " // i_outer, PRINT_ALWAYS)
	     RAISE_ERROR("Potential_Simple_calc out of space to store cut_bonds information", error)
	   endif
	 end do
	 call finalise(cut_bonds)
	 if (current_verbosity() >= PRINT_ANAL) then
	   ! prefix should be "UNCARVED_CLUSTER"
	   call write(at, 'stdout')
	 endif
	 call calc(this, at, f=f, args_str=new_args_str, error=error)
	 PASS_ERROR_WITH_INFO('potential_calc after calc with carve_cluster=F', error)
	 if (do_rescale_r)  f = f*r_scale
       endif ! do_carve_cluster
    else ! little_clusters and single_cluster are false..

       ! For IP, call setup_atoms() hook now in case any properties must be added.
       ! This must be done *before* we assign pointers to force, local_e etc.
       if (associated(this%ip)) then
          call setup_atoms(this%ip, at)
       end if

       do_calc_force = calc_force .or. present(f)
       do_calc_energy = calc_energy .or. present(e)
       do_calc_local_e = calc_local_e .or. present(local_e)
       do_calc_df = calc_df .or. present(df)
       do_calc_virial = calc_virial .or. present(virial)

       if(do_calc_virial .or. do_calc_energy .or. do_calc_force .or. do_calc_local_e) then

          if (do_calc_force) then
             if (.not. present(f)) then
                if (.not. has_property(at, 'force')) call add_property(at, 'force', 0.0_dp, n_cols=3)
             else
                force_ptr => f
             end if
          end if

          if (do_calc_local_e) then
             if (.not. present(local_e)) then
                if (.not. has_property(at, 'local_e')) call add_property(at, 'local_e', 0.0_dp)
             else
                local_e_ptr => local_e
             end if
          end if

          if (do_calc_energy) then
             if (.not. present(e)) then
                e_ptr => my_e
             else
                e_ptr => e
             end if
          end if

          if (do_calc_virial) then
             if (.not. present(virial)) then
                virial_ptr => my_virial
             else
                virial_ptr => virial
             end if
          end if

          ! Do pointer assignments after all properties have been added
          if (do_calc_force .and. .not. present(f)) then
             if (.not. assign_pointer(at, 'force', force_ptr)) call system_abort('Potential_Simple_calc: cannot assign force_ptr')
          end if
          if (do_calc_local_e .and. .not. present(local_e)) then
             if (.not. assign_pointer(at, 'local_e', local_e_ptr)) call system_abort('Potential_Simple_calc: cannot assign local_e_ptr')
          end if
          
          if(associated(this%ip)) then

             if (do_calc_virial) then
                if (.not. do_calc_energy .and. .not. do_calc_local_e .and. .not. do_calc_force) then
                   call Calc(this%ip, at, virial=virial_ptr, args_str=args_str, error=error)
                else if (.not. do_calc_energy .and. .not. do_calc_local_e .and. do_calc_force) then
                   call Calc(this%ip, at, f=force_ptr, virial=virial_ptr, args_str=args_str, error=error)
                else if (.not. do_calc_energy .and. do_calc_local_e .and. .not. do_calc_force) then
                   call Calc(this%ip, at, local_e=local_e_ptr, virial=virial_ptr, args_str=args_str, error=error)
                else if (.not. do_calc_energy .and. do_calc_local_e .and. do_calc_force) then
                   call Calc(this%ip, at, local_e=local_e_ptr, f=force_ptr, virial=virial_ptr, args_str=args_str, error=error)
                else if (do_calc_energy .and. .not. do_calc_local_e .and. .not. do_calc_force) then
                   call Calc(this%ip, at, energy=e_ptr, virial=virial_ptr, args_str=args_str, error=error)
                else if (do_calc_energy .and. .not. do_calc_local_e .and. do_calc_force) then
                   call Calc(this%ip, at, energy=e_ptr, f=force_ptr, virial=virial_ptr, args_str=args_str, error=error)
                else if (do_calc_energy .and. do_calc_local_e .and. .not. do_calc_force) then
                   call Calc(this%ip, at, energy=e_ptr, local_e=local_e_ptr, virial=virial_ptr, args_str=args_str, error=error)
                else if (do_calc_energy .and. do_calc_local_e .and. do_calc_force) then
                   call Calc(this%ip, at, energy=e_ptr, local_e=local_e_ptr, f=force_ptr, virial=virial_ptr, args_str=args_str, error=error)
                end if
             else
                if (.not. do_calc_energy .and. .not. do_calc_local_e .and. .not. do_calc_force) then
                   call Calc(this%ip, at, args_str=args_str, error=error)
                else if (.not. do_calc_energy .and. .not. do_calc_local_e .and. do_calc_force) then
                   call Calc(this%ip, at, f=force_ptr, args_str=args_str, error=error)
                else if (.not. do_calc_energy .and. do_calc_local_e .and. .not. do_calc_force) then
                   call Calc(this%ip, at, local_e=local_e_ptr, args_str=args_str, error=error)
                else if (.not. do_calc_energy .and. do_calc_local_e .and. do_calc_force) then
                   call Calc(this%ip, at, local_e=local_e_ptr, f=force_ptr, args_str=args_str, error=error)
                else if (do_calc_energy .and. .not. do_calc_local_e .and. .not. do_calc_force) then
                   call Calc(this%ip, at, energy=e_ptr, args_str=args_str, error=error)
                else if (do_calc_energy .and. .not. do_calc_local_e .and. do_calc_force) then
                   call Calc(this%ip, at, energy=e_ptr, f=force_ptr, args_str=args_str, error=error)
                else if (do_calc_energy .and. do_calc_local_e .and. .not. do_calc_force) then
                   call Calc(this%ip, at, energy=e_ptr, local_e=local_e_ptr, args_str=args_str, error=error)
                else if (do_calc_energy .and. do_calc_local_e .and. do_calc_force) then
                   call Calc(this%ip, at, energy=e_ptr, local_e=local_e_ptr, f=force_ptr, args_str=args_str, error=error)
                end if
             end if
#ifdef HAVE_TB
          else if(associated(this%tb)) then

             if (do_calc_virial) then
                if (.not. do_calc_energy .and. .not. do_calc_local_e .and. .not. do_calc_force) then
                   call Calc(this%tb, at, virial=virial_ptr, args_str=args_str, error=error)
                else if (.not. do_calc_energy .and. .not. do_calc_local_e .and. do_calc_force) then
                   call Calc(this%tb, at, forces=force_ptr, virial=virial_ptr, args_str=args_str, error=error)
                else if (.not. do_calc_energy .and. do_calc_local_e .and. .not. do_calc_force) then
                   call Calc(this%tb, at, local_e=local_e_ptr, virial=virial_ptr, args_str=args_str, error=error)
                else if (.not. do_calc_energy .and. do_calc_local_e .and. do_calc_force) then
                   call Calc(this%tb, at, local_e=local_e_ptr, forces=force_ptr, virial=virial_ptr, args_str=args_str, error=error)
                else if (do_calc_energy .and. .not. do_calc_local_e .and. .not. do_calc_force) then
                   call Calc(this%tb, at, energy=e_ptr, virial=virial_ptr, args_str=args_str, error=error)
                else if (do_calc_energy .and. .not. do_calc_local_e .and. do_calc_force) then
                   call Calc(this%tb, at, energy=e_ptr, forces=force_ptr, virial=virial_ptr, args_str=args_str, error=error)
                else if (do_calc_energy .and. do_calc_local_e .and. .not. do_calc_force) then
                   call Calc(this%tb, at, energy=e_ptr, local_e=local_e_ptr, virial=virial_ptr, args_str=args_str, error=error)
                else if (do_calc_energy .and. do_calc_local_e .and. do_calc_force) then
                   call Calc(this%tb, at, energy=e_ptr, local_e=local_e_ptr, forces=force_ptr, virial=virial_ptr, args_str=args_str, error=error)
                end if
             else
                if (.not. do_calc_energy .and. .not. do_calc_local_e .and. .not. do_calc_force) then
                   call Calc(this%tb, at, args_str=args_str, error=error)
                else if (.not. do_calc_energy .and. .not. do_calc_local_e .and. do_calc_force) then
                   call Calc(this%tb, at, forces=force_ptr, args_str=args_str, error=error)
                else if (.not. do_calc_energy .and. do_calc_local_e .and. .not. do_calc_force) then
                   call Calc(this%tb, at, local_e=local_e_ptr, args_str=args_str, error=error)
                else if (.not. do_calc_energy .and. do_calc_local_e .and. do_calc_force) then
                   call Calc(this%tb, at, local_e=local_e_ptr, forces=force_ptr, args_str=args_str, error=error)
                else if (do_calc_energy .and. .not. do_calc_local_e .and. .not. do_calc_force) then
                   call Calc(this%tb, at, energy=e_ptr, args_str=args_str, error=error)
                else if (do_calc_energy .and. .not. do_calc_local_e .and. do_calc_force) then
                   call Calc(this%tb, at, energy=e_ptr, forces=force_ptr, args_str=args_str, error=error)
                else if (do_calc_energy .and. do_calc_local_e .and. .not. do_calc_force) then
                   call Calc(this%tb, at, energy=e_ptr, local_e=local_e_ptr, args_str=args_str, error=error)
                else if (do_calc_energy .and. do_calc_local_e .and. do_calc_force) then
                   call Calc(this%tb, at, energy=e_ptr, local_e=local_e_ptr, forces=force_ptr, args_str=args_str, error=error)
                end if                
             end if

#endif
          elseif(associated(this%filepot)) then

             if (do_calc_virial) then
                if (.not. do_calc_energy .and. .not. do_calc_local_e .and. .not. do_calc_force) then
                   call Calc(this%filepot, at, virial=virial_ptr, args_str=args_str, error=error)
                else if (.not. do_calc_energy .and. .not. do_calc_local_e .and. do_calc_force) then
                   call Calc(this%filepot, at, forces=force_ptr, virial=virial_ptr, args_str=args_str, error=error)
                else if (.not. do_calc_energy .and. do_calc_local_e .and. .not. do_calc_force) then
                   call Calc(this%filepot, at, local_e=local_e_ptr, virial=virial_ptr, args_str=args_str, error=error)
                else if (.not. do_calc_energy .and. do_calc_local_e .and. do_calc_force) then
                   call Calc(this%filepot, at, local_e=local_e_ptr, forces=force_ptr, virial=virial_ptr, args_str=args_str, error=error)
                else if (do_calc_energy .and. .not. do_calc_local_e .and. .not. do_calc_force) then
                   call Calc(this%filepot, at, energy=e_ptr, virial=virial_ptr, args_str=args_str, error=error)
                else if (do_calc_energy .and. .not. do_calc_local_e .and. do_calc_force) then
                   call Calc(this%filepot, at, energy=e_ptr, forces=force_ptr, virial=virial_ptr, args_str=args_str, error=error)
                else if (do_calc_energy .and. do_calc_local_e .and. .not. do_calc_force) then
                   call Calc(this%filepot, at, energy=e_ptr, local_e=local_e_ptr, virial=virial_ptr, args_str=args_str, error=error)
                else if (do_calc_energy .and. do_calc_local_e .and. do_calc_force) then
                   call Calc(this%filepot, at, energy=e_ptr, local_e=local_e_ptr, forces=force_ptr, virial=virial_ptr, args_str=args_str, error=error)
                end if
             else
                if (.not. do_calc_energy .and. .not. do_calc_local_e .and. .not. do_calc_force) then
                   call Calc(this%filepot, at, args_str=args_str, error=error)
                else if (.not. do_calc_energy .and. .not. do_calc_local_e .and. do_calc_force) then
                   call Calc(this%filepot, at, forces=force_ptr, args_str=args_str, error=error)
                else if (.not. do_calc_energy .and. do_calc_local_e .and. .not. do_calc_force) then
                   call Calc(this%filepot, at, local_e=local_e_ptr, args_str=args_str, error=error)
                else if (.not. do_calc_energy .and. do_calc_local_e .and. do_calc_force) then
                   call Calc(this%filepot, at, local_e=local_e_ptr, forces=force_ptr, args_str=args_str, error=error)
                else if (do_calc_energy .and. .not. do_calc_local_e .and. .not. do_calc_force) then
                   call Calc(this%filepot, at, energy=e_ptr, args_str=args_str, error=error)
                else if (do_calc_energy .and. .not. do_calc_local_e .and. do_calc_force) then
                   call Calc(this%filepot, at, energy=e_ptr, forces=force_ptr, args_str=args_str, error=error)
                else if (do_calc_energy .and. do_calc_local_e .and. .not. do_calc_force) then
                   call Calc(this%filepot, at, energy=e_ptr, local_e=local_e_ptr, args_str=args_str, error=error)
                else if (do_calc_energy .and. do_calc_local_e .and. do_calc_force) then
                   call Calc(this%filepot, at, energy=e_ptr, local_e=local_e_ptr, forces=force_ptr, args_str=args_str, error=error)
                end if
             end if

          elseif(associated(this%callbackpot)) then

             if (do_calc_virial) then
                if (.not. do_calc_energy .and. .not. do_calc_local_e .and. .not. do_calc_force) then
                   call Calc(this%callbackpot, at, virial=virial_ptr, args_str=args_str, error=error)
                else if (.not. do_calc_energy .and. .not. do_calc_local_e .and. do_calc_force) then
                   call Calc(this%callbackpot, at, forces=force_ptr, virial=virial_ptr, args_str=args_str, error=error)
                else if (.not. do_calc_energy .and. do_calc_local_e .and. .not. do_calc_force) then
                   call Calc(this%callbackpot, at, local_e=local_e_ptr, virial=virial_ptr, args_str=args_str, error=error)
                else if (.not. do_calc_energy .and. do_calc_local_e .and. do_calc_force) then
                   call Calc(this%callbackpot, at, local_e=local_e_ptr, forces=force_ptr, virial=virial_ptr, args_str=args_str, error=error)
                else if (do_calc_energy .and. .not. do_calc_local_e .and. .not. do_calc_force) then
                   call Calc(this%callbackpot, at, energy=e_ptr, virial=virial_ptr, args_str=args_str, error=error)
                else if (do_calc_energy .and. .not. do_calc_local_e .and. do_calc_force) then
                   call Calc(this%callbackpot, at, energy=e_ptr, forces=force_ptr, virial=virial_ptr, args_str=args_str, error=error)
                else if (do_calc_energy .and. do_calc_local_e .and. .not. do_calc_force) then
                   call Calc(this%callbackpot, at, energy=e_ptr, local_e=local_e_ptr, virial=virial_ptr, args_str=args_str, error=error)
                else if (do_calc_energy .and. do_calc_local_e .and. do_calc_force) then
                   call Calc(this%callbackpot, at, energy=e_ptr, local_e=local_e_ptr, forces=force_ptr, virial=virial_ptr, args_str=args_str, error=error)
                end if
             else
                if (.not. do_calc_energy .and. .not. do_calc_local_e .and. .not. do_calc_force) then
                   call Calc(this%callbackpot, at, args_str=args_str, error=error)
                else if (.not. do_calc_energy .and. .not. do_calc_local_e .and. do_calc_force) then
                   call Calc(this%callbackpot, at, forces=force_ptr, args_str=args_str, error=error)
                else if (.not. do_calc_energy .and. do_calc_local_e .and. .not. do_calc_force) then
                   call Calc(this%callbackpot, at, local_e=local_e_ptr, args_str=args_str, error=error)
                else if (.not. do_calc_energy .and. do_calc_local_e .and. do_calc_force) then
                   call Calc(this%callbackpot, at, local_e=local_e_ptr, forces=force_ptr, args_str=args_str, error=error)
                else if (do_calc_energy .and. .not. do_calc_local_e .and. .not. do_calc_force) then
                   call Calc(this%callbackpot, at, energy=e_ptr, args_str=args_str, error=error)
                else if (do_calc_energy .and. .not. do_calc_local_e .and. do_calc_force) then
                   call Calc(this%callbackpot, at, energy=e_ptr, forces=force_ptr, args_str=args_str, error=error)
                else if (do_calc_energy .and. do_calc_local_e .and. .not. do_calc_force) then
                   call Calc(this%callbackpot, at, energy=e_ptr, local_e=local_e_ptr, args_str=args_str, error=error)
                else if (do_calc_energy .and. do_calc_local_e .and. do_calc_force) then
                   call Calc(this%callbackpot, at, energy=e_ptr, local_e=local_e_ptr, forces=force_ptr, args_str=args_str, error=error)
                end if
             end if

          elseif(this%is_wrapper) then
             !
             ! put here hardcoded energy and force functions
             !
             !if(present(e)) e = wrapper_energy(at)
             !if(present(f)) call wrapper_force(at, f)
             RAISE_ERROR("Potential_Simple_Calc: hardcoded wrapper functions are not defined", error)
          else
             RAISE_ERROR ("Potential_Simple_Calc: Potential_Simple is not initialised", error)
          end if
	 PASS_ERROR_WITH_INFO('potential_calc after actual calc', error)

          if (calc_energy) then
             call set_value(at%params, 'energy', e_ptr)
          end if

          if (calc_virial) then
             call set_value(at%params, 'virial', virial_ptr)
          end if

          ! Copy force and local_e to properties as well if necessary
          if (calc_force .and. present(f)) then
             if (.not. has_property(at, 'force')) call add_property(at, 'force', 0.0_dp, n_cols=3)
             if (.not. assign_pointer(at, 'force', force_ptr)) then
	       RAISE_ERROR('Potential_Simple_calc: cannot assign force_ptr', error)
	     endif
             force_ptr(:,:) = f(:,:)
          end if

          if (calc_local_e .and. present(local_e)) then
             if (.not. has_property(at, 'local_e')) call add_property(at, 'local_e', 0.0_dp)
             if (.not. assign_pointer(at, 'local_e', local_e_ptr)) then
	        RAISE_ERROR('Potential_Simple_calc: cannot assign local_e_ptr', error)
	     endif
             local_e_ptr(:) = local_e(:)
          end if

       end if

       if (do_calc_df) then ! do forces by finite difference

          ! must remove 'calc_df' from args_str if it's there (and the rest, or properties get overwritten)
          call initialise(params)
          call read_string(params, my_args_str)
          call remove_value(params, 'calc_df')
          call remove_value(params, 'calc_force')
          call remove_value(params, 'calc_energy')
          call remove_value(params, 'calc_local_e')
          call remove_value(params, 'calc_virial')
          new_args_str = write_string(params)
          call finalise(params)

          if (.not. present(df)) then
             if (.not. has_property(at, 'df')) call add_property(at, 'df', 0.0_dp, n_cols=3)
             if (.not. assign_pointer(at, 'df', df_ptr)) then
	       RAISE_ERROR('Potential_Simple_calc: cannot assign df_ptr', error)
	     endif
          else
             df_ptr => df
          end if

          do i=1,at%N
             do k=1,3
                pos_save = at%pos(k,i)
                at%pos(k,i) = pos_save + delta
                call calc_dists(at)
                call calc(this, at, e=e_plus, args_str=new_args_str, err=err)
                at%pos(k,i) = pos_save - delta
                call calc_dists(at)
                call calc(this, at, e=e_minus, args_str=new_args_str, err=err)
                at%pos(k,i) = pos_save
                call calc_dists(at)
                df_ptr(k,i) = (e_minus-e_plus)/(2.0_dp*delta) ! force is -ve gradient
             end do
          end do

          ! Copy to property as well if necessary
          if (calc_df .and. present(df)) then
             if (.not. has_property(at, 'df')) call add_property(at, 'df', 0.0_dp, n_cols=3)
             if (.not. assign_pointer(at, 'df', df_ptr)) then
	       RAISE_ERROR('Potential_Simple_calc: cannot assign df_ptr', error)
	     endif
             df_ptr(:,:) = df(:,:)
          end if
       end if
    end if

  end subroutine Potential_Simple_Calc

!  function Potential_Simple_Brief_Description(this) result desc
!    type(Potential_Simple), intent(inout) :: this
!    character(30), intent(out) :: desc
!    if(associated(this%tb)) then
!       desc = ''
!    elseif(associated(this%ip)) then
!       desc = Brief_Description(this%ip)
!    elseif(associated(this%filepot)) then
!       desc = ''
!    elseif(this%is_wrapper) then
!       desc = ''
!    else
!       call system_abort ("Potential_Simple_Brief_Description: Potential_Simple is not initialised")
!    end if
!
!  end function Potential_Simple_Brief_Description


  subroutine Potential_Simple_Print(this, file, error)
    type(Potential_Simple), intent(inout) :: this
    type(Inoutput), intent(inout),optional:: file
    integer, intent(out), optional :: error

    INIT_ERROR(error)

    if(associated(this%ip)) then
       call Print(this%ip, file=file)
#ifdef HAVE_TB
    else if(associated(this%tb)) then
       call Print(this%tb, file=file)
#endif
    elseif(associated(this%filepot)) then
       call Print(this%filepot, file=file)
    elseif(associated(this%callbackpot)) then
       call Print(this%callbackpot, file=file)
    elseif(this%is_wrapper) then
       call print("Potential_Simple: wrapper Potential_Simple")
    else
       RAISE_ERROR ("Potential_Simple_Print: no potential type is set", error)
    end if

  end subroutine Potential_Simple_Print

  subroutine Potential_Simple_setup_parallel(this, at, e, local_e, f, virial, args_str, error)
    type(Potential_Simple), intent(inout) :: this
    type(Atoms), intent(inout) :: at     !% The atoms structure to compute energy and forces
    real(dp), intent(out), optional :: e                   !% Total energy
    real(dp), intent(out), optional :: local_e(:)          !% Energy per atom    
    real(dp), intent(out), optional :: f(:,:)              !% Forces, dimensioned \texttt{(3,at%N)}
    real(dp), intent(out), optional :: virial(3,3)         !% Virial
    character(len=*), intent(in), optional :: args_str
    integer, intent(out), optional :: error

    INIT_ERROR(error)

    if(associated(this%ip)) then
       call setup_parallel(this%ip, at, e, local_e, f, virial)
#ifdef HAVE_TB
    else if(associated(this%tb)) then
       return
#endif
    elseif(associated(this%filepot)) then
       return
    elseif(associated(this%callbackpot)) then
       return
    elseif(this%is_wrapper) then
       return
    else
       RAISE_ERROR ("Potential_Simple_Print: Potential_Simple is not initialised", error)
    end if

  end subroutine Potential_Simple_setup_parallel

  subroutine Potential_Simple_set_callback(this, callback, error)
    type(Potential_Simple), intent(inout) :: this
    interface
       subroutine callback(at)
#ifdef HAVE_QUIPPY
         integer, intent(in) :: at(12)
#else
         use Atoms_module, only: Atoms
         type(Atoms), intent(inout) :: at
#endif
       end subroutine callback
    end interface
    integer, intent(out), optional :: error

    INIT_ERROR(error)
    
    if (.not. associated(this%callbackpot)) then
      RAISE_ERROR('Potential_Simple_set_callback: this Potential_Simple is not a CallbackPot', error)
    endif
    call set_callback(this%callbackpot, callback)

  end subroutine Potential_Simple_set_callback

end module Potential_Simple_module
