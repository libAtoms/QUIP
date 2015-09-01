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

  use error_module
  use system_module, only : dp, inoutput, PRINT_VERBOSE, PRINT_NERD, PRINT_ALWAYS, PRINT_ANAL, current_verbosity, INPUT, optional_default, parse_string, verbosity_push, verbosity_pop
  use extendable_str_module
  use dictionary_module
  use paramreader_module
  use linearalgebra_module
  use table_module
  use atoms_types_module
  use atoms_module
  use connection_module, only: is_in_subregion
  use cinoutput_module
  use clusters_module
  use partition_module

  use QUIP_Common_module
  use MPI_context_module
  use IP_module
#ifdef HAVE_TB
  use TB_module
#endif
  use FilePot_module
  use CallbackPot_module
  use SocketPot_module

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
     type(SocketPot_type), pointer :: socketpot => null()
     logical :: is_wrapper
     logical :: little_clusters
     logical :: force_using_fd
     logical :: virial_using_fd
     
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

    logical is_TB, is_IP, is_FilePot, is_wrapper, is_callbackpot, is_socketpot
    type(Dictionary) :: params

    INIT_ERROR(error)

    call Finalise(this)

    call initialise(params)
    call param_register(params, 'TB', 'false', is_TB, help_string="If true, a tight-binding model")
    call param_register(params, 'IP', 'false', is_IP, help_string="If true, an interatomic potential model")
    call param_register(params, 'FilePot', 'false', is_FilePot, help_string="If true, a potential that interacts with another executable by reading/writing files")
    call param_register(params, 'wrapper', 'false', is_wrapper, help_string="If true, a hardcoded wrapper function")
    call param_register(params, 'CallbackPot', 'false', is_CallbackPot, help_string="If true, a callback potential (calls arbitrary passed by user)")
    call param_register(params, 'SocketPot', 'false', is_SocketPot, help_string="If true, a socket potential that communicates via TCP/IP sockets")
    call param_register(params, 'little_clusters', 'false', this%little_clusters, help_string="If true, uses little cluster, calculate forces only")
    call param_register(params, 'force_using_fd', 'F', this%force_using_fd, &
      help_string="If true, and if 'force' is also present in the calc argument list, calculate forces using finite difference.")
    call param_register(params, 'virial_using_fd', 'F', this%virial_using_fd, &
      help_string="If true, and if 'virial' is also present in the calc argument list, calculate virial using finite difference.")


    if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='Potential_Simple_Initialise_str args_str')) then
      call system_abort("Potential_Simple_Initialise_str failed to parse args_str='"//trim(args_str)//"'")
    endif
    call finalise(params)

    if (count( (/is_TB, is_IP, is_FilePot, is_wrapper, is_callbackpot, is_socketpot /) ) /= 1) then
      call system_abort("Potential_Simple_Initialise_str found too few or too many Potential_Simple types args_str='"//trim(args_str)//"'")
    endif

    if (is_IP) then
       allocate(this%ip)
       if(present(param_str)) then
          if (this%little_clusters) then
             call Initialise(this%ip, args_str, param_str, error=error)
          else
             call Initialise(this%ip, args_str, param_str, mpi_obj, error=error)
          end if
	  PASS_ERROR(error)
       else
          RAISE_ERROR('Potential_Simple_initialise: no param_str present during IP init', error)
       endif
#ifdef HAVE_TB
    else if (is_TB) then
       allocate(this%tb)
       if(present(param_str)) then
          if (this%little_clusters) then
             call Initialise(this%tb, args_str, param_str, error=error)
          else
             call Initialise(this%tb, args_str, param_str, mpi_obj=mpi_obj, error=error)
          end if
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
       if (this%little_clusters) then
          call Initialise(this%filepot, args_str, error=error)
       else
          call Initialise(this%filepot, args_str, mpi_obj, error=error)
       end if

       PASS_ERROR(error)
    else if (is_CallbackPot) then
       allocate(this%CallbackPot)
       if (this%little_clusters) then
          call Initialise(this%callbackpot, args_str, error=error)
       else
          call Initialise(this%callbackpot, args_str, mpi=mpi_obj, error=error)
       end if
       PASS_ERROR(error)
    else if (is_SocketPot) then
       allocate(this%SocketPot)
       if (this%little_clusters) then
          call Initialise(this%socketpot, args_str, error=error)
       else
          call Initialise(this%socketpot, args_str, mpi=mpi_obj, error=error)
       end if
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
    elseif(associated(this%socketpot)) then
       call Finalise(this%socketpot)
       deallocate(this%socketpot)
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
    elseif(associated(this%socketpot)) then
       Potential_Simple_cutoff = cutoff(this%socketpot)
    else
       Potential_Simple_cutoff = 0.0_dp
    end if
  end function Potential_Simple_cutoff

  recursive subroutine Potential_Simple_Calc(this, at, args_str, error)
    type(Potential_Simple), intent(inout) :: this
    type(Atoms), intent(inout) :: at     !% The atoms structure to compute energy and forces
    character(len=*), intent(in), optional :: args_str
    integer, intent(out), optional :: error

    real(dp) :: energy, virial(3,3), deform(3,3), lat_save(3,3)
    real(dp), allocatable :: allpos_save(:,:)
    real(dp), pointer :: at_force_ptr(:,:), at_local_energy_ptr(:), at_local_virial_ptr(:,:)

    integer:: i,j,k,n, zero_loc(1)
    real(dp):: e_plus, e_minus, pos_save, r_scale, E_scale, cluster_box_buffer
    type(Dictionary) :: params
    logical :: single_cluster, little_clusters, dummy, do_rescale_r, do_rescale_E, has_cluster_box_buffer
    integer :: partition_k_clusters
    real(dp) :: partition_nneightol
    logical :: partition_do_ripening
    character(len=STRING_LENGTH) :: my_args_str, cluster_args_str, new_args_str
    integer, pointer, dimension(:) :: hybrid_mark, cluster_index, termindex, modified_hybrid_mark
    real(dp), pointer, dimension(:) :: weight_region1, at_prop_ptr_r, cluster_prop_ptr_r
    integer, allocatable, dimension(:) :: hybrid_mark_saved
    logical, allocatable, dimension(:) :: cluster_mask
    real(dp), allocatable, dimension(:) :: weight_region1_saved
    real(dp), pointer, dimension(:,:) :: f_cluster
    logical :: do_carve_cluster
    type(Table) :: cluster_info, cut_bonds, t
    integer, pointer :: cut_bonds_p(:,:), old_cut_bonds_p(:,:)
    integer :: i_inner, i_outer, n_non_term
    type(Atoms) :: cluster
    character(len=STRING_LENGTH), target :: calc_force, calc_energy, calc_local_energy, calc_virial, calc_local_virial
    logical :: do_calc_force, do_calc_energy, do_calc_local_energy, do_calc_virial, do_calc_local_virial

    integer, pointer :: cluster_mark_p(:), at_prop_ptr_i(:), cluster_prop_ptr_i(:)
    integer, pointer :: old_cluster_mark_p(:)
    character(len=STRING_LENGTH) :: run_suffix
    logical :: force_using_fd, use_ridders
    real(dp) :: force_fd_delta
    logical :: virial_using_fd
    real(dp) :: virial_fd_delta
    real(dp) :: origin(3), extent(3,3), extent_inv(3,3), subregion_center(3)

    character(len=STRING_LENGTH) :: read_extra_param_list, read_extra_property_list
    character(STRING_LENGTH) :: tmp_params_array(100), copy_keys(100)
    integer :: n_copy, n_params, prop_j

    integer, parameter :: ridders_ntab = 10
    real(dp), parameter :: ridders_con = 1.4_dp
    real(dp), parameter :: ridders_con2 = ridders_con**2
    real(dp), parameter :: ridders_safe = 2.0_dp
    integer :: ridders_itab, ridders_jtab
    real(dp) :: ridders_hh, ridders_error, ridders_errort, ridders_fac, ridders_a(ridders_ntab,ridders_ntab)


    INIT_ERROR(error)

    if (at%N <= 0) then
      RAISE_ERROR("Potential_Simple_Calc called with at%N <= 0", error)
    endif

    if (present(args_str)) then
      call print('Potential_Simple_calc got args_str "'//trim(args_str)//'"', PRINT_VERBOSE)
      my_args_str = args_str
    else
      call print('Potential_Simple_calc got no args_str', PRINT_VERBOSE)
      my_args_str = ""
    endif

    call initialise(params)
    call param_register(params, 'single_cluster', 'F', single_cluster, &
      help_string="If true, calculate all active/transition atoms with a single big cluster")
    call param_register(params, 'run_suffix', '', run_suffix, &
      help_string="suffix to append to hybrid_mark field used")
    call param_register(params, 'carve_cluster', 'T', do_carve_cluster, &
      help_string="If true, calculate active region atoms by carving out a cluster")
    call param_register(params, 'little_clusters', 'F', little_clusters, &
      help_string="If true, calculate forces (only) by doing each atom separately surrounded by a little buffer cluster")
    call param_register(params, 'partition_k_clusters', '0', partition_k_clusters, &
      help_string="If given and K > 1, partition QM core region into K sub-clusters using partition_qm_list() routine")
    call param_register(params, 'partition_nneightol', '0.0_dp', partition_nneightol, &
      help_string="If given override at%nneightol used to generate connectivity matrix for partitioning")
    call param_register(params, 'partition_do_ripening', 'F', partition_do_ripening, &
      help_string="If true, enable digestive ripening to refine the METIS-generated K-way partitioning (not yet implemented)")
    call param_register(params, 'cluster_box_buffer', '0.0_dp', cluster_box_buffer, has_value_target=has_cluster_box_buffer, &
      help_string="If present, quickly cut out atoms in box within cluster_box_radius of any active atoms, rather than doing it properly")
    call param_register(params, 'r_scale', '1.0_dp', r_scale, has_value_target=do_rescale_r, &
      help_string="rescale calculated positions (and correspondingly forces) by this factor")
    call param_register(params, 'E_scale', '1.0_dp', E_scale, has_value_target=do_rescale_E, &
      help_string="rescale calculate energies (and correspondingly forces) by this factor")
    call param_register(params, 'use_ridders', 'F', use_ridders, &
      help_string="If true and using numerical derivatives, use the Ridders method.")
    call param_register(params, 'force_using_fd', 'F', force_using_fd, &
      help_string="If true, and if 'force' is also present in the argument list, calculate forces using finite difference.")
    call param_register(params, 'force_fd_delta', '1.0e-4', force_fd_delta, &
      help_string="Displacement to use with finite difference force calculation")
    call param_register(params, 'virial_using_fd', 'F', virial_using_fd, &
      help_string="If true, and if 'virial' is also present in the argument list, calculate virial using finite difference.")
    call param_register(params, 'virial_fd_delta', '1.0e-4', virial_fd_delta, &
      help_string="Displacement to use with finite difference virial calculation")
    call param_register(params, 'energy', '', calc_energy, &
      help_string="If present, calculate energy and put it in field with this string as name")
    call param_register(params, 'force', '', calc_force, &
      help_string="If present, calculate force and put it in field with this string as name")
    call param_register(params, 'local_energy', '', calc_local_energy, &
      help_string="If present, calculate local_energy and put it in field with this string as name")
    call param_register(params, 'virial', '', calc_virial, &
      help_string="If present, calculate virial and put it in field with this string as name")
    call param_register(params, 'local_virial', '', calc_local_virial, &
      help_string="If present, calculate local_virial and put it in field with this string as name")
    call param_register(params, "read_extra_param_list", '', read_extra_param_list, &
         help_string="if single_cluster=T and carve_cluster=T, extra params to copy back from cluster")
    call param_register(params, "read_extra_property_list", '', read_extra_property_list, &
         help_string="if single_cluster=T and carve_cluster=T, extra properties to copy back from cluster")

    if (.not. param_read_line(params, my_args_str, ignore_unknown=.true.,task='Potential_Simple_Calc_str args_str') ) then
      RAISE_ERROR("Potential_Simple_calc failed to parse args_str='"//trim(my_args_str)//"'", error)
    endif
    call finalise(params)

    if(this%force_using_fd) then
       force_using_fd = .true.
    end if
    if(this%virial_using_fd) then
       virial_using_fd = .true.
    end if

    if (count((/single_cluster, little_clusters, partition_k_clusters > 1/)) > 1) then
         RAISE_ERROR('Potential_Simple_calc: single_cluster and little_clusters options are mutually exclusive', error)
    endif

    if (len_trim(calc_force) > 0) call assign_property_pointer(at, trim(calc_force), at_force_ptr)
    if (len_trim(calc_local_energy) > 0) call assign_property_pointer(at, trim(calc_local_energy), at_local_energy_ptr)
    if (len_trim(calc_local_virial) > 0) call assign_property_pointer(at, trim(calc_local_virial), at_local_virial_ptr)

    if (little_clusters) then

       if (.not. this%little_clusters) then
          RAISE_ERROR('Potential_Simple_calc: little_clusters=T in calc() args_str but not in initialise() args_str.', error)
       end if

       if (len_trim(calc_energy) > 0 .or. len_trim(calc_local_energy) > 0 .or.  len_trim(calc_virial) > 0 .or. len_trim(calc_local_virial) > 0 .or. &
           len_trim(calc_force) <= 0) then
          RAISE_ERROR('Potential_Simple_calc: little_clusters option only supports calculation of forces, not energies, local energies or virials', error)
       endif

       ! modified args_str options for create_cluster_info() and carve_cluster()
       call initialise(params)
       call read_string(params, my_args_str)
       call set_value(params, 'run_suffix', run_suffix) ! ensure that run_suffix doesn't have value T_NONE
       if (do_rescale_r) call set_value(params, 'do_rescale_r', .true.)
       cluster_args_str = write_string(params, real_format='f16.8')
       call finalise(params)

       ! modified args_str options for potential calc() on clusters
       ! must remove "little_clusters" from args_str so that recursion terminates
       ! We rescale cluster explicity, so remove from args_str so that underlying potential doesn't do it
       call initialise(params)
       call read_string(params, my_args_str)
       call set_value(params, 'run_suffix', run_suffix) ! ensure that run_suffix doesn't have value T_NONE
       call remove_value(params, 'little_clusters')
       call remove_value(params, 'r_scale')
       call remove_value(params, 'E_scale')
       new_args_str = write_string(params, real_format='f16.8')
       call finalise(params)
 
       if (.not. assign_pointer(at, 'hybrid_mark'//trim(run_suffix), hybrid_mark)) then
            RAISE_ERROR('Potential_Simple_calc: cannot assign pointer to hybrid_mark property ', error)
       endif

       if (.not. any(hybrid_mark == HYBRID_ACTIVE_MARK)) then
          at_force_ptr = 0.0_dp
          return
       end if

       ! Save hybrid_mark and weight_region1, because we're going to overwrite them
       allocate(hybrid_mark_saved(at%N))
       hybrid_mark_saved = hybrid_mark

       if (has_property(at, 'weight_region1'//trim(run_suffix))) then
          dummy = assign_pointer(at, 'weight_region1'//trim(run_suffix), weight_region1)
          allocate(weight_region1_saved(at%N))
          weight_region1_saved = weight_region1
       end if

       ! call ourselves once for each active or transition atom 
       at_force_ptr = 0.0_dp
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
          call create_hybrid_weights(at, cluster_args_str)
          cluster_info = create_cluster_info_from_mark(at, cluster_args_str,mark_name='hybrid_mark'//trim(run_suffix),error=error)
	  PASS_ERROR_WITH_INFO("potential_calc: creating little cluster ="//i//" from hybrid_mark", error)
	  call carve_cluster(at, cluster_args_str, cluster_info, cluster)
	  call finalise(cluster_info)

          ! Reassign pointers - create_cluster_info_from_mark() might have broken them
          if (has_property(at, 'hybrid_mark'//trim(run_suffix))) &
               dummy = assign_pointer(at, 'hybrid_mark'//trim(run_suffix), hybrid_mark)
          if (has_property(at, 'weight_region1'//trim(run_suffix))) &
               dummy = assign_pointer(at, 'weight_region1'//trim(run_suffix), weight_region1)

	  if (current_verbosity() >= PRINT_NERD) then
	    call write(cluster, 'stdout', prefix='LITTLE_CLUSTER')
	  endif
          call print('ARGS0 | '//cluster_args_str,PRINT_VERBOSE)

          ! Disable MPI for duration of calc() call, since we're doing parallelisation at level of clusters
          call calc(this, cluster, args_str=new_args_str)
	  if (.not. assign_pointer(cluster, trim(calc_force), f_cluster)) then
	     RAISE_ERROR('Potential_Simple_calc: small_clusters failed to get a valid '//trim(calc_force)//' property in cluster from calc',error)
	  endif

          if (do_rescale_r)  f_cluster = f_cluster*r_scale
          if (do_rescale_E)  f_cluster = f_cluster*E_scale
          at_force_ptr(:,i) = f_cluster(:,1)
	   call finalise(cluster)
       end do

       if (this%mpi%active) then
          call sum_in_place(this%mpi, at_force_ptr)
       end if
       hybrid_mark = hybrid_mark_saved
       deallocate(hybrid_mark_saved)

       if (allocated(weight_region1_saved)) then
          weight_region1 = weight_region1_saved
          deallocate(weight_region1_saved)
       end if

    else if (partition_k_clusters > 1) then

       ! Split QM cluster into K pieces using METIS graph partitioning library call

       if (partition_nneightol == 0.0_dp) partition_nneightol = at%nneightol
       call partition_qm_list(at, partition_k_clusters, partition_nneightol, partition_do_ripening)
       call write(at, 'partitioned.xyz')

       ! ... FIXME need to loop over the sub-clusters and call ourselves recursively on each

    else if (single_cluster) then

       if (len_trim(calc_energy) > 0 .or. len_trim(calc_local_energy) > 0 .or.  len_trim(calc_virial) > 0 .or. len_trim(calc_local_virial) > 0 .or. &
           len_trim(calc_force) <= 0) then
          RAISE_ERROR('Potential_Simple_calc: single_cluster option only supports calculation of forces, not energies, local energies or virials', error)
       endif

       if (.not. assign_pointer(at, 'hybrid_mark'//trim(run_suffix), hybrid_mark)) then
            RAISE_ERROR('Potential_Simple_calc: single_cluster cannot assign pointer to hybrid_mark'//trim(run_suffix)//' property ', error)
       endif

       if (.not. any(hybrid_mark == HYBRID_ACTIVE_MARK)) then
          at_force_ptr = 0.0_dp
          return
       end if

       ! modified args_str options for create_cluster_info() and carve_cluster()
       call initialise(params)
       call read_string(params, my_args_str)
       call set_value(params, 'run_suffix', run_suffix) ! ensure that run_suffix doesn't have value T_NONE
       if (do_rescale_r) call set_value(params, 'do_rescale_r', .true.)
       cluster_args_str = write_string(params, real_format='f16.8')
       call finalise(params)

       ! modified args_str options for potential calc() on clusters
       ! must remove "single_cluster" from args_str so that recursion terminates
       ! We rescale cluster explicity, so remove from args_str so that underlying potential doesn't do it
       call initialise(params)
       call read_string(params, my_args_str)
       call set_value(params, 'run_suffix', run_suffix) ! ensure that run_suffix doesn't have value T_NONE
       call remove_value(params, 'single_cluster')
       call remove_value(params, 'r_scale')
       call remove_value(params, 'E_scale')
       new_args_str = write_string(params, real_format='f16.8')
       call finalise(params)

       ! call ourselves on a cluster formed from marked atoms
       call print('Potential_Simple_calc: constructing single_cluster', PRINT_VERBOSE)

       if (do_carve_cluster) then
	 call print('Potential_Simple_calc: carving cluster', PRINT_VERBOSE)

         if (has_cluster_box_buffer) then
            call print('Doing quick cluster carving using cluster box buffer '//cluster_box_buffer//' A')
            
            call estimate_origin_extent(at, hybrid_mark == HYBRID_ACTIVE_MARK, cluster_box_buffer, origin, extent)
            call matrix3x3_inverse(extent,extent_inv)
            subregion_center = origin + 0.5_dp*sum(extent,2)
            allocate(cluster_mask(at%n))
            do i=1, at%N
               cluster_mask(i) = is_in_subregion(at%pos(:,i), subregion_center, at%lattice, at%g, extent_inv)
            end do
            call select(cluster, at, mask=cluster_mask, orig_index=.true.)
            deallocate(cluster_mask)
            
            ! move orig_index property to index//run_suffix
            call assign_property_pointer(cluster, 'orig_index', cluster_index, error=error)
            PASS_ERROR(error)
            call add_property(cluster, 'index'//trim(run_suffix), cluster_index)
            call remove_property(cluster, 'orig_index', error=error)
            PASS_ERROR(error)
            ! there are no terminating atoms
            call add_property(cluster, 'termindex'//trim(run_suffix), 0)
         else
            cluster_info = create_cluster_info_from_mark(at, cluster_args_str, mark_name='hybrid_mark'//trim(run_suffix), error=error)
            PASS_ERROR_WITH_INFO("potential_calc: creating cluster info from hybrid_mark", error)

            ! Check there are no repeated indices among the non-termination atoms in the cluster
            n_non_term = count(cluster_info%int(6,1:cluster_info%n) == 0)
            t = int_subtable(cluster_info,(/ (i,i=1,n_non_term) /),(/1/))
            if (multiple_images(t)) then
               RAISE_ERROR('Potential_Simple_calc: single_cluster=T not yet implemented when cluster contains repeated periodic images', error)
            endif
            call finalise(t)
            
            call carve_cluster(at, cluster_args_str, cluster_info, cluster, mark_name='hybrid_mark'//trim(run_suffix), error=error)
            PASS_ERROR_WITH_INFO("potential_calc: carving cluster", error)
            call finalise(cluster_info)
         end if
         if (current_verbosity() >= PRINT_NERD) then
            call write(cluster, 'stdout', prefix='CLUSTER')
         endif
	 if (.not. assign_pointer(cluster, 'index'//trim(run_suffix), cluster_index)) then
	      RAISE_ERROR('Potential_Simple_calc: cluster is missing index property', error)
	 endif
	 if (.not. assign_pointer(cluster, 'termindex'//trim(run_suffix), termindex)) then
	      RAISE_ERROR('Potential_Simple_calc: cluster is missing termindex property', error)
	 endif
         call print('ARGS1 | '//cluster_args_str,PRINT_VERBOSE)

	 call calc(this, cluster, args_str=new_args_str, error=error)
	 PASS_ERROR_WITH_INFO('potential_calc after calc in carve_cluster', error)

	 call assign_property_pointer(cluster, trim(calc_force), f_cluster, error=error)
	 PASS_ERROR_WITH_INFO('Potential_Simple_calc: single_cluster failed to get a valid '//trim(calc_force)//' property in cluster from calc',error)
	 if (do_rescale_r)  f_cluster = f_cluster*r_scale
         if (do_rescale_E)  f_cluster = f_cluster*E_scale

         ! Reassign pointers - create_cluster_info_from_mark() might have broken them
         if (has_property(at, 'hybrid_mark'//trim(run_suffix))) &
              dummy = assign_pointer(at, 'hybrid_mark'//trim(run_suffix), hybrid_mark)

	 ! copy forces for all active and transition atoms
	 at_force_ptr = 0.0_dp
	 do i=1,cluster%N
	    if (termindex(i) /= 0) cycle ! skip termination atoms
	    if (hybrid_mark(cluster_index(i)) == HYBRID_ACTIVE_MARK .or. &
		hybrid_mark(cluster_index(i)) == HYBRID_TRANS_MARK) &
		  at_force_ptr(:,cluster_index(i)) = f_cluster(:,i)
	 end do
         nullify(f_cluster)

         ! might need to copy some params from cluster to main Atoms object, for compatibility with filepot
         if (len_trim(read_extra_param_list) > 0) then
            call parse_string(read_extra_param_list, ':', tmp_params_array, n_params, error=error)
            PASS_ERROR(error)

            n_copy = 0
            do i=1,n_params
               if (has_key(cluster%params, trim(tmp_params_array(i)))) then
                  n_copy = n_copy + 1
                  copy_keys(n_copy) = tmp_params_array(i)
                  call print("Potential_simple calc() copying param key "//trim(copy_keys(n_copy)), PRINT_VERBOSE)
               else if  (has_key(cluster%params, trim(tmp_params_array(i))//trim(run_suffix))) then
                  n_copy = n_copy + 1
                  copy_keys(n_copy) =  trim(tmp_params_array(i))//trim(run_suffix)
                  call print("Potential_simple calc() copying param key "//trim(copy_keys(n_copy)), PRINT_VERBOSE)
               end if
            end do
            call subset(cluster%params, copy_keys(1:n_copy), at%params, out_no_initialise=.true.)
         end if

         ! do the same for any extra properties to be copied back
         if (len_trim(read_extra_property_list) > 0) then
            call parse_string(read_extra_property_list, ':', tmp_params_array, n_params, error=error)
            PASS_ERROR(error)

            n_copy = 0
            do i=1,n_params
               if (has_key(cluster%properties, trim(tmp_params_array(i)))) then
                  n_copy = n_copy + 1
                  copy_keys(n_copy) = tmp_params_array(i)
                  call print("Potential_simple calc() copying property key "//trim(copy_keys(n_copy)), PRINT_VERBOSE)
               else if  (has_key(cluster%properties, trim(tmp_params_array(i))//trim(run_suffix))) then
                  n_copy = n_copy + 1
                  copy_keys(n_copy) =  trim(tmp_params_array(i))//trim(run_suffix)
                  call print("Potential_simple calc() copying property key "//trim(copy_keys(n_copy)), PRINT_VERBOSE)
               end if
            end do

            do j=1, n_copy
               prop_j = lookup_entry_i(cluster%properties, copy_keys(j))
               if (cluster%properties%entries(prop_j)%type == T_INTEGER_A) then
                  if (.not. has_property(at, copy_keys(j))) call add_property(at, copy_keys(j), 0)
                  call assign_property_pointer(at, copy_keys(j), at_prop_ptr_i, error=error)
                  PASS_ERROR(error)
                  call assign_property_pointer(cluster, copy_keys(j), cluster_prop_ptr_i, error=error)
                  PASS_ERROR(error)
                  do i=1,cluster%N
                     at_prop_ptr_i(cluster_index(i)) = cluster_prop_ptr_i(i)
                  end do
               else if (cluster%properties%entries(prop_j)%type == T_REAL_A) then
                  if (.not. has_property(at, copy_keys(j))) call add_property(at, copy_keys(j), 0.0_dp)
                  call assign_property_pointer(at, copy_keys(j), at_prop_ptr_r, error=error)
                  PASS_ERROR(error)
                  call assign_property_pointer(cluster, copy_keys(j), cluster_prop_ptr_r, error=error)
                  PASS_ERROR(error)
                  do i=1,cluster%N
                     at_prop_ptr_r(cluster_index(i)) = cluster_prop_ptr_r(i)
                  end do
               else
                  RAISE_ERROR('unsupported property type '//cluster%properties%entries(prop_j)%type, error)
               end if
            end do
         end if

	 call finalise(cluster)
       else ! not do_carve_cluster
	 call print('Potential_Simple_calc: not carving cluster', PRINT_VERBOSE)
	 cluster_info = create_cluster_info_from_mark(at, trim(cluster_args_str) // " cluster_same_lattice", cut_bonds, mark_name='hybrid_mark'//trim(run_suffix), error=error)
	 PASS_ERROR_WITH_INFO('potential_calc creating cluster info from hybrid mark with carve_cluster=F', error)

	 call add_property(at, 'cluster_mark'//trim(run_suffix), HYBRID_NO_MARK)
	 call add_property(at, 'old_cluster_mark'//trim(run_suffix), HYBRID_NO_MARK)
	 if (.not. assign_pointer(at, 'cluster_mark'//trim(run_suffix), cluster_mark_p)) then
	   RAISE_ERROR("Potential_Simple_calc failed to assing pointer for cluster_mark"//trim(run_suffix)//" pointer", error)
	 endif
	 if (.not. assign_pointer(at, 'old_cluster_mark'//trim(run_suffix), old_cluster_mark_p)) then
	   RAISE_ERROR("Potential_Simple_calc failed to assing pointer for old_cluster_mark"//trim(run_suffix)//" pointer", error)
	 endif
	 old_cluster_mark_p = cluster_mark_p
	 cluster_mark_p = HYBRID_NO_MARK
	 if (.not. assign_pointer(at, 'modified_hybrid_mark'//trim(run_suffix), modified_hybrid_mark)) then
	   cluster_mark_p = hybrid_mark
	 else
	   cluster_mark_p = modified_hybrid_mark
	 endif

         !save cut bonds in cut_bonds property
         call add_property(at, 'old_cut_bonds'//trim(run_suffix), 0, n_cols=MAX_CUT_BONDS)
	 if (.not. assign_pointer(at, 'old_cut_bonds'//trim(run_suffix), old_cut_bonds_p)) then
	   RAISE_ERROR("Potential_Simple_calc failed to assing pointer for cut_bonds pointer", error)
	 endif
	 call add_property(at, 'cut_bonds'//trim(run_suffix), 0, n_cols=MAX_CUT_BONDS)
	 if (.not. assign_pointer(at, 'cut_bonds'//trim(run_suffix), cut_bonds_p)) then
	   RAISE_ERROR("Potential_Simple_calc failed to assing pointer for cut_bonds pointer", error)
	 endif
         old_cut_bonds_p = cut_bonds_p
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
	 call calc(this, at, args_str=new_args_str, error=error)
	 PASS_ERROR_WITH_INFO('potential_calc after calc with carve_cluster=F', error)
	 if (do_rescale_r)  at_force_ptr = at_force_ptr*r_scale
	 if (do_rescale_E)  at_force_ptr = at_force_ptr*E_scale
       endif ! do_carve_cluster
    else ! little_clusters and single_cluster are false..

       ! For IP, call setup_atoms() hook now in case any properties must be added.
       ! This must be done *before* we assign pointers to force, local_e etc.
       if (associated(this%ip)) then
          call setup_atoms(this%ip, at)
       end if

       do_calc_force = (len_trim(calc_force) > 0) .and. .not. force_using_fd
       do_calc_energy = len_trim(calc_energy) > 0
       do_calc_local_energy = len_trim(calc_local_energy) > 0
       do_calc_virial = (len_trim(calc_virial) > 0) .and. .not. virial_using_fd
       do_calc_local_virial = len_trim(calc_local_virial) > 0 
       call print("do_calc_force=        "//do_calc_force//" calc_force="//trim(calc_force), PRINT_VERBOSE)
       call print("force_using_fd=       "//force_using_fd, PRINT_VERBOSE)
       call print("virial_using_fd=      "//virial_using_fd, PRINT_VERBOSE)
       call print("do_calc_energy=       "//do_calc_energy//" calc_energy="//trim(calc_energy), PRINT_VERBOSE)
       call print("do_calc_local_energy= "//do_calc_local_energy//" calc_local_energy="//trim(calc_local_energy), PRINT_VERBOSE)
       call print("do_calc_virial=       "//do_calc_virial//" calc_virial="//trim(calc_virial), PRINT_VERBOSE)
       call print("do_calc_local_virial= "//do_calc_local_virial//" calc_local_virial="//trim(calc_local_virial), PRINT_VERBOSE)

       if(do_calc_virial .or. do_calc_energy .or. do_calc_force .or. do_calc_local_energy .or. do_calc_local_virial) then
          if(associated(this%ip)) then

             if (do_calc_local_virial) then
                if (do_calc_virial) then
                   if (.not. do_calc_energy .and. .not. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%ip, at, virial=virial, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. .not. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%ip, at, f=at_force_ptr, virial=virial, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%ip, at, local_e=at_local_energy_ptr, virial=virial, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%ip, at, local_e=at_local_energy_ptr, f=at_force_ptr, virial=virial, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (do_calc_energy .and. .not. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%ip, at, energy=energy, virial=virial, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (do_calc_energy .and. .not. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%ip, at, energy=energy, f=at_force_ptr, virial=virial, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (do_calc_energy .and. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%ip, at, energy=energy, local_e=at_local_energy_ptr, virial=virial, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (do_calc_energy .and. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%ip, at, energy=energy, local_e=at_local_energy_ptr, f=at_force_ptr, virial=virial, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   end if
                else
                   if (.not. do_calc_energy .and. .not. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%ip, at, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. .not. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%ip, at, f=at_force_ptr, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%ip, at, local_e=at_local_energy_ptr, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%ip, at, local_e=at_local_energy_ptr, f=at_force_ptr, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (do_calc_energy .and. .not. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%ip, at, energy=energy, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (do_calc_energy .and. .not. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%ip, at, energy=energy, f=at_force_ptr, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (do_calc_energy .and. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%ip, at, energy=energy, local_e=at_local_energy_ptr, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (do_calc_energy .and. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%ip, at, energy=energy, local_e=at_local_energy_ptr, f=at_force_ptr, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   end if
                end if
             else ! no local virials
                if (do_calc_virial) then
                   if (.not. do_calc_energy .and. .not. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%ip, at, virial=virial, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. .not. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%ip, at, f=at_force_ptr, virial=virial, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%ip, at, local_e=at_local_energy_ptr, virial=virial, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%ip, at, local_e=at_local_energy_ptr, f=at_force_ptr, virial=virial, args_str=args_str, error=error)
                   else if (do_calc_energy .and. .not. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%ip, at, energy=energy, virial=virial, args_str=args_str, error=error)
                   else if (do_calc_energy .and. .not. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%ip, at, energy=energy, f=at_force_ptr, virial=virial, args_str=args_str, error=error)
                   else if (do_calc_energy .and. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%ip, at, energy=energy, local_e=at_local_energy_ptr, virial=virial, args_str=args_str, error=error)
                   else if (do_calc_energy .and. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%ip, at, energy=energy, local_e=at_local_energy_ptr, f=at_force_ptr, virial=virial, args_str=args_str, error=error)
                   end if
                else
                   if (.not. do_calc_energy .and. .not. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%ip, at, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. .not. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%ip, at, f=at_force_ptr, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%ip, at, local_e=at_local_energy_ptr, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%ip, at, local_e=at_local_energy_ptr, f=at_force_ptr, args_str=args_str, error=error)
                   else if (do_calc_energy .and. .not. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%ip, at, energy=energy, args_str=args_str, error=error)
                   else if (do_calc_energy .and. .not. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%ip, at, energy=energy, f=at_force_ptr, args_str=args_str, error=error)
                   else if (do_calc_energy .and. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%ip, at, energy=energy, local_e=at_local_energy_ptr, args_str=args_str, error=error)
                   else if (do_calc_energy .and. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%ip, at, energy=energy, local_e=at_local_energy_ptr, f=at_force_ptr, args_str=args_str, error=error)
                   end if
                end if
             endif
#ifdef HAVE_TB
          else if(associated(this%tb)) then

             if (do_calc_local_virial) then
                if (do_calc_virial) then
                   if (.not. do_calc_energy .and. .not. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%tb, at, virial=virial, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. .not. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%tb, at, forces=at_force_ptr, virial=virial, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%tb, at, local_e=at_local_energy_ptr, virial=virial, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%tb, at, local_e=at_local_energy_ptr, forces=at_force_ptr, virial=virial, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (do_calc_energy .and. .not. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%tb, at, energy=energy, virial=virial, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (do_calc_energy .and. .not. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%tb, at, energy=energy, forces=at_force_ptr, virial=virial, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (do_calc_energy .and. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%tb, at, energy=energy, local_e=at_local_energy_ptr, virial=virial, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (do_calc_energy .and. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%tb, at, energy=energy, local_e=at_local_energy_ptr, forces=at_force_ptr, virial=virial, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   end if
                else
                   if (.not. do_calc_energy .and. .not. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%tb, at, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. .not. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%tb, at, forces=at_force_ptr, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%tb, at, local_e=at_local_energy_ptr, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%tb, at, local_e=at_local_energy_ptr, forces=at_force_ptr, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (do_calc_energy .and. .not. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%tb, at, energy=energy, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (do_calc_energy .and. .not. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%tb, at, energy=energy, forces=at_force_ptr, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (do_calc_energy .and. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%tb, at, energy=energy, local_e=at_local_energy_ptr, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (do_calc_energy .and. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%tb, at, energy=energy, local_e=at_local_energy_ptr, forces=at_force_ptr, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   end if                
                end if
             else ! no local virials
                if (do_calc_virial) then
                   if (.not. do_calc_energy .and. .not. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%tb, at, virial=virial, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. .not. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%tb, at, forces=at_force_ptr, virial=virial, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%tb, at, local_e=at_local_energy_ptr, virial=virial, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%tb, at, local_e=at_local_energy_ptr, forces=at_force_ptr, virial=virial, args_str=args_str, error=error)
                   else if (do_calc_energy .and. .not. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%tb, at, energy=energy, virial=virial, args_str=args_str, error=error)
                   else if (do_calc_energy .and. .not. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%tb, at, energy=energy, forces=at_force_ptr, virial=virial, args_str=args_str, error=error)
                   else if (do_calc_energy .and. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%tb, at, energy=energy, local_e=at_local_energy_ptr, virial=virial, args_str=args_str, error=error)
                   else if (do_calc_energy .and. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%tb, at, energy=energy, local_e=at_local_energy_ptr, forces=at_force_ptr, virial=virial, args_str=args_str, error=error)
                   end if
                else
                   if (.not. do_calc_energy .and. .not. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%tb, at, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. .not. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%tb, at, forces=at_force_ptr, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%tb, at, local_e=at_local_energy_ptr, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%tb, at, local_e=at_local_energy_ptr, forces=at_force_ptr, args_str=args_str, error=error)
                   else if (do_calc_energy .and. .not. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%tb, at, energy=energy, args_str=args_str, error=error)
                   else if (do_calc_energy .and. .not. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%tb, at, energy=energy, forces=at_force_ptr, args_str=args_str, error=error)
                   else if (do_calc_energy .and. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%tb, at, energy=energy, local_e=at_local_energy_ptr, args_str=args_str, error=error)
                   else if (do_calc_energy .and. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%tb, at, energy=energy, local_e=at_local_energy_ptr, forces=at_force_ptr, args_str=args_str, error=error)
                   end if                
                end if
             endif

#endif
          elseif(associated(this%filepot)) then

             if (do_calc_local_virial) then
                if (do_calc_virial) then
                   if (.not. do_calc_energy .and. .not. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%filepot, at, virial=virial, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. .not. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%filepot, at, forces=at_force_ptr, virial=virial, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%filepot, at, local_e=at_local_energy_ptr, virial=virial, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%filepot, at, local_e=at_local_energy_ptr, forces=at_force_ptr, virial=virial, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (do_calc_energy .and. .not. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%filepot, at, energy=energy, virial=virial, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (do_calc_energy .and. .not. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%filepot, at, energy=energy, forces=at_force_ptr, virial=virial, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (do_calc_energy .and. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%filepot, at, energy=energy, local_e=at_local_energy_ptr, virial=virial, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (do_calc_energy .and. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%filepot, at, energy=energy, local_e=at_local_energy_ptr, forces=at_force_ptr, virial=virial, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   end if
                else
                   if (.not. do_calc_energy .and. .not. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%filepot, at, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. .not. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%filepot, at, forces=at_force_ptr, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%filepot, at, local_e=at_local_energy_ptr, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%filepot, at, local_e=at_local_energy_ptr, forces=at_force_ptr, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (do_calc_energy .and. .not. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%filepot, at, energy=energy, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (do_calc_energy .and. .not. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%filepot, at, energy=energy, forces=at_force_ptr, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (do_calc_energy .and. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%filepot, at, energy=energy, local_e=at_local_energy_ptr, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (do_calc_energy .and. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%filepot, at, energy=energy, local_e=at_local_energy_ptr, forces=at_force_ptr, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   end if                
                end if
             else ! no local virials
                if (do_calc_virial) then
                   if (.not. do_calc_energy .and. .not. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%filepot, at, virial=virial, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. .not. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%filepot, at, forces=at_force_ptr, virial=virial, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%filepot, at, local_e=at_local_energy_ptr, virial=virial, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%filepot, at, local_e=at_local_energy_ptr, forces=at_force_ptr, virial=virial, args_str=args_str, error=error)
                   else if (do_calc_energy .and. .not. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%filepot, at, energy=energy, virial=virial, args_str=args_str, error=error)
                   else if (do_calc_energy .and. .not. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%filepot, at, energy=energy, forces=at_force_ptr, virial=virial, args_str=args_str, error=error)
                   else if (do_calc_energy .and. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%filepot, at, energy=energy, local_e=at_local_energy_ptr, virial=virial, args_str=args_str, error=error)
                   else if (do_calc_energy .and. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%filepot, at, energy=energy, local_e=at_local_energy_ptr, forces=at_force_ptr, virial=virial, args_str=args_str, error=error)
                   end if
                else
                   if (.not. do_calc_energy .and. .not. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%filepot, at, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. .not. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%filepot, at, forces=at_force_ptr, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%filepot, at, local_e=at_local_energy_ptr, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%filepot, at, local_e=at_local_energy_ptr, forces=at_force_ptr, args_str=args_str, error=error)
                   else if (do_calc_energy .and. .not. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%filepot, at, energy=energy, args_str=args_str, error=error)
                   else if (do_calc_energy .and. .not. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%filepot, at, energy=energy, forces=at_force_ptr, args_str=args_str, error=error)
                   else if (do_calc_energy .and. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%filepot, at, energy=energy, local_e=at_local_energy_ptr, args_str=args_str, error=error)
                   else if (do_calc_energy .and. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%filepot, at, energy=energy, local_e=at_local_energy_ptr, forces=at_force_ptr, args_str=args_str, error=error)
                   end if                
                end if
             endif

          elseif(associated(this%callbackpot)) then

             if (do_calc_local_virial) then
                if (do_calc_virial) then
                   if (.not. do_calc_energy .and. .not. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%callbackpot, at, virial=virial, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. .not. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%callbackpot, at, forces=at_force_ptr, virial=virial, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%callbackpot, at, local_e=at_local_energy_ptr, virial=virial, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%callbackpot, at, local_e=at_local_energy_ptr, forces=at_force_ptr, virial=virial, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (do_calc_energy .and. .not. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%callbackpot, at, energy=energy, virial=virial, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (do_calc_energy .and. .not. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%callbackpot, at, energy=energy, forces=at_force_ptr, virial=virial, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (do_calc_energy .and. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%callbackpot, at, energy=energy, local_e=at_local_energy_ptr, virial=virial, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (do_calc_energy .and. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%callbackpot, at, energy=energy, local_e=at_local_energy_ptr, forces=at_force_ptr, virial=virial, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   end if
                else
                   if (.not. do_calc_energy .and. .not. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%callbackpot, at, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. .not. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%callbackpot, at, forces=at_force_ptr, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%callbackpot, at, local_e=at_local_energy_ptr, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%callbackpot, at, local_e=at_local_energy_ptr, forces=at_force_ptr, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (do_calc_energy .and. .not. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%callbackpot, at, energy=energy, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (do_calc_energy .and. .not. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%callbackpot, at, energy=energy, forces=at_force_ptr, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (do_calc_energy .and. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%callbackpot, at, energy=energy, local_e=at_local_energy_ptr, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (do_calc_energy .and. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%callbackpot, at, energy=energy, local_e=at_local_energy_ptr, forces=at_force_ptr, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   end if                
                end if
             else ! no local virials
                if (do_calc_virial) then
                   if (.not. do_calc_energy .and. .not. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%callbackpot, at, virial=virial, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. .not. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%callbackpot, at, forces=at_force_ptr, virial=virial, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%callbackpot, at, local_e=at_local_energy_ptr, virial=virial, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%callbackpot, at, local_e=at_local_energy_ptr, forces=at_force_ptr, virial=virial, args_str=args_str, error=error)
                   else if (do_calc_energy .and. .not. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%callbackpot, at, energy=energy, virial=virial, args_str=args_str, error=error)
                   else if (do_calc_energy .and. .not. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%callbackpot, at, energy=energy, forces=at_force_ptr, virial=virial, args_str=args_str, error=error)
                   else if (do_calc_energy .and. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%callbackpot, at, energy=energy, local_e=at_local_energy_ptr, virial=virial, args_str=args_str, error=error)
                   else if (do_calc_energy .and. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%callbackpot, at, energy=energy, local_e=at_local_energy_ptr, forces=at_force_ptr, virial=virial, args_str=args_str, error=error)
                   end if
                else
                   if (.not. do_calc_energy .and. .not. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%callbackpot, at, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. .not. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%callbackpot, at, forces=at_force_ptr, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%callbackpot, at, local_e=at_local_energy_ptr, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%callbackpot, at, local_e=at_local_energy_ptr, forces=at_force_ptr, args_str=args_str, error=error)
                   else if (do_calc_energy .and. .not. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%callbackpot, at, energy=energy, args_str=args_str, error=error)
                   else if (do_calc_energy .and. .not. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%callbackpot, at, energy=energy, forces=at_force_ptr, args_str=args_str, error=error)
                   else if (do_calc_energy .and. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%callbackpot, at, energy=energy, local_e=at_local_energy_ptr, args_str=args_str, error=error)
                   else if (do_calc_energy .and. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%callbackpot, at, energy=energy, local_e=at_local_energy_ptr, forces=at_force_ptr, args_str=args_str, error=error)
                   end if                
                end if
             endif


          elseif(associated(this%socketpot)) then

             if (do_calc_local_virial) then
                if (do_calc_virial) then
                   if (.not. do_calc_energy .and. .not. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%socketpot, at, virial=virial, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. .not. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%socketpot, at, forces=at_force_ptr, virial=virial, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%socketpot, at, local_e=at_local_energy_ptr, virial=virial, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%socketpot, at, local_e=at_local_energy_ptr, forces=at_force_ptr, virial=virial, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (do_calc_energy .and. .not. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%socketpot, at, energy=energy, virial=virial, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (do_calc_energy .and. .not. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%socketpot, at, energy=energy, forces=at_force_ptr, virial=virial, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (do_calc_energy .and. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%socketpot, at, energy=energy, local_e=at_local_energy_ptr, virial=virial, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (do_calc_energy .and. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%socketpot, at, energy=energy, local_e=at_local_energy_ptr, forces=at_force_ptr, virial=virial, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   end if
                else
                   if (.not. do_calc_energy .and. .not. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%socketpot, at, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. .not. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%socketpot, at, forces=at_force_ptr, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%socketpot, at, local_e=at_local_energy_ptr, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%socketpot, at, local_e=at_local_energy_ptr, forces=at_force_ptr, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (do_calc_energy .and. .not. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%socketpot, at, energy=energy, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (do_calc_energy .and. .not. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%socketpot, at, energy=energy, forces=at_force_ptr, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (do_calc_energy .and. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%socketpot, at, energy=energy, local_e=at_local_energy_ptr, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   else if (do_calc_energy .and. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%socketpot, at, energy=energy, local_e=at_local_energy_ptr, forces=at_force_ptr, local_virial = at_local_virial_ptr, args_str=args_str, error=error)
                   end if                
                end if
             else ! no local virials
                if (do_calc_virial) then
                   if (.not. do_calc_energy .and. .not. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%socketpot, at, virial=virial, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. .not. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%socketpot, at, forces=at_force_ptr, virial=virial, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%socketpot, at, local_e=at_local_energy_ptr, virial=virial, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%socketpot, at, local_e=at_local_energy_ptr, forces=at_force_ptr, virial=virial, args_str=args_str, error=error)
                   else if (do_calc_energy .and. .not. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%socketpot, at, energy=energy, virial=virial, args_str=args_str, error=error)
                   else if (do_calc_energy .and. .not. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%socketpot, at, energy=energy, forces=at_force_ptr, virial=virial, args_str=args_str, error=error)
                   else if (do_calc_energy .and. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%socketpot, at, energy=energy, local_e=at_local_energy_ptr, virial=virial, args_str=args_str, error=error)
                   else if (do_calc_energy .and. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%socketpot, at, energy=energy, local_e=at_local_energy_ptr, forces=at_force_ptr, virial=virial, args_str=args_str, error=error)
                   end if
                else
                   if (.not. do_calc_energy .and. .not. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%socketpot, at, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. .not. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%socketpot, at, forces=at_force_ptr, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%socketpot, at, local_e=at_local_energy_ptr, args_str=args_str, error=error)
                   else if (.not. do_calc_energy .and. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%socketpot, at, local_e=at_local_energy_ptr, forces=at_force_ptr, args_str=args_str, error=error)
                   else if (do_calc_energy .and. .not. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%socketpot, at, energy=energy, args_str=args_str, error=error)
                   else if (do_calc_energy .and. .not. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%socketpot, at, energy=energy, forces=at_force_ptr, args_str=args_str, error=error)
                   else if (do_calc_energy .and. do_calc_local_energy .and. .not. do_calc_force) then
                      call Calc(this%socketpot, at, energy=energy, local_e=at_local_energy_ptr, args_str=args_str, error=error)
                   else if (do_calc_energy .and. do_calc_local_energy .and. do_calc_force) then
                      call Calc(this%socketpot, at, energy=energy, local_e=at_local_energy_ptr, forces=at_force_ptr, args_str=args_str, error=error)
                   end if                
                end if
             endif

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

          if (len_trim(calc_energy) /= 0) then
             call set_value(at%params, trim(calc_energy), energy)
             call print("Setting parameter `"//trim(calc_energy)//"' to "//energy, PRINT_VERBOSE)
          end if
          if (len_trim(calc_virial) /= 0) then
             call set_value(at%params, trim(calc_virial), virial)
          end if

       end if

       if (force_using_fd .and. len_trim(calc_force) > 0) then ! do forces by finite difference

          call print("Calculating force by finite differences with displacement="//force_fd_delta, PRINT_VERBOSE)

          ! must remove 'force_using_fd' from args_str if it's there (and the rest, or properties get overwritten)
          call initialise(params)
          call read_string(params, my_args_str)
          call set_value(params, 'run_suffix', run_suffix) ! ensure that run_suffix doesn't have value T_NONE
          call remove_value(params, 'force_using_fd')
          call remove_value(params, 'force')
          call remove_value(params, 'energy')
          call remove_value(params, 'local_energy')
          call remove_value(params, 'virial')
	  call set_value(params, "energy", "fd_energy")
          new_args_str = write_string(params, real_format='f16.8')
          call finalise(params)
          call print("calc_force="//trim(calc_force), PRINT_ANAL)
          call print("New_args_str: "//new_args_str, PRINT_ANAL)


          if(use_ridders) then
             do i = 1, at%N
                do k = 1, 3
                   pos_save = at%pos(k,i)
                   ridders_hh = force_fd_delta

                   at%pos(k,i) = pos_save + ridders_hh
                   if(at%cutoff > 0) call calc_dists(at)
                   call calc(this, at, args_str=new_args_str, error=error)
                   call get_param_value(at, "fd_energy", e_plus, error=error)
                   PASS_ERROR_WITH_INFO("Potential_Simple_calc doing fd forces failed to get energy property fd_energy", error)

                   at%pos(k,i) = pos_save - ridders_hh
                   if(at%cutoff > 0) call calc_dists(at)
                   call calc(this, at, args_str=new_args_str, error=error)
                   call get_param_value(at, "fd_energy", e_minus, error=error)
                   PASS_ERROR_WITH_INFO("Potential_Simple_calc doing fd forces failed to get energy property fd_energy", error)

                   ridders_a(1,1) = (e_plus - e_minus) / (2.0_dp*ridders_hh)
                   ridders_error = huge(1.0_dp)

                   do ridders_itab = 2, ridders_ntab
                      ridders_hh = ridders_hh / ridders_con

                      at%pos(k,i) = pos_save + ridders_hh
                      if(at%cutoff > 0) call calc_dists(at)
                      call calc(this, at, args_str=new_args_str, error=error)
                      call get_param_value(at, "fd_energy", e_plus, error=error)
                      PASS_ERROR_WITH_INFO("Potential_Simple_calc doing fd forces failed to get energy property fd_energy", error)

                      at%pos(k,i) = pos_save - ridders_hh
                      if(at%cutoff > 0) call calc_dists(at)
                      call calc(this, at, args_str=new_args_str, error=error)
                      call get_param_value(at, "fd_energy", e_minus, error=error)
                      PASS_ERROR_WITH_INFO("Potential_Simple_calc doing fd forces failed to get energy property fd_energy", error)

                      ridders_a(1,ridders_itab) = (e_plus - e_minus) / (2.0_dp*ridders_hh)
                      ridders_fac = ridders_con2

                      do ridders_jtab = 2, ridders_itab
                         ridders_a(ridders_jtab,ridders_itab) = (ridders_a(ridders_jtab-1,ridders_itab)*ridders_fac - &
                            ridders_a(ridders_jtab-1,ridders_itab-1))/(ridders_fac-1.0_dp)
                         ridders_fac = ridders_con2*ridders_fac
                         ridders_errort = max(abs(ridders_a(ridders_jtab,ridders_itab)-ridders_a(ridders_jtab-1,ridders_itab)), &
                            abs(ridders_a(ridders_jtab,ridders_itab)-ridders_a(ridders_jtab-1,ridders_itab-1)))

                         if( ridders_errort <= ridders_error ) then
                            ridders_error = ridders_errort
                            at_force_ptr(k,i) = -ridders_a(ridders_jtab,ridders_itab)
                         endif
                      enddo ! jtab

                      if( abs(ridders_a(ridders_itab,ridders_itab)-ridders_a(ridders_itab-1,ridders_itab-1) ) >= ridders_SAFE*ridders_error) exit
                   enddo ! itab
                   at%pos(k,i) = pos_save
                   !error_force(alpha,i) = abs(error)
                enddo ! k
             enddo ! i
          else
             do i=1,at%N
                do k=1,3
                   pos_save = at%pos(k,i)
                   at%pos(k,i) = pos_save + force_fd_delta
                   if(at%cutoff > 0) call calc_dists(at)
                   call calc(this, at, args_str=new_args_str, error=error)
                   call get_param_value(at, "fd_energy", e_plus, error=error)
                   PASS_ERROR_WITH_INFO("Potential_Simple_calc doing fd forces failed to get energy property fd_energy", error)
                   at%pos(k,i) = pos_save - force_fd_delta
                   if(at%cutoff > 0) call calc_dists(at)
                   call calc(this, at, args_str=new_args_str, error=error)
                   call get_param_value(at, "fd_energy", e_minus, error=error)
                   PASS_ERROR_WITH_INFO("Potential_Simple_calc doing fd forces failed to get energy property fd_energy", error)
                   at%pos(k,i) = pos_save
                   if(at%cutoff > 0) call calc_dists(at)
                   at_force_ptr(k,i) = (e_minus-e_plus)/(2.0_dp*force_fd_delta) ! force is -ve gradient
                end do
             end do
          endif

          call remove_value(at%params, "fd_energy")
          call print("Done with finite difference force calculation", PRINT_VERBOSE)
       end if
       if(virial_using_fd .and. len_trim(calc_virial) > 0) then ! do virial by finite difference
          call print("Calculating virial by finite differences with displacement="//virial_fd_delta, PRINT_VERBOSE)

          ! must remove 'virial_using_fd' from args_str if it's there (and the rest, or properties get overwritten)
          call initialise(params)
          call read_string(params, my_args_str)
          call set_value(params, 'run_suffix', run_suffix) ! ensure that run_suffix doesn't have value T_NONE
          call remove_value(params, 'virial_using_fd')
          call remove_value(params, 'force')
          call remove_value(params, 'energy')
          call remove_value(params, 'local_energy')
          call remove_value(params, 'virial')
          call set_value(params, "energy", "fd_energy")
          new_args_str = write_string(params, real_format='f16.8')
          call finalise(params)

          call print("Virial finite difference calculation: new_args_str=`"//new_args_str//"'", PRINT_VERBOSE)
          
          allocate(allpos_save(3,at%N))
          allpos_save = at%pos
          lat_save = at%lattice
          do i=1,3
             do j=i,3
                
                deform = 0.0_dp
                call add_identity(deform)
                deform(i,j) = deform(i,j) + virial_fd_delta
                if(i /= j) deform(j,i) = deform(j,i) + virial_fd_delta
                call set_lattice(at, matmul(deform,lat_save), scale_positions=.true.)
                if(at%cutoff > 0) call calc_dists(at)
!                call verbosity_push(PRINT_ANAL)
!                call print(at)
!                call verbosity_pop()
                 call calc(this, at, args_str=new_args_str, error=error)
!                call verbosity_push(PRINT_ANAL)
!                call print(at)
!                call verbosity_pop()
                call get_param_value(at, "fd_energy", e_plus, error=error)
                PASS_ERROR_WITH_INFO("Potential_Simple_calc doing fd virial failed to get energy property fd_energy", error)
                
                deform = 0.0_dp
                call add_identity(deform)
                deform(i,j) = deform(i,j) - virial_fd_delta
                if(i /= j) deform(j,i) = deform(j,i) - virial_fd_delta
                call set_lattice(at, matmul(deform,lat_save), scale_positions=.true.)
                if(at%cutoff > 0) call calc_dists(at)
                call calc(this, at, args_str=new_args_str, error=error)
                call get_param_value(at, "fd_energy", e_minus, error=error)
                PASS_ERROR_WITH_INFO("Potential_Simple_calc doing fd virial failed to get energy property fd_energy", error)
                virial(i,j) = -(e_plus-e_minus)/(2.0_dp*virial_fd_delta)
                virial(j,i) = virial(i,j)
             end do
          end do
          call set_value(at%params, trim(calc_virial), virial)
          at%pos = allpos_save
          call set_lattice(at, lat_save, scale_positions=.false.)
          call remove_value(at%params, "fd_energy")
          if(at%cutoff > 0) call calc_dists(at)
          deallocate(allpos_save)
          call print("Done with finite difference virial calculation", PRINT_VERBOSE)
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
    elseif(associated(this%socketpot)) then
       call Print(this%socketpot, file=file)
    elseif(this%is_wrapper) then
       call print("Potential_Simple: wrapper Potential_Simple")
    else
       RAISE_ERROR ("Potential_Simple_Print: no potential type is set", error)
    end if

  end subroutine Potential_Simple_Print

  subroutine Potential_Simple_setup_parallel(this, at, args_str, error)
    type(Potential_Simple), intent(inout) :: this
    type(Atoms), intent(inout) :: at     !% The atoms structure to compute energy and forces
    character(len=*), intent(in) :: args_str
    integer, intent(out), optional :: error

    real(dp) :: e                   
    real(dp), allocatable :: f(:,:)              
    type(Dictionary) :: params
    character(STRING_LENGTH) :: calc_energy, calc_force

    INIT_ERROR(error)

    call initialise(params)
    call param_register(params, "energy", "", calc_energy, help_string="If present, calculate energy.  Also name of parameter to put energy into")
    call param_register(params, "force", "", calc_force, help_string="If present, calculate forces.  Also name of property to put forces into")
    if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='Potential_Simple_setup_parallel args_str')) then
       RAISE_ERROR("Potential_Simple_setup_parallel failed to parse args_str='"//trim(args_str)//"'", error)
    endif
    call finalise(params)


    if(associated(this%ip)) then
       if (len_trim(calc_energy) > 0) then
	  if (len_trim(calc_force) > 0) then
	     allocate(f(3,at%n))
	     call setup_parallel(this%ip, at, energy=e, f=f, args_str=args_str)
	     deallocate(f)
	  else
	     call setup_parallel(this%ip, at, energy=e, args_str=args_str)
	  endif
       else
	  if (len_trim(calc_force) > 0) then
	     allocate(f(3,at%n))
	     call setup_parallel(this%ip, at, f=f, args_str=args_str)
	     deallocate(f)
	  else
	     call setup_parallel(this%ip, at, args_str=args_str)
	  endif
       endif
#ifdef HAVE_TB
    else if(associated(this%tb)) then
       return
#endif
    elseif(associated(this%filepot)) then
       return
    elseif(associated(this%callbackpot)) then
       return
    elseif(associated(this%socketpot)) then
       return
    elseif(this%is_wrapper) then
       return
    else
       RAISE_ERROR ("Potential_Simple_Print: Potential_Simple is not initialised", error)
    end if

    if (allocated(f)) deallocate(f)

  end subroutine Potential_Simple_setup_parallel

  subroutine Potential_Simple_set_callback(this, callback, error)
    type(Potential_Simple), intent(inout) :: this
    interface
       subroutine callback(at)
         integer, intent(in) :: at(12)
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
