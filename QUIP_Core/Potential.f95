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

!%    This module encapsulates all the interatomic potentials implemented in QUIP
!%    
!%    A Potential object represents an interatomic potential, a
!%    tight binding model or an interface to an external code used to
!%    perform calculations. It is initialised from an `args_str`
!%    describing the type of potential, and an XML formatted string
!%    `param_str` giving the parameters.
!%
!%    Types of Potential:
!%    
!%      ====================  ==========================================================
!%      `args_str` prefix      Description
!%      ====================  ==========================================================   
!%      ``IP``                Interatomic Potential
!%      ``TB``                Tight Binding Model
!%      ``FilePot``           File potential, used to communicate with external program
!%      ``CallbackPot``       Callback potential, computation done by Python function
!%      ``Sum``               Sum of two other potentials
!%      ``ForceMixing``       Combination of forces from two other potentials
!%      ====================  ==========================================================
!%    
!%    
!%    Types of interatomic potential available:
!%    
!%      ======================== ==========================================================
!%      `args_str` prefix        Description
!%      ======================== ==========================================================   
!%      ``IP BOP``               Bond order potential for metals
!%      ``IP BornMayer``         Born-Mayer potential for oxides
!%                               (e.g. BKS potential for silica)
!%      ``IP Brenner``           Brenner (1990) potential for carbon
!%      ``IP Brenner_2002``      Brenner (2002) reactive potential for carbon
!%      ``IP Brenner_Screened``  Interface to Pastewka et al. screened Brenner reactive
!%                               potential for carbon
!%      ``IP Coulomb``           Coulomb interaction: support direct summation,
!%                               Ewald and damped shifted force Coulomb potential
!%      ``IP Einstein``          Einstein crystal potential
!%      ``IP EAM_ErcolAd``       Embedded atom potential of Ercolessi and Adams
!%      ``IP FB``                Flikkema and Bromley potential
!%      ``IP FS``                Finnis-Sinclair potential for metals
!%      ``IP FX``                Wrapper around ttm3f water potential of
!%                               Fanourgakis-Xantheas
!%      ``IP GAP``               Gaussian approximation potential
!%      ``IP Glue``              Generic implementation of 'glue' potential
!%      ``IP HFdimer``           Simple interatomic potential for an HF dimer, from
!%                               MP2 calculations
!%      ``IP KIM``               Interface to KIM, the Knowledgebase of Interatomic
!%                               potential Models (www.openkim.org)
!%      ``IP LJ``                Lennard-Jones potential
!%      ``IP Morse``             Morse potential
!%      ``IP PartridgeSchwenke`` Partridge-Schwenke model for a water monomer
!%      ``IP SW``                Stillinger-Weber potential for silicon
!%      ``IP SW_VP``             Combined Stillinger-Weber and Vashista potential
!%                               for Si and :mol:`SiO_2`.
!%      ``IP Si_MEAM``           Silicon modified embedded attom potential
!%      ``IP Sutton_Chen``       Sutton-Chen potential
!%      ``IP TS``                Tangney-Scandolo polarisable potential for oxides
!%      ``IP Tersoff``           Tersoff potential for silicon
!%      ``IP WaterDimer_Gillan`` 2-body potential for water dimer
!%      ======================== ==========================================================   
!%    
!%    Types of tight binding potential available:
!%    
!%      ======================= ==========================================================
!%      `args_str` prefix       Description
!%      ======================= ==========================================================   
!%      ``TB Bowler``           Bowler tight binding model 
!%      ``TB DFTB``             Density functional tight binding
!%      ``TB GSP``              Goodwin-Skinner-Pettifor tight binding model
!%      ``TB NRL_TB``           Naval Research Laboratory tight binding model
!%      ======================= ==========================================================   
!%    
!%    Examples of the XML parameters for each of these potential can be
!%    found in the `QUIP_Core/parameters <http://src.tcm.phy.cam.ac.uk/viewvc/jrk33/repo/tags/QUIP_release/QUIP_Core/parameters/>`_
!%    directory of the QUIP svn repository.
!%
!%!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#ifdef HAVE_LOCAL_E_MIX
#define HAVE_HYBRID
#endif
#ifdef HAVE_ONIOM
#define HAVE_HYBRID
#endif

#include "error.inc"
module Potential_module

  use error_module
  use system_module, only : dp, inoutput, print, PRINT_ALWAYS, PRINT_NORMAL, PRINT_VERBOSE, PRINT_NERD, initialise, finalise, INPUT, &
   optional_default, current_verbosity, mainlog, round, verbosity_push_decrement, verbosity_push, verbosity_pop, print_warning, system_timer, system_abort, operator(//)
  use units_module, only : GPA
  use periodictable_module, only :  ElementCovRad, ElementMass
  use extendable_str_module, only : extendable_str, initialise, read, string, finalise
  use linearalgebra_module , only : norm, trace, matrix3x3_det, normsq, least_squares, add_identity, inverse, diagonalise, symmetric_linear_solve, operator(.fne.), operator(.mult.), operator(.feq.), print
  use dictionary_module, only : dictionary, STRING_LENGTH, lookup_entry_i, write_string, get_value, has_key, read_string, set_value, remove_value, initialise, finalise
  use paramreader_module, only : param_register, param_read_line
  use mpi_context_module, only : mpi_context
  use table_module, only : table, find, int_part, wipe, append, finalise
  use minimization_module , only : minim, n_minim, fire_minim, test_gradient, n_test_gradient
  use connection_module, only : connection
  use atoms_types_module, only : atoms, assign_pointer, add_property, assign_property_pointer, add_property_from_pointer, diff_min_image, distance_min_image
  use atoms_module, only : has_property, cell_volume, neighbour, n_neighbours, set_lattice, is_nearest_neighbour, &
   get_param_value, remove_property, calc_connect, set_cutoff, set_param_value, calc_dists, atoms_repoint, finalise, assignment(=)
  use cinoutput_module, only : cinoutput, write
  use dynamicalsystem_module, only : dynamicalsystem, ds_print_status, advance_verlet1, advance_verlet2
  use clusters_module, only : HYBRID_ACTIVE_MARK, HYBRID_NO_MARK, HYBRID_BUFFER_MARK, create_embed_and_fit_lists_from_cluster_mark, create_embed_and_fit_lists, &
   create_hybrid_weights, add_cut_hydrogens, create_cluster_info_from_mark, bfs_grow, carve_cluster

  use QUIP_Common_module
#ifdef HAVE_TB
  use TB_module
#endif
  use Potential_simple_module

  use adjustablepotential_module, only: adjustable_potential_init, adjustable_potential_optimise, &
       adjustable_potential_force, adjustable_potential_finalise

  implicit none

#ifndef POTENTIAL_NO_DEFAULT_PRIVATE
  private
#endif

  !*************************************************************************
  !*
  !*  Potential header
  !*
  !*************************************************************************

  integer :: hack_restraint_i(2)
  real(dp) :: hack_restraint_r, hack_restraint_k
  public :: hack_restraint_i, hack_restraint_r, hack_restraint_k

  public :: Potential
  !%  Potential type which abstracts all QUIP interatomic potentials
  !%
  !%  Provides interface to all energy/force/virial calculating schemes,
  !%  including actual calculations, as well as abstract hybrid schemes
  !%  such as LOTF, Force Mixing, ONIOM, and the local energy scheme.
  !%
  !%  Typically a Potential is constructed from an initialisation
  !%  args_str and an XML parameter file, e.g. in Fortran::
  !%
  !%      type(InOutput) :: xml_file
  !%      type(Potential) :: pot
  !%      ...
  !%      call initialise(xml_file, 'SW.xml', INPUT)
  !%      call initialise(pot, 'IP SW', param_file=xml_file)
  !%
  !%  Or, equivaently in Python::
  !%
  !%      pot = Potential('IP SW', param_filename='SW.xml')
  !%
  !%  creates a Stillinger-Weber potential using the parameters from
  !%  the file 'SW.xml'. The XML parameters can also be given directly
  !%  as a string, via the `param_str` argument.
  !%
  !%  The main workhorse is the :meth:`calc` routine, which is used
  !%  internally to perform all calculations, e.g. to calculate forces::
  !%
  !%      type(Atoms) :: at
  !%      real(dp) :: force(3,8)
  !%      ...
  !%      call diamond(at, 5.44, 14)
  !%      call randomise(at%pos, 0.01)
  !%      call set_cutoff(at, cutoff(pot))
  !%      call calc_connect(at)
  !%      call calc(pot, at, force=force)
  !%
  !%  Note how the 'Atoms%cutoff' attribute is set to the cutoff of
  !%  this Potential, and then the neighbour lists are built with the
  !%  :meth:`~quippy.atoms.Atoms.calc_connect` routine.
  !%
  !%  A Potential can be used to optimise the geometry of an
  !%  :class:`~quippy.atoms.Atoms` structure, using the :meth:`minim` routine,
  !%  (or, in Python, via  the :class:`Minim` wrapper class).
  type Potential
     type(MPI_context) :: mpi
     character(len=STRING_LENGTH) :: init_args_pot1, init_args_pot2, xml_label, xml_init_args, calc_args = ""

     logical :: is_simple = .false.
     type(Potential_simple) :: simple

     type(Potential), pointer :: l_mpot1 => null() ! local copies of potentials, if they are created by potential_initialise from an args string
     type(Potential), pointer :: l_mpot2 => null() ! local copies of potentials, if they are created by potential_initialise from an args string

     logical :: is_sum = .false.
     type(Potential_Sum), pointer :: sum => null()

     logical :: is_forcemixing = .false.
     type(Potential_FM), pointer :: forcemixing => null()

     logical :: is_evb = .false.
     type(Potential_EVB), pointer :: evb => null()

#ifdef HAVE_LOCAL_E_MIX
     logical :: is_local_e_mix = .false.
     type(Potential_Local_E_Mix), pointer :: local_e_mix => null()
#endif

#ifdef HAVE_ONIOM
     logical :: is_oniom = .false.
     type(Potential_ONIOM), pointer :: oniom => null()
#endif

     logical :: is_cluster = .false.
     type(Potential_Cluster), pointer :: cluster => null()

     logical :: do_rescale_r, do_rescale_E
     real(dp) :: r_scale, E_scale

  end type Potential

  public :: Initialise, Potential_Filename_Initialise
  interface Initialise
     module procedure Potential_Initialise, Potential_Initialise_inoutput
  end interface

  public :: Finalise
  interface Finalise
     module procedure Potential_Finalise
  end interface

  public :: Setup_Parallel
  interface Setup_Parallel
     module procedure Potential_Setup_Parallel
  end interface

  public :: Print
  interface Print
     module procedure Potential_Print
  end interface


  !% Return the cutoff of this 'Potential', in Angstrom. This is the
  !% minimum neighbour connectivity cutoff that should be used: if
  !% you're doing MD you'll want to use a slightly larger cutoff so
  !% that new neighbours don't drift in to range between connectivity
  !% updates
  public :: Cutoff
  interface Cutoff
     module procedure Potential_Cutoff
  end interface

  public :: Calc

  !% Apply this Potential to the Atoms object
  !% 'at', which must have connectivity information
  !% (i.e. 'Atoms%calc_connect' should have been called). The
  !% optional arguments determine what should be calculated and how
  !% it will be returned. Each physical quantity has a
  !% corresponding optional argument, which can either be an 'True'
  !% to store the result inside the Atoms object (i.e. in
  !% Atoms%params' or in 'Atoms%properties' with the
  !% default name, a string to specify a different property or
  !% parameter name, or an array of the the correct shape to
  !% receive the quantity in question, as set out in the table
  !% below.
  !%
  !%   ================  ============= ================ =========================
  !%   Array argument    Quantity      Shape            Default storage location
  !%   ================  ============= ================ =========================
  !%   ``energy``        Energy        ``()``  	        ``energy`` param
  !%   ``local_energy``  Local energy  ``(at.n,)``      ``local_energy`` property
  !%   ``force``         Force         ``(3,at.n)``     ``force`` property
  !%   ``virial``        Virial tensor ``(3,3)``        ``virial`` param
  !%   ``local_virial``  Local virial  ``(3,3,at.n)``   ``local_virial`` property
  !%   ================  ============= ================ =========================
  !%
  !% The 'args_str' argument is an optional string  containing
  !% additional arguments which depend on the particular Potential
  !% being used. 
  !%
  !% Not all Potentials support all of these quantities: an error
  !% will be raised if you ask for something that is not supported.
  interface Calc
     module procedure Potential_Calc
  end interface

  !% Minimise the configuration 'at' under the action of this
  !% Potential.  Returns number of minimisation steps taken. If
  !% an error occurs or convergence is not reached within 'max_steps'
  !% steps, 'status' will be set to 1 on exit.
  !%
  !% Example usage (in Python, Fortran code is similar. See
  !% :ref:`geomopt` in the quippy tutorial for full
  !% explanation)::
  !%
  !%>      at0 = diamond(5.44, 14)
  !%>      at0.calc_connect()
  !%>      pot = Potential('IP SW', param_str='''<SW_params n_types="1">
  !%>              <comment> Stillinger and Weber, Phys. Rev. B  31 p 5262 (1984)</comment>
  !%>              <per_type_data type="1" atomic_num="14" />
  !%>
  !%>              <per_pair_data atnum_i="14" atnum_j="14" AA="7.049556277" BB="0.6022245584"
  !%>                p="4" q="0" a="1.80" sigma="2.0951" eps="2.1675" />
  !%>
  !%>              <per_triplet_data atnum_c="14" atnum_j="14" atnum_k="14"
  !%>                lambda="21.0" gamma="1.20" eps="2.1675" />
  !%>             </SW_params>''')
  !%>      pot.minim(at0, 'cg', 1e-7, 100, do_pos=True, do_lat=True)
  public :: Minim
  interface Minim
     module procedure Potential_Minim
  end interface

  public :: test_gradient
  interface test_gradient
     module procedure pot_test_gradient
  end interface

  public :: n_test_gradient
  interface n_test_gradient
     module procedure pot_n_test_gradient
  end interface

  public :: set_callback
  interface set_callback
     module procedure potential_set_callback
  end interface

  ! Allow Potential_Precon_Minim access to some Potential Functions
  public :: energy_func
  public :: gradient_func
  public :: print_hook
  public :: pack_pos_dg
  public :: unpack_pos_dg
  public :: fix_atoms_deform_grad
  public :: prep_atoms_deform_grad
  public :: max_rij_change
  public :: constrain_virial_post


  public :: potential_minimise
  type Potential_minimise
     real(dp) :: minim_pos_lat_preconditioner = 1.0_dp

     real(dp) :: minim_save_lat(3,3)
     character(len=STRING_LENGTH) :: minim_args_str
     logical :: minim_do_pos, minim_do_lat
     integer :: minim_n_eval_e, minim_n_eval_f, minim_n_eval_ef
     real(dp) :: pos_lat_preconditioner_factor
     type(Atoms), pointer :: minim_at => null()
     real(dp), pointer :: last_connect_x(:) => null()
     type(Potential), pointer :: minim_pot => null()
     type(CInoutput), pointer :: minim_cinoutput_movie => null()
     real(dp), dimension(3,3) :: external_pressure = 0.0_dp
    
    logical :: connectivity_rebuilt = .false.

  end type Potential_minimise

  type(Potential), pointer :: parse_pot
  logical, save :: parse_in_pot, parse_in_pot_done, parse_matched_label

  public :: run

  interface run
     module procedure dynamicalsystem_run
  end interface run

#include "Potential_Sum_header.f95"
#include "Potential_ForceMixing_header.f95"
#include "Potential_EVB_header.f95"

#ifdef HAVE_LOCAL_E_MIX
#include "Potential_Local_E_Mix_header.f95"
#endif
#ifdef HAVE_ONIOM
#include "Potential_ONIOM_header.f95"
#endif
#include "Potential_Cluster_header.f95"

  ! Public interfaces from Potential_Hybrid_utils.f95

  public :: bulk_modulus
  interface bulk_modulus
     module procedure potential_bulk_modulus
  end interface bulk_modulus

  public :: test_local_virial
  interface test_local_virial
     module procedure potential_test_local_virial
  end interface test_local_virial

  contains

  !*************************************************************************
  !*
  !*  Potential routines
  !*
  !*************************************************************************

recursive subroutine potential_Filename_Initialise(this, args_str, param_filename, bulk_scale, mpi_obj, error)
  type(Potential), intent(inout) :: this
  character(len=*), intent(in) :: args_str !% Valid arguments are 'Sum', 'ForceMixing', 'EVB', 'Local_E_Mix' and 'ONIOM', and any type of simple_potential
  character(len=*), intent(in) :: param_filename !% name of xml parameter file for potential initializers
  type(Atoms), optional, intent(inout) :: bulk_scale !% optional bulk structure for calculating space and E rescaling
  type(MPI_Context), intent(in), optional :: mpi_obj
  integer, intent(out), optional :: error

  type(inoutput) :: io

  INIT_ERROR(error)

  if (len_trim(param_filename) > 0) call initialise(io, param_filename, INPUT, master_only=.true.)
  call initialise(this, args_str, io, bulk_scale=bulk_scale, mpi_obj=mpi_obj, error=error)
  PASS_ERROR(error)
  call finalise(io)

end subroutine potential_Filename_Initialise

subroutine potential_initialise_inoutput(this, args_str, io_obj, bulk_scale, mpi_obj, error)
  type(Potential), intent(inout) :: this
  character(len=*), intent(in) :: args_str !% Valid arguments are 'Sum', 'ForceMixing', 'EVB', 'Local_E_Mix' and 'ONIOM', and any type of simple_potential
  type(InOutput), intent(in) :: io_obj !% name of xml parameter inoutput for potential initializers
  type(Atoms), optional, intent(inout) :: bulk_scale !% optional bulk structure for calculating space and E rescaling
  type(MPI_Context), intent(in), optional :: mpi_obj
  integer, intent(out), optional :: error

  type(extendable_str) :: es

  INIT_ERROR(error)

  call initialise(es)
  if (io_obj%initialised) then
     if (present(mpi_obj)) then
       call read(es, io_obj%unit, convert_to_string=.true., mpi_comm=mpi_obj%communicator)
     else
       call read(es, io_obj%unit, convert_to_string=.true.)
     endif
  else
     call initialise(es)
  endif

  call initialise(this, args_str=args_str, param_str=string(es), bulk_scale=bulk_scale, mpi_obj=mpi_obj, error=error)
  PASS_ERROR(error)
  call finalise(es)

end subroutine potential_initialise_inoutput

recursive subroutine potential_initialise(this, args_str, pot1, pot2, param_str, bulk_scale, mpi_obj, error)
  type(Potential), intent(inout) :: this
  character(len=*), optional, intent(in) :: args_str !% Valid arguments are 'Sum', 'ForceMixing', 'EVB', 'Local_E_Mix' and 'ONIOM', and any type of simple_potential
  type(Potential), optional, intent(in), target :: pot1 !% Optional first Potential upon which this Potential is based
  type(Potential), optional, intent(in), target :: pot2 !% Optional second potential
  character(len=*), optional, intent(in) :: param_str !% contents of xml parameter file for potential initializers, if needed
  type(Atoms), optional, intent(inout) :: bulk_scale !% optional bulk structure for calculating space and E rescaling
  type(MPI_Context), optional, intent(in) :: mpi_obj
  integer, intent(out), optional :: error

  type(Potential), pointer :: u_pot1, u_pot2
  type(Dictionary) :: params

  logical :: has_xml_label, minimise_bulk, has_target_vol, has_target_B, has_r_scale, has_E_scale
  real(dp) :: target_vol, target_B, r_scale, vol, B
  character(len=STRING_LENGTH) :: my_args_str, init_args_str, first_tag, xml_label
  type(Atoms) :: bulk
  integer :: it, tag_start, tag_end

  INIT_ERROR(error)

  call finalise(this)

  my_args_str = ''
  if (present(args_str)) my_args_str = args_str

  ! Default init_args are xml_label=LABEL, extracting LABEL from first XML tag
  if (len_trim(my_args_str) == 0) then
     if (.not. present(param_str)) then
        RAISE_ERROR('No init_args given and param_str not present', error)
     end if
     if (len_trim(param_str) == 0) then
        RAISE_ERROR('No init_args given and param_str is empty', error)
     end if
     tag_start = index(param_str, '<')
     tag_end = index(param_str, '>')
     xml_label = param_str(tag_start+1:tag_end-1)
     call print_warning('Potential_initialise using default init_args "Potential xml_label='//trim(xml_label)//'"')
     my_args_str = 'Potential xml_label='//xml_label
  end if

  call initialise(params)
  call param_register(params, 'xml_label', '', this%xml_label, has_value_target=has_xml_label, help_string="Label in xml file Potential stanza to match")
  call param_register(params, 'calc_args', '', this%calc_args, help_string="Default calc_args that are passed each time calc() is called")
  if(.not. param_read_line(params, my_args_str, ignore_unknown=.true.,task='Potential_Initialise args_str')) then
    RAISE_ERROR("Potential_initialise failed to parse args_str='"//trim(my_args_str)//"'", error)
  endif
  call finalise(params)

  my_args_str = ""
  if(has_xml_label) then
     call Potential_read_params_xml(this, param_str)
     my_args_str = trim(this%xml_init_args)
  endif
  my_args_str = trim(my_args_str) // " " // trim(args_str)

  call initialise(params)
  call param_register(params, 'init_args_pot1', '', this%init_args_pot1, help_string="Argument string for initializing pot1 (for non-simple potentials")
  call param_register(params, 'init_args_pot2', '', this%init_args_pot2, help_string="Argument string for initializing pot2 (for non-simple potentials")
  call param_register(params, 'Sum', 'false', this%is_sum, help_string="Potential that's a sum of 2 other potentials")
  call param_register(params, 'ForceMixing', 'false', this%is_forcemixing, help_string="Potential that's force-mixing of 2 other potentials")
  call param_register(params, 'EVB', 'false', this%is_evb, help_string="Potential using empirical-valence bond to mix 2 other potentials")
#ifdef HAVE_LOCAL_E_MIX
  call param_register(params, 'Local_E_Mix', 'false', this%is_local_e_mix, help_string="Potential that's local energy mixing of 2 other potentials")
#endif /* HAVE_LOCAL_E_MIX */
#ifdef HAVE_ONIOM
  call param_register(params, 'ONIOM', 'false', this%is_oniom, help_string="Potential from ONIOM mixing of two other potential energies")
#endif /* HAVE_ONIOM */
  call param_register(params, 'Cluster', 'false', this%is_cluster, help_string="Potential evaluated using clusters")
  call param_register(params, 'do_rescale_r', 'F', this%do_rescale_r, help_string="If true, rescale distances by factor r_scale.")
  call param_register(params, 'r_scale', '1.0', this%r_scale, has_value_target=has_r_scale, help_string="Recaling factor for distances. Default 1.0.")
  call param_register(params, 'do_rescale_E', 'F', this%do_rescale_E, help_string="If true, rescale energy by factor E_scale.")
  call param_register(params, 'E_scale', '1.0', this%E_scale, has_value_target=has_E_scale, help_string="Recaling factor for energy. Default 1.0.")
  call param_register(params, "minimise_bulk", "F", minimise_bulk, help_string="If true, minimise bulk_scale structure before measuring eqm. volume and bulk modulus for rescaling")
  call param_register(params, "target_vol", "0.0", target_vol, has_value_target=has_target_vol, &
       help_string="Target volume per cell used if do_rescale_r=T Unit is A^3.")
  call param_register(params, "target_B", "0.0", target_B, has_value_target=has_target_B, &
       help_string="Target bulk modulus used if do_rescale_E=T. Unit is GPa.")

  if(.not. param_read_line(params, my_args_str, ignore_unknown=.true.,task='Potential_Initialise args_str')) then
    RAISE_ERROR("Potential_initialise failed to parse args_str='"//trim(my_args_str)//"'", error)
  endif
  call finalise(params)

  if (present(pot1) .and. len_trim(this%init_args_pot1) > 0) then
    RAISE_ERROR("Potential_initialise got both pot and args_str with init_args_pot1 passed in, conflict", error)
   endif
  if (present(pot2) .and. len_trim(this%init_args_pot2) > 0) then
    RAISE_ERROR("Potential_initialise got both pot2 and args_str with init_args_pot2 passed in, conflict", error)
  endif

  if (len_trim(this%init_args_pot1) > 0) then
    allocate(this%l_mpot1)
    call initialise(this%l_mpot1, args_str=this%init_args_pot1, param_str=param_str, bulk_scale=bulk_scale, mpi_obj=mpi_obj, error=error)
    PASS_ERROR_WITH_INFO("Initializing pot1", error)
    u_pot1 => this%l_mpot1
  else
    u_pot1 => pot1
  endif
  if (len_trim(this%init_args_pot2) > 0) then
    allocate(this%l_mpot2)
    call initialise(this%l_mpot2, args_str=this%init_args_pot2, param_str=param_str, bulk_scale=bulk_scale, mpi_obj=mpi_obj, error=error)
    PASS_ERROR_WITH_INFO("Initializing pot2", error)
    u_pot2 => this%l_mpot2
  else
    if (present(pot2)) then
      u_pot2 => pot2
    else
      nullify(u_pot2)
    endif
  endif

  if (present(mpi_obj)) this%mpi = mpi_obj

  this%is_simple = .not. any( (/ this%is_sum, this%is_forcemixing, this%is_evb &
#ifdef HAVE_LOCAL_E_MIX
  , this%is_local_e_mix &
#endif
#ifdef HAVE_ONIOM
  , this%is_oniom  &
#endif
  , this%is_cluster &
  /) )

  if (count( (/this%is_simple, this%is_sum, this%is_forcemixing, this%is_evb &
#ifdef HAVE_LOCAL_E_MIX
          , this%is_local_e_mix &
#endif /* HAVE_LOCAL_E_MIX */
#ifdef HAVE_ONIOM
          , this%is_oniom &
#endif /* HAVE_ONIOM */
          , this%is_cluster &
          /) ) /= 1) then
     RAISE_ERROR("Potential_initialise found too few or two many Potential types args_str='"//trim(my_args_str)//"'", error)
  end if

  if (this%is_simple) then
    if (present(bulk_scale) .and. .not. has_target_vol .and. .not. has_target_B) call print("Potential_initialise Simple ignoring bulk_scale passed in", PRINT_ALWAYS)

    call initialise(this%simple, my_args_str, param_str, mpi_obj, error=error)
    PASS_ERROR_WITH_INFO("Initializing pot", error)
  else if (this%is_sum) then
    if (present(bulk_scale)) call print("Potential_initialise Sum ignoring bulk_scale passed in", PRINT_ALWAYS)

    if (.not. associated(u_pot2)) then
      RAISE_ERROR('Potential_initialise: two potentials needs for sum potential', error)
    endif

    allocate(this%sum)
    call initialise(this%sum, my_args_str, u_pot1, u_pot2, mpi_obj, error=error)
    PASS_ERROR_WITH_INFO("Initializing sum", error)

  else if (this%is_forcemixing) then

    allocate(this%forcemixing)
    if(associated(u_pot2)) then
       call initialise(this%forcemixing, my_args_str, u_pot1, u_pot2, bulk_scale, mpi_obj, error=error)
    else
       ! if only one pot is given, assign it to QM. this is useful for time-embedding LOTF
       call initialise(this%forcemixing, my_args_str, qmpot=u_pot1, reference_bulk=bulk_scale, mpi=mpi_obj, error=error)
    endif
    PASS_ERROR_WITH_INFO("Initializing forcemixing", error)

  else if (this%is_evb) then
    if (present(bulk_scale)) call print("Potential_initialise EVB ignoring bulk_scale passed in", PRINT_ALWAYS)

    !1 is enough, everything is the same apart from the passed topology
    if (associated(u_pot2)) then
      call print('WARNING! Potential_initialise EVB ignoring u_pot2 (it will use u_pot1 twice)', PRINT_ALWAYS)
    endif

    allocate(this%evb)
    call initialise(this%evb, my_args_str, u_pot1, mpi_obj, error=error)
    PASS_ERROR_WITH_INFO("Initializing evb", error)

#ifdef HAVE_LOCAL_E_MIX
  else if (this%is_local_e_mix) then
    if (.not. associated(u_pot2)) then
      RAISE_ERROR('Potential_initialise: two potentials needs for local_e_mix potential', error)
   endif

    allocate(this%local_e_mix)
    call initialise(this%local_e_mix, my_args_str, u_pot1, u_pot2, bulk_scale, mpi_obj, error=error)
    PASS_ERROR_WITH_INFO("Initializing local_e_mix", error)
#endif /* HAVE_LOCAL_E_MIX */
#ifdef HAVE_ONIOM
  else if (this%is_oniom) then
    if (.not. associated(u_pot2)) then
      RAISE_ERROR('Potential_initialise: two potentials needs for oniom potential', error)
    endif

    allocate(this%oniom)
    call initialise(this%oniom, my_args_str, u_pot1, u_pot2, bulk_scale, mpi_obj, error=error)
    PASS_ERROR_WITH_INFO("Initializing oniom", error)
#endif /* HAVE_ONIOM */

  else if (this%is_cluster) then
     allocate(this%cluster)
     call initialise(this%cluster, my_args_str, u_pot1, mpi_obj, error=error)
     PASS_ERROR_WITH_INFO("Initializing Cluster", error)
  end if

  if (this%is_simple .and. (this%do_rescale_r .or. this%do_rescale_E)) then
     if (this%do_rescale_r) then
        if (.not. has_r_scale .and. .not. has_target_vol) then
           RAISE_ERROR("potential_initialise: do_rescale_r=T but neither r_scale nor target_vol present", error)
        end if
     end if
     if (this%do_rescale_E) then
        if (.not. has_E_scale .and. .not. has_target_B) then
           RAISE_ERROR("potential_initialise: do_rescale_E=T but neither E_scale nor target_B present", error)
        end if
     end if
     
     ! Automatic calculation of rescaling factors r_scale and E_scale to match target_vol and/or target_B
     if (has_target_vol .or. has_target_B) then
        if (.not. present(bulk_scale)) then
           RAISE_ERROR("potential_initialise: target_vol or target_B present but no bulk_scale provided", error)
        end if
        
        bulk = bulk_scale
        call set_cutoff(bulk, cutoff(this)+1.0_dp)
        call calc_connect(bulk)
        if (minimise_bulk) then
           call print("MINIMISING bulk in potential", PRINT_VERBOSE)
           call verbosity_push_decrement()
           it = minim(this, at=bulk, method='cg', convergence_tol=1e-6_dp, max_steps=100, linminroutine='NR_LINMIN', &
                do_print = .false., do_lat = .true., args_str=args_str)
           call verbosity_pop()
        endif

        r_scale = 1.0_dp
        if (this%do_rescale_r) then
           vol = cell_volume(bulk)
           r_scale = (vol/target_vol)**(1.0_dp/3.0_dp)
           call print("do_rescale_r: original cell volume (from minim) is "//vol//" A^3")
           call print("do_rescale_r: setting potential r_scale to "//r_scale//" to match target cell volume of "//target_vol//" A^3")
        end if

        if (this%do_rescale_E) then
           call bulk_modulus(this, bulk, B, vol, minimise_bulk, args_str=args_str)
           this%E_scale = target_B/(B*r_scale**3)
           call print("do_rescale_E: original cell volume (from quadratic fit) is "//vol//" A^3")
           call print("do_rescale_E: original bulk modulus (from quadratic fit) is "//B//" GPa")
           call print("do_rescale_E: setting potential E_scale to "//this%E_scale//" to match target bulk modulus of "//target_B//" GPa")
        end if

        if (this%do_rescale_r) this%r_scale = r_scale
     end if
  end if

  end subroutine potential_initialise

  recursive subroutine potential_finalise(this, error)
    type(Potential), intent(inout) :: this
    integer, intent(out), optional :: error

  INIT_ERROR(error)

    if (this%is_simple) then
       call finalise(this%simple)
    else if (this%is_sum) then
       call finalise(this%sum)
       deallocate(this%sum)
    else if (this%is_forcemixing) then
       call finalise(this%forcemixing)
       deallocate(this%forcemixing)
    else if (this%is_evb) then
       call finalise(this%evb)
       deallocate(this%evb)
#ifdef HAVE_LOCAL_E_MIX
    else if (this%is_local_e_mix) then
       call finalise(this%local_e_mix)
       deallocate(this%local_e_mix)
#endif /* HAVE_LOCAL_E_MIX */
#ifdef HAVE_ONIOM
    else if (this%is_oniom) then
       call finalise(this%oniom)
       deallocate(this%oniom)
#endif /* HAVE_ONIOM */
    else if (this%is_cluster) then
       call finalise(this%cluster)
       deallocate(this%cluster)
    end if

    if (associated(this%l_mpot1)) call finalise(this%l_mpot1)
    if (associated(this%l_mpot2)) call finalise(this%l_mpot2)
    nullify(this%l_mpot1)
    nullify(this%l_mpot2)

    this%is_simple = .false.
    this%is_sum = .false.
    this%is_forcemixing = .false.
    this%is_evb = .false.
#ifdef HAVE_LOCAL_E_MIX
    this%is_local_e_mix = .false.
#endif
#ifdef HAVE_ONIOM
    this%is_oniom = .false.
#endif
    this%is_cluster = .false.

  end subroutine potential_finalise

  recursive subroutine potential_calc(this, at, energy, force, virial, local_energy, local_virial, args_str, error)
    type(Potential), intent(inout) :: this
    type(Atoms), intent(inout) :: at
    real(dp), intent(out), optional :: energy
    real(dp), intent(out), optional :: force(:,:)
    real(dp), intent(out), optional :: virial(3,3)
    real(dp), intent(out), optional :: local_virial(:,:)
    real(dp), intent(out), optional :: local_energy(:)
    character(len=*), intent(in), optional :: args_str
    integer, intent(out), optional :: error

    real(dp), pointer :: at_force_ptr(:,:), at_local_energy_ptr(:), at_local_virial_ptr(:,:)
    type(Dictionary) :: params
    character(len=STRING_LENGTH) :: calc_energy, calc_force, calc_virial, calc_local_energy, calc_local_virial, extra_args_str
    character(len=STRING_LENGTH) :: use_calc_energy, use_calc_force, use_calc_virial, use_calc_local_energy, use_calc_local_virial
    real(dp) :: r_scale, E_scale
    logical :: has_r_scale, has_E_scale
    integer i

    INIT_ERROR(error)

    if (cutoff(this) > 0.0_dp) then
       ! For Potentials which need connectivity information, ensure Atoms cutoff is >= Potential cutoff

       if (.not. at%connect%initialised) then
          call print('Potential_calc: setting Atoms cutoff to Potential cutoff ('//cutoff(this)//')', PRINT_VERBOSE)
          call set_cutoff(at, cutoff(this))
          call calc_connect(at)
       end if

       if (at%cutoff < cutoff(this)) then
          RAISE_ERROR('Potential_calc: cutoff of Atoms object ('//at%cutoff//') < Potential cutoff ('//cutoff(this)//')', error)
       end if
    end if
    
    calc_energy = ""
    calc_virial = ""
    calc_force = ""
    calc_local_energy = ""
    calc_local_virial = ""

    call initialise(params)
    call param_register(params, "energy", "", calc_energy, help_string="If present, calculate energy and put it in field with this string as name")
    call param_register(params, "virial", "", calc_virial, help_string="If present, calculate virial and put it in field with this string as name")
    call param_register(params, "force", "", calc_force, help_string="If present, calculate force and put it in field with this string as name")
    call param_register(params, "local_energy", "", calc_local_energy, help_string="If present, calculate local energy and put it in field with this string as name")
    call param_register(params, "local_virial", "", calc_local_virial, help_string="If present, calculate local virial and put it in field with this string as name")
    call param_register(params, "r_scale", "0.0", r_scale, has_value_target=has_r_scale, help_string="Distance rescale factor. Overrides r_scale init arg")
    call param_register(params, "E_scale", "0.0", E_scale, has_value_target=has_E_scale, help_string="Energy rescale factor. Overrides E_scale init arg")
    if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='Potential_Calc args_str')) then
       RAISE_ERROR('Potential_Calc failed to parse args_str="'//trim(args_str)//'"', error)
    endif
    call finalise(params)

    extra_args_str = ""

    ! create property/param names and possibly storage
    use_calc_energy = trim(calc_energy)
    if (present(energy) .and. len_trim(calc_energy) == 0) then ! have optional and no args_str - make new name for param
       use_calc_energy = "energy"
       do while (lookup_entry_i(at%params, trim(use_calc_energy)) > 0)
	 use_calc_energy = "T"//trim(use_calc_energy)
       end do
       extra_args_str = trim(extra_args_str) // " energy="//trim(use_calc_energy)
    endif

    use_calc_force = trim(calc_force)

    if (present(force)) then
       if (len_trim(calc_force) == 0) then ! have force optional but not args_str - add property from pointer in new property name
	  use_calc_force = "force"
	  do while (has_property(at, trim(use_calc_force)))
	    use_calc_force = "T"//trim(use_calc_force)
	  end do
	  extra_args_str = trim(extra_args_str) // " force="//trim(use_calc_force)
	  call add_property_from_pointer(at, trim(use_calc_force), force)
       else ! has force optional _and_ args_str - add property with its own storage
	  call add_property(at, trim(calc_force), 0.0_dp, n_cols=3, ptr2=at_force_ptr)
       endif
    else if (len_trim(calc_force) > 0) then ! no optional, have args_str - add property with its own storage
       call add_property(at, trim(calc_force), 0.0_dp, n_cols=3, ptr2=at_force_ptr)
    endif

    use_calc_virial = trim(calc_virial) ! have optional and no args_str - make new name for param
    if (present(virial) .and. len_trim(calc_virial) == 0) then
       use_calc_virial = "virial"
       do while (lookup_entry_i(at%params, trim(use_calc_virial)) > 0)
	 use_calc_virial = "T"//trim(use_calc_virial)
       end do
       extra_args_str = trim(extra_args_str) // " virial="//trim(use_calc_virial)
    endif

    use_calc_local_energy = trim(calc_local_energy)
    if (present(local_energy)) then
       if (len_trim(calc_local_energy) == 0) then ! have local_energy optional but not args_str - add property from pointer in new property name
	  use_calc_local_energy = "local_energy"
	  do while (has_property(at, trim(use_calc_local_energy)))
	    use_calc_local_energy = "T"//trim(use_calc_local_energy)
	  end do
	  extra_args_str = trim(extra_args_str) // " local_energy="//trim(use_calc_local_energy)
	  call add_property_from_pointer(at, trim(use_calc_local_energy), local_energy)
       else ! has local_energy optional _and_ args_str - add property with its own storage
	  call add_property(at, trim(calc_local_energy), 0.0_dp, n_cols=1, ptr=at_local_energy_ptr)
       endif
    else if (len_trim(calc_local_energy) > 0) then ! no optional, have args_str - add property with its own storage
       call add_property(at, trim(calc_local_energy), 0.0_dp, n_cols=1, ptr=at_local_energy_ptr)
    endif

    use_calc_local_virial = trim(calc_local_virial)
    if (present(local_virial)) then
       if (len_trim(calc_local_virial) == 0) then ! have local_virial optional but not args_str - add property from pointer in new property name
	  use_calc_local_virial = "local_virial"
	  do while (has_property(at, trim(use_calc_local_virial)))
	    use_calc_local_virial = "T"//trim(use_calc_local_virial)
	  end do
	  extra_args_str = trim(extra_args_str) // " local_virial="//trim(use_calc_local_virial)
	  call add_property_from_pointer(at, trim(use_calc_local_virial), local_virial)
       else ! has local_virial optional _and_ args_str - add property with its own storage
	  call add_property(at, trim(calc_local_virial), 0.0_dp, n_cols=9, ptr2=at_local_virial_ptr)
       endif
    else if (len_trim(calc_local_virial) > 0) then ! no optional, have args_str - add property with its own storage
       call add_property(at, trim(calc_local_virial), 0.0_dp, n_cols=9, ptr2=at_local_virial_ptr)
    endif

    ! Set r_scale, E_scale in extra_args_str if not present in args_str
    if (this%do_rescale_r .and. .not. has_r_scale) extra_args_str = trim(extra_args_str)//" r_scale="//this%r_scale
    if (this%do_rescale_E .and. .not. has_E_scale) extra_args_str = trim(extra_args_str)//" E_scale="//this%E_scale

    ! do actual calculation using args_str
    if (this%is_simple) then
       call Calc(this%simple, at, trim(this%calc_args)//" "//trim(args_str)//" "//trim(extra_args_str), error=error)
       PASS_ERROR(error)
    else if (this%is_sum) then
       call Calc(this%sum, at, trim(this%calc_args)//" "//trim(args_str)//" "//trim(extra_args_str), error=error)
       PASS_ERROR(error)
    else if (this%is_forcemixing) then
       call Calc(this%forcemixing, at, trim(this%calc_args)//" "//trim(args_str)//" "//trim(extra_args_str), error=error)
       PASS_ERROR(error)
    else if (this%is_evb) then
       call Calc(this%evb, at, trim(this%calc_args)//" "//trim(args_str)//" "//trim(extra_args_str), error=error)
       PASS_ERROR(error)
#ifdef HAVE_LOCAL_E_MIX
    else if (this%is_local_e_mix) then
      call Calc(this%local_e_mix, at, trim(this%calc_args)//" "//trim(args_str)//" "//trim(extra_args_str), error=error)
      PASS_ERROR(error)
#endif /* HAVE_local_e_mix */
#ifdef HAVE_ONIOM
    else if (this%is_oniom) then
      call Calc(this%oniom, at, trim(this%calc_args)//" "//trim(args_str)//" "//trim(extra_args_str), error=error)
      PASS_ERROR(error)
#endif /* HAVE_ONIOM */
    else if (this%is_cluster) then
       call Calc(this%cluster, at, trim(this%calc_args)//" "//trim(args_str)//" "//trim(extra_args_str), error=error)
       PASS_ERROR(error)
    else
      RAISE_ERROR('Potential_Calc: no potential type is set',error)
    endif

    ! get values from property/param if necessary, and remove temporary if only optional was passed in
    if (present(energy)) then
      call get_param_value(at, trim(use_calc_energy), energy)
      if (len_trim(calc_energy) == 0) call remove_value(at%params, trim(use_calc_energy))
    endif
    if (present(force)) then
      if (len_trim(calc_force) == 0) then ! property should be a pointer to force optional
	call remove_property(at, trim(use_calc_force))
      else
        force = at_force_ptr
      endif
    endif
    if (present(virial)) then
      call get_param_value(at, trim(use_calc_virial), virial)
      if (len_trim(calc_virial) == 0) call remove_value(at%params, trim(use_calc_virial))
    endif
    if (present(local_energy)) then
      if (len_trim(calc_local_energy) == 0) then ! property should be pointer to local_energy optional
	call remove_property(at, trim(use_calc_local_energy))
      else
	local_energy = at_local_energy_ptr
      endif
    endif
    if (present(local_virial)) then
      if (len_trim(calc_local_virial) == 0) then ! property should be pointer to local_virial optional
	call remove_property(at, trim(use_calc_local_virial))
      else
	local_virial = at_local_virial_ptr
      endif
    endif

  end subroutine potential_calc

  subroutine Potential_setup_parallel(this, at, args_str, error)
    type(Potential), intent(inout) :: this
    type(Atoms), intent(inout) :: at     !% The atoms structure to compute energy and forces
    character(len=*), intent(in) :: args_str
    integer, intent(out), optional :: error

    INIT_ERROR(error)

    if(this%is_simple) then
       call setup_parallel(this%simple, at, args_str)
    endif
  end subroutine Potential_setup_parallel

  recursive subroutine potential_print(this, file, error)
    type(Potential), intent(inout) :: this
    type(Inoutput), intent(inout), optional :: file
    integer, intent(out), optional :: error

    INIT_ERROR(error)

    if (this%is_simple) then
       call Print('Potential containing potential')
       call Print(this%simple, file=file)
    else if (this%is_sum) then
       call Print(this%sum, file=file)
    else if (this%is_forcemixing) then
       call Print(this%forcemixing, file=file)
    else if (this%is_evb) then
       call Print(this%evb, file=file)
#ifdef HAVE_LOCAL_E_MIX
    else if (this%is_local_e_mix) then
       call Print(this%local_e_mix, file=file)
#endif /* HAVE_LOCAL_E_MIX */
#ifdef HAVE_ONIOM
    else if (this%is_oniom) then
       call Print(this%oniom, file=file)
#endif /* HAVE_oniom */
    else if (this%is_cluster) then
       call Print(this%cluster, file=file)
    else
       RAISE_ERROR('Potential_Print: no potential type is set', error)
    end if
  end subroutine potential_print


  recursive function potential_cutoff(this, error)
    type(Potential), intent(in) :: this
    integer, intent(out), optional :: error
    real(dp) :: potential_cutoff

    INIT_ERROR(error)

    if (this%is_simple) then
       potential_cutoff = cutoff(this%simple)
    else if (this%is_sum) then
       potential_cutoff = cutoff(this%sum)
    else if (this%is_forcemixing) then
       potential_cutoff = cutoff(this%forcemixing)
    else if (this%is_evb) then
       potential_cutoff = cutoff(this%evb)
#ifdef HAVE_LOCAL_E_MIX
    else if (this%is_local_e_mix) then
       potential_cutoff = cutoff(this%local_e_mix)
#endif /* HAVE_local_e_mix */
#ifdef HAVE_ONIOM
    else if (this%is_oniom) then
       potential_cutoff = cutoff(this%oniom)
#endif /* HAVE_ONIOM */
    else if (this%is_cluster) then
       potential_cutoff = cutoff(this%cluster)
    else
       RAISE_ERROR('Potential_Cutoff: no potential type is set', error)
    end if
  end function potential_cutoff


  subroutine potential_set_calc_args_str(this, args_str, error)
     type(Potential), intent(inout) :: this
     character(len=*), intent(in) :: args_str
     integer, intent(out), optional :: error

     INIT_ERROR(error)

     if(len_trim(args_str) > len(this%calc_args)) then
        RAISE_ERROR('potential_set_calc_args_str: args_str longer than the avaiable calc_args field', error)
     endif

     this%calc_args = trim(args_str)

  endsubroutine potential_set_calc_args_str

  function potential_get_calc_args_str(this, error)
     type(Potential), intent(inout) :: this
     integer, intent(out), optional :: error
     character(len=STRING_LENGTH) :: potential_get_calc_args_str

     INIT_ERROR(error)

     if(len_trim(this%calc_args) > len(potential_get_calc_args_str)) then
        RAISE_ERROR('potential_get_calc_args_str: calc_args field too long', error)
     endif

     potential_get_calc_args_str = trim(this%calc_args)

  endfunction potential_get_calc_args_str

  !*************************************************************************
  !*
  !*  Minimisation routines
  !*
  !*************************************************************************


  function potential_minim(this, at, method, convergence_tol, max_steps, linminroutine, do_print, print_inoutput, print_cinoutput, &
       do_pos, do_lat, args_str, eps_guess, fire_minim_dt0, fire_minim_dt_max, external_pressure, use_precond, &
       hook_print_interval, error)
    type(Potential), intent(inout), target :: this !% potential to evaluate energy/forces with
    type(Atoms), intent(inout), target :: at !% starting configuration
    character(*), intent(in)    :: method !% passed to minim()
    real(dp),     intent(in)    :: convergence_tol !% Minimisation is treated as converged once $|\mathbf{\nabla}f|^2 <$
                                                    !% 'convergence_tol'. 
    integer,      intent(in)    :: max_steps  !% Maximum number of steps
    character(*), intent(in),optional    :: linminroutine !% Name of the line minisation routine to use, passed to base minim()
    logical, optional :: do_print !% if true, print configurations using minim's hook()
    type(inoutput), intent(inout), optional, target :: print_inoutput !% inoutput object to print configs to, needed if do_print is true
    type(cinoutput), intent(inout), optional, target :: print_cinoutput !% cinoutput object to print configs to, needed if do_print is true
    logical, optional :: do_pos, do_lat !% do relaxation w.r.t. positions and/or lattice (if neither is included, do both)
    character(len=*), intent(in), optional :: args_str !% arguments to pass to calc()
    real(dp), intent(in), optional :: eps_guess !% eps_guess argument to pass to minim
    real(dp), intent(in), optional :: fire_minim_dt0 !% if using fire minim, initial value for time step
    real(dp), intent(in), optional :: fire_minim_dt_max !% if using fire minim, max value for time step
    real(dp), dimension(3,3), optional :: external_pressure
    logical, intent(in), optional :: use_precond
    integer, intent(in), optional :: hook_print_interval !% how often to print xyz from hook function
    integer, intent(out), optional :: error !% set to 1 if an error occurred during minimisation
    integer::potential_minim

    character(len=100) :: use_method

    logical :: my_use_precond
    integer n_iter, n_iter_tot
    real(dp), allocatable :: x(:)
    real(dp) :: deform_grad(3,3)
    logical my_do_print
    logical done
    real(dp) :: my_eps_guess
    real(dp) :: my_fire_minim_dt0, my_fire_minim_dt_max

    real(dp) :: initial_E, final_E, mass
    type(potential_minimise) am
    integer :: am_data_size
    character, allocatable :: am_data(:)
    character, dimension(1) :: am_mold
    integer :: status

    INIT_ERROR(error)

    my_use_precond = optional_default(.false., use_precond)

    am_data_size = size(transfer(am, am_mold))
    allocate(am_data(am_data_size))

    use_method = trim(method)

    my_fire_minim_dt0 = optional_default(1.0_dp, fire_minim_dt0)
    my_fire_minim_dt_max = optional_default(20.0_dp, fire_minim_dt_max)

    my_eps_guess = optional_default(1.0e-2_dp/at%N, eps_guess)
    if (my_eps_guess .feq. 0.0_dp) my_eps_guess = 1.0e-2_dp/at%N

    call calc_connect(at)
    am%minim_at => at
    am%pos_lat_preconditioner_factor = am%minim_pos_lat_preconditioner*am%minim_at%N

    if (present(args_str)) then
      am%minim_args_str = args_str
    else
      am%minim_args_str = ""
    endif
    am%minim_pot => this

    if (.not.present(do_pos) .and. .not. present(do_lat)) then
      am%minim_do_pos = .true.
      am%minim_do_lat = .true.
    else
      am%minim_do_pos = optional_default(.false., do_pos)
      am%minim_do_lat = optional_default(.false., do_lat)
    endif

    am%external_pressure = 0.0_dp
    if (present(external_pressure)) then
       am%external_pressure = external_pressure
       if( (am%external_pressure(1,1) .fne. am%external_pressure(2,2)) .or. &
       & (am%external_pressure(1,1) .fne. am%external_pressure(3,3)) .or. &
       & (am%external_pressure(1,2) .fne. 0.0_dp) .or. &
       & (am%external_pressure(1,3) .fne. 0.0_dp) .or. &
       & (am%external_pressure(2,1) .fne. 0.0_dp) .or. &
       & (am%external_pressure(2,3) .fne. 0.0_dp) .or. &
       & (am%external_pressure(3,1) .fne. 0.0_dp) .or. &
       & (am%external_pressure(3,2) .fne. 0.0_dp) ) then
          if(trim(use_method) /= 'fire') then
             call print_warning('Anisotrpic pressure is being used. Switching to fire_minim.')
             use_method = 'fire'
          endif
       endif
    endif

    my_do_print = optional_default(.false., do_print)

    if (my_do_print .and. .not. present(print_inoutput) .and. .not. present(print_cinoutput)) &
         call system_abort("potential_minim: do_print is true, but no print_inoutput or print_cinoutput present")

    if (my_do_print) then
       if (present(print_cinoutput)) then
          am%minim_cinoutput_movie => print_cinoutput
          if (this%is_forcemixing) &
               this%forcemixing%minim_cinoutput_movie => print_cinoutput
#ifdef HAVE_LOCAL_E_MIX
          if (this%is_local_e_mix) &
               this%local_e_mix%minim_cinoutput_movie => print_cinoutput
#endif /* HAVE_LOCAL_E_MIX */
#ifdef HAVE_ONIOM
          if (this%is_oniom) &
               this%oniom%minim_cinoutput_movie => print_cinoutput
#endif /* HAVE_ONIOM */
       end if
    else
      nullify(am%minim_cinoutput_movie)
    endif

    am%minim_n_eval_e = 0
    am%minim_n_eval_f = 0

    allocate(x(9+am%minim_at%N*3))
    allocate(am%last_connect_x(size(x)))
    am%last_connect_x=1.0e38_dp
    deform_grad = 0.0_dp; call add_identity(deform_grad)
    call pack_pos_dg(am%minim_at%pos, deform_grad, x, am%pos_lat_preconditioner_factor)

    am_data = transfer(am, am_data)
    if (trim(use_method) == 'cg_n') then
       n_iter = n_minim(x, both_func, my_use_precond, apply_precond_func, initial_E, final_E, my_eps_guess, max_steps, convergence_tol, print_hook, &
            hook_print_interval=hook_print_interval, data=am_data, error=error)
       PASS_ERROR(error)
    else if (trim(use_method) == 'fire') then
       if (has_property(at, 'mass')) then
          mass = at%mass(1)
       else
          mass = ElementMass(at%Z(1))
       end if
       n_iter = fire_minim(x, mass, dummy_energy_func, gradient_func, my_fire_minim_dt0, convergence_tol, max_steps, &
            print_hook, hook_print_interval=hook_print_interval, data=am_data, dt_max=my_fire_minim_dt_max, status=status)
    else
       n_iter = minim(x, energy_func, gradient_func, use_method, convergence_tol, max_steps, linminroutine, &
            print_hook, hook_print_interval=hook_print_interval, eps_guess=my_eps_guess, data=am_data, status=status)
    endif
    call print("minim relax w.r.t. both n_iter " // n_iter, PRINT_VERBOSE)
    call unpack_pos_dg(x, am%minim_at%N, am%minim_at%pos, deform_grad, 1.0_dp/am%pos_lat_preconditioner_factor)
    call prep_atoms_deform_grad(deform_grad, am%minim_at, am)
    call calc_connect(am%minim_at)
    n_iter_tot = n_iter
    done = .true.
    deallocate(am%last_connect_x)
    deallocate(x)
    call print("MINIM_N_EVAL E " // am%minim_n_eval_e // " F " // am%minim_n_eval_f // &
      " EF " // am%minim_n_eval_ef, PRINT_VERBOSE)

    potential_minim = n_iter_tot

    deallocate(am_data)

  end function potential_minim

  subroutine potential_test_local_virial(this, at, args_str)
    type(Atoms), intent(inout), target :: at
    type(Potential), target, intent(inout) :: this
    character(len=*), intent(in), optional :: args_str

    real(dp) :: dx, r_ik
    real(dp), dimension(3) :: diff_ik
    real(dp), dimension(3,3) :: virial
    real(dp), dimension(:), allocatable :: local_energy_p, local_energy_m
    real(dp), dimension(:,:), allocatable :: tmp_local_virial, pos0
    real(dp), dimension(:,:,:), allocatable :: local_virial, local_virial_fd

    integer :: d, i, k, n, alpha

    call print_warning("TLV: potential_test_local_virial: your cell has to be big enough so no multiple images are used for any atom.")

    if(at%cutoff < cutoff(this)) call set_cutoff(at,cutoff(this))

    call calc_connect(at)

    allocate(tmp_local_virial(9,at%N), local_virial(3,3,at%N), local_virial_fd(3,3,at%N), &
       local_energy_p(at%N), local_energy_m(at%N), pos0(3,at%N))

    call calc(this,at,local_virial = tmp_local_virial, virial = virial)
    local_virial = reshape(tmp_local_virial,(/3,3,at%N/))

    call print("TLV: RMS difference between virial and local virial summed over atoms = "//sqrt( sum( (virial - sum(local_virial,dim=3))**2 ) ))

    pos0 = at%pos

    dx = 1.0_dp

    print"(a,a16,a20)", "TLV", "dx", "RMS"
    do d = 1, 10
       dx = dx * 0.1_dp
       local_virial_fd = 0.0_dp

       do i = 1, at%N

          do alpha = 1, 3
             at%pos(alpha,i) = pos0(alpha,i) + dx
             call calc_connect(at)
             call calc(this,at,local_energy=local_energy_p,args_str=args_str)

             at%pos(alpha,i) = pos0(alpha,i) - dx
             call calc_connect(at)
             call calc(this,at,local_energy=local_energy_m,args_str=args_str)

             at%pos(alpha,i) = pos0(alpha,i) 
             call calc_connect(at)

             do n = 1, n_neighbours(at,i)
                k = neighbour(at,i,n,distance = r_ik, diff = diff_ik)

                if( r_ik > cutoff(this) ) cycle

                local_virial_fd(:,alpha,i) = local_virial_fd(:,alpha,i) + diff_ik * ( local_energy_p(k) - local_energy_m(k) ) / 2.0_dp / dx
             enddo
          enddo
       enddo

       print"(a,f16.10,e20.10)", "TLV", dx, sqrt(sum((local_virial_fd - local_virial)**2))
    enddo

  endsubroutine potential_test_local_virial

  subroutine undo_travel(at)
    type(Atoms), intent(inout) :: at

!    if (any(at%travel /= 0)) then
!      at%pos = at%pos + (at%lattice .mult. at%travel)
!      at%travel = 0
!    endif
end subroutine undo_travel

  ! test the gradient given a potential, string to pass to calc, and atoms structure
  function pot_test_gradient(pot, at, do_pos, do_lat, args_str, dir_field)
    type(Atoms), intent(inout), target :: at
    type(Potential), target, intent(inout) :: pot
    logical, intent(in), optional :: do_pos, do_lat
    character(len=*), intent(in), optional :: args_str
    character(len=*), intent(in), optional :: dir_field
    logical pot_test_gradient

    real(dp) :: dg(3,3)
    real(dp), allocatable :: x(:), dir_x(:)
    real(dp), pointer :: dir_p(:,:)
    integer :: am_data_size
    type(potential_minimise) :: am
    character, allocatable :: am_data(:)
    character, dimension(1) :: am_mold

    am_data_size = size(transfer(am, am_mold))
    allocate(am_data(am_data_size))

    pot_test_gradient = .true.

    if (present(args_str)) then
      am%minim_args_str = args_str
    else
      am%minim_args_str = ""
    endif
    am%minim_pot => pot

    ! neither do_pos not do_lat is specified, so do both by default
    if (.not. present(do_pos) .and. .not. present(do_lat)) then
      am%minim_do_pos= .true.
      am%minim_do_lat= .true.
    else ! something is specified, so do exactly that
      am%minim_do_pos = optional_default(.false., do_pos)
      am%minim_do_lat = optional_default(.false., do_lat)
    endif

    am%pos_lat_preconditioner_factor = am%minim_pos_lat_preconditioner*at%N

    if (present(dir_field)) then
      if (.not. assign_pointer(at, trim(dir_field), dir_p)) &
	call system_abort("pot_test_gradient got dir_field '"//trim(dir_field)//"' but can't assign pointer to it")
      call print("pot_test_gradient using '"//trim(dir_field)//"' field for gradient test direction")
    else
      nullify(dir_p)
    endif

    am%minim_at => at
    allocate(x(at%N*3+9))
    allocate(am%last_connect_x(at%N*3+9))
    am%last_connect_x=1.0e38_dp
    dg = 0.0_dp; call add_identity(dg)
    call pack_pos_dg(at%pos, dg, x, am%pos_lat_preconditioner_factor)
    if (associated(dir_p)) then
      allocate(dir_x(at%N*3+9))
      dg = 0.0_dp
      call pack_pos_dg(dir_p, dg, dir_x, am%pos_lat_preconditioner_factor)
    endif
    am_data = transfer(am, am_data)
    if (associated(dir_p)) then
      pot_test_gradient = test_gradient(x, energy_func, gradient_func, dir=dir_x, data=am_data)
    else
      pot_test_gradient = test_gradient(x, energy_func, gradient_func, data=am_data)
    endif
    deallocate(x)
    if (allocated(dir_x)) deallocate(dir_x)
    deallocate(am%last_connect_x)

    deallocate(am_data)

  end function pot_test_gradient

  ! test the gradient given a potential, string to pass to calc, and atoms structure
  subroutine pot_n_test_gradient(pot, at, do_pos, do_lat, args_str, dir_field)
    type(Potential), target, intent(inout) :: pot
    type(Atoms), intent(inout), target :: at
    logical, intent(in), optional :: do_pos, do_lat
    character(len=*), intent(in), optional :: args_str
    character(len=*), intent(in), optional :: dir_field

    real(dp), allocatable :: x(:), dir_x(:)
    real(dp), pointer :: dir_p(:,:)
    real(dp) :: dg(3,3)
    integer :: am_data_size
    type(potential_minimise) :: am
    character, allocatable :: am_data(:)
    character, dimension(1) :: am_mold

    am_data_size = size(transfer(am, am_mold))
    allocate(am_data(am_data_size))

    if (present(args_str)) then
      am%minim_args_str = args_str
    else
      am%minim_args_str = ""
    endif
    am%minim_pot => pot

    ! neither do_pos not do_lat is specified, so do both by default
    if (.not. present(do_pos) .and. .not. present(do_lat)) then
      am%minim_do_pos= .true.
      am%minim_do_lat= .true.
    else ! something is specified, so do exactly that
      am%minim_do_pos = optional_default(.false., do_pos)
      am%minim_do_lat = optional_default(.false., do_lat)
    endif

    am%pos_lat_preconditioner_factor = am%minim_pos_lat_preconditioner*at%N

    if (present(dir_field)) then
      if (.not. assign_pointer(at, trim(dir_field), dir_p)) &
	call system_abort("pot_test_gradient got dir_field '"//trim(dir_field)//"' but can't assign pointer to it")
      call print("pot_test_gradient using '"//trim(dir_field)//"' field for gradient test direction")
    else
      nullify(dir_p)
    endif

    am%minim_at => at
    allocate(x(at%N*3+9))
    allocate(am%last_connect_x(at%N*3+9))
    am%last_connect_x=1.0e38_dp
    dg = 0.0_dp; call add_identity(dg)
    call pack_pos_dg(at%pos, dg, x, am%pos_lat_preconditioner_factor)
    if (associated(dir_p)) then
      allocate(dir_x(at%N*3+9))
      dg = 0.0_dp
      call pack_pos_dg(dir_p, dg, dir_x, am%pos_lat_preconditioner_factor)
    endif
    am_data = transfer(am, am_data)
    if (associated(dir_p)) then
      call n_test_gradient(x, energy_func, gradient_func, data=am_data, dir=dir_x)
    else
      call n_test_gradient(x, energy_func, gradient_func, data=am_data)
    endif
    deallocate(x)
    if (allocated(dir_x)) deallocate(dir_x)
    deallocate(am%last_connect_x)

    deallocate(am_data)

  end subroutine pot_n_test_gradient

  ! hook for minim to print configuration during position relaxation
  subroutine print_hook(x,dx,E,done,do_print,am_data)
    real(dp), intent(in) :: x(:), dx(:)
    real(dp), intent(in) :: E
    logical, intent(out) :: done
    logical, optional, intent(in) :: do_print
    character, intent(in), optional :: am_data(:)

    real(dp), allocatable :: f(:,:)
    real(dp), pointer :: f_p(:,:)
    integer, pointer, dimension(:) :: hybrid_mark
    real(dp) :: virial(3,3), deform_grad(3,3), deform_grad_inv(3,3)
    type(potential_minimise) :: am
    logical :: my_do_print

    if (.not. present(am_data)) call system_abort("potential_minimise print_hook must have am_data")
    my_do_print = optional_default(.true., do_print)

    am = transfer(am_data, am)
    ! beware of transfer and pointers !!!
    call atoms_repoint(am%minim_at)

    if (associated(am%minim_cinoutput_movie)) then
      if (size(x) /= am%minim_at%N*3+9) call system_abort("Called gradient_func() with size mismatch " // &
        size(x) // " /= " // am%minim_at%N // "*3+9")

      call unpack_pos_dg(x, am%minim_at%N, am%minim_at%pos, deform_grad, 1.0_dp/am%pos_lat_preconditioner_factor)
      call prep_atoms_deform_grad(deform_grad, am%minim_at, am)

      allocate(f(3,am%minim_at%N))
      f = reshape(dx(10:), (/3,am%minim_at%N/))
      f = transpose(deform_grad) .mult. f

      virial = reshape(dx(1:9), (/3,3/) )
      call inverse(deform_grad, deform_grad_inv)
      virial = virial .mult. transpose(deform_grad_inv)

      ! Set atom parameters to print energy, maxforce and maxvirial
      call set_value(am%minim_at%params, 'E', E)
      if (am%minim_do_pos) then
         call set_value(am%minim_at%params, 'MaxForce', maxval(norm(f,1)))
         call set_value(am%minim_at%params, 'df2', normsq(reshape(f,(/3*am%minim_at%N/))))


         if (am%minim_pot%is_forcemixing) then
	    !NB This should really be an arbitrary field name, but no obvious way to pass print_hook info on what 
	    !NB the name is (no run_suffix args_str argument, for example)
            if (assign_pointer(am%minim_at, 'hybrid_mark', hybrid_mark)) then 
               if (.not. am%minim_pot%forcemixing%minimise_mm) then
                  call set_value(am%minim_at%params, 'QM_MaxForce', &
                       maxval(abs(f(:,find(hybrid_mark == HYBRID_ACTIVE_MARK)))))
            
                  call set_value(am%minim_at%params, 'QM_df2', &
                       normsq(reshape(f(:,find(hybrid_mark == HYBRID_ACTIVE_MARK)),&
                       (/3*count(hybrid_mark == HYBRID_ACTIVE_MARK)/))))

               end if
            end if
         end if

      end if

      if (am%minim_do_lat) &
	call set_value(am%minim_at%params, 'MaxVirial', maxval(abs(virial)))

      if (assign_pointer(am%minim_at, "forces", f_p)) f_p = -f
      if (assign_pointer(am%minim_at, "force", f_p)) f_p = -f ! compatibility with crack program

      if (my_do_print) then
         if (associated(am%minim_cinoutput_movie)) then
	    if (trim(am%minim_cinoutput_movie%filename) == "stdout") then
	       if (.not. am%minim_pot%mpi%active .or. (am%minim_pot%mpi%active .and. am%minim_pot%mpi%my_proc == 0)) &
		  call write(am%minim_cinoutput_movie, am%minim_at, prefix="MINIM_TRAJ")
	    else
	       if (.not. am%minim_pot%mpi%active .or. (am%minim_pot%mpi%active .and. am%minim_pot%mpi%my_proc == 0)) &
		  call write(am%minim_cinoutput_movie, am%minim_at)
	    end if
         end if
      end if

      deallocate(f)
      call fix_atoms_deform_grad(deform_grad, am%minim_at,am)
    endif
    done = .false.

  end subroutine print_hook

  function dummy_energy_func(xx, am_data)
    real(dp) :: xx(:)
    character(len=1), optional :: am_data(:)
    real(dp) :: dummy_energy_func

    dummy_energy_func = 0.0_dp

  end function dummy_energy_func


  ! compute energy
  function energy_func(x, am_data)
    real(dp) :: x(:)
    character(len=1), optional :: am_data(:)
    real(dp) :: energy_func

    real(dp) :: max_atom_rij_change
    real(dp) :: deform_grad(3,3)
    type(potential_minimise)  :: am
    real(dp), pointer :: minim_applied_force(:,:)

    call system_timer("energy_func")

    if (.not. present(am_data)) call system_abort("potential_minimise energy_func must have am_data")
    am = transfer(am_data, am)
    call atoms_repoint(am%minim_at)

    am%minim_n_eval_e = am%minim_n_eval_e + 1

    if (size(x) /= am%minim_at%N*3+9) call system_abort("Called energy_func() with size mismatch " // &
      size(x) // " /= " // am%minim_at%N // "*3+9")

    ! Note: TB will return 0 for cutoff(am%minim_pot), but TB does its own calc_connect, so doesn't matter
    max_atom_rij_change = max_rij_change(am%last_connect_x, x, cutoff(am%minim_pot), &
      1.0_dp/am%pos_lat_preconditioner_factor)

    call print("energy_func got x " // x, PRINT_NERD)

    call unpack_pos_dg(x, am%minim_at%N, am%minim_at%pos, deform_grad, 1.0_dp/am%pos_lat_preconditioner_factor)
    call prep_atoms_deform_grad(deform_grad, am%minim_at, am)

    ! Safety factor of 1.1, just in case
    ! Note: TB will return 0 for cutoff(am%minim_pot), but TB does its own calc_connect, so doesn't matter
    if (1.1*max_atom_rij_change >= am%minim_at%cutoff - cutoff(am%minim_pot)) then
      call print("gradient_func: Do calc_connect, atoms moved " // max_atom_rij_change // &
        "*1.1 >= buffer " // (am%minim_at%cutoff - cutoff(am%minim_pot)), PRINT_NERD)
      call calc_connect(am%minim_at)
      am%last_connect_x = x
    else
      call print("gradient_func: Do calc_dists, atoms moved " // max_atom_rij_change // &
        " *1.1 < buffer " // (am%minim_at%cutoff - cutoff(am%minim_pot)), PRINT_NERD)
      call calc_dists(am%minim_at)
    end if

    call print("energy_func using am%minim_at", PRINT_NERD)
    if (current_verbosity() >= PRINT_NERD) call write(am%minim_at, 'stdout')

    call calc(am%minim_pot, am%minim_at, energy = energy_func, args_str = am%minim_args_str)
    call print ("energy_func got energy " // energy_func, PRINT_NERD)

    energy_func = energy_func + cell_volume(am%minim_at)*trace(am%external_pressure) / 3.0_dp
    call print ("energy_func got enthalpy " // energy_func, PRINT_NERD)

    ! add extra forces
    if (assign_pointer(am%minim_at, "minim_applied_force", minim_applied_force)) then
      energy_func = energy_func - sum(minim_applied_force*am%minim_at%pos)
    end if

    call fix_atoms_deform_grad(deform_grad, am%minim_at, am)
    call pack_pos_dg(am%minim_at%pos, deform_grad, x, am%pos_lat_preconditioner_factor)

    am_data = transfer(am, am_data)

    call system_timer("energy_func")

  end function energy_func

  ! compute gradient (forces and virial)
  ! x is vectorized version of atomic positions
  ! result is vectorized version of forces
  function gradient_func(x, am_data)
    real(dp) :: x(:)
    character(len=1), optional :: am_data(:)
    real(dp) :: gradient_func(size(x))

    real(dp) :: deform_grad(3,3), virial(3,3), deform_grad_inv(3,3)
    real(dp), allocatable :: f(:,:)
    real(dp) :: max_atom_rij_change
    integer :: t_i(2), i
    integer, pointer, dimension(:) :: move_mask, fixed_pot
    real(dp), pointer :: minim_applied_force(:,:)

    type(potential_minimise) :: am

    call system_timer("gradient_func")

    if (.not. present(am_data)) call system_abort("potential_minimise gradient_func must have am_data")
    am = transfer(am_data, am)
    ! beware of transfer and pointers !!!
    call atoms_repoint(am%minim_at)

    am%minim_n_eval_e = am%minim_n_eval_e + 1

    am%minim_n_eval_f = am%minim_n_eval_f + 1

    if (size(x) /= am%minim_at%N*3+9) call system_abort("Called gradient_func() with size mismatch " // &
      size(x) // " /= " // am%minim_at%N // "*3+9")

    ! Note: TB will return 0 for cutoff(am%minim_pot), but TB does its own calc_connect, so doesn't matter
    max_atom_rij_change = max_rij_change(am%last_connect_x, x, cutoff(am%minim_pot), &
      1.0_dp/am%pos_lat_preconditioner_factor)

    if (current_verbosity() >= PRINT_NERD) call print("gradient_func got x " // x, PRINT_NERD)

    call unpack_pos_dg(x, am%minim_at%N, am%minim_at%pos, deform_grad, 1.0_dp/am%pos_lat_preconditioner_factor)
    call prep_atoms_deform_grad(deform_grad, am%minim_at, am)

    ! Safety factor of 1.1, just in case
    ! Note: TB will return 0 for cutoff(am%minim_pot), but TB does its own calc_connect, so doesn't matter
    if (1.1*max_atom_rij_change >= am%minim_at%cutoff - cutoff(am%minim_pot)) then
      call print("gradient_func: Do calc_connect, atoms moved " // max_atom_rij_change // &
        "*1.1 >= buffer " // (am%minim_at%cutoff - cutoff(am%minim_pot)), PRINT_NERD)
      call calc_connect(am%minim_at)
      am%last_connect_x = x
    else
      call print("gradient_func: Do calc_dists, atoms moved " // max_atom_rij_change // &
        " *1.1 < buffer " // (am%minim_at%cutoff - cutoff(am%minim_pot)), PRINT_NERD)
      call calc_dists(am%minim_at)
    end if

    if (current_verbosity() >= PRINT_NERD) call add_property(am%minim_at, "force", 0.0_dp, 3)

    allocate(f(3,am%minim_at%N))
    f = 0.0_dp
    virial = 0.0_dp
    if (am%minim_do_pos .and. am%minim_do_lat) then
      call calc(am%minim_pot, am%minim_at, force = f, virial = virial, args_str = am%minim_args_str)
    else if (am%minim_do_pos) then
      call calc(am%minim_pot, am%minim_at, force = f, args_str = am%minim_args_str)
    else
      call calc(am%minim_pot, am%minim_at, virial = virial, args_str = am%minim_args_str)
    endif

    call print("gradient_func used am%minim_at, got forces", PRINT_NERD)
    if (current_verbosity() >= PRINT_NERD) call write(am%minim_at,'stdout')

    ! zero forces if fixed by potential
    if (am%minim_do_pos .and. assign_pointer(am%minim_at, "fixed_pot", fixed_pot)) then
      do i=1, am%minim_at%N
	if (fixed_pot(i) /= 0) f(:,i) = 0.0_dp
      end do
    endif

    ! Zero force on any fixed atoms
    if (assign_pointer(am%minim_at, 'move_mask', move_mask)) then
       do i=1,am%minim_at%N
          if (move_mask(i) == 0) f(:,i) = 0.0_dp
       end do
    end if

    ! add extra forces
    if (assign_pointer(am%minim_at, "minim_applied_force", minim_applied_force)) then
      f = f + minim_applied_force
      ! val = val - sum(minim_applied_force*am%minim_at%pos)
    endif

    if (current_verbosity() >= PRINT_NERD) then
       call print ("gradient_func got f", PRINT_NERD)
       call print(f, PRINT_NERD)
       call print ("gradient_func got virial", PRINT_NERD)
       call print(virial, PRINT_NERD)
    end if
    virial = virial - am%external_pressure*cell_volume(am%minim_at)
    
    if (current_verbosity() >= PRINT_NERD) then
       call print ("gradient_func got virial, external pressure subtracted", PRINT_NERD)
       call print(virial, PRINT_NERD)
    end if

    f = transpose(deform_grad) .mult. f

    call inverse(deform_grad, deform_grad_inv)
    virial = virial .mult. transpose(deform_grad_inv)

    call constrain_virial_post(am%minim_at, virial)

    call pack_pos_dg(-f, -virial, gradient_func, 1.0_dp/am%pos_lat_preconditioner_factor)
    if (current_verbosity() >= PRINT_NERD) then
      call print ("gradient_func packed as", PRINT_NERD)
      call print (gradient_func, PRINT_NERD)
    endif

    t_i = maxloc(abs(f))
    call print ("gradient_func got max force component " // f(t_i(1),t_i(2)) // " on atom " // t_i(2) // &
      " max virial component " // maxval(abs(virial)))

    deallocate(f)

    call fix_atoms_deform_grad(deform_grad, am%minim_at,am)
    call pack_pos_dg(am%minim_at%pos, deform_grad, x, am%pos_lat_preconditioner_factor)
    am_data = transfer(am, am_data)

    call system_timer("gradient_func")

  end function gradient_func

! Cristoph Ortner's simple Hessian preconditioner, dense matrix inverse right now
  subroutine apply_precond_func(x,g,P_g,am_data,error)
    real(dp) :: x(:), g(:), P_g(:)
    character(len=1), optional :: am_data(:)
    integer, optional :: error

    type(potential_minimise) :: am
    real(dp) :: deform_grad(3,3)
    real(dp), allocatable :: P(:,:)
    real(dp) :: c0, c1
    integer :: i, jj, j
    real(dp) :: max_atom_rij_change
    integer :: n_nearest
    integer, pointer :: move_mask(:)

    INIT_ERROR(error)

    am = transfer(am_data, am)
    call atoms_repoint(am%minim_at)

    max_atom_rij_change = max_rij_change(am%last_connect_x, x, cutoff(am%minim_pot), &
      1.0_dp/am%pos_lat_preconditioner_factor)
    !max_atom_rij_change = 1.038_dp !FIXME this line shouldn't be here either

    call print("apply_precond_func got x " // x, PRINT_NERD)

    call unpack_pos_dg(x, am%minim_at%N, am%minim_at%pos, deform_grad, 1.0_dp/am%pos_lat_preconditioner_factor)
    call prep_atoms_deform_grad(deform_grad, am%minim_at, am)

    ! Safety factor of 1.1, just in case
    ! Note: TB will return 0 for cutoff(am%minim_pot), but TB does its own calc_connect, so doesn't matter
    if (1.1*max_atom_rij_change >= am%minim_at%cutoff - cutoff(am%minim_pot)) then
      call print("apply_precond_func: Do calc_connect, atoms moved " // max_atom_rij_change // "*1.1 >= buffer " // &
        (am%minim_at%cutoff - cutoff(am%minim_pot)), PRINT_NERD)
      call calc_connect(am%minim_at)
      am%last_connect_x = x
    else
      call print("apply_precond_func: Do calc_dists, atoms moved " // max_atom_rij_change // " *1.1 < buffer " // &
        (am%minim_at%cutoff - cutoff(am%minim_pot)), PRINT_NERD)
      call calc_dists(am%minim_at)
    end if

    ! identity
    ! c0 = 1.0_dp
    ! c1 = 0.0_dp
    ! simplest expression of connectivity
    c0 = 0.01_dp
    c1 = 1.0_dp

    allocate(P(3*am%minim_at%N,3*am%minim_at%N))
    P = 0.0_dp
    do i=1, am%minim_at%N
      n_nearest = 0
      do jj=1, n_neighbours(am%minim_at, i)
	if (is_nearest_neighbour(am%minim_at, i, jj)) then
	  n_nearest = n_nearest + 1
	  j = neighbour(am%minim_at, i, jj)
	  P(3*(i-1)+1,3*(j-1)+1) = -c1
	  P(3*(i-1)+2,3*(j-1)+2) = -c1
	  P(3*(i-1)+3,3*(j-1)+3) = -c1
	endif
      end do
      P(3*(i-1)+1,3*(i-1)+1) = c0+c1*n_nearest
      P(3*(i-1)+2,3*(i-1)+2) = c0+c1*n_nearest
      P(3*(i-1)+3,3*(i-1)+3) = c0+c1*n_nearest
    end do

    if (assign_pointer(am%minim_at, 'move_mask', move_mask)) then
       do i=1,am%minim_at%N
          if (move_mask(i) == 0) then
	    P(3*(i-1)+1:3*(i-1)+3,:) = 0.0_dp
	    P(:,3*(i-1)+1:3*(i-1)+3) = 0.0_dp
	    P(3*(i-1)+1,3*(i-1)+1) = 1.0_dp
	    P(3*(i-1)+2,3*(i-1)+2) = 1.0_dp
	    P(3*(i-1)+3,3*(i-1)+3) = 1.0_dp
	  endif
       end do
    endif

    P_g(1:9) = g(1:9)
    call symmetric_linear_solve(P, g(10:10+am%minim_at%N*3-1), P_g(10:10+am%minim_at%N*3-1))

  end subroutine apply_precond_func
  ! compute energy and gradient (forces and virial)
  ! x is vectorized version of atomic positions
  ! result is vectorized version of forces
  subroutine both_func(x,val,grad,am_data,error)
    real(dp) :: x(:)
    real(dp) :: val
    real(dp) :: grad(:)
    character(len=1), optional :: am_data(:)
    integer, optional :: error

    real(dp) :: deform_grad(3,3), virial(3,3), deform_grad_inv(3,3)
    real(dp), allocatable :: f(:,:)
    real(dp) :: max_atom_rij_change
    integer :: t_i(2), i
    integer, pointer, dimension(:) :: move_mask, fixed_pot
    real(dp), pointer :: minim_applied_force(:,:)

    type(potential_minimise) :: am

    real(dp) :: hack_restraint_E, hack_restraint_F(3), dr

    INIT_ERROR(error)

    call system_timer("both_func")

    if (.not. present(am_data)) call system_abort("potential_minimise both_func must have am_data")
    am = transfer(am_data, am)
    ! beware of transfer and pointers !!!
    call atoms_repoint(am%minim_at)

    am%minim_n_eval_ef = am%minim_n_eval_ef + 1

    if (size(x) /= am%minim_at%N*3+9) call system_abort("Called both_func() with size mismatch " // &
      size(x) // " /= " // am%minim_at%N // "*3+9")

    ! Note: TB will return 0 for cutoff(am%minim_pot), but TB does its own calc_connect, so doesn't matter
    max_atom_rij_change = max_rij_change(am%last_connect_x, x, cutoff(am%minim_pot), &
      1.0_dp/am%pos_lat_preconditioner_factor)
    !max_atom_rij_change = 1.038_dp ! FIXME this line shouldn't be here, right?

    call print("both_func got x " // x, PRINT_NERD)

    call unpack_pos_dg(x, am%minim_at%N, am%minim_at%pos, deform_grad, 1.0_dp/am%pos_lat_preconditioner_factor)
    call prep_atoms_deform_grad(deform_grad, am%minim_at, am)

    ! Safety factor of 1.1, just in case
    ! Note: TB will return 0 for cutoff(am%minim_pot), but TB does its own calc_connect, so doesn't matter
    if (1.1*max_atom_rij_change >= am%minim_at%cutoff - cutoff(am%minim_pot)) then
      call print("both_func: Do calc_connect, atoms moved " // max_atom_rij_change // "*1.1 >= buffer " // &
        (am%minim_at%cutoff - cutoff(am%minim_pot)), PRINT_NERD)
      call calc_connect(am%minim_at)
      am%last_connect_x = x
    else
      call print("both_func: Do calc_dists, atoms moved " // max_atom_rij_change // " *1.1 < buffer " // &
        (am%minim_at%cutoff - cutoff(am%minim_pot)), PRINT_NERD)
      call calc_dists(am%minim_at)
    end if

    call print("both_func using am%minim_at", PRINT_NERD)
    if (current_verbosity() >= PRINT_NERD) call write(am%minim_at,'stdout')

    allocate(f(3,am%minim_at%N))
    f = 0.0_dp
    virial = 0.0_dp
    if (am%minim_do_pos .and. am%minim_do_lat) then
      call calc(am%minim_pot, am%minim_at, energy = val, force = f, virial = virial, args_str = am%minim_args_str)
    else if (am%minim_do_pos) then
      call calc(am%minim_pot, am%minim_at, energy = val, force = f, args_str = am%minim_args_str)
    else
      call calc(am%minim_pot, am%minim_at, energy = val, virial = virial, args_str = am%minim_args_str)
    endif

if (am%minim_do_pos) then
   if (all(hack_restraint_i > 0) .and. all(hack_restraint_i <= am%minim_at%n)) then
      call print("hack_restraint i "//hack_restraint_i// " k " // hack_restraint_k // " r " // hack_restraint_r, PRINT_ALWAYS)
      dr = distance_min_image(am%minim_at, hack_restraint_i(1), hack_restraint_i(2)) 
      hack_restraint_E = 0.5_dp * hack_restraint_k*(dr - hack_restraint_r)**2
      val = val + hack_restraint_E
      hack_restraint_F = hack_restraint_k*(dr - hack_restraint_r)*diff_min_image(am%minim_at, hack_restraint_i(1), hack_restraint_i(2))/dr
      call print("hack_restraint dr-r "//(dr-hack_restraint_r)//"E "//hack_restraint_E// " F " // hack_restraint_F, PRINT_ALWAYS)
      f(:,hack_restraint_i(1)) = f(:,hack_restraint_i(1)) + hack_restraint_F
      f(:,hack_restraint_i(2)) = f(:,hack_restraint_i(2)) - hack_restraint_F
   endif
endif

    if (assign_pointer(am%minim_at, "minim_applied_force", minim_applied_force)) then
      f = f + minim_applied_force
      val = val - sum(minim_applied_force*am%minim_at%pos)
    endif

    call print ("both_func got f", PRINT_NERD)
    call print(f, PRINT_NERD)
    call print ("both_func got virial", PRINT_NERD)
    call print(virial, PRINT_NERD)

    virial = virial - am%external_pressure*cell_volume(am%minim_at)
    call print ("both_func got virial, external pressure subtracted", PRINT_NERD)
    call print(virial, PRINT_NERD)

    val = val + cell_volume(am%minim_at)*trace(am%external_pressure) / 3.0_dp
    call print ("both_func got enthalpy " // val, PRINT_NERD)
    ! zero forces if fixed by potential
    if (am%minim_do_pos .and. assign_pointer(am%minim_at, "fixed_pot", fixed_pot)) then
      do i=1, am%minim_at%N
	if (fixed_pot(i) /= 0) f(:,i) = 0.0_dp
      end do
    endif

    ! Zero force on any fixed atoms
    if (assign_pointer(am%minim_at, 'move_mask', move_mask)) then
       do i=1,am%minim_at%N
          if (move_mask(i) == 0) f(:,i) = 0.0_dp
       end do
    end if

    f = transpose(deform_grad) .mult. f

    call inverse(deform_grad, deform_grad_inv)
    virial = virial .mult. transpose(deform_grad_inv)

    call constrain_virial_post(am%minim_at, virial)

    call pack_pos_dg(-f, -virial, grad, 1.0_dp/am%pos_lat_preconditioner_factor)
    call print ("both_func gradient packed as", PRINT_NERD)
    call print (grad, PRINT_NERD)

    t_i = maxloc(abs(f))
    call print ("both_func got max force component " // f(t_i(1),t_i(2)) // " on atom " // t_i(2) // &
      " max virial component " // maxval(abs(virial)))

    deallocate(f)

    call fix_atoms_deform_grad(deform_grad, am%minim_at,am)
    call pack_pos_dg(am%minim_at%pos, deform_grad, x, am%pos_lat_preconditioner_factor)

    am_data = transfer(am, am_data)

    call system_timer("both_func")

  end subroutine both_func

  function max_rij_change(last_connect_x, x, r_cut, lat_factor)
    real(dp), intent(inout) :: last_connect_x(:)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in) :: r_cut, lat_factor
    real(dp) :: max_rij_change

    real(dp) :: F0(3,3), F0_evals(3), dF(3,3), dF_evals(3)

    call print("max_rij_change, size(last_connect_x) " // size(last_connect_x) // &
      " size(x) " // size(x), PRINT_NERD)
    if (size(last_connect_x) == size(x)) then
      F0 = lat_factor*reshape(last_connect_x(1:9), (/ 3,3 /))
      dF = lat_factor*reshape(x(1:9), (/ 3,3 /)) - F0
      call diagonalise(matmul(F0,transpose(F0)), F0_evals)
      call diagonalise(matmul(dF,transpose(dF)), dF_evals)
      F0_evals = sqrt(F0_evals)
      dF_evals = sqrt(dF_evals)
      call print("max_rij_change = " // maxval(abs(F0_evals)) //"*2.0_dp*"// &
        maxval(abs(last_connect_x-x))*sqrt(3.0_dp)  // &
        " + " // maxval(abs(dF_evals))//"*"//r_cut, PRINT_VERBOSE)
      max_rij_change = maxval(abs(F0_evals))*2.0_dp*maxval(abs(last_connect_x-x))*sqrt(3.0_dp) + &
                 maxval(abs(dF_evals))*r_cut
    else
      call system_abort("max_rij_change: size(last_connect_x)="//size(last_connect_x)// &
        " /= size(x)="//size(x))
    endif

  end function max_rij_change

  ! prepare atoms structure given a deformation gradient matrix
  subroutine prep_atoms_deform_grad(deform_grad, at, am)
    real(dp), intent(in) :: deform_grad(3,3)
    type(Atoms), intent(inout) :: at
    type(potential_minimise), intent(inout) :: am

    am%minim_save_lat = at%lattice
    call set_lattice(at, deform_grad .mult. am%minim_save_lat, scale_positions=.true.)
  end subroutine prep_atoms_deform_grad

  ! remove effect of deformation gradient applied by prep_atoms_deform_grad
  subroutine fix_atoms_deform_grad(deform_grad, at,am)
    real(dp), intent(in) :: deform_grad(3,3)
    type(Atoms), intent(inout) :: at
    type(potential_minimise), intent(inout) :: am

    real(dp) :: deform_grad_inv(3,3)

    call inverse(deform_grad, deform_grad_inv)
    at%pos = deform_grad_inv .mult. at%pos
    call set_lattice(at, am%minim_save_lat, scale_positions=.false.)

  end subroutine fix_atoms_deform_grad

  ! unpack a deformation grad and a pos array from 1-D array into a atoms%pos 2-D array
  subroutine unpack_pos_dg(xx, at_N, at_pos, dg, lat_factor)
    real(dp), intent(in) :: xx(:)
    integer :: at_N
    real(dp), intent(inout) :: at_pos(:,:)
    real(dp) :: dg(3,3)
    real(dp), intent(in) :: lat_factor

    if (3*at_N+9 /= size(xx)) call system_abort("Called unpack_pos with mismatching sizes x " // size(xx) // " at " // at_N)

    dg = lat_factor*reshape(xx(1:9), (/ 3,3 /) )
    at_pos = reshape(xx(10:), (/ 3, at_N /) )

end subroutine unpack_pos_dg

  ! pack a 3xN 2-D array into a 1-D array
  subroutine pack_pos_dg(x2d, dg2d, x, lat_factor)
    real(dp), intent(in) :: x2d(:,:)
    real(dp), intent(in)  :: dg2d(3,3)
    real(dp), intent(out) :: x(:)
    real(dp), intent(in)  :: lat_factor

    if (size(x2d,1) /= 3) call system_abort("Called pack with mismatching size(x2d,1) " // &
      size(x2d,1) // " != 3")
    
    if (size(dg2d,1) /= 3 .or. size(dg2d,2) /= 3) &
      call system_abort("Called pack with mismatching size(dg2d,1) " // shape(dg2d) // " != 3x3")
    if (size(x) /= (size(x2d)+size(dg2d))) &
      call system_abort("Called pack with mismatching size x " // size(x) // " x2d " // size(x2d) // &
      " dg2d " // size(dg2d))

    x(1:9) = lat_factor*reshape(dg2d, (/ 9 /) )

      x(10:) = reshape(x2d, (/ size(x2d) /) )

end subroutine pack_pos_dg

  subroutine Potential_read_params_xml(this, param_str)
     type(Potential), intent(inout), target :: this
     character(len=*), intent(in) :: param_str

     type(xml_t) :: fxml

     if (len(trim(param_str)) <= 0) &
     call system_abort('Potential_read_params_xml: invalid param_str length '//len(trim(param_str)) )

     parse_in_pot = .false.
     parse_in_pot_done = .false.
     parse_matched_label = .false.
     parse_pot => this

     call open_xml_string(fxml, param_str)

     call parse(fxml,  &
       startElement_handler = Potential_startElement_handler, &
       endElement_handler = Potential_endElement_handler)
     call close_xml_t(fxml)

     if(.not. parse_in_pot_done) then
	if (len_trim(param_str) > 10000) then
	   call system_abort('Potential_read_params_xml: could not initialise potential from xml_label. param_str(1:10000)='//trim(param_str(1:10000)))
	else
	   call system_abort('Potential_read_params_xml: could not initialise potential from xml_label. param_str='//trim(param_str))
	endif
     endif

  endsubroutine Potential_read_params_xml

  subroutine Potential_startElement_handler(URI, localname, name, attributes)
     character(len=*), intent(in)   :: URI
     character(len=*), intent(in)   :: localname
     character(len=*), intent(in)   :: name
     type(dictionary_t), intent(in) :: attributes
   
     integer :: status
     character(len=STRING_LENGTH) :: value

     if(name == 'Potential') then ! new Potential stanza

        if(parse_in_pot) &
           call system_abort("Potential_startElement_handler entered <Potential> with parse_in true. Probably a forgotten /> at the end of a tag.")

        if(parse_matched_label) return ! we already found an exact match for this label

        call QUIP_FoX_get_value(attributes, 'label', value, status)
        if(status /= 0) value = ''

        if(len(trim(parse_pot%xml_label)) > 0) then ! we were passed in a label
           if(trim(value) == trim(parse_pot%xml_label)) then ! exact match
              parse_matched_label = .true.
              parse_in_pot = .true.
           else ! no match
              parse_in_pot = .false.
           endif
        else ! no label passed in
           call system_abort("Potential_startElement_handler: no label passed in")
        endif

        call QUIP_FoX_get_value(attributes, 'init_args', value, status)
        if(status == 0) then
           read (value, '(a)') parse_pot%xml_init_args
        else
           call system_abort("Potential_startElement_handler: no init_args attribute found")
        endif
        
     endif

  endsubroutine Potential_startElement_handler

  subroutine Potential_endElement_handler(URI, localname, name)
     character(len=*), intent(in)   :: URI
     character(len=*), intent(in)   :: localname
     character(len=*), intent(in)   :: name

     if(parse_in_pot) then
        if(name == 'Potential') then
           parse_in_pot = .false.
           parse_in_pot_done = .true.
        endif
     endif

  endsubroutine Potential_endElement_handler

  subroutine potential_set_callback(this, callback)
    type(Potential), intent(inout) :: this
    interface
       subroutine callback(at)
         integer, intent(in) :: at(12)
       end subroutine callback
    end interface
    
    if (this%is_simple) then
       call set_callback(this%simple, callback)
    else
       call system_abort('potential_set_callback() only implemented for simple Potentials.')
    end if

  end subroutine potential_set_callback


#include "Potential_Sum_routines.f95"
#include "Potential_ForceMixing_routines.f95"
#include "Potential_EVB_routines.f95"

#ifdef HAVE_LOCAL_E_MIX
#include "Potential_Local_E_Mix_routines.f95"
#endif
#ifdef HAVE_ONIOM
#include "Potential_ONIOM_routines.f95"
#endif

#include "Potential_Cluster_routines.f95"
#include "Potential_Hybrid_utils.f95"


  !% Run 'n_steps' of dynamics using forces from Potential 'pot'.
  !%
  !% For each step, forces are evaluated using the Potential
  !% 'pot' and the DynamicalSystem is advanced by a time 'dt'
  !% (default 1 fs).  'n_steps' (default 10 steps) are carried out in
  !% total, with snapshots saved every 'save_interval' steps. The
  !% connectivity is recalculated every 'connect_interval' steps.
  !% 'args_str' can be used to supply extra arguments to 'Potential%calc'.
  subroutine DynamicalSystem_run(this, pot, dt, n_steps, hook, hook_interval, summary_interval, write_interval, connect_interval, trajectory, args_str, error)
    type(DynamicalSystem), intent(inout), target :: this
    type(Potential), intent(inout) :: pot
    real(dp), intent(in) :: dt
    integer, intent(in) :: n_steps
    integer, intent(in), optional :: summary_interval, hook_interval, write_interval, connect_interval
    type(CInOutput), intent(inout), optional :: trajectory
    character(len=*), intent(in), optional :: args_str
    integer, intent(out), optional :: error
    interface
       subroutine hook()
       end subroutine hook
    end interface    
    
    integer :: n, my_summary_interval, my_hook_interval, my_write_interval, my_connect_interval
    real(dp) :: e
    real(dp), pointer, dimension(:,:) :: f
    character(len=STRING_LENGTH) :: my_args_str
    type(Dictionary) :: params

    INIT_ERROR(error)

    my_summary_interval = optional_default(1, summary_interval)
    my_hook_interval = optional_default(1, hook_interval)
    my_write_interval = optional_default(1, write_interval)
    my_connect_interval = optional_default(1, connect_interval)
    my_args_str = optional_default("", args_str)
    call initialise(params)
    call read_string(params, my_args_str)
    if (.not. has_key(params, 'energy')) call set_value(params, 'energy')
    if (.not. has_key(params, 'force'))  call set_value(params, 'force')
    my_args_str = write_string(params)
    call finalise(params)

    call calc_connect(this%atoms)
    call calc(pot, this%atoms, args_str=my_args_str, error=error)
    PASS_ERROR(error)
    call set_value(this%atoms%params, 'time', this%t)
    if (.not. get_value(this%atoms%params, 'energy', e)) &
         call system_abort("dynamicalsystem_run failed to get energy")
    if (.not. assign_pointer(this%atoms, 'force', f)) &
         call system_abort("dynamicalsystem_run failed to get forces")
    if (my_summary_interval > 0) call ds_print_status(this, epot=e)
    call hook()
    if (present(trajectory)) call write(trajectory, this%atoms)

    ! initialize accelerations from forces, so first call to verlet1 will be correct
    this%atoms%acc(1,:) = f(1,:)/this%atoms%mass
    this%atoms%acc(2,:) = f(2,:)/this%atoms%mass
    this%atoms%acc(3,:) = f(3,:)/this%atoms%mass

    do n=1,n_steps
       call advance_verlet1(this, dt)
       call calc(pot, this%atoms, args_str=my_args_str, error=error)
       PASS_ERROR(error)
       call advance_verlet2(this, dt, f)
       if (.not. get_value(this%atoms%params, 'energy', e)) &
            call system_abort("dynamicalsystem_run failed to get energy")
       if (.not. assign_pointer(this%atoms, 'force', f)) &
            call system_abort("dynamicalsystem_run failed to get forces")
       if (my_summary_interval > 0 .and. mod(n, my_summary_interval) == 0) call ds_print_status(this, epot=e)
       call set_value(this%atoms%params, 'time', this%t)

       if (my_hook_interval > 0 .and. mod(n,my_hook_interval) == 0) call hook()
       if (present(trajectory) .and. my_write_interval > 0 .and. mod(n,my_write_interval) == 0) call write(trajectory, this%atoms)
       if (my_connect_interval > 0 .and. mod(n,my_connect_interval) == 0) call calc_connect(this%atoms)
    end do

  end subroutine DynamicalSystem_run

  subroutine constrain_DG(at, deform_grad)
    type(Atoms), intent(in) :: at
    real(dp), intent(inout) :: deform_grad(3,3)

    logical :: minim_constant_volume
    real(dp) :: scaled_ident(3,3), DG_det
    integer :: error

    minim_constant_volume = .false.
    call get_param_value(at, "Minim_Constant_Volume", minim_constant_volume, error=error)

    ! zero out trace component, which changes volume
    if (minim_constant_volume) then
      DG_det = matrix3x3_det(deform_grad)
      scaled_ident = 0.0_dp
      call add_identity(scaled_ident)
      scaled_ident = scaled_ident * (DG_det**(-1.0/3.0))
      deform_grad = matmul(deform_grad, scaled_ident)
    endif
  end subroutine constrain_DG

  subroutine constrain_virial_post(at, virial)
    type(Atoms), intent(in) :: at
    real(dp), intent(inout) :: virial(3,3)

    integer :: error
    logical :: minim_hydrostatic_strain
    logical :: minim_lattice_fix_mask(3,3)
    real(dp) :: minim_lattice_fix(3,3)
    real(dp) :: scaled_ident

    real(dp) :: virial_trace

    minim_hydrostatic_strain = .false.
    call get_param_value(at, "Minim_Hydrostatic_Strain", minim_hydrostatic_strain, error=error)
    CLEAR_ERROR(error)
    minim_lattice_fix = 0.0_dp
    call get_param_value(at, "Minim_Lattice_Fix", minim_lattice_fix, error=error)
    CLEAR_ERROR(error)
    minim_lattice_fix_mask = (minim_lattice_fix /= 0.0_dp)

    ! project onto identity
    if (minim_hydrostatic_strain) then
      virial_trace = trace(virial)
      virial = 0.0_dp
      virial(1,1) = virial_trace/3.0_dp
      virial(2,2) = virial_trace/3.0_dp
      virial(3,3) = virial_trace/3.0_dp
    endif
    ! Zero out components corresponding to fixed lattice elements
    if (any(minim_lattice_fix_mask)) then
       virial = merge(0.0_dp, virial, minim_lattice_fix_mask)
    end if

  end subroutine constrain_virial_post



end module Potential_module
