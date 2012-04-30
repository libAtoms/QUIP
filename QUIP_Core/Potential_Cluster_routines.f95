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
!X Cluster routines to be included in Potential.f95
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  recursive subroutine potential_Cluster_initialise(this, args_str, inner_pot, mpi_obj, error)
    type(Potential_Cluster), intent(inout) :: this
    character(len=*), intent(in) :: args_str
    type(Potential), intent(inout), target :: inner_pot
    type(MPI_Context), intent(in), optional :: mpi_obj
    integer, intent(out), optional :: error

    type(Dictionary) :: params
    logical :: minimise_bulk
    logical :: do_rescale_r, do_rescale_E, do_tb_defaults

    INIT_ERROR(error)

    call initialise(params)
    call param_register(params, "run_suffix", "", this%run_suffix, help_string="Suffix to apply to hybrid mark properties$")
    call param_register(params, "r_scale", "1.0", this%r_scale_pot1, help_string="Rescaling factor for cluster positions")
    call param_register(params, "E_scale", "1.0", this%E_scale_pot1, help_string="Rescaling factor for cluster energies (and hence also forces)")

    if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='Potential_Cluster_initialise args_str') ) then
      RAISE_ERROR("Potential_Cluster_initialise failed to parse args_str='"//trim(args_str)//"'",error)
    endif
    call finalise(params)

    if (this%r_scale_pot1 <= 0.0_dp) this%r_scale_pot1 = 1.0_dp
    if (this%E_scale_pot1 <= 0.0_dp) this%E_scale_pot1 = 1.0_dp

    call print ("Rescaling positions in region1 potential by " // this%r_scale_pot1 // " to match lattice constants")
    call print ("Rescaling energies in region1 potential by " // this%E_scale_pot1 // " to match bulk modulus")

    this%inner_pot => inner_pot

  end subroutine

  recursive subroutine potential_Cluster_finalise(this)
    type(Potential_Cluster), intent(inout) :: this

    nullify(this%inner_pot)

  end subroutine potential_Cluster_finalise


  recursive function potential_Cluster_cutoff(this)
    type(Potential_Cluster), intent(in) :: this
    real(dp) :: potential_Cluster_cutoff

    potential_Cluster_cutoff = cutoff(this%inner_pot)

  end function potential_Cluster_cutoff

  recursive subroutine potential_Cluster_print(this, file)
    type(Potential_Cluster),          intent(inout) :: this
    type(Inoutput),intent(inout),optional:: file
    
    call print('Cluster potential:', file=file)
    call print('Inner Potential:', file=file)
    call print('================', file=file)
    call print(this%inner_pot, file=file)
    call print('')
    call print('r_scale_pot1=' // this%r_scale_pot1 // ' E_scale_pot1=' // this%E_scale_pot1, file=file)
    call print('')
  end subroutine potential_Cluster_print


  recursive subroutine potential_Cluster_calc(this, at, args_str, error)
    type(Potential_Cluster), intent(inout) :: this
    type(Atoms), intent(inout) :: at
    character(len=*), intent(in), optional :: args_str
    integer, intent(out), optional :: error

    type(Dictionary) :: params
    character(STRING_LENGTH) :: calc_energy, calc_force, calc_virial, calc_local_energy, calc_local_virial, run_suffix
    logical single_cluster, little_clusters

    INIT_ERROR(error)

    if (.not. associated(this%inner_pot)) then
      RAISE_ERROR("Potential_Cluster_calc: this%inner_pot", error)
    endif

    call initialise(params)
    call param_register(params, "run_suffix", trim(this%run_suffix), run_suffix, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, "energy", "", calc_energy, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, "force", "", calc_force, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, "virial", "", calc_virial, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, "local_energy", "", calc_local_energy, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, "local_virial", "", calc_local_virial, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params, 'single_cluster', 'F', single_cluster, &
      help_string="If true, calculate all active/transition atoms with a single big cluster")
    call param_register(params, 'little_clusters', 'F', little_clusters, &
      help_string="If true, calculate forces (only) by doing each atom separately surrounded by a little buffer cluster")
    if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='Potential_Cluster_Calc args_str') ) then
      RAISE_ERROR("Potential_Cluster_calc_energy failed to parse args_str='"//trim(args_str)//"'",error)
    endif
    call finalise(params)
    
    if (len_trim(calc_energy) /= 0 .or. len_trim(calc_local_energy) /= 0 .or. &
        len_trim(calc_virial) /= 0 .or. len_trim(calc_local_virial) /= 0) then
       RAISE_ERROR("Potential_Cluster_calc: can only calculate forces, not energies or virials", error)
    end if

    if (single_cluster .and. little_clusters) then
         RAISE_ERROR('Potential_Cluster_calc: single_cluster and little_clusters options are mutually exclusive', error)
    end if

  end subroutine potential_Cluster_calc


  subroutine calc_single_cluster(this, at, args_str, error)
    type(Potential_Cluster), intent(inout) :: this
    type(Atoms), intent(inout) :: at
    character(len=*), intent(in), optional :: args_str
    integer, intent(out), optional :: error

    INIT_ERROR(error)

  end subroutine calc_single_cluster

  subroutine calc_little_clusters(this, at, args_str, error)
    type(Potential_Cluster), intent(inout) :: this
    type(Atoms), intent(inout) :: at
    character(len=*), intent(in), optional :: args_str
    integer, intent(out), optional :: error

    INIT_ERROR(error)

  end subroutine calc_little_clusters
