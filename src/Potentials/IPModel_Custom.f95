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
!X IPModel_Custom
!X
!% Customised interatomic potential: 
!%
!% Energy and Force routines are hardwired
!% Cutoff is hardwired
!% 
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#include "error.inc"

module IPModel_Custom_module

use error_module
use system_module, only : dp, inoutput, print, PRINT_NERD, operator(//)
use dictionary_module
use paramreader_module
use linearalgebra_module
use atoms_types_module
use atoms_module
use topology_module

use mpi_context_module
use QUIP_Common_module

implicit none
private

include 'IPModel_interface.h'

public :: IPModel_Custom
type IPModel_Custom
  real(dp) :: cutoff = 0.0_dp
  real(dp) :: kbond = 0.0_dp
  real(dp) :: kangle = 0.0_dp
  real(dp) :: bond_r0 = 0.0_dp
  real(dp) :: angle_cos0 = 0.0_dp
end type IPModel_Custom

logical, private :: parse_in_ip, parse_matched_label
type(IPModel_Custom), private, pointer :: parse_ip

interface Initialise
  module procedure IPModel_Custom_Initialise_str
end interface Initialise

interface Finalise
  module procedure IPModel_Custom_Finalise
end interface Finalise

interface Print
  module procedure IPModel_Custom_Print
end interface Print

interface Calc
  module procedure IPModel_Custom_Calc
end interface Calc

contains

subroutine IPModel_Custom_Initialise_str(this, args_str, param_str, error)
  type(IPModel_Custom), intent(inout) :: this
  character(len=*), intent(in) :: args_str, param_str
  type(Dictionary) :: params
  integer, optional, intent(out):: error


  INIT_ERROR(error)
  call Finalise(this)

  call initialise(params)
  call param_register(params, 'kbond', '0.0', this%kbond, help_string='Strength of quadratic restraint on C-H bonds.  Potential is kconf*(r-r0)^2')
  call param_register(params, 'kangle', '0.0', this%kangle, help_string='Strength of quadratic restraint on H-C-H cosines.  Potential is kconf*(cos(theta)-cos(theta0))^2')
  call param_register(params, 'bond_r0', '0.0', this%bond_r0, help_string='Equilibrium bond length for C-H bonds.')
  call param_register(params, 'angle_cos0', '0.0', this%angle_cos0, help_string='Cosine of equilibrium bond angle for H-C-H triplets.')
  call param_register(params, 'cutoff', '0.0', this%cutoff, help_string='Cutoff for finding methane monomers')
  if(.not. param_read_line(params, args_str, ignore_unknown=.true., task='IPModel_Custom_Initialise args_str')) then
     RAISE_ERROR("IPModel_Custom_Init failed to parse args_str='"//trim(args_str)//"'", error)
  end if
  call finalise(params)

end subroutine IPModel_Custom_Initialise_str

subroutine IPModel_Custom_Finalise(this)
  type(IPModel_Custom), intent(inout) :: this

  ! Add finalisation code here

end subroutine IPModel_Custom_Finalise


subroutine IPModel_Custom_Calc(this, at, e, local_e, f, virial, local_virial, args_str, mpi, error)
   type(IPModel_Custom), intent(inout):: this
   type(Atoms), intent(inout)      :: at
   real(dp), intent(out), optional :: e, local_e(:)
   real(dp), intent(out), optional :: f(:,:), local_virial(:,:)   !% Forces, dimensioned as \texttt{f(3,at%N)}, local virials, dimensioned as \texttt{local_virial(9,at%N)} 
   real(dp), intent(out), optional :: virial(3,3)
   character(len=*), optional      :: args_str
   type(MPI_Context), intent(in), optional :: mpi
   integer, intent(out), optional :: error

   ! Add calc() code here

   ! Confining potential for methanes - first find general monomers, then place
   ! harmonic restraints on bonds and angles.

   ! NB be sure to use the neighbour list passed from lammps

   integer, dimension(:,:), allocatable :: monomer_index
   logical, dimension(:), allocatable :: is_associated
   real(dp) :: energy, force(3,at%N)

   type(Dictionary) :: params
   character(STRING_LENGTH) :: atom_mask_name
   logical :: has_atom_mask_name
   logical, dimension(:), pointer :: atom_mask_pointer

   integer :: mon_i, atom_i, rank_j, atom_j, rank_k, atom_k
   real(dp) :: e_pair, d_epair_dr, rij(3), rik(3), rij_mag, rik_mag, rij_norm(3), rik_norm(3)
   real(dp) :: e_trip, d_etrip_dcos, cos_ijk, fij(3), fik(3)
   real(dp) :: virial_i(3,3), virial_j(3,3), virial_k(3,3)


   INIT_ERROR(error)

   if(present(args_str)) then
      call initialise(params)

      call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, &
      help_string="Name of a logical property in the atoms object. For monomers where this property is true the Potential is " // &
      "calculated.")

      if (.not. param_read_line(params,args_str,ignore_unknown=.true.,task='general_dimer_calc args_str')) then
         RAISE_ERROR("general_dimer_calc failed to parse args_str='"//trim(args_str)//"'", error)
      endif

      call finalise(params)

      if( has_atom_mask_name ) then
         if (.not. assign_pointer(at, trim(atom_mask_name), atom_mask_pointer)) then
            RAISE_ERROR("general_dimer_calc did not find "//trim(atom_mask_name)//" property in the atoms object.", error)
         endif
      else
         atom_mask_pointer => null()
      endif
   endif

   if (present(e)) e = 0.0_dp
   if (present(local_e)) then
      call check_size('Local_E',local_e,(/at%N/),'IPModel_Custom_Calc', error)
      local_e = 0.0_dp
   endif
   if (present(f)) then
      call check_size('Force',f,(/3,at%Nbuffer/),'IPModel_Custom_Calc', error)
      f = 0.0_dp
   end if
   if (present(virial)) virial = 0.0_dp
   if (present(local_virial)) then
      call check_size('Local_virial',local_virial,(/9,at%Nbuffer/),'IPModel_Custom_Calc', error)
      local_virial = 0.0_dp
   endif

   allocate(is_associated(at%N))
   is_associated = .false.
   call find_general_monomer(at, monomer_index, (/6, 1, 1, 1, 1/), is_associated, this%cutoff, general_ordercheck=.true., error=error)

   if(.not. all(is_associated)) then
      call print("WARNING: IP Custom: not all atoms assigned to a methane monomer. If you have partial monomers this is OK.", PRINT_NERD)
   end if

   ! First, loop over monomers.  These are also the centres of angle triplets.
   do mon_i = 1, size(monomer_index, 2)
      atom_i = monomer_index(1, mon_i)
      ! Only evaluate the monomer if the central C atom is local - should check
      ! how lammps accounts for forces on non-local atoms
      if (associated(atom_mask_pointer)) then
         if (.not. atom_mask_pointer(atom_i)) cycle
      end if

      do rank_j = 1, n_neighbours(at, atom_i, max_dist=this%cutoff)
         atom_j = neighbour(at, atom_i, rank_j, distance=rij_mag, diff=rij, cosines=rij_norm, max_dist=this%cutoff)
         if (rij_mag .feq. 0.0_dp) cycle
         ! This really shouldn't happen in any normal cases, but account for it anyways
         if (.not. (any(monomer_index(2:5, mon_i) .eq. atom_j))) cycle

         ! Pairwise forces and energies
         e_pair = this%kbond * (rij_mag - this%bond_r0)**2
         d_epair_dr = 2.0_dp * this%kbond * (rij_mag - this%bond_r0)

         if (present(e)) e = e + e_pair
         if (present(local_e)) then
            ! Eh, let's just concentrate all the 'local' quantities on the monomer centres.
            local_e(atom_i) = local_e(atom_i) + e_pair
         end if
         if (present(f)) then
            f(:, atom_j) = f(:, atom_j) - 1.0_dp * d_epair_dr * rij_norm
            f(:, atom_i) = f(:, atom_i) + 1.0_dp * d_epair_dr * rij_norm
         end if
         if (present(virial) .or. present(local_virial)) then
            virial_i = -1.0 * d_epair_dr * rij_mag * (rij_norm .outer. rij_norm)
         end if
         if (present(virial)) virial = virial + virial_i
         if (present(local_virial)) then
            local_virial(:, atom_i) = local_virial(:, atom_i) + reshape(virial_i, (/9/))
         end if

         do rank_k = 1, n_neighbours(at, atom_i, max_dist=this%cutoff)
            atom_k = neighbour(at, atom_i, rank_k, distance=rik_mag, diff=rik, cosines=rik_norm, max_dist=this%cutoff)
            ! Again, shouldn't happen, but better to be safe
            if (.not. (any(monomer_index(2:5, mon_i) .eq. atom_k))) cycle
            if (atom_k <= atom_j) cycle

            cos_ijk = sum(rij_norm*rik_norm)
            e_trip = this%kangle * (cos_ijk - this%angle_cos0)**2
            d_etrip_dcos = 2.0_dp * this%kangle * (cos_ijk - this%angle_cos0)
            if (present(e)) e = e + e_trip
            if (present(local_e)) then
               ! Hm, the triplet local quantities can be assigned to the outer atoms
               local_e(atom_j) = local_e(atom_j) + 0.5_dp * e_trip
               local_e(atom_k) = local_e(atom_k) + 0.5_dp * e_trip
            end if
            if (present(f) .or. present(virial) .or. present(local_virial)) then
               fij = -1.0_dp * d_etrip_dcos * (rik_norm - rij_norm*cos_ijk) / rij_mag
               fik = -1.0_dp * d_etrip_dcos * (rij_norm - rik_norm*cos_ijk) / rik_mag
            end if
            if (present(f)) then
               f(:, atom_i) = f(:, atom_i) + fij + fik
               f(:, atom_j) = f(:, atom_j) - fij
               f(:, atom_k) = f(:, atom_k) - fik
            end if
            if (present(virial) .or. present(local_virial)) then
               virial_j = -1.0_dp * d_etrip_dcos * ((rik_norm .outer. rij_norm) - (rij_norm .outer. rij_norm)*cos_ijk)
               virial_k = -1.0_dp * d_etrip_dcos * ((rij_norm .outer. rik_norm) - (rik_norm .outer. rik_norm)*cos_ijk)
            end if
            if (present(virial)) virial = virial + virial_j + virial_k
            if (present(local_virial)) then
               local_virial(:, atom_i) = local_virial(:, atom_i) + reshape(virial_j, (/9/)) + reshape(virial_k, (/9/))
            end if
         end do
      end do
   end do

   deallocate(monomer_index)
   deallocate(is_associated)

end subroutine IPModel_Custom_Calc


subroutine IPModel_Custom_Print(this, file)
  type(IPModel_Custom), intent(in) :: this
  type(Inoutput), intent(inout),optional :: file

  call Print("IPModel_Custom : Custom Potential", file=file)
  call Print("IPModel_Custom : cutoff = " // this%cutoff, file=file)
  call Print("IPModel_Custom : kbond = " // this%kbond, file=file)
  call Print("IPModel_Custom : kangle = " // this%kangle, file=file)
  call Print("IPModel_Custom : bond_r0 = " // this%bond_r0, file=file)
  call Print("IPModel_Custom : angle_cos0 = " // this%angle_cos0, file=file)

end subroutine IPModel_Custom_Print


end module IPModel_Custom_module
