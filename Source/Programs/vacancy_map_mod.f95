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

module vacancy_map_module
use libatoms_module
use potential_module
implicit none
private

#if defined(HAVE_LOCAL_E_MIX) || defined(HAVE_ONIOM)
public :: init_hybrid
#endif
#ifdef HAVE_LOTF
public :: init_lotf_forcemix
#endif

#if defined(HAVE_LOCAL_E_MIX) || defined(HAVE_ONIOM) || defined(HAVE_LOTF)
contains
#endif

#if defined(HAVE_LOCAL_E_MIX) || defined(HAVE_ONIOM)
subroutine init_hybrid(pot1, pot2, hybridpot, hybrid_args_str)
  type(Potential), intent(inout), target :: pot1, pot2
  type(Potential), intent(out) :: hybridpot
  character(len=*), optional, intent(in) :: hybrid_args_str

  type(Atoms) :: bulk, bulk1

  ! potential 1
  call fcc(bulk, 4.02406_dp)
  call supercell(bulk1, bulk, 4, 4, 4)
  bulk1%Z(:) = 13
  call set_cutoff(bulk1, cutoff(pot1)+0.5_dp)
  call calc_connect(bulk1)

  ! should do minimise_bulk instead, but no mechanism for doing it just for pot1
  call initialise(hybridpot, "Hybrid terminate=F r_scale=1.0005154009139747 " // trim(hybrid_args_str), pot2, pot1, bulk_scale=bulk1)

  call print(hybridpot)

end subroutine
#endif

#ifdef HAVE_LOTF
subroutine init_lotf_forcemix(pot1, pot2, lotfpot, buffer_hops, lotf_args_str)
  type(Potential), intent(inout), target :: pot1, pot2
  type(Potential), intent(out) :: lotfpot
  integer, intent(in) :: buffer_hops
  character(len=*), optional, intent(in) :: lotf_args_str

  type(Atoms) :: bulk, bulk1

  ! potential 1
  call fcc(bulk, 4.02406_dp)
  call supercell(bulk1, bulk, 4, 4, 4)
  bulk1%Z(:) = 13
  call set_cutoff(bulk1, cutoff(pot1)+0.5_dp)
  call calc_connect(bulk1)

  call print("WARNING: No way to r_scale in init_lotf_forcemix", PRINT_ALWAYS)

  ! should do minimise_bulk instead, but no mechanism for doing it just for pot1
  call initialise(lotfpot, args_str="LOTF fit_method=force_mixing_abrupt small_clusters=F periodic_y terminate=F buffer_hops="//buffer_hops//" "// trim(lotf_args_str), pot1=pot2, pot2=pot1, bulk_scale=bulk1)

  call print(lotfpot)

end subroutine
#endif

end module vacancy_map_module
