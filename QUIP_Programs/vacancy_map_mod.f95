module vacancy_map_module
use libatoms_module
use quip_module
use structures_module
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
  type(MetaPotential), intent(out) :: hybridpot
  character(len=*), optional, intent(in) :: hybrid_args_str

  type(Atoms) :: bulk, bulk1

  ! potential 1
  call fcc(bulk, 4.02406_dp)
  call supercell(bulk1, bulk, 4, 4, 4)
  bulk1%Z(:) = 13
  call set_cutoff(bulk1, cutoff(pot1)+0.5_dp)
  call calc_connect(bulk1)

  ! should do minimise_bulk instead, but no mechanism for doing it just for pot1
  call initialise(hybridpot, "Hybrid terminate=F r_scale=1.0005154009139747 " // trim(hybrid_args_str), pot2, pot1, bulk1)

  call print(hybridpot)

end subroutine
#endif

#ifdef HAVE_LOTF
subroutine init_lotf_forcemix(pot1, pot2, lotfpot, buffer_hops, lotf_args_str)
  type(Potential), intent(inout), target :: pot1, pot2
  type(MetaPotential), intent(out) :: lotfpot
  integer, intent(in) :: buffer_hops
  character(len=*), optional, intent(in) :: lotf_args_str

  type(Atoms) :: bulk, bulk1

  ! potential 1
  call fcc(bulk, 4.02406_dp)
  call supercell(bulk1, bulk, 4, 4, 4)
  bulk1%Z(:) = 13
  call set_cutoff(bulk1, cutoff(pot1)+0.5_dp)
  call calc_connect(bulk1)

  call print("WARNING: No way to r_scale in init_lotf_forcemix", ERROR)

  ! should do minimise_bulk instead, but no mechanism for doing it just for pot1
  call initialise(lotfpot, "LOTF fit_method=force_mixing_abrupt small_clusters=F periodic_y terminate=F buffer_hops="//buffer_hops//" "// trim(lotf_args_str), pot2, pot1, bulk1)

  call print(lotfpot)

end subroutine
#endif

end module vacancy_map_module
