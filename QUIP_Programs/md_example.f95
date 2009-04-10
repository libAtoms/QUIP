program example_LJ
use libAtoms_module
use QUIP_module
implicit none

  character(len=10240) lj_str
  type(Potential) pot
  type(Atoms) at
  type(DynamicalSystem) ds

  real(dp) :: e
  real(dp), allocatable :: f(:,:)
  integer it

  call system_initialise()

  lj_str = '<LJ_params n_types="1">' // &
   '<per_type_data type="1" atomic_num="29" />' // &
   '<per_pair_data type1="1" type2="1" sigma="2.0" eps6="2.0" eps12="2.0" cutoff="4.0" shifted="T" />' // &
   '</LJ_params>'

  call Initialise(pot, "IP LJ", lj_str)

  call print(pot)

  call atoms_read_xyz_filename(at, "md.xyz")
  call ds_initialise(ds, at)

  call set_cutoff(ds%atoms, cutoff(pot))


  allocate(f(3,ds%atoms%N))
  call rescale_velo(ds, 300.0_dp)

  call print(at)
  call print(ds%atoms)

  do it=1, 3
    if (mod(it,5) == 1) then
      call calc_connect(ds%atoms)
    end if

    call calc(pot, ds%atoms, e = e, f = f)

    call advance_verlet(ds, 1.0_dp, f)
    call ds_print_status(ds, 'D', e)
  end do

  call print(ds%atoms)

end program
