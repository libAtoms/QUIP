! HEADER_CHECK

program lotf_metapot

  use libAtoms_module
  use QUIP_Module

  implicit None

  type(Atoms) :: at
  type(Potential) :: pot
  real(dp), allocatable :: f(:,:)
  integer, pointer :: hybrid(:)
  integer :: i, j, n_cycle, n_extrap
  type(Table) :: embedlist
  type(MetaPotential) :: lotf
  type(DynamicalSystem) :: ds, ds_saved
  type(CInOutput) :: movie

  call system_initialise(NORMAL)
  call initialise(pot, 'FilePot command=./castep_driver.sh')

  call read(at, "start.xyz")
  call initialise(movie, "movie.xyz", action=OUTPUT)
  

  n_cycle = 10
  n_extrap = 10

  allocate(f(3,at%N))

  call set_cutoff(at, 4.0_dp)
  call calc_connect(at)

  ! Mark QM region
  call add_property(at, 'hybrid', 0)
  if (.not. assign_pointer(at, 'hybrid', hybrid)) &
       call system_abort('Cannot assign hybrid pointer')
  hybrid = 1

  ! initialise Dynamical System
  call initialise(ds, at)
  call rescale_velo(ds, 300.0_dp)

  ! Set up metapotential for LOTF
  call initialise(lotf, 'ForceMixing method=lotf_adj_pot_svd fit_hops=0 buffer_hops=0 '//&
       'randomise_buffer=F terminate=F min_images_only=T qm_args_str={cluster_same_lattice=T  single_cluster=T cluster_calc_connect=T}', pot)

  call print_title('LOTF Metapotential')
  call print(lotf)

  ! Predictor-corrector dynamics
  call Print_Title('Bootstrap')

  call calc(lotf, ds%atoms, f=f) ! bootstrap the adjustable potential

  do i=1,n_cycle

     ! Extrapolation
     call print_title('Extrapolation')
     ds_saved = ds

     do j=1,n_extrap
        if (j == 1) then
           call calc(lotf, ds%atoms, f=f, args_str="lotf_do_qm=F lotf_do_init=T lotf_do_map=T")
        else
           call calc(lotf, ds%atoms, f=f, args_str="lotf_do_qm=F lotf_do_init=F")
        end if
        call advance_verlet(ds, 1.0_dp, f)
        call ds_print_status(ds, 'E', instantaneous=.true.)
     end do

     call print_title('QM Force Evaluation')
     ! QM force calculation and optimisation of adj pot
     call calc(lotf, ds%atoms, f=f, args_str="lotf_do_qm=T lotf_do_init=F lotf_do_fit=T")

     ! Interpolation
     call print_title('Interpolation')
     ds = ds_saved
     do j=1,n_extrap
        call calc(lotf, ds%atoms, f=f, args_str="lotf_do_qm=F lotf_do_init=F lotf_do_interp=T lotf_interp="&
             //(real(j-1,dp)/real(n_extrap,dp)))
        call advance_verlet(ds, 1.0_dp, f)
        call ds_print_status(ds, 'I', instantaneous=.true.)
     end do

     call print_title('Connectivity update')
     call calc_connect(ds%atoms)

     call write(movie, ds%atoms)

  end do

  deallocate(f)
  call finalise(pot)
  call finalise(at)
  call finalise(ds)
  call finalise(movie)
  call system_finalise
  
end program lotf_metapot
