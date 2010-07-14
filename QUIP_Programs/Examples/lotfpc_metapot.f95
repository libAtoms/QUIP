! HJ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HJ X
! HJ X   libAtoms+QUIP: atomistic simulation library
! HJ X
! HJ X   Portions of this code were written by
! HJ X     Albert Bartok-Partay, Silvia Cereda, Gabor Csanyi, James Kermode,
! HJ X     Ivan Solt, Wojciech Szlachta, Csilla Varnai, Steven Winfield.
! HJ X
! HJ X   Copyright 2006-2010.
! HJ X
! HJ X   These portions of the source code are released under the GNU General
! HJ X   Public License, version 2, http://www.gnu.org/copyleft/gpl.html
! HJ X
! HJ X   If you would like to license the source code under different terms,
! HJ X   please contact Gabor Csanyi, gabor@csanyi.net
! HJ X
! HJ X   Portions of this code were written by Noam Bernstein as part of
! HJ X   his employment for the U.S. Government, and are not subject
! HJ X   to copyright in the USA.
! HJ X
! HJ X
! HJ X   When using this software, please cite the following reference:
! HJ X
! HJ X   http://www.libatoms.org
! HJ X
! HJ X  Additional contributions by
! HJ X    Alessio Comisso, Chiara Gattinoni, and Gianpietro Moras
! HJ X
! HJ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

program lotf_metapot

  use libAtoms_module
  use QUIP_Module

  implicit None

  type(Atoms) :: dia, at
  type(InOutput) :: xml
  type(Potential) :: pot1, pot2
  real(dp), allocatable :: f(:,:), f_hyb(:,:)
  integer, pointer :: hybrid(:)
  integer :: i, j, n_cycle, n_extrap
  type(Table) :: embedlist
  type(MetaPotential) :: lotf, forcemix
  type(DynamicalSystem) :: ds, ds_saved

  call system_initialise(PRINT_NORMAL,seed=2)
  call initialise(xml, 'lotf_params.xml')
  call initialise(pot1, 'IP SW label="PRB_31_plus_H"', xml)
  call rewind(xml)
  call initialise(pot2, 'IP SW label="eps_2.3"', xml)
  
  call diamond(dia, 5.44_dp, (/14/))
  call supercell(at, dia, 2, 2, 2)
  call randomise(at%pos, 0.1_dp)

  n_cycle = 1
  n_extrap = 10

  allocate(f(3,at%N),f_hyb(3,at%N))

  call set_cutoff(at, cutoff(pot1)+2.0_dp)
  call calc_connect(at)

  ! Mark QM region
  call add_property(at, 'hybrid', 0)
  if (.not. assign_pointer(at, 'hybrid', hybrid)) &
       call system_abort('Cannot assign hybrid pointer')

  call append(embedlist, (/1,0,0,0/))
  call bfs_grow(at, embedlist, 2, nneighb_only=.true.)

  hybrid(int_part(embedlist,1)) = 1

  call print('embed list')
  call print(embedlist)

  ! Set up metapotentials for LOTF and for force mixing with the same buffer size
  call initialise(lotf, 'ForceMixing method=lotf_adj_pot_svd fit_hops=3 buffer_hops=3 '//&
       'randomise_buffer=T qm_args_str={single_cluster=T cluster_calc_connect=T}', pot1, pot2, dia)
  call initialise(forcemix, 'ForceMixing method=force_mixing buffer_hops=3 '//&
       'qm_args_str={single_cluster=T cluster_calc_connect=T}', pot1, pot2)

  call initialise(ds, at)
  
  call print_title('LOTF Metapotential')
  call print(lotf)

  call print_title('Force Mixing Metapotential')
  call print(forcemix)

  ! Predictor-corrector dynamics
  call Print_Title('Bootstrap')
  call calc(lotf, ds%atoms, f=f) ! bootstrap the adjustable potential
  do i=1,n_cycle

     ! Update QM region
     call print_title('Quantum Selection')
     ! ...
     

     ! Extrapolation
     call print_title('Extrapolation')
     call ds_save_state(ds_saved, ds)

     do j=1,n_extrap
        if (j == 1) then
           call calc(lotf, ds%atoms, f=f, args_str="lotf_do_qm=F lotf_do_init=T lotf_do_map=T")
        else
           call calc(lotf, ds%atoms, f=f, args_str="lotf_do_qm=F lotf_do_init=F")
        end if
        call calc(forcemix, ds%atoms, f=f_hyb)
        call advance_verlet(ds, 1.0_dp, f)
        call ds_print_status(ds, 'E', instantaneous=.true.)
        call print('E err '//ds%t//' '//rms_diff(f_hyb, f)//' '//maxval(abs(f_hyb-f)))
     end do

     call print_title('QM Force Evaluation')
     ! QM force calculation and optimisation of adj pot
     call calc(lotf, ds%atoms, f=f, args_str="lotf_do_qm=T lotf_do_init=F lotf_do_fit=T")

     ! Interpolation
     call print_title('Interpolation')
     call ds_restore_state(ds, ds_saved)
     do j=1,n_extrap
        call calc(lotf, ds%atoms, f=f, args_str="lotf_do_qm=F lotf_do_init=F lotf_do_interp=T lotf_interp="&
             //(real(j-1,dp)/real(n_extrap,dp)))
        call calc(forcemix, ds%atoms, f=f_hyb)
        call advance_verlet(ds, 1.0_dp, f)
        call ds_print_status(ds, 'I', instantaneous=.true.)
        call print('I err '//ds%t//' '//rms_diff(f_hyb, f)//' '//maxval(abs(f_hyb-f)))
     end do

     call print_title('Connectivity update')
     call calc_connect(ds%atoms)

  end do

  deallocate(f, f_hyb)
  call finalise(xml)
  call finalise(pot1)
  call finalise(pot2)
  call finalise(at)
  call finalise(ds)
  call finalise(dia)
  call system_finalise
  
end program lotf_metapot
