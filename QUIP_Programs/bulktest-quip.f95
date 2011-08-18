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

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X     When using this software, the following should be referenced:
!X
!X     Gabor Csanyi, Tristan Albaret, Mike C. Payne and Alessandro De Vita
!X     "Learn on the fly": a hybrid classical and quantum-mechanical
!X         molecular dynamics simulation
!X     Physical Review Letters 93 p. 175503 (2004) >>PDF [626 KB]
!X
!X     Gabor Csanyi, T. Albaret, G. Moras, M. C. Payne, A. De Vita
!X     Multiscale hybrid simulation methods for material systems
!X     J. Phys. Cond. Mat. 17 R691-R703 Topical Review (2005)
!X
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

program bulktest

  use libAtoms_module
  use potential_module

  implicit none

  type(Dictionary) :: params
 
  real(dp) :: time_step, init_temp, sim_temp
  integer :: seed, fit_hops, embed_hops
  character(FIELD_LENGTH) :: classicalpot_args, qmpot_args
  character(FIELD_LENGTH) :: xml, pot_args

  type(DynamicalSystem) :: ds
  type(Atoms) :: at, dia
  type(inoutput) :: movie, xmlfile
  type(table) :: embedlist, fitlist
  real(dp), allocatable :: f(:,:)
  integer :: i, step

  type(Potential) :: classicalpot, qmpot
  type(Potential) :: pot
  

  ! initialise program
  call system_initialise(PRINT_NORMAL)
  call initialise(movie, "movie.xyz")

  ! Setup parameters
  call param_register(params, 'time_step', '1.0', time_step, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'seed', '0', seed, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'init_temp', '300.0', init_temp, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'sim_temp', '300.0', sim_temp, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'classicalpot', 'IP SW', classicalpot_args, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'qmpot', 'TB Bowler',  qmpot_args, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'xml', 'lotf.xml', xml, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'pot', 'LOTF buffer_hops=3 small_clusters=F', pot_args, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'embed_hops', '2', embed_hops, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'fit_hops', '3', fit_hops, help_string="No help yet.  This source file was $LastChangedBy$")

  if (.not. param_read_args(params, (/ (i, i=1,cmd_arg_count()) /), .true.)) &
       call system_abort('Error reading command line arguemnts')

  call initialise(xmlfile, xml)
  call initialise(classicalpot, classicalpot_args, xmlfile)
  call rewind(xmlfile)
  call initialise(qmpot, qmpot_args, xmlfile)
  call finalise(xmlfile)

  call initialise(pot, pot_args, classicalpot, qmpot)
  call print(pot)

  ! create some atoms
  call diamond(dia, 5.44_dp)
  call supercell(at, dia, 5,5,5)
  call set_atoms(at, 14)
  call atoms_set_cutoff(at, 4.0_dp)
  call randomise(at%pos, 0.01_dp)
  call calc_connect(at)

  allocate(f(3,at%N))

  ! initialise dynamics
  call initialise(ds, at)
  call rescale_velo(ds, init_temp)
  call zero_momentum(ds)

  ds%thermostat = THERMOSTAT_LANGEVIN_THERM
  ds%thermal_tau = 1000.0_dp
  ds%sim_temp = sim_temp

  ! create list of embedded atoms
  call append(embedlist, (/1,0,0,0/))
  call BFS_grow(ds%atoms, embedlist, embed_hops)

  ! grow the embed list to include a fit zone
  fitlist = embedlist
  call BFS_grow(ds%atoms, fitlist, fit_hops)

  call list_to_property(ds%atoms, embedlist, 'embed')
  call list_to_property(ds%atoms, fitlist, 'fit')

  call set_embed(pot, embedlist)
  call set_fit(pot, fitlist)

  ! Do the dynamics
  step = 0
  do while (.true.)
     call calc(pot, ds%atoms, f=f)
     call advance_verlet(ds, time_step, f)
     call ds_print_status(ds, 'D')
     if (mod(step,100) == 0) then
        call print_xyz(ds%atoms, movie, real_format='f12.8', all_properties=.true.)
     end if
     step = step + 1
  end do
  
  call adjustable_potential_finalise()
  call finalise(at)
  call finalise(ds)
  call finalise(movie)
  call finalise(embedlist)
  call finalise(fitlist)
  call finalise(pot)
  call finalise(classicalpot)
  call finalise(qmpot)

  call system_finalise()

end program bulktest
