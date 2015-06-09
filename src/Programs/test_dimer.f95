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
!X     Learn-on-the-fly (LOTF) hybrid molecular dynamics code
!X
!X     This source code is confidential, all distribution is 
!X     prohibited. Making unauthorized copies is also prohibited
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

program test_dimer

  use libAtoms_module
  use LOTF_module

  implicit none

  type(Atoms)::at
  real(dp)::f(3,2)
  type(table)::fitlist
  integer::i

  call system_initialise()

  call atoms_initialise(at, 2, reshape((/4.0_dp,0.0_dp,0.0_dp,0.0_dp,4.0_dp,0.0_dp,0.0_dp,0.0_dp,4.0_dp/), (/3,3/)))
  at%Z = (/14,14/)
  at%pos = 0.0_dp
  at%pos(1,1) = 1.0_dp
  at%pos(1,2) = -1.0_dp
  call set_cutoff(at, 6.0_dp)
  call calc_connect(at)
  call print(at)

  write(line, *) "Dimer energies and forces:" ; call print(line)
  do i=1,150
     at%pos(1,2) = 2.0_dp+i*0.01_dp
     call calc_dists(at)
     call SW_force_noopt(at, f)
     write(line, *) at%pos(1,2), SW_energy(at), f(1,1) ; call print(line)
  end do
  write(line, *) ; call print(line)

  at%pos(1,2) = 2.3_dp
  call calc_dists(at)
  call print(at)
  
  write(line, *) "SW energy: ", SW_energy(at);     call print(line)
  write(line, *) "SW force:"; call print(line)
  call SW_force_noopt(at, f)
  call print(transpose(f))
  write(line, *) ; call print(line)

  call append(fitlist, 1,reshape(f(:,1), (/3/)))
  call append(fitlist, 2,reshape(f(:,2), (/3/)))

  call print('fitlist')
  call print(fitlist)

  call adjustable_potential_init(at, fitlist, 0)
  call adjustable_potential_optimise(at, real_part(fitlist))
  call adjustable_potential_print_params()

  at%pos(1,2) = 2.35_dp
  call calc_dists(at)  
  write(line, *) "Adjustable potential at d = ", distance(at, 1, 2, (/0,0,0/)) ; call print(line)
  call adjustable_potential_force(at, f)
  call print(transpose(f))
  call finalise(fitlist)

  call finalise(at)
  call system_finalise()
end program test_dimer
