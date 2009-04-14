!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X     Learn-on-the-fly (LOTF) hybrid molecular dynamics code
!X
!X    
!X     Authors: Gabor Csanyi, Alessio Commisso, Steven Winfield
!X     James Kermode, Gianpietro Moras, Michael Payne, Alessandro De Vita
!X
!X     
!X     Copyright 2005, All Rights Reserved 
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

! $Id: test_dimer.f95,v 1.3 2008-04-09 11:26:56 jrk33 Exp $

! $Log: not supported by cvs2svn $
! Revision 1.2  2008/04/09 11:16:39  jrk33
! Made compile with changes to code. Not tested.
!
! Revision 1.1.1.1  2008/03/13 17:36:35  gc121
! this module contains the LOTF top level programs
!
! Revision 1.11  2007/05/08 15:13:53  jrk33
! Uncommented code - Adjustable potential is now updated for new multiple image convention
!
! Revision 1.10  2007/05/02 23:45:27  gc121
! updated to reflect changes in libAtoms; AdjPot not done yet, so most of code is commented out for now
!
! Revision 1.9  2007/05/02 22:46:35  gc121
! started to update after the changes in libAtoms
!
! Revision 1.8  2006/06/20 17:23:19  gc121
! added new copyright notice to include James, Gian, Mike and Alessandro
!

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
  call atoms_set_cutoff(at, 6.0_dp)
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
