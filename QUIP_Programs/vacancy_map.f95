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

program test
use libatoms_module
use potential_module
use atoms_minimisation_module

implicit none
  type(Atoms) at, at2, at3, bulk
  type(Potential) pot
  type(inoutput) params, out, in
  integer:: it, i
  logical:: status
  real(dp) :: e, ebulk, maxforce, center(3), d
  real(dp), allocatable :: f(:,:)
  real(dp) :: lat(3,3), virial(3,3)
  integer n_groups, i_group, group_skip, ii, i_do
  character(len=256) arg

  call system_initialise()

  call Initialise(params, "quip_params.xml")
  call Initialise(pot, 'IP EAM_Ercolessi_Adams', params)
  !call print(pot)
  call initialise(out, 'out.xyz')

  call initialise(in, 'stdin')
  call read_xyz(at, in)
  call set_cutoff(at, cutoff(pot)+0.5_dp)

  allocate(f(3,at%N))
  center(:) = 0.0_dp

  if (cmd_arg_count() == 0) then
    n_groups = 0
  else if (cmd_arg_count() == 2 .or. cmd_arg_count() == 3) then
    call get_cmd_arg(1, arg)
    read (arg, *) n_groups
    call get_cmd_arg(2, arg)
    read (arg, *) i_group
    if (i_group < 1 .or. i_group > n_groups) call system_abort("i_group("//i_group//") out of range of 1..n_groups ("//n_groups//")")
    call print ("dividing atoms into " // n_groups // " groups, doing group " // i_group)
    if (cmd_arg_count() == 3) then
      call get_cmd_arg(3, arg)
      read (arg, *) group_skip
    else
      group_skip = 0
    endif
  else
    call system_abort("Usage: vacancy_map [ n_groups i_group [ group_skip ] ]")
  endif

  call system_timer("initial_relaxation")

  ! compute bulk energy
  lat(:,1) = (/2_dp, 0_dp, 2_dp/)
  lat(:,2) = (/2_dp, 2_dp, 0_dp/)
  lat(:,3) = (/0_dp, 2_dp, 2_dp/)
  call initialise(bulk,1, lat)
  bulk%pos(:,1) = (/ 0.0_dp, 0.0_dp, 0.0_dp /)
  bulk%Z(1) = 13
  call set_cutoff(bulk, cutoff(pot)+0.5_dp)
  call calc_connect(bulk)
  it = minim(bulk, pot, 'cg', 'NR_LINMIN', 0.01_dp, 100, do_print = .false., print_inoutput = out, do_lat = .true.)
  call calc(pot, bulk, e=ebulk)
  call print("bulk")
  call print(bulk)
  call print("Bulk energy " // ebulk)

  call calc_connect(at)

  ! first relax the initial starting config
  call calc(pot, at, e=e,f=f,virial=virial)
  call print ("Start unrelaxed E " // e // " max f component " // maxval(abs(f)) // " max virial component " // maxval(abs(virial)))
  call print_xyz(at, out, "Dislocation unrelaxed")

  it = minim(at, pot, 'cg', 'NR_LINMIN', 0.01_dp, 100, do_print = .false., print_inoutput = out, do_lat = .true.)
  it = minim(at, pot, 'cg', 'NR_LINMIN', 0.01_dp, 100, do_print = .false., print_inoutput = out, do_pos = .true.)
  it = minim(at, pot, 'cg', 'NR_LINMIN', 0.01_dp, 100, do_print = .false., print_inoutput = out, do_lat = .true.)
  it = minim(at, pot, 'cg', 'NR_LINMIN', 0.01_dp, 100, do_print = .false., print_inoutput = out, do_pos = .true.)

  call calc(pot, at, e=e,f=f,virial=virial)
  call print ("Start relaxed E " // e // " max f component " // maxval(abs(f)) // " max virial component " // maxval(abs(virial)))

  call print_xyz(at, out, "Dislocation relaxed")

  deallocate(f)

  call system_timer("initial_relaxation")

  ! compute vacancy calculations

  ii = 0
  i_do = 0
  do i=1, at%N
     d = distance_min_image(at, i,center) 
     if(d > 40.0_dp) cycle

     ii = ii + 1

     if (n_groups > 0) then
       if (mod(ii-1, n_groups)+1 /= i_group) cycle
     endif

     i_do = i_do + 1
     if (n_groups > 0 .and. i_do <= group_skip) cycle

     call print("")
     call print("Starting vacancy calulation corresponding to atom "//i)

     call system_timer("vacancy_relaxation")
     call supercell(at2, at, 1,10,1)
     call remove_atoms(at2, i)
     call calc_connect(at2)
     call print_xyz(at2, out, "Vacancy unrelaxed")
     !call verbosity_push_increment()
     it = minim(at2, pot, 'cg', 'NR_LINMIN', 0.01_dp, 100, do_print = .false., print_inoutput = out, do_pos = .true.)
     !call verbosity_pop()
     call print_xyz(at2, out, "Vacancy relaxed")
     call calc_dists(at2)
     allocate(f(3,at2%N-1))
     call calc(pot, at2, e=e,f=f)
     maxforce = maxval(norm(f,1))
     deallocate(f)
     call print('Atom '//i//' vacancy position: '//at%pos(1,i)//' '//at%pos(2,i)//' '//at%pos(3,i)//' Energy: '//e//' Max force after '//it//' iterations: '//maxforce)
    call system_timer("vacancy_relaxation")
    call print("")
  end do

  call system_finalise()
end program
