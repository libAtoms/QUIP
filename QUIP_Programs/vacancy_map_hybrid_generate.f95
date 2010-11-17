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
use vacancy_map_module

implicit none
  type(Atoms) at, at2, bulk, bulk1, bulk2
  type(Potential), target :: pot1, pot2
  type(Potential) :: hybridpot
  type(inoutput) params, out, in
  type(Table)    :: qmlist, core_list
  integer:: it, i, pass
  logical:: status
  integer, pointer :: hybrid_mark(:)
  integer center_i
  real(dp) :: e, ebulk, center(3), d
  real(dp), allocatable :: f(:,:)
  real(dp) :: lat(3,3), virial(3,3)
  character(len=256) arg
  type(MPI_Context) :: mpi
  integer :: imin, imax
  real(dp) :: dmax
  logical :: core_qm, vac_qm
  real(dp) :: core_qm_x, core_qm_z, r2
  type(Dictionary) :: cli_params
  character(len=1024) :: file_suffix
  integer :: n_supercell

  integer :: vac_qm_hops, buffer_hops

  call system_initialise(enable_timing=.true.,verbosity=PRINT_NORMAL)

  call initialise(mpi)

  call initialise(cli_params)
  call param_register(cli_params, "core_qm", "T", core_qm, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params, "buffer_hops", "1", buffer_hops, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params, "core_qm_x", "10.0", core_qm_x, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params, "core_qm_z", "5.0", core_qm_z, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params, "vac_qm", "F", vac_qm, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params, "vac_qm_hops", "1", vac_qm_hops, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params, "n_supercell", "4", n_supercell, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params, "imin", "0", imin, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params, "imax", "0", imax, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params, "dmax", "1.0e38", dmax, help_string="No help yet.  This source file was $LastChangedBy$")
  if (.not. param_read_args(cli_params)) then
    call print("Usage: vacancy_map_hybrid_generate [buffer_hops=i(1)]", PRINT_ALWAYS)
    call print("  [core_qm=l(T)] [core_qm_x=r(10.0)] core_qm_z=r(5.0)]", PRINT_ALWAYS)
    call print("  [vac_qm=l(F)] vac_qm_hops=i(1)] [n_supercell=i(4)]", PRINT_ALWAYS)
    call print("  [imin=i(0)] [imax=i(0)] [dmax=r(1.0e38)]", PRINT_ALWAYS)
    call system_abort("Confused by CLI arguments")
  end if
  call finalise(cli_params)

  call print("Using core_qm " // core_qm // " vac_qm " // vac_qm)
  if (core_qm .or. vac_qm) then
    call print("buffer_hops " // buffer_hops)
    if (core_qm) then
      call print("Using core_qm_x " // core_qm_x // " core_qm_z " // core_qm_z)
    endif
    if (vac_qm) then
      call print("Using vac_qm_hops " // vac_qm_hops)
    endif
  endif

  call Initialise(params, "quip_params_hybrid_generate.xml")
  call Initialise(pot1, 'TB NRL-TB label=Aluminum SCF=GLOBAL_U GLOBAL_U=20', params, mpi)
  call rewind(params)
  call Initialise(pot2, 'IP EAM_Ercolessi_Adams', params)
  call finalise(params)

  call print("Doing init_hybrid")
  call init_hybrid(pot1, pot2, hybridpot)

  ! read initial config
  call initialise(out, 'out.xyz', action=OUTPUT)
  call initialise(in, 'start.xyz')
  call read_xyz(at, in)
  call finalise(in)

  call add_property(at, 'hybrid_mark', 0)
  call add_property(at, 'weight_region1', 0.0_dp)

  !call print("minimizing initial config using mm")
  !call set_cutoff(at, cutoff(pot2)+0.5_dp)
  !call calc_connect(at)
  !it = minim(hybridpot, at, 'cg', 0.01_dp, 100, 'NR_LINMIN', do_print = .false., do_lat = .true., do_pos = .true.)

  allocate(f(3,at%N))
  call set_cutoff(at, cutoff(pot2))
  call calc_connect(at)
  call calc(hybridpot, at, e=e,f=f,virial=virial)
  call print ("Start relaxed E " // e // " max f component " // maxval(abs(f)) // " max virial component " // maxval(abs(virial)))
  deallocate(f)

  call print_xyz(at, out, "Dislocation_relaxed")

  ! compute vacancy starting configs with weights

  call supercell(at2, at, 1,n_supercell,1) ! length along the dislocation
  call set_cutoff(at2, 3.2_dp)
  call calc_connect(at2)

  center(:) = 0.0_dp

  if (.not. assign_pointer(at2, 'hybrid_mark', hybrid_mark)) &
    call system_abort("Couldn't find hybrid_mark in at2")

  hybrid_mark = HYBRID_NO_MARK

  mainlog%default_real_precision=2
  file_suffix = ""
  if (core_qm .or. vac_qm) then
    file_suffix = trim(file_suffix) // "buffer_hops_" // buffer_hops
    if (core_qm) then
      file_suffix = trim(file_suffix) // ".core_qm_x_"//core_qm_x//".core_qm_z_"//core_qm_z
    endif
  endif
  mainlog%default_real_precision=16

  ! core qm region
  if(core_qm) then

     call wipe(core_list)
     do i=1,at2%N
        r2 = (at2%pos(1,i)/core_qm_x)**2+(at2%pos(3,i)/core_qm_z)**2
        if(r2 < 1.0_dp) call append(core_list, (/i, 0, 0 ,0/))
     end do

     call print(core_list, PRINT_VERBOSE)
     
     hybrid_mark(int_part(core_list, 1)) = HYBRID_ACTIVE_MARK

     call create_local_energy_weights(at2, 1, buffer_hops) ! trans_width, buffer_width

     call print_xyz(at2, "start."//trim(file_suffix)//".xyz", properties="pos:Z:hybrid_mark:weight_region1")
  end if


  if(vac_qm) then
  
     ! add vacancy qm region
     if (imin <= 0) imin = 1
     if (imax > at2%N) imax = at%N
     do i=imin, imax
        
        d = distance_min_image(at, i, center) 
        if (d > dmax) cycle
        
        call print("Making config with hybrid centered around atom " // i)
        
        hybrid_mark = HYBRID_NO_MARK
        
        call wipe(qmlist)

        call append(qmlist, (/i, 0,0, 0/))
        call bfs_grow(at2, qmlist, vac_qm_hops, nneighb_only = .false.) ! qm region

        hybrid_mark(int_part(qmlist, 1)) = HYBRID_ACTIVE_MARK
        if(core_qm) hybrid_mark(int_part(core_list, 1)) = HYBRID_ACTIVE_MARK
        
        call create_local_energy_weights(at2, 1, buffer_hops) ! trans_width, buffer_width
        
        call print_xyz(at2, "start.vac_i_"//i//".vac_qm_hops_"//vac_qm_hops//"."//trim(file_suffix)//".xyz", &
             comment="hybrid_center_i="//i//" hybrid_center='"//at2%pos(:,i)//"'", &
             properties="pos:Z:hybrid_mark:weight_region1")
        
     end do
  end if

  call system_finalise()
end program
