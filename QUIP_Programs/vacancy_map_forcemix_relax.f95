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

#ifdef HAVE_HYBRID
#ifdef HAVE_LOTF
program vacancy_map_forcemix_relax
use libatoms_module
use potential_module
use vacancy_map_module

implicit none
  type(Atoms) at
  type(Potential), target ::pot1, pot2
  type(Potential)         ::lotfpot
  type(inoutput) params, out
  integer :: it
  real(dp) :: e, maxforce, d
  real(dp) :: hybrid_center(3)
  integer :: hybrid_center_i
  real(dp), allocatable :: f(:,:)
  real(dp) :: lat(3,3), virial(3,3)
  character(len=256) arg
  type(MPI_Context) :: mpi
  integer vac_i
  character(len=1024) :: comment, in_file
  type(Dictionary) :: cli_params, comment_params
  logical :: doing_restart, use_n_minim
  real(dp) :: eps_guess
  character(len=1024) :: lotf_args_str
  integer :: buffer_hops
  integer :: n_embed
  type(Table) :: embedlist, fitlist
  integer, pointer :: hybrid_mark_p(:)
  integer :: i, ii

  call system_initialise(enable_timing=.true.)
  ! trace_memory=.true.

  call initialise(cli_params)
  in_file=''
  lotf_args_str=''
  call param_register(cli_params, 'restart', 'F', doing_restart, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params, 'use_n_minim', 'F', use_n_minim, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params, 'in_file', 'stdin', in_file, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params, 'vac_i', '0', vac_i, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params, 'buffer_hops', '3', buffer_hops, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params, 'lotf_args_str', '', lotf_args_str, help_string="No help yet.  This source file was $LastChangedBy$")

  if (.not. param_read_args(cli_params, do_check=.true.)) then
    call system_abort("Usage: vacancy_map_forcemix_relax [ restart=L(F) ] [ use_n_minim=L(F) ] [ vac_i=n(0) ] [ in_file=filename(stdin) ] [ buffer_hops=n(3) ] [ lotf_args_str=args ]")
    call system_abort("Confused by CLI argument")
  endif
  call finalise(cli_params)

  call print("in_file = "//trim(in_file))
  call print("restart = "//doing_restart)
  call print("use_n_minim = "//use_n_minim)
  call print("vac_i = "//vac_i)
  call print("lotf_args_str = "//trim(lotf_args_str))

  if (doing_restart) then
    eps_guess = 1.0e-2_dp
  else
    eps_guess = 0.0_dp
  endif

  call initialise(mpi)
  call Potential_Initialise_filename(pot1, 'TB NRL-TB label=Aluminum', 'quip_params_hybrid_relax.xml', mpi_obj=mpi)
  call Potential_Initialise_filename(pot2, 'IP EAM_Ercolessi_Adams', 'quip_params_hybrid_relax.xml', mpi_obj=mpi)

  call init_lotf_forcemix(pot1, pot2, lotfpot, buffer_hops, trim(lotf_args_str))

  call initialise(out, 'out.xyz', action=OUTPUT)
  call read_xyz(at, in_file, comment=comment, mpi_comm=mpi%communicator)

  if (vac_i < 0 .or. vac_i > at%N) call system_abort("vac_i " // vac_i // " out of range : 0 or 1 .. N=" // at%N)

  call initialise(comment_params)
  call param_register(comment_params,"hybrid_center_i", "0", hybrid_center_i, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(comment_params,"hybrid_center", "1e38 1e38 1e38", hybrid_center, help_string="No help yet.  This source file was $LastChangedBy$")
  if (.not. param_read_line(comment_params, comment, ignore_unknown=.true.)) then
    call print("Failed to parse hybrid_center_i and hybrid_center from comment line '"//trim(comment)//"'")
  endif
  call finalise(comment_params)

  if (vac_i > 0) then
    call print("Removing atom " // vac_i // " from pos " // at%pos(:,vac_i))
    call remove_atoms(at, vac_i)
  endif

  call print_xyz(at, out, comment=trim(comment)//" Unrelaxed_hybrid_disloc", all_properties=.true.)

  allocate(f(3,at%N))

  if (.not. assign_pointer(at, "hybrid_mark", hybrid_mark_p)) &
    call system_abort("at had no hybrid_mark field")
  n_embed = count(hybrid_mark_p == HYBRID_REGION1_MARK)
  call table_allocate(embedlist,4,0,0,0)
  call table_allocate(fitlist,4,0,0,0)
  ii = 0
  do i=1, at%N
    if (hybrid_mark_p(i) == HYBRID_REGION1_MARK) then
      ii = ii + 1
      call append(embedlist, (/ i, 0, 0, 0 /) )
      call append(fitlist, (/ i, 0, 0, 0 /) )
    endif
  end do
  call set_embed(lotfpot, embedlist)
  call set_fit(lotfpot, fitlist)

  call set_cutoff(at, cutoff(pot2)+0.5_dp)
  call calc_connect(at)
  call verbosity_push_increment(5)
  call setup_parallel(pot2, at, e=e, f=f)
  it = minim(lotfpot, at, 'cg', 0.001_dp, 500, 'LINMIN_DERIV', do_print = .true., print_inoutput = out, &
    do_pos = .true., eps_guess = eps_guess, use_n_minim=use_n_minim)
  call verbosity_pop()
  call print_xyz(at, out, comment=trim(comment)//" Relaxed_hybrid_disloc", all_properties=.true.)
  call calc_dists(at)

  call calc(lotfpot, at, e=e,f=f)
  maxforce = maxval(norm(f,1))
  deallocate(f)

  call print('Atom '//hybrid_center_i//' position: '//hybrid_center//' Energy: '//e// &
    ' Max force after '//it//' iterations: '//maxforce)

  call system_finalise()
end program
#else
vacancy_map_forcemix_relax without HAVE_HYBRID and HAVE_LOTF makes no sense
#endif /* HAVE_LOTF */
#else
vacancy_map_forcemix_relax without HAVE_HYBRID and HAVE_LOTF makes no sense
#endif /* HAVE_HYBRID */
