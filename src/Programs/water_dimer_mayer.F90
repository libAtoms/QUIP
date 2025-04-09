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

#include "error.inc"

program mc
  use libatoms_module
  use potential_module
  implicit none

  integer :: Nsobol
  character(len=FIELD_LENGTH) :: sobolfilename
  real(dp), allocatable :: sobol(:,:)
  type(Potential) :: pot
  type(MPI_context) :: mpi_glob
  type(Cinoutput) :: atoms_out_file, atoms_in_file
  type(Atoms) :: at
  type(Dictionary) :: params
  type(extendable_str) :: params_es
  character(len=FIELD_LENGTH) :: atoms_in_filename, atoms_out_filename, params_in_filename
  character(len=FIELD_LENGTH) :: pot_init_args, pot_calc_args
  integer ::  i,j
  integer :: error = ERROR_NONE
  real(dp) :: E
  logical :: do_timing, quiet_calc, have_Roo
  real(dp), pointer :: pos(:,:), force(:,:)
  real(dp) :: Roo, drOH, dHOH, dr(3)
  real(dp), allocatable :: Roo_array(:), E_Array(:)
 

  call system_initialise()
  call initialise(mpi_glob)
  
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  ! get program parameters
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  call initialise(params)
  call param_register(params, 'atoms_in_file', 'stdin', atoms_in_filename, help_string="Initial atomic data filename")
  call param_register(params, 'atoms_out_file', 'mc.xyz', atoms_out_filename, help_string="Output filename for MC samples")
  call param_register(params, 'params_in_file', 'quip_params.xml', params_in_filename, help_string="QUIP XML parameter filename")
  call param_register(params, 'pot_init_args', PARAM_MANDATORY, pot_init_args, help_string="string to initialise potential")
  call param_register(params, 'pot_calc_args', '', pot_calc_args, help_string="args string for potential calc")
  call param_register(params, 'quiet_calc', 'T', quiet_calc, help_string="do calc() quietly")
  call param_register(params, 'do_timing', 'F', do_timing, help_string="if true, do timing")
  call param_register(params, 'Roo', '3.0', Roo, has_value_target=have_Roo, help_string="set the O-O distance to this")
  call param_register(params, 'drOH', PARAM_MANDATORY, drOH, help_string="set the max deviation of rOH (angstroms)")
  call param_register(params, 'dHOH', PARAM_MANDATORY, dHOH, help_string="set the max deviation of HOH (degrees)")
  call param_register(params, 'sobolfile', 'sobol12', sobolfilename, help_string="name of file with sobol integration grid points")
  call param_register(params, 'Nsobol', PARAM_MANDATORY, Nsobol, help_string="Number of sobol integration grid points")
  

  ! parse params
  if (.not. param_read_args(params)) then
    call param_print_help(params)
    call system_abort("Error reading params from command line")
  endif

  if (len(trim(pot_init_args)) == 0) then
     call param_print_help(params)
     call system_abort("get_params() got empty pot_init_args")
  end if
  call finalise(params)



  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  ! Initialise simulation
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  ! Turn on timing
  if (do_timing) call enable_timing()

  ! Read in atoms
  call initialise(atoms_in_file, atoms_in_filename, INPUT, mpi=mpi_glob)
  call read(at, atoms_in_file, error=error)
  HANDLE_ERROR(error)

  ! Read in QUIP parameters and initialise potential
  if (len_trim(params_in_filename) > 0) then
     call read(params_es, params_in_filename, convert_to_string=.true., mpi_comm=mpi_glob%communicator, mpi_id=mpi_glob%my_proc)
  else
     call initialise(params_es)
  endif
  call initialise(pot, args_str=pot_init_args, param_str=string(params_es), mpi_obj=mpi_glob)
  call print(pot)

  call finalise(params_es)

  call initialise(atoms_out_file, atoms_out_filename, OUTPUT, mpi=mpi_glob)

  call set_cutoff(at, cutoff(pot))

  call print('Sobol integration using '//Nsobol//' points')
  call print('drOH = '//drOH)
  call print('dHOH = '//dHOH)

  ! read in Sobol sequence
  allocate(sobol(Nsobol,12))
  open(unit=21, file=trim(sobolfilename))
  do i=1,Nsobol
     read(21,*) sobol(i,1), sobol(i,2), sobol(i,3), sobol(i,4), sobol(i,5), sobol(i,6), sobol(i,7), sobol(i,8), sobol(i,9), sobol(i,10), sobol(i,11), sobol(i,12) 
  end do


  ! calculate connection
  call calc_connect(at)

  !call assign_pointer(at, 'pos', pos)
  !call assign_pointer(at, 'force', force)



  ! set up Roo grid

  allocate(Roo_array(65))
  allocate(E_array(65))
  Roo_array(1) = 2.1_dp
  do i=2,40
     Roo_array(i) = Roo_array(i-1)*1.02_dp
  end do
  do i=41,size(Roo_array)
     Roo_array(i) = Roo_array(i-1)*1.05_dp
  end do

  call print('Roo array:')
  call print(Roo_array)

  call system_timer("sobol_loop")


  do i=1,Nsobol
     at%pos(:,1) = (/0.0_dp, 0.0_dp, 0.0_dp/)
     at%pos(:,4) = (/0.0_dp, 0.0_dp, 0.0_dp/)
     call rotate_dimer(at, sobol(i,:))

     Roo=0.0_dp
     do j=1,size(Roo_array)

        ! shift to the next Roo distance
        at%pos(:,4) = at%pos(:,4) + (/0.0_dp, 0.0_dp, Roo_array(j)-Roo/)
        at%pos(:,5) = at%pos(:,5) + (/0.0_dp, 0.0_dp, Roo_array(j)-Roo/)
        at%pos(:,6) = at%pos(:,6) + (/0.0_dp, 0.0_dp, Roo_array(j)-Roo/)
        Roo = Roo_array(j)

        ! calculate new energy and force
        if (quiet_calc) call verbosity_push_decrement()
        call calc(pot, at, args_str=trim(pot_calc_args)//" energy")
        if (quiet_calc) call verbosity_pop()
        
        call get_param_value(at, "energy", E_array(j))

        !call print(pos)
        !call print(force)
        
        !call write(atoms_out_file, at, real_format='%18.10f')
     end do
     call print("E= "//E_array)
  end do
  call system_timer("sobol_loop")

  deallocate(sobol)
  deallocate(Roo_array)
  deallocate(E_array)

  call system_finalise()

contains

  subroutine rotate_dimer(at, vec)
    type(Atoms), intent(inout) :: at
    real(dp), intent(in) :: vec(12)
    type(Quaternion) :: quat
    real(dp), parameter :: rOH = 0.958_dp, HOH = 104.5_dp

    ! Atoms are in the order OHHOHH
       
    ! change OH distances
    at%pos(:,2) = at%pos(:,1) + (/0.0_dp, 0.0_dp, rOH + (vec(1)-0.5_dp)*2.0_dp*drOH/)
    at%pos(:,3) = at%pos(:,1) + (/0.0_dp, 0.0_dp, rOH + (vec(2)-0.5_dp)*2.0_dp*drOH/)
    at%pos(:,5) = at%pos(:,4) + (/0.0_dp, 0.0_dp, rOH + (vec(3)-0.5_dp)*2.0_dp*drOH/)
    at%pos(:,6) = at%pos(:,4) + (/0.0_dp, 0.0_dp, rOH + (vec(4)-0.5_dp)*2.0_dp*drOH/)
    
    ! change OH angles
    quat = rotation( (/1.0_dp,0.0_dp,0.0_dp/), (HOH + (vec(5)-0.5_dp)*2.0_dp*dHOH)/180.0_dp*PI)
    at%pos(:,2) = at%pos(:,1)+ rotate(at%pos(:,2)-at%pos(:,1), quat)
    
    quat = rotation( (/1.0_dp,0.0_dp,0.0_dp/), (HOH + (vec(6)-0.5_dp)*2.0_dp*dHOH)/180.0_dp*PI)
    at%pos(:,5) = at%pos(:,4)+ rotate(at%pos(:,5)-at%pos(:,4), quat)

    ! uniform random rotation. formula from:
    ! K. Shoemake. Uniform random rotations. In D. Kirk, editor, Graphics Gems III, pages 124-132. Academic, New York, 1992.
    ! via  http://planning.cs.uiuc.edu/node198.html
    
    ! rotate first water monomer around its oxygen
    call initialise(quat, sqrt(1.0_dp-vec(7))*sin(2.0_dp*PI*vec(8)), sqrt(1.0_dp-vec(7))*cos(2.0_dp*PI*vec(8)), sqrt(vec(7))*sin(2.0_dp*PI*vec(9)), sqrt(vec(7))*cos(2.0_dp*PI*vec(9))) ! quaternion representing random rotation
    at%pos(:,2) = at%pos(:,1)+rotate(at%pos(:,2)-at%pos(:,1),quat)
    at%pos(:,3) = at%pos(:,1)+rotate(at%pos(:,3)-at%pos(:,1),quat)
       
    ! rotate second water monomer around its oxygen
    call initialise(quat, sqrt(1.0_dp-vec(10))*sin(2.0_dp*PI*vec(11)), sqrt(1.0_dp-vec(10))*cos(2.0_dp*PI*vec(11)), sqrt(vec(10))*sin(2.0_dp*PI*vec(12)), sqrt(vec(10))*cos(2.0_dp*PI*vec(12))) ! quaternion representing random rotation
    at%pos(:,5) = at%pos(:,4)+rotate(at%pos(:,5)-at%pos(:,4),quat)
    at%pos(:,6) = at%pos(:,4)+rotate(at%pos(:,6)-at%pos(:,4),quat)

    call finalise(quat)

  end subroutine rotate_dimer

end program mc
