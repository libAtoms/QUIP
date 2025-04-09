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

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X cp2k_filepot program
!X
!% outdated filepot program for CP2K driver using hard-coded input file 
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

!FilePotential for the CP2K driver
!FilePot uses go_cp2k from the cp2k_driver_module
!reads filepot.0.xyz  >>  runs CP2K  >>  writes filepot.0.out with forces

!Initialisation:
!For QS or MM calculation:
!   call initialise(CP2K_pot,'FilePot command=/Users/csilla/QUIP/build.darwin_x86_64_g95/cp2k_filepot property_list=pos min_cutoff=0.0')
!For QM/MM calculation:
!   call initialise(CP2K_pot,'FilePot command=/Users/csilla/QUIP/build.darwin_x86_64_g95/cp2k_filepot property_list=pos:QM_flag min_cutoff=0.0')

!Calc (args_str is the same as for go_cp2k):
!   call calc(CP2K_pot,ds%atoms,energy=energy,forces=f1,args_str=trim(args_str))


program cp2k_filepot

  use libatoms_module
  use cp2k_driver_module
  use potential_module
!  use atoms_module
!  use topology_module

  implicit none

    type(Atoms)                   :: my_atoms
    real(dp)                      :: energy
    real(dp), allocatable         :: f0(:,:)
    real(dp), pointer             :: forces_p(:,:)
    type(InOutput)                :: xyz
 
    !Input
    type(Dictionary)              :: params_in
    character(len=STRING_LENGTH)  :: args_str
    character(len=STRING_LENGTH)   :: Run_Type
    character(len=STRING_LENGTH)   :: Print_PSF
    character(len=STRING_LENGTH)   :: coord_file
    character(len=STRING_LENGTH)   :: new_coord_file
    character(len=STRING_LENGTH)   :: cp2k_program
    character(len=STRING_LENGTH)   :: fileroot_str
    character(len=STRING_LENGTH)   :: basis_set_file, potential_file, dft_file, global_file, cell_file
  logical                               :: clean_up_files
!  logical                     :: Delete_Metal_Connections
    type(Table) :: intrares_impropers
logical :: have_silica_potential
logical :: QM_list_changed
logical :: dummy
integer, pointer :: old_cluster_mark_p(:), cluster_mark_p(:)
integer :: i
  real(dp),dimension(3)          :: old_QM_cell
!    call system_initialise(verbosity=PRINT_SILENT,enable_timing=.true.)
    call system_initialise(verbosity=PRINT_NORMAL,enable_timing=.true.)
    call system_timer('program')

    !read parameters
    call initialise(params_in)
    basis_set_file=''
    cell_file=''
    dft_file=''
    global_file=''
    call param_register(params_in, 'Run_Type', 'MM', Run_Type, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params_in, 'PSF_Print', 'NO_PSF', Print_PSF, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params_in, 'coord_file', 'filepot.0.xyz',coord_file, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params_in, 'new_coord_file', 'filepot.0.out',new_coord_file, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params_in, 'cp2k_program', 'cp2k_serial',cp2k_program, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params_in, 'root', '', fileroot_str, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params_in, 'basis_set_file', '', basis_set_file, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params_in, 'potential_file', '', potential_file, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params_in, 'dft_file', '', dft_file, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params_in, 'global_file', '', global_file, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params_in, 'cell_file', '', cell_file, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params_in, 'clean_up_files', 'T', clean_up_files, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params_in, 'have_silica_potential', 'F', have_silica_potential, help_string="No help yet.  This source file was $LastChangedBy$")
!    call param_register(params_in, 'Delete_Metal_Connections', 'T', Delete_Metal_Connections, help_string="No help yet.  This source file was $LastChangedBy$")

    if (.not.param_read_args(params_in,ignore_unknown=.true.,task='cp2k_filepot args_str')) then
      call system_abort('could not parse argument line')
    end if

    call finalise(params_in)

    !print parameters
    call print('Run parameters:')
    call print('  Run_Type '//trim(Run_Type))
    call print('  PSF_Print '//trim(Print_PSF))
    call print('  coord_file '//trim(coord_file))
    call print('  new_coord_file '//trim(new_coord_file))
    call print('  cp2k_program '//trim(cp2k_program))
    call print('  root '//trim(fileroot_str))
    call print('  basis_set_file '//trim(basis_set_file))
    call print('  potential_file '//trim(potential_file))
    call print('  dft_file '//trim(dft_file))
    call print('  global_file '//trim(global_file))
    call print('  cell_file '//trim(cell_file))
    call print('  have_silica_potential '//have_silica_potential)
    call print('---------------------------------------')
    call print('')

! starts here

    call read_xyz(my_atoms,coord_file)

    !create connectivity & topology
    if (have_silica_potential) then
       call set_cutoff(my_atoms,SILICA_2body_CUTOFF)
    else
       call set_cutoff(my_atoms,0._dp)
    endif
    call calc_connect(my_atoms)
!    if (Delete_Metal_Connections) call delete_metal_connects(my_atoms)
    call map_into_cell(my_atoms)
    call calc_dists(my_atoms)
    if (trim(Run_Type).ne.'QS') then
       call create_CHARMM(my_atoms,do_CHARMM=.true.,intrares_impropers=intrares_impropers)
    endif

if (trim(Run_Type).eq.'QMMM_CORE') call system_abort('not yet tested.')
if ((trim(Run_Type).eq.'QMMM_EXTENDED').or.&
    (trim(Run_Type).eq.'QMMM_CORE')) then
   if (.not.has_property(my_atoms, 'old_cluster_mark')) call system_abort('no old_cluster_mark')
   if (.not.has_property(my_atoms, 'cluster_mark')) call system_abort('no cluster_mark')
   
   dummy = assign_pointer(my_atoms, 'old_cluster_mark', old_cluster_mark_p)
   dummy = assign_pointer(my_atoms, 'cluster_mark', cluster_mark_p)

QM_list_changed = .false.
   do i=1,my_atoms%N
      if (old_cluster_mark_p(i).ne.cluster_mark_p(i)) then
         if (any((/old_cluster_mark_p(i),cluster_mark_p(i)/).eq.HYBRID_NO_MARK) .or. &
             any((/old_cluster_mark_p(i),cluster_mark_p(i)/).eq.HYBRID_TERM_MARK) ) then
            QM_list_changed = .true.
         endif
      endif
   enddo
   call set_value(my_atoms%params,'QM_list_changed',QM_list_changed)
call print('set_value QM_list_changed '//QM_list_changed)

endif
    !generate args_str for go_cp2k
    write (args_str,'(a)') 'Run_Type='//trim(Run_Type)//' PSF_Print='//trim(Print_PSF)
    if (trim(cp2k_program).ne.'') args_str = trim(args_str)//' cp2k_program='//trim(cp2k_program)
    if (trim(fileroot_str).ne.'') args_str = trim(args_str)//' root='//trim(fileroot_str)
    if (trim(basis_set_file).ne.'') args_str = trim(args_str)//' basis_set_file='//trim(basis_set_file)
    if (trim(potential_file).ne.'') args_str = trim(args_str)//' potential_file='//trim(potential_file)
    if (trim(dft_file).ne.'') args_str = trim(args_str)//' dft_file='//trim(dft_file)
    if (trim(global_file).ne.'') args_str = trim(args_str)//' global_file='//trim(global_file)
    if (trim(cell_file).ne.'') args_str = trim(args_str)//' cell_file='//trim(cell_file)
    if (clean_up_files) then
       args_str = trim(args_str)//' clean_up_files=T'
    else
       args_str = trim(args_str)//' clean_up_files=F'
    endif
    if (have_silica_potential) then
       args_str = trim(args_str)//' have_silica_potential=T'
    else
       args_str = trim(args_str)//' have_silica_potential=F'
    endif

    call print('FILEPOT |'//args_str,verbosity=PRINT_SILENT)

    !call CP2K
    allocate(f0(3,my_atoms%N))
    call go_cp2k(my_atoms=my_atoms, forces=f0, energy=energy, args_str=args_str,intrares_impropers=intrares_impropers)
    !momentum conservation
    f0 = sum0(f0, my_atoms)

if ((trim(Run_Type).eq.'QMMM_EXTENDED')) then
    if (.not.get_value(my_atoms%params,'QM_cell',old_QM_cell)) call system_abort('Something has gone wrong. No QM_cell property found.')
endif

    !write energy and forces
    call set_value(my_atoms%params,'energy',energy)
    call add_property(my_atoms,'force', 0.0_dp, n_cols=3)
    if (.not. assign_pointer(my_atoms, 'force', forces_p)) then
       call system_abort("filepot_read_output needed forces, but couldn't find force in my_atoms ")
    endif
    forces_p = f0

    call initialise(xyz,new_coord_file,action=OUTPUT)
    call print_xyz(my_atoms,xyz,all_properties=.true.,real_format='f20.13')
    call finalise(xyz)

    deallocate(f0)
    call finalise(my_atoms)
    call verbosity_push(PRINT_NORMAL)
    call system_timer('FILEPOT |')
    call verbosity_push(PRINT_SILENT)
    call system_finalise

contains

! momentum conservation over all atoms, mass weighted
  function sum0(force,at) result(force0)

    real(dp), dimension(:,:), intent(in) :: force
    type(Atoms),              intent(in) :: at
    real(dp) :: force0(size(force,1),size(force,2))
    integer :: i
    real(dp) :: sumF(3), sum_weight

    sumF = sum(force,2)
    call print('Sum of the forces is '//sumF(1:3))

    if ((sumF(1) .feq. 0.0_dp) .and.  (sumF(2) .feq. 0.0_dp) .and.  (sumF(3) .feq. 0.0_dp)) then
       call print('Sum of the forces is zero.')
       force0 = force
    else
      !F_corr weighted by element mass
      sum_weight = sum(ElementMass(at%Z(1:at%N)))

      force0(1,:) = force(1,:) - sumF(1) * ElementMass(at%Z(:)) / sum_weight
      force0(2,:) = force(2,:) - sumF(2) * ElementMass(at%Z(:)) / sum_weight
      force0(3,:) = force(3,:) - sumF(3) * ElementMass(at%Z(:)) / sum_weight
    endif

    sumF = sum(force0,2)
    call print('Sum of the forces after mom.cons.: '//sumF(1:3))

  end function sum0

end program cp2k_filepot
