! last modified -- 2009-02-19 -- Csilla
! takes an xyz file, identifies the residues and outputs a PSF and a PDB file.

program xyz2pdb

!  use libatoms_module

  use atoms_module,            only: atoms, finalise, &
                                     read_xyz, &
                                     set_cutoff, &
                                     calc_connect, print, &
                                     distance_min_image, bond_length
  use dictionary_module,       only: dictionary, initialise, finalise, &
                                     set_value
  use linearalgebra_module,    only: find_in_array
  use paramreader_module,      only: param_register, param_read_args, &
                                     FIELD_LENGTH, PARAM_MANDATORY
  use periodictable_module,    only: ElementCovRad
  use system_module,           only: dp, inoutput, &
                                     system_initialise, system_finalise, &
                                     system_timer, system_abort, &
                                     print, verbosity_push, &
                                     operator(//), &
                                     INPUT, OUTPUT, &
                                     SILENT, NORMAL, ANAL, NERD
  use table_module,            only: table, finalise, int_part, delete
  use topology_module,         only: create_CHARMM, delete_bond, &
                                     write_brookhaven_pdb_file, &
                                     write_psf_file, &
                                     MM_RUN

  implicit none

    type(Atoms)                 :: my_atoms
    type(Table)                 :: intrares_impropers

    type(Dictionary)            :: params_in
    character(len=FIELD_LENGTH) :: Library, &
                                   xyz_file
    real(dp)                    :: Neighbour_Tolerance
    logical                     :: Delete_Metal_Connections
    logical                     :: print_xsc
    character(80)               :: Root, &
                                   pdb_name, &
                                   psf_name, &
                                   xsc_name
    integer                     :: len_name
    logical                     :: ex

    call system_initialise(verbosity=silent,enable_timing=.true.)
    call verbosity_push(NORMAL)
    call system_timer('program')

   ! reading in run parameters
    call initialise(params_in)
    call param_register(params_in, 'File', PARAM_MANDATORY, xyz_file)
    call param_register(params_in, 'Residue_Library', 'protein_res.CHARMM.lib',Library)
    call param_register(params_in, 'Neighbour_Tolerance', '1.2', Neighbour_Tolerance)
    call param_register(params_in, 'Delete_Metal_Connections', 'T', Delete_Metal_Connections)
    call param_register(params_in, 'Print_XSC', 'F', print_xsc)
    if (.not. param_read_args(params_in, do_check = .true.)) then
      call print_usage
      call system_abort('could not parse argument line')
    end if

   ! filenames
    len_name = len_trim(xyz_file)
    if (any(xyz_file(len_name-3:len_name).eq.(/'.xyz','.XYZ'/))) then
       Root = xyz_file(1:len_name-4)
    else
       Root = xyz_file(1:len_name)
    endif
    pdb_name = trim(Root)//'.pdb'
    psf_name = trim(Root)//'.psf'
    xsc_name = trim(Root)//'.xsc'

   ! print run parameters
    call print('Input:')
    call print('  XYZ file: '//trim(xyz_file))
    call print('  Residue Library: '//trim(Library))
    call print('  Print_XSC (NAMD cell file): '//print_xsc)
    call print('  Delete_Metal_Connections'//Delete_Metal_Connections)
    call print('  Neighbour Tolerance: '//Neighbour_Tolerance)
    call print('Output:')
    call print('  PDB file: '//trim(pdb_name))
    call print('  PSF file: '//trim(psf_name))
    if (print_xsc) call print('  XSC file: '//trim(xsc_name))
    call print('CHARMM atomic types and Brookhaven PDB format are used')
    call print('')

    call finalise(params_in)

   ! check if input files exist
    inquire(file=trim(xyz_file),exist=ex)
    if (.not.ex) then
       call print_usage
       call system_abort('Missing xyz File '//xyz_file)
    endif
    inquire(file=trim(Library),exist=ex)
    if (.not.ex) then
       call print_usage
       call system_abort('Missing Residue_Library '//Library)
    endif

   ! read in XYZ
    call print('Reading in coordinates...')
    call read_xyz(my_atoms,xyz_file)

   ! calculating connectivities
    call print('Calculating connectivities...')
   ! use nonuniform connect_cutoff, include only nearest neighbours, otherwise the adjacency matrix will contain disjoint regions
    call set_cutoff(my_atoms,0.0_dp)
    my_atoms%nneightol = Neighbour_Tolerance
    call calc_connect(my_atoms)
   ! remove bonds for metal ions - everything but H, C, N, O, Si, P, S
    if (Delete_Metal_Connections) call delete_metal_connects(my_atoms)
!    call print(my_atoms%connect)


   ! identify residues
    call print('Identifying residues...')
    call set_value(my_atoms%params,'Library',trim(Library))
    call create_CHARMM(my_atoms,do_CHARMM=.true.,intrares_impropers=intrares_impropers)

   ! print output PDB and PSF files
    call print('Writing files with CHARMM format...')
    call write_psf_file(my_atoms,psf_file=trim(psf_name),run_type=MM_RUN,intrares_impropers=intrares_impropers)
    call write_brookhaven_pdb_file(my_atoms,trim(pdb_name),run_type=MM_RUN)
    if (Print_XSC) call write_xsc_file(my_atoms,xsc_file=trim(xsc_name))

    call finalise(intrares_impropers)
    call finalise(my_atoms)
    call print ('Finished.')

    call system_timer('program')
    call system_finalise

contains

   ! remove bonds for metal ions - everything but H, C, N, O, Si, P, S, Cl
  subroutine delete_metal_connects(my_atoms)

    type(Atoms), intent(inout) :: my_atoms
    integer :: i, j

    do i = 1, my_atoms%N
       if (any(my_atoms%Z(i).eq.(/1,6,7,8,14,15,16,17/))) cycle
!       call print('found metal atom '//i)
       do j = 1, my_atoms%N
          if (i.ne.j) call delete_bond(my_atoms, i, j)
       enddo
    enddo

  end subroutine delete_metal_connects

  subroutine print_usage

    call print('Usage: xyz2pdb File=filename.xyz [Residue_Library=library] [Print_XSC] [Delete_Metal_Connections=T] [Neighbour_tolerance=1.2]')
    call print('')
    call print('  File=filename,        where your input file has extension .xyz')
    call print('  [Residue_Library=library],  optional, default is protein_res.CHARMM.lib')
    call print('  [Print_XSC],          optional, whether to print NAMD cell file, default is false')
    call print('  [Delete_Metal_Connections=T],  optional, default is true, only calculates connection for H,C,N,O,Si,P,S,Cl')
    call print('                           should work fine - only modify if you want bonds with other elements in your structure')
    call print('  [Neighbour_Tolerance=1.2],  optional, default is 1.2, should work fine - do not poke it ')
    call print('                          unless you know you have unusally long bonds in your structure')

    call print('')

  end subroutine

  subroutine write_xsc_file(at, xsc_file)

    character(len=*),  intent(in) :: xsc_file
    type(Atoms),       intent(in) :: at

    type(InOutput) :: xsc

    call initialise(xsc,trim(xsc_file),action=OUTPUT)
    call print('   XSC file: '//trim(xsc%filename))

    call print('#NAMD cell file',file=xsc)
    call print('#$LABELS step a_x a_y a_z b_x b_y b_z c_x c_y c_z o_x o_y o_z',file=xsc)
    call print('1 '//reshape(at%lattice,(/9/))//' '//(/0,0,0/),file=xsc)
!write_string(this%params,real_format='f18.6')
    call finalise(xsc)

  end subroutine write_xsc_file

end program xyz2pdb
