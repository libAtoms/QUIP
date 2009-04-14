! A program that takes an XYZ file and outputs a PSF and a PDB file.
! Also needs a residue library.
! written by Csilla Varnai
! last modified -- 2009-02-19

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
  use system_module,           only: dp, system_initialise, system_finalise, &
                                     system_timer, system_abort, &
                                     print, verbosity_push, &
                                     operator(//), &
                                     SILENT, NORMAL, ANAL, NERD
  use table_module,            only: table, finalise, int_part, delete
  use topology_module,         only: create_CHARMM, &
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
    character(80)               :: Root, &
                                   pdb_name, &
                                   psf_name
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

   ! print run parameters
    call print('Input:')
    call print('  XYZ file: '//trim(xyz_file))
    call print('  Residue Library: '//trim(Library))
    call print('  Delete_Metal_Connections'//Delete_Metal_Connections)
    call print('  Neighbour Tolerance: '//Neighbour_Tolerance)
    call print('Output:')
    call print('  PDB file: '//trim(pdb_name))
    call print('  PSF file: '//trim(psf_name))
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

   ! removes a bond between i and j, if the bond is present
  subroutine delete_bond(my_atoms, i, j)

    type(Atoms), intent(inout) :: my_atoms
    integer, intent(in) :: i, j

    integer :: ii, jj, kk, k, ll, change
    integer, allocatable :: bond_table(:,:)

    if (i.eq.j) return

    if (i.lt.j) then
       ii = i
       jj = j
    else
       ii = j
       jj = i
    endif

    allocate (bond_table(4,my_atoms%connect%neighbour1(ii)%N))
    bond_table = int_part(my_atoms%connect%neighbour1(ii))
    kk = find_in_array(bond_table(1,:),jj)
    if (kk.gt.0) then
!       call print('found bond to delete for atoms '//ii//' '//jj)
       call delete(my_atoms%connect%neighbour1(ii),kk)
    endif
    deallocate(bond_table)

   ! if I delete a bond from neighbour1, I should update in neighbour2 that
   ! it's the $(n-1)$-th (or the last is now takes the place of the deleted) neighbour from now on
    if (kk.gt.0) then
       do k = kk, my_atoms%connect%neighbour1(ii)%N     ! k-th neighbour of ii
          ll = my_atoms%connect%neighbour1(ii)%int(1,k)     ! is ll
          allocate (bond_table(2,my_atoms%connect%neighbour2(ll)%N))
          bond_table = int_part(my_atoms%connect%neighbour2(ll))
          change = find_in_array(bond_table(1,:),ii)    ! ll has an account for ii
          if (change.eq.0) call system_abort('Found nonsymmetrical connectivity.')
          my_atoms%connect%neighbour2(ll)%int(2,change) = k ! this account is updated now
          deallocate (bond_table)
       enddo
    endif

    allocate (bond_table(2,my_atoms%connect%neighbour2(jj)%N))
    bond_table = int_part(my_atoms%connect%neighbour2(jj))
    kk = find_in_array(bond_table(1,:),ii)
    if (kk.gt.0) then
!       call print('found bond to delete for atoms '//jj//' '//ii)
       call delete(my_atoms%connect%neighbour2(jj),kk)
    endif
    deallocate(bond_table)

  end subroutine delete_bond

  subroutine print_usage

    call print('Usage: xyz2pdb File=filename.xyz [Residue_Library=library] [Delete_Metal_Connections=T] [Neighbour_tolerance=1.2]')
    call print('')
    call print('  File=filename,        where your input file has extension .xyz')
    call print('  [Residue_Library=library],  optional, default is protein_res.CHARMM.lib')
    call print('  [Delete_Metal_Connections=T],  optional, default is true, only calculates connection for H,C,N,O,Si,P,S,Cl')
    call print('                           should work fine - only modify if you want bonds with other elements in your structure')
    call print('  [Neighbour_Tolerance=1.2],  optional, default is 1.2, should work fine - do not poke it ')
    call print('                          unless you know you have unusally long bonds in your structure')

    call print('')

  end subroutine print_usage

end program xyz2pdb
