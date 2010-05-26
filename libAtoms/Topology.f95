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

module topology_module

  use atoms_module,            only: atoms, connection, print, finalise, &
                                     map_into_cell, calc_dists, &
                                     assignment(=), &
                                     add_property, has_property, &
                                     read_line, parse_line, &
                                     atoms_n_neighbours, atoms_neighbour, bond_length, &
                                     distance_min_image, &
                                     DEFAULT_NNEIGHTOL, set_cutoff, remove_bond, calc_connect, &
				     assign_pointer, is_nearest_neighbour_abs_index
  use clusters_module,         only: bfs_step, add_cut_hydrogens
  use dictionary_module,       only: get_value, value_len
  use linearalgebra_module,    only: find_in_array, find, &
                                     print
  use periodictable_module,    only: ElementName, ElementMass, ElementCovRad
  use structures_module,       only: find_motif
  use system_module,           only: dp, inoutput, initialise, finalise, &
                                     INPUT, OUTPUT, INOUT, &
                                     system_timer, &
                                     optional_default, &
                                     print, print_title, &
                                     string_to_int, string_to_real, round, &
                                     parse_string, read_line, &
                                     ERROR, SILENT, NORMAL, VERBOSE, NERD, ANAL, &
                                     operator(//), allocatable_array_pointers, &
				     verbosity_push, verbosity_pop
#ifndef HAVE_QUIPPY
  use system_module,           only: system_abort
#endif
  use table_module,            only: table, initialise, finalise, &
                                     append, allocate, delete, &
                                     int_part, TABLE_STRING_LENGTH
  use units_module,            only: MASSCONVERT, PI


  implicit none

  private :: next_motif, write_psf_section, create_bond_list, write_pdb_file
  private :: create_angle_list, create_dihedral_list
  private :: create_improper_list, get_property
  private :: create_pos_dep_charges, calc_fc

  public  :: delete_metal_connects, &
             write_brookhaven_pdb_file, &
             write_cp2k_pdb_file, &
             write_psf_file, &
             write_psf_file_arb_pos, &
             create_residue_labels, &
             create_residue_labels_arb_pos, &
             NONE_RUN, &
             QS_RUN, &
             MM_RUN, &
             QMMM_RUN_CORE, &
             QMMM_RUN_EXTENDED


!parameters for Run_Type
  integer, parameter :: NONE_RUN = -100
  integer, parameter :: QS_RUN = -1
  integer, parameter :: MM_RUN = 0
  integer, parameter :: QMMM_RUN_CORE = 1
  integer, parameter :: QMMM_RUN_EXTENDED = 2

!parameters for the Residue Library
  integer,  parameter, private :: MAX_KNOWN_RESIDUES = 200 !Maximum number of residues in the library
                                                           !(really about 116 at the mo)
  integer,  parameter, private :: MAX_ATOMS_PER_RES = 50   !Maximum number of atoms in one residue
  integer,  parameter, private :: MAX_IMPROPERS_PER_RES = 12   !Maximum number of impropers in one residue

  real(dp), parameter, public :: SILICON_2BODY_CUTOFF = 3.8_dp  !Si-Si 2body interaction
  real(dp), parameter, public :: SILICA_2BODY_CUTOFF = 7.5_dp !Si-O, O-O 2body interaction
  real(dp), parameter, public :: SILICON_3BODY_CUTOFF = 3.8_dp !Si-Si-Si, Si-Si-O 3body cutoff
  real(dp), parameter, public :: SILICA_3BODY_CUTOFF = 3.6_dp !Si-O-Si, O-Si-O, Si-O-H 3body cutoff
!  real(dp), parameter, public :: SILICON_2BODY_CUTOFF = 2.8_dp  !Si-Si 2body interaction
!  real(dp), parameter, public :: SILICA_2BODY_CUTOFF = 5.5_dp !Si-O, O-O 2body interaction
!  real(dp), parameter, public :: SILICON_3BODY_CUTOFF = 2.8_dp !Si-Si-Si, Si-Si-O 3body cutoff
!  real(dp), parameter, public :: SILICA_3BODY_CUTOFF = 2.6_dp !Si-O-Si, O-Si-O, Si-O-H 3body cutoff

  real(dp), parameter, private :: Danny_R_q = 2.0_dp
  real(dp), parameter, private :: Danny_Delta_q = 0.1_dp

  real(dp), parameter, private :: Danny_Si_Si_cutoff = 2.8_dp
  real(dp), parameter, private :: Danny_Si_Si_O_cutoff = 2.8_dp
  real(dp), parameter, private :: Danny_Si_O_cutoff = 2.6_dp
  real(dp), parameter, private :: Danny_Si_H_cutoff = 0._dp
!  real(dp), parameter, private :: Danny_O_O_cutoff = 2.6_dp
  real(dp), parameter, private :: Danny_O_O_cutoff = 0._dp
  real(dp), parameter, private :: Danny_O_H_cutoff = 1.0_dp
  real(dp), parameter, private :: Danny_H_H_cutoff = 0._dp

contains


  !% Topology calculation using arbitrary (usually avgpos) coordinates, as a wrapper to find_residue_labels
  !%
  subroutine create_residue_labels_arb_pos(at,do_CHARMM,intrares_impropers,have_silica_potential,pos_field_for_connectivity)
    type(Atoms),           intent(inout),target :: at
    logical,     optional, intent(in)    :: do_CHARMM
    type(Table), optional, intent(out)   :: intrares_impropers
    logical,     optional, intent(in)    :: have_silica_potential
    character(len=*), optional, intent(in) :: pos_field_for_connectivity

    real(dp), pointer :: use_pos(:,:)
    type(Connection) :: t_connect
    type(Atoms) :: at_copy
    logical :: do_have_silica_potential
    logical :: use_pos_is_pos

    ! save a copy
    at_copy = at

    use_pos_is_pos = .false.
    ! find desired position field (pos, avgpos, whatever)
    if (present(pos_field_for_connectivity)) then
      if (.not. assign_pointer(at, trim(pos_field_for_connectivity), use_pos)) &
	call system_abort("calc_topology can't find pos field '"//trim(pos_field_for_connectivity)//"'")
      if (trim(pos_field_for_connectivity) == 'pos') use_pos_is_pos = .true.
    else
      if (.not. assign_pointer(at, 'avgpos', use_pos)) &
	call system_abort("calc_topology can't find default pos field avgpos")
    endif

    do_have_silica_potential = optional_default(.false.,have_silica_potential)

    ! copy desired pos to pos, and new connectivity
    !NB don't do if use_pos => pos
    if (.not. use_pos_is_pos) at_copy%pos = use_pos
    if (do_have_silica_potential) then
      call set_cutoff(at_copy,SILICA_2body_CUTOFF)
    else
      call set_cutoff(at_copy,0.0_dp)
    endif
    call calc_connect(at_copy, alt_connect=t_connect)

    ! now create labels using this connectivity object
    ! if cutoff is set to 0, nneighb_only doesn't matter
    ! if cutoff is large, then nneighb_only must be true
    !    so might as well pass true
    call create_residue_labels(at,do_CHARMM,intrares_impropers,nneighb_only=.true.,alt_connect=t_connect,have_silica_potential=do_have_silica_potential)
    call finalise(t_connect)

  end subroutine create_residue_labels_arb_pos

  !% Topology calculation. Recognize residues and save atomic names, residue names etc. in the atoms object.
  !% Generally useable for AMBER and CHARMM force fields, in case of CHARMM, the MM charges are assigned here, too.
  !% do_CHARMM=.true. is the default
  !% Optionally outputs the intraresidual impropers.
  !% The algorithm will fail if the Residue_Library is not present in atoms%params or in the directory.
  !% Optionally use hysteretic_connect if passed as alt_connect.
  !% Optionally could use hysteretic neighbours instead of nearest neighbours, if the cutoff of the
  !% alt_connect were the same as at%cutoff(_break).
  !%
  subroutine create_residue_labels(at,do_CHARMM,intrares_impropers, nneighb_only,alt_connect,have_silica_potential) !, hysteretic_neighbours)

    type(Atoms),           intent(inout),target :: at
    logical,     optional, intent(in)    :: do_CHARMM
    type(Table), optional, intent(out)   :: intrares_impropers
    logical, intent(in), optional :: nneighb_only
    type(Connection), intent(in), optional, target :: alt_connect
    logical,     optional, intent(in)    :: have_silica_potential
!    logical, optional, intent(in) :: hysteretic_neighbours

    character(*), parameter  :: me = 'create_residue_labels_pos: '
    logical :: remove_Si_H_silica_bonds = .true.

    type(Inoutput)                       :: lib
    character(4)                         :: cha_res_name(MAX_KNOWN_RESIDUES), Cres_name
    character(3)                         :: pdb_res_name(MAX_KNOWN_RESIDUES), pres_name
    integer                              :: residue_number(at%N)
    character(4)                         :: atom_name(at%N), atom_name_PDB(at%N)
    real(dp)                             :: atom_charge(at%N)
    type(Table)                          :: residue_type, list
    logical                              :: unidentified(at%N)
    integer, allocatable, dimension(:,:) :: motif
    integer                              :: i, m, n, nres
    character(4), allocatable, dimension(:) :: at_names, at_names_PDB
    real(dp),     allocatable, dimension(:) :: at_charges
    logical                              :: my_do_charmm
    integer                              :: atom_type_index, &
                                            atom_type_PDB_index, &
                                            atom_res_name_index, &
                                            atom_mol_name_index, &
                                            atom_res_number_index, &
                                            atom_charge_index
    logical                              :: ex
    character(len=value_len)             :: residue_library
    integer                              :: i_impr, n_impr
    integer, allocatable                 :: imp_atoms(:,:)
    real(dp)                             :: mol_charge_sum
    logical                              :: found_residues
    type(Table)                          :: atom_Si, atom_SiO, SiOH_list
    real(dp), dimension(:), allocatable  :: charge
integer :: j,atom_i, ji
!logical :: silanol
    type(Table) :: bondH,bondSi
    integer :: bond_H,bond_Si
!    integer                             :: qm_flag_index, pos_indices(3)
!    logical                             :: do_qmmm
    type(Connection), pointer :: use_connect
!    logical :: use_hysteretic_neighbours
type(Table) :: O_atom, O_neighb
integer :: hydrogen
logical :: silica_potential

!    integer, pointer :: mol_id(:)
    type(allocatable_array_pointers), allocatable :: molecules(:)

    call system_timer('create_residue_labels_pos')

    silica_potential = optional_default(.false.,have_silica_potential)

    if (present(alt_connect)) then
      use_connect => alt_connect
    else
      use_connect => at%connect
    endif
    if (.not.use_connect%initialised) call system_abort(me//'No connectivity data present in atoms structure')    
!    use_hysteretic_neighbours = optional_default(.false.,hysteretic_neighbours)

    my_do_charmm = optional_default(.true.,do_CHARMM)

    residue_library = ''
    call print_title('Creating CHARMM format')
    ex = .false.
    ex = get_value(at%params,'Library',residue_library)
    if (ex) call print('Library: '//trim(residue_library))
    if (.not.ex) call system_abort('create_residue_labels_pos: no residue library specified, but topology generation requested')

    !Open the residue library
    !call print('Opening library...')
    call initialise(lib,trim(residue_library),action=INPUT)

    !Set all atoms as initially unidentified
    unidentified = .true.

    !Read each of the residue motifs from the library
    n = 0
    nres = 0
    mol_charge_sum = 0._dp
    found_residues = .false.
    if (present(intrares_impropers)) call initialise(intrares_impropers,4,0,0,0,0)
    call allocate(residue_type,1,0,0,0,1000)
    call print('Identifying atoms...')

!!!!!!!!!!!!!!! DANNY POTENTIAL !!!!!!!!!!!!!!!!
    if (silica_potential) then

   ! SIO residue for Danny potential if Si atom is present in the atoms structure
       if (any(at%Z(1:at%N).eq.14)) then
          call print('|-Looking for SIO residue, not from the library...')
          call print('| |-Found... will be treated as 1 molecule, 1 residue...')
          !all this bulk will be 1 residue
          n = n + 1
          cha_res_name(n) = 'SIO2'
          pdb_res_name(n) = '' !not used
          call append(residue_type,(/n/))
          nres = nres + 1

          !Add Si atoms
          call initialise(atom_Si,4,0,0,0,0)
          do i = 1,at%N
             if (at%Z(i).eq.14) then
                call append(atom_Si,(/i,0,0,0/))
             endif
          enddo
          call print(atom_Si%N//' Si atoms found in total')
          !Add O atoms
          call bfs_step(at,atom_Si,atom_SiO,nneighb_only=.true.,min_images_only=.true.,alt_connect=use_connect)
          call print(atom_SiO%N//' O atoms found in total')
!          if (any(at%Z(atom_SiO%int(1,1:atom_SiO%N)).eq.1)) call system_abort('Si-H bond')
          if (remove_Si_H_silica_bonds) then
             do i=1,atom_SiO%N !check Hs bonded to Si. There shouldn't be any,removing the bond.
                 if (at%Z(atom_SiO%int(1,i)).eq.1) then
                    call print('WARNING! Si and H are very close',verbosity=ERROR)
                    bond_H = atom_SiO%int(1,i)
                    call initialise(bondH,4,0,0,0,0)
                    call append(bondH,(/bond_H,0,0,0/))
                    call bfs_step(at,bondH,bondSi,nneighb_only=.true.,min_images_only=.true.,alt_connect=use_connect)
                    do j = 1,bondSi%N
                       if (at%Z(bondSi%int(1,i)).eq.14) then
                          bond_Si = bondSi%int(1,i)
                          call print('WARNING! Remove Si '//bond_Si//' and H '//bond_H//' bond ('//distance_min_image(at,bond_H,bond_Si)//')',verbosity=ERROR)
                          call remove_bond(use_connect,bond_H,bond_Si)
                       endif
                    enddo
                 endif
!                call print('atom_SiO has '//at%Z(atom_SiO%int(1,i)))
             enddo
          endif
          !Add H atoms -- if .not.remove_Si_H_bonds, we might include whole water molecules at this stage, adding the remaining -OH.
      !    call bfs_step(at,atom_SiO,atom_SIOH,nneighb_only=.true.,min_images_only=.true.)
          call add_cut_hydrogens(at,atom_SiO,alt_connect=use_connect)
          call print(atom_SiO%N//' O/H atoms found in total')

          !check if none of these atoms are identified yet
          if (any(.not.unidentified(atom_Si%int(1,1:atom_Si%N)))) then
             call system_abort('already identified atoms found again.')
          endif
          if (any(.not.unidentified(atom_SiO%int(1,1:atom_SiO%N)))) then
!             call system_abort('already identified atoms found again.')
             do i = 1,atom_SiO%N,-1
                if (.not.unidentified(atom_SiO%int(1,i))) then
                   call print('delete from SiO2 list already identified atom '//atom_SiO%int(1,1:atom_SiO%N))
                   call delete(atom_SiO,i)
                endif
             enddo
          endif
          unidentified(atom_Si%int(1,1:atom_Si%N)) = .false.
          unidentified(atom_SiO%int(1,1:atom_SiO%N)) = .false.

          !add atom, residue and molecule names
          !SIO
          do i = 1, atom_Si%N                              !atom_Si  only has Si atoms
             atom_i = atom_Si%int(1,i)
             atom_name(atom_i) = 'SIO'
             atom_name_PDB(atom_i) = 'SIO'
          enddo
          !OSB, OSI & HSI
          do i = 1, atom_SiO%N                             !atom_SiO only has O,H atoms
             atom_i = atom_SiO%int(1,i)
             !OSB & OSI
             if (at%Z(atom_i).eq.8) then
                call initialise(O_neighb,4,0,0,0)
                call initialise(O_atom,4,0,0,0)
                call append(O_atom,(/atom_i,0,0,0/))
                call bfs_step(at,O_atom,O_neighb,nneighb_only=.true.,min_images_only=.true.,alt_connect=use_connect)
                !check nearest neighbour number = 2
                if (O_neighb%N.ne.2) then
                   call print('WARNING! silica O '//atom_i//'has '//O_neighb%N//'/=2 nearest neighbours',ERROR)
                   call print('neighbours: '//O_neighb%int(1,1:O_neighb%N))
                endif
                !check if it has a H nearest neighbour
                hydrogen = find_in_array(at%Z(O_neighb%int(1,1:O_neighb%N)),1)
                if (hydrogen.ne.0) then
                   atom_name(atom_i) = 'OSI' !silanol O
                   atom_name_PDB(atom_i) = 'OSI' !silanol O
!                   call print('Found OH silanol oxygen.'//atom_SiO%int(1,i)//' hydrogen: '//O_neighb%int(1,hydrogen))
                   !check if it has only 1 H nearest neighbour
                   if (hydrogen.lt.O_neighb%N) then
                      if(find_in_array(at%Z(O_neighb%int(1,hydrogen+1:O_neighb%N)),1).gt.0) &
                         call system_abort('More than 1 H neighbours of O '//atom_i)
                   endif
                else
                   atom_name(atom_i) = 'OSB' !bridging O
                   atom_name_PDB(atom_i) = 'OSB' !bridging O
!                   call print('Found OB bridging oxygen.'//atom_SiO%int(1,i))
                endif
                call finalise(O_atom)
                call finalise(O_neighb)
             !HSI
             elseif (at%Z(atom_SiO%int(1,i)).eq.1) then
                atom_name(atom_SiO%int(1,i)) = 'HSI'
                atom_name_PDB(atom_SiO%int(1,i)) = 'HSI'
             else
                call system_abort('Non O/H atom '//atom_i//'!?')
             endif
          enddo

          !Add all the silica atoms together
          call initialise(SiOH_list,4,0,0,0,0)
          call append (SiOH_list,atom_Si)
          call append (SiOH_list,atom_SiO)

          !Residue numbers
          residue_number(SiOH_list%int(1,1:SiOH_list%N)) = nres

          !Charges
          call create_pos_dep_charges(at,SiOH_list,charge) !,residue_names=cha_res_name(residue_type%int(1,residue_number(1:at%N))))
          atom_charge(SiOH_list%int(1,1:SiOH_list%N)) = 0._dp
          atom_charge(SiOH_list%int(1,1:SiOH_list%N)) = charge(SiOH_list%int(1,1:SiOH_list%N))
call print("overall silica charge: "//sum(atom_charge(SiOH_list%int(1,1:SiOH_list%N))))
          call print('Atomic charges: ',ANAL)
          call print('   ATOM     CHARGE',ANAL)
          do i=1,at%N
             call print('   '//i//'   '//atom_charge(i),verbosity=ANAL)
          enddo

          call finalise(atom_Si)
          call finalise(atom_SiO)
          call finalise(SiOH_list)

       else
          call print('WARNING! have_silica_potential is true, but found no silicon atoms in the atoms object!',ERROR)
       endif

    endif
!!!!!!!!!!!!!!! END DANNY POTENTIAL !!!!!!!!!!!!!!!!

    do 

       ! Pull the next residue template from the library
       if (my_do_charmm) then
          call next_motif(lib,cres_name,pres_name,motif,atom_names=at_names,atom_charges=at_charges,atom_names_PDB=at_names_PDB, n_impr=n_impr,imp_atoms=imp_atoms,do_CHARMM=.true.)
       else
          call next_motif(lib,cres_name,pres_name,motif,atom_names=at_names,atom_names_PDB=at_names_PDB, do_CHARMM=.false.)
       endif

       if (cres_name=='NONE') then
          found_residues = .true.
          exit
       endif

       ! Store its CHARMM (3 char) and PDB (3 char) names
       ! e.g. cha/pdb_res_name(5) corresponds to the 5th residue found in the library
       n = n + 1
       cha_res_name(n) = cres_name
       pdb_res_name(n) = pres_name

       ! Search the atom structure for this residue
       call print('|-Looking for '//cres_name//'...',verbosity=ANAL)
       call find_motif(at,motif,list,mask=unidentified,nneighb_only=nneighb_only,alt_connect=use_connect) !,hysteretic_neighbours=use_hysteretic_neighbours)

       if (list%N > 0) then
          
          call print('| |-Found '//list%N//' occurrences of '//cres_name//' with charge '//(sum(at_charges(1:size(at_charges)))))
          mol_charge_sum = mol_charge_sum + list%N * sum(at_charges(1:size(at_charges)))

          ! Loop over all found instances of the residue

          do m = 1, list%N

             !Mark the newly identified atoms
!1 row is 1 residue
             unidentified(list%int(:,m)) = .false.

             !Store the residue info
             nres = nres + 1
             call append(residue_type,(/n/))
             ! if residue_type%int(1,2) == 3 then residue no. 2 matches the 3rd residue in the library
             residue_number(list%int(:,m)) = nres
                  !e.g. if residue_number(i) = j        j-th residue in the atoms object (1 to 8000 in case of 8000 H2O)
                  !        residue_type(1,j)   = k        k-th residue in the library file, in order
                  !        cha_res_name(k)   = 'ALA'    name of the k-th residue in the library file
                  !then atom 'i' is in a residue 'ALA'
             atom_name(list%int(:,m)) = at_names
             atom_name_PDB(list%int(:,m)) = at_names_PDB
             if (my_do_charmm) then
                atom_charge(list%int(:,m)) = at_charges
               ! intraresidual IMPROPERs
                if (present(intrares_impropers)) then
                   do i_impr = 1, n_impr
                      call append(intrares_impropers,list%int(imp_atoms(1:4,i_impr),m))
!                      call print('Added intraresidual improper '//intrares_impropers%int(1:4,intrares_impropers%N))
                   enddo
                endif
             endif
          end do

       end if

       call print('|',verbosity=ANAL)

    end do

   ! check if residue library is empty
    if (.not.found_residues) call system_abort('Residue library '//trim(lib%filename)//' does not contain any residues!')

    call print('Finished.')
    call print(nres//' residues found in total')

    if (any(unidentified)) then
       call print(count(unidentified)//' unidentified atoms',verbosity=ERROR)
       call print(find(unidentified))
       do i=1,at%N
	  if (unidentified(i)) then
	    call print(ElementName(at%Z(i))//' atom '//i//' has avgpos: '//round(at%pos(1,i),5)//&
	      ' '//round(at%pos(2,i),5)//' '//round(at%pos(3,i),5),verbosity=ERROR)
	    call print(ElementName(at%Z(i))//' atom '//i//' has number of neighbours: '//atoms_n_neighbours(at,i),verbosity=ERROR)
	    do ji=1, atoms_n_neighbours(at, i)
	      j = atoms_neighbour(at, i, ji)
	      call print("  neighbour " // j // " is of type " // ElementName(at%Z(j)), verbosity=ERROR)
	    end do
	  endif
       enddo

       ! THIS IS WHERE THE CALCULATION OF NEW PARAMETERS SHOULD GO
      call system_abort('create_residue_labels_pos: Unidentified atoms')

    else
       call print('All atoms identified')
       call print('Total charge of the molecule: '//round(mol_charge_sum,5))
    end if

   ! add data to store CHARMM topology
    call add_property(at,'atom_type',repeat(' ',TABLE_STRING_LENGTH))
    call add_property(at,'atom_type_PDB',repeat(' ',TABLE_STRING_LENGTH))
    call add_property(at,'atom_res_name',repeat(' ',TABLE_STRING_LENGTH))
    call add_property(at,'atom_mol_name',repeat(' ',TABLE_STRING_LENGTH))
    call add_property(at,'atom_res_number',0)
    call add_property(at,'atom_charge',0._dp)

    atom_type_index = get_property(at,'atom_type')
    atom_type_PDB_index = get_property(at,'atom_type_PDB')
    atom_res_name_index = get_property(at,'atom_res_name')
    atom_mol_name_index = get_property(at,'atom_mol_name')
    atom_res_number_index = get_property(at,'atom_res_number')
    atom_charge_index = get_property(at,'atom_charge')

    at%data%str(atom_type_index,1:at%N) = 'X'
    at%data%str(atom_type_PDB_index,1:at%N) = 'X'
    at%data%str(atom_res_name_index,1:at%N) = 'X'
    at%data%str(atom_mol_name_index,1:at%N) = 'X'
    at%data%int(atom_res_number_index,1:at%N) = 0
    at%data%real(atom_charge_index,1:at%N) = 0._dp

    at%data%str(atom_res_name_index,1:at%N) = cha_res_name(residue_type%int(1,residue_number(1:at%N)))
    at%data%int(atom_res_number_index,1:at%N) = residue_number(1:at%N)
    !NB workaround for pgf90 bug (as of 9.0-1)
    ! at%data%str(atom_type_index,1:at%N) = adjustl(atom_name(1:at%N))
    do i=1, at%N
      at%data%str(atom_type_index,i) = adjustl(atom_name(i))
      at%data%str(atom_type_PDB_index,i) = adjustl(atom_name_PDB(i))
    end do
    !NB end of workaround for pgf90 bug (as of 9.0-1)

    if (my_do_charmm) then
       allocate(molecules(at%N))
       call find_molecule_ids(at,molecules,nneighb_only=nneighb_only,alt_connect=alt_connect)
       do i=1, size(molecules)
	 if (allocated(molecules(i)%i_a)) then
           ! special case for silica molecule
           if (silica_potential .and. count(at%Z(molecules(i)%i_a) == 14) /=0) then
	       at%data%str(atom_mol_name_index,molecules(i)%i_a) = at%data%str(atom_res_name_index,molecules(i)%i_a)
	   ! special case for single atoms
	   elseif (size(molecules(i)%i_a) == 1) then
	     at%data%str(atom_mol_name_index,molecules(i)%i_a) = at%data%str(atom_res_name_index,molecules(i)%i_a)
	   ! special case for H2O
	   else if (size(molecules(i)%i_a) == 3) then
	     if (count(at%Z(molecules(i)%i_a) == 8) == 1 .and. &
	         count(at%Z(molecules(i)%i_a) == 1) == 2) then
	       at%data%str(atom_mol_name_index,molecules(i)%i_a) = at%data%str(atom_res_name_index,molecules(i)%i_a)
	     else ! default
call print("Found molecule containing "//size(molecules(i)%i_a)//" atoms and not water, single atom or silica")
	       at%data%str(atom_mol_name_index,molecules(i)%i_a) = "M"//i
	     endif
	   else ! default
call print("Found molecule containing "//size(molecules(i)%i_a)//" atoms and not water, single atom or silica")
	     at%data%str(atom_mol_name_index,molecules(i)%i_a) = "M"//i
	   endif
	 end if ! allocated(molecules)
       end do ! i=1,size(molecules)
       deallocate(molecules)
       at%data%real(atom_charge_index,1:at%N) = atom_charge(1:at%N)
    endif

    if (any(at%data%int(atom_res_number_index,1:at%N).le.0)) &
       call system_abort('create_residue_labels_pos: atom_res_number is not >0 for every atom')
    if (any(at%data%str(atom_type_index,1:at%N).eq.'X')) &
       call system_abort('create_residue_labels_pos: atom_type is not saved for at least one atom')
    if (any(at%data%str(atom_res_name_index,1:at%N).eq.'X')) &
       call system_abort('create_residue_labels_pos: atom_res_name is not saved for at least one atom')

    !Free up allocations
    call finalise(residue_type)
    call finalise(list)
    if (allocated(motif)) deallocate(motif)

    !Close the library
    call finalise(lib)

    call system_timer('create_residue_labels_pos')

  end subroutine create_residue_labels

  subroutine find_molecule_ids(at,molecules,nneighb_only,alt_connect)
    type(Atoms), intent(inout) :: at
    type(allocatable_array_pointers), optional :: molecules(:)
    logical, intent(in), optional :: nneighb_only
    type(Connection), intent(in), optional, target :: alt_connect

    integer, pointer :: mol_id(:)
    integer :: last_seed, last_mol_id, seed_at
    logical :: added_something
    type(Table) :: cur_molec, next_atoms

    !if property not present
    call add_property(at,'mol_id',0)
    if (.not. assign_pointer(at,'mol_id',mol_id)) &
      call system_abort("find_molecule_ids can't assign mol_id")
    !if present
    mol_id=0

    last_seed = 0
    last_mol_id = 0
    do while (any(mol_id == 0))
      do seed_at=last_seed+1, at%N
	if (mol_id(seed_at) /= 0) cycle

	! find molecule from seed
	call initialise(cur_molec)
	call append(cur_molec, (/ seed_at, 0, 0, 0 /) )
	added_something = .true.
	! look for neighbours
	do while (added_something)
	   call initialise(next_atoms)
	   call bfs_step(at,cur_molec,next_atoms,nneighb_only = nneighb_only, min_images_only = .true., alt_connect=alt_connect)
	   if (next_atoms%N > 0) then
	    added_something = .true.
	    call append(cur_molec, next_atoms)
	   else
	    added_something = .false.
	   endif
	end do

	! set mol_id
	mol_id(cur_molec%int(1,1:cur_molec%N)) = last_mol_id+1
	! store in molecules
	if (present(molecules)) then
	  allocate(molecules(last_mol_id+1)%i_a(cur_molec%N))
	  molecules(last_mol_id+1)%i_a = cur_molec%int(1,1:cur_molec%N)
	endif
	call finalise(cur_molec)

	last_mol_id = last_mol_id + 1

	last_seed = seed_at
      end do ! seed_at
    end do ! while any(mol_id == 0)

    call finalise(cur_molec)
    call finalise(next_atoms)
  end subroutine find_molecule_ids


  !% This is the subroutine that reads residues from a library.
  !% Used by create_CHARMM, can read from CHARMM and AMBER residue libraries.
  !% do_CHARMM=.true. is the default
  !
  subroutine next_motif(library,res_name,pdb_name,motif,atom_names,atom_charges,atom_names_PDB, n_impr,imp_atoms,do_CHARMM)
    
    type(Inoutput),                   intent(in)  :: library
    character(4),                     intent(out) :: res_name
    character(3),                     intent(out) :: pdb_name
    integer,             allocatable, intent(out) :: motif(:,:)
    character(4),        allocatable, intent(out) :: atom_names(:), atom_names_PDB(:)
    real(dp), optional,  allocatable, intent(out) :: atom_charges(:)
    logical,  optional,               intent(in)  :: do_CHARMM
    integer,  optional,  allocatable, intent(out) :: imp_atoms(:,:)
    integer,  optional,               intent(out) :: n_impr

    character(20), dimension(10) :: fields
    integer                      :: status, num_fields, data(7), i, n_at, max_num_fields
    type(Table)                  :: motif_table
    character(4)                 :: tmp_at_names(MAX_ATOMS_PER_RES), tmp_at_names_PDB(MAX_ATOMS_PER_RES)
    real(dp)                     :: tmp_at_charges(MAX_ATOMS_PER_RES),check_charge
    logical                      :: my_do_charmm
    character(len=1024)          :: line
   ! for improper generation
    integer                      :: imp_fields, tmp_imp_atoms(4,MAX_IMPROPERS_PER_RES)

    my_do_charmm = .true.
    if (present(do_CHARMM)) my_do_charmm = do_CHARMM

    if (my_do_charmm) then
       max_num_fields = 10
       imp_fields = 5
    else !do AMBER
       max_num_fields = 8
    endif

    status = 0

    do while(status==0)
       line = read_line(library,status)
       if (line(1:8)=='%residue') exit
    end do
    
    if (status/=0) then
       res_name = 'NONE'
       pdb_name = 'NON'
       return
    end if

    call parse_string(line,' ',fields,num_fields)
    res_name = trim(adjustl(fields(3)))
    pdb_name = trim(adjustl(fields(4)))

    call allocate(motif_table,7,0,0,0,20)
    n_at = 0
    check_charge=0._dp
   ! residue structure [& charges]
    do
       call parse_line(library,' ',fields,num_fields)
       if (num_fields < max_num_fields-1) exit ! last column of protein library (atom_name_PDB) is optional
       do i = 1, 7
          data(i) = string_to_int(fields(i))
       end do
       call append(motif_table,data)
       n_at = n_at + 1
       tmp_at_names(n_at) = fields(8)
       if (my_do_charmm) then
          tmp_at_charges(n_at) = string_to_real(fields(9))
          check_charge=check_charge+tmp_at_charges(n_at)
          if(num_fields == 10) then
             tmp_at_names_PDB(n_at) = fields(10)
          else
             tmp_at_names_PDB(n_at) = tmp_at_names(n_at)
          end if
       endif
    end do
    if (my_do_charmm) then
      ! intra amino acid IMPROPER generation here
       n_impr = 0
       do
          if (num_fields < imp_fields) exit
          if (trim(fields(1)).ne.'IMP') call system_abort('wrong improper format, should be: "IMP 1 4 7 10" with the 1st in the middle')
          n_impr = n_impr + 1
          do i = 1,4
             tmp_imp_atoms(i,n_impr) = string_to_int(fields(i+1))
          enddo
          call parse_line(library,' ',fields,num_fields)
       enddo
    endif

!    if (abs(mod(check_charge,1.0_dp)).ge.0.0001_dp .and. &
!        abs(mod(check_charge,1.0_dp)+1._dp).ge.0.0001_dp .and. &    !for -0.9999...
!        abs(mod(check_charge,1.0_dp)-1._dp).ge.0.0001_dp) then    !for +0.9999...
!       call print('WARNING next_motif: Charge of '//res_name//' residue is :'//round(check_charge,4),verbosity=ERROR)
!    endif

    allocate(motif(motif_table%N,7))

    motif = transpose(int_part(motif_table))

    allocate(atom_names(n_at))
    atom_names = tmp_at_names(1:n_at)
    allocate(atom_names_PDB(n_at))
    if (my_do_charmm) then
       allocate(atom_charges(n_at))
       atom_charges = tmp_at_charges(1:n_at)
       allocate(imp_atoms(4,n_impr))
       imp_atoms(1:4,1:n_impr) = tmp_imp_atoms(1:4,1:n_impr)
       atom_names_PDB = tmp_at_names_PDB(1:n_at)
    else
       atom_names_PDB = atom_names
    endif

    call finalise(motif_table)
    
  end subroutine next_motif

  !% Writes PDB format using the pdb_format passed as an input
  !% printing charges into the last column of PDB file 
  !% Sample lines:
  !%ATOM      1  CT3 ALA A   1       0.767   0.801  13.311  0.00  0.00     ALA   C  -0.2700
  !%ATOM      2   HA ALA A   1       0.074  -0.060  13.188  0.00  0.00     ALA   H   0.0900
  !%ATOM      3   HA ALA A   1       0.176   1.741  13.298  0.00  0.00     ALA   H   0.0900
  !
  subroutine write_pdb_file(at,pdb_file,pdb_format,run_type_string)

    character(len=*),           intent(in) :: pdb_file
    type(Atoms),                intent(in) :: at
    character(len=*),           intent(in) :: pdb_format
    character(len=*), optional, intent(in) :: run_type_string

    type(Inoutput)           :: pdb
    character(103)           :: sor !Brookhaven only uses 88
    integer                  :: mm
    character(4)             :: QM_prefix_atom_mol_name
    integer                  :: atom_type_index, &
                                atom_res_name_index, &
                                atom_mol_name_index, &
                                atom_res_number_index, &
                                atom_charge_index
    character(len=1024)      :: my_run_type_string
    real(dp)                 :: cell_lengths(3), cell_angles(3)

    call system_timer('write_pdb_file')
    my_run_type_string = optional_default('',run_type_string)

    call initialise(pdb,trim(pdb_file),action=OUTPUT)
    call print('   PDB file: '//trim(pdb%filename))
!    call print('REMARK'//at%N,file=pdb)
!lattice information could be added in a line like this:
!CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1          
    call lattice_xyz_to_abc(at%lattice, cell_lengths, cell_angles)
    sor=''
    write(sor, '(a6,3f9.3,3f7.2,a16)') 'CRYST1', cell_lengths(:), cell_angles(:), ' P 1           1'
    call print(sor, file=pdb)

!    if ((trim(my_run_type_string).eq.'QMMM_CORE') .or. &
!        (trim(my_run_type_string).eq.'QMMM_EXTENDED')) then
!       qm_flag_index = get_property(at,'cluster_mark')
!       if (trim(my_run_type_string).eq.'QMMM_EXTENDED') run_type=QMMM_RUN_EXTENDED
!       if (trim(my_run_type_string).eq.'QMMM_CORE') run_type=QMMM_RUN_CORE
!    endif
    atom_type_index = get_property(at,'atom_type_PDB')
    atom_res_name_index = get_property(at,'atom_res_name')
    atom_mol_name_index = get_property(at,'atom_mol_name')
    atom_res_number_index = get_property(at,'atom_res_number')
    atom_charge_index = get_property(at,'atom_charge')


    do mm=1,at%N
      ! e.g. CP2K needs different name for QM molecules, if use isolated atoms
       sor = ''
       QM_prefix_atom_mol_name = ''
       QM_prefix_atom_mol_name = trim(at%data%str(atom_mol_name_index,mm))
!does not work for contiguous molecules, e.g. silica
!!!       if ((trim(run_type_string).eq.'QMMM_CORE') .or. &
!!!           (trim(run_type_string).eq.'QMMM_EXTENDED')) then
!!!          if (at%data%int(qm_flag_index,mm).ge.QMMM_RUN_CORE .and. at%data%int(qm_flag_index,mm).le.run_type) then
!!!             QM_prefix_atom_mol_name = 'QM'//trim(at%data%str(atom_mol_name_index,mm))
!!!!             call print('QM molecule '//QM_prefix_atom_mol_name)
!!!          endif
!!!       endif
!             call print('molecule '//QM_prefix_atom_mol_name)
!       call print('writing PDB file: atom type '//at%data%str(atom_type_index,mm))
       write(sor,trim(pdb_format)) 'ATOM  ',mm,at%data%str(atom_type_index,mm),at%data%str(atom_res_name_index,mm),at%data%int(atom_res_number_index,mm), &
                             at%pos(1:3,mm),0._dp,0._dp,QM_prefix_atom_mol_name,ElementName(at%Z(mm)),at%data%real(atom_charge_index,mm)
!                             at%pos(1:3,mm),0._dp,0._dp,this%atom_res_name(mm),ElementName(at%Z(mm)),this%atom_charge(mm)
       call print(sor,file=pdb)
    enddo
    call print('END',file=pdb)

    call finalise(pdb)

    call system_timer('write_pdb_file')

  end subroutine write_pdb_file

  !% Writes Brookhaven PDB format
  !% charges are printed into the PSF file, too
  !% use CHARGE_EXTENDED keyword i.e. reads charges from the last column of PDB file
  !% Sample line:
  !%ATOM      1  CT3 ALA A   1       0.767   0.801  13.311  0.00  0.00     ALA   C  -0.2700
  !%ATOM      2   HA ALA A   1       0.074  -0.060  13.188  0.00  0.00     ALA   H   0.0900
  !%ATOM      3   HA ALA A   1       0.176   1.741  13.298  0.00  0.00     ALA   H   0.0900
  !
  subroutine write_brookhaven_pdb_file(at,pdb_file,run_type_string)

    character(len=*),           intent(in) :: pdb_file
    type(atoms),                intent(in) :: at
    character(len=*), optional, intent(in) :: run_type_string

!    character(*), parameter  :: pdb_format = '(a6,i5,1x,a4,1x,a4,1x,i4,1x,3x,3f8.3,2f6.2,10x,a2,2x,f7.4)'
    character(*), parameter  :: pdb_format = '(a6,i5,1x,a4,1x,a4,i5,1x,3x,3f8.3,2f6.2,6x,a4,1x,a2,2x,f7.4)'

  !Brookhaven PDB format
  !       sor(1:6)   = 'ATOM  '
  !       sor(7:11)  = mm
  !       sor(13:16) = this%atom_type(mm)
  !       sor(18:21) = this%res_name(mm)
  !!      sor(22:22) = ' A'                   !these two are now
  !       sor(23:26) = residue_number(mm)     !   merged to handle >9999 residues
  !       sor(31:38) = at%pos(1,mm)
  !       sor(39:46) = at%pos(2,mm)
  !       sor(47:54) = at%pos(3,mm)
  !       sor(55:60) = '  0.00'
  !       sor(61:66) = '  0.00'
  !!      sor(72:75) = 'MOL1'
  !       sor(77:78) = ElementName(at%Z(mm))
  !       sor(79:86) = this%atom_charge(mm)

    call system_timer('write_brookhaven_pdb_file')

    call write_pdb_file(at,pdb_file,pdb_format,run_type_string)

    call system_timer('write_brookhaven_pdb_file')

  end subroutine write_brookhaven_pdb_file

  !% Writes modified PDB format for accurate coordinates and charges, for CP2K
  !% Use CHARGE_EXTENDED keyword i.e. reads charges from the last column of PDB file 
  !% Sample line:
  !%ATOM      1  CT3 ALA A   1       0.76700000   0.80100000  13.31100000  0.00  0.00     ALA   C  -0.27
  !%ATOM      2   HA ALA A   1       0.07400000  -0.06000000  13.18800000  0.00  0.00     ALA   H   0.09
  !%ATOM      3   HA ALA A   1       0.17600000   1.74100000  13.29800000  0.00  0.00     ALA   H   0.09
  !
  subroutine write_cp2k_pdb_file(at,pdb_file,run_type_string)

    character(len=*),  intent(in)  :: pdb_file
    type(atoms),       intent(in)  :: at
    character(len=*), optional, intent(in) :: run_type_string

!    character(*), parameter  :: pdb_format = '(a6,i5,1x,a4,1x,a4,1x,i4,1x,3x,3f13.8,2f6.2,10x,a2,2x,f6.4)'
    character(*), parameter  :: pdb_format = '(a6,i5,1x,a4,1x,a4,i5,1x,3x,3f13.8,2f6.2,5x,a4,1x,a2,2x,f7.4)'

  !CP2K modified PDB format 
  !       sor(1:6)   = 'ATOM  '
  !       sor(7:11)  = mm
  !       sor(13:16) = this%atom_type(mm)
  !       sor(18:21) = this%res_name(mm)
  !!      sor(22:22) = ' A'                   !these two are now
  !       sor(23:26) = residue_number(mm)     !   merged to handle >9999 residues
  !       sor(31:43) = at%pos(1,mm)
  !       sor(44:56) = at%pos(2,mm)
  !       sor(57:69) = at%pos(3,mm)
  !       sor(70:75) = '  0.00'
  !       sor(76:81) = '  0.00'
  !!      sor(87:90) = 'MOL1'
  !       sor(92:93) = ElementName(at%Z(mm))
  !       sor(96:)   = this%atom_charge(mm)

    call system_timer('write_cp2k_pdb_file')

    call write_pdb_file(at,pdb_file,pdb_format,run_type_string)

    call system_timer('write_cp2k_pdb_file')

  end subroutine write_cp2k_pdb_file

  subroutine write_psf_file_arb_pos(at,psf_file,run_type_string,intrares_impropers,imp_filename,add_silica_23body,pos_field_for_connectivity)
    character(len=*),           intent(in) :: psf_file
    type(atoms),                intent(inout) :: at
    character(len=*), optional, intent(in) :: run_type_string
    type(Table),      optional, intent(in) :: intrares_impropers
    character(80),    optional, intent(in) :: imp_filename
    logical,          optional, intent(in) :: add_silica_23body
    character(len=*), optional, intent(in) :: pos_field_for_connectivity

    real(dp), pointer :: use_pos(:,:)
    type(Connection) :: t_connect
    type(Atoms) :: at_copy
    !character(len=TABLE_STRING_LENGTH), pointer :: atom_res_name_p(:)
    logical :: use_pos_is_pos, do_add_silica_23body

    ! save a copy
    at_copy = at

    do_add_silica_23body = optional_default(.false., add_silica_23body)

    ! find desired position field (pos, avgpos, whatever)
    use_pos_is_pos = .false.
    if (present(pos_field_for_connectivity)) then
      if (.not. assign_pointer(at, trim(pos_field_for_connectivity), use_pos)) &
	call system_abort("calc_topology can't find pos field '"//trim(pos_field_for_connectivity)//"'")
      if (trim(pos_field_for_connectivity) == 'pos') use_pos_is_pos = .true.
    else
      if (.not. assign_pointer(at, 'avgpos', use_pos)) &
	call system_abort("calc_topology can't find default pos field avgpos")
    endif

    ! copy desired pos to pos, and new connectivity
    !NB don't do if use_pos => pos
    if (.not. use_pos_is_pos) at_copy%pos = use_pos
    if (do_add_silica_23body) then
      call set_cutoff(at_copy,SILICA_2body_CUTOFF)
    else
      call set_cutoff(at_copy,0.0_dp)
    endif
    call calc_connect(at_copy, alt_connect=t_connect)

    ! now create labels using this connectivity object
    ! if cutoff is set to 0, nneighb_only doesn't matter
    ! if cutoff is large, then nneighb_only must be true
    !    so might as well pass true
    call write_psf_file (at,psf_file,run_type_string,intrares_impropers,imp_filename,add_silica_23body,nneighb_only=.true.,alt_connect=t_connect)
    call finalise(t_connect)
  end subroutine write_psf_file_arb_pos

  !% Writes PSF topology file, to be used with the PDB coordinate file.
  !% PSF contains the list of atoms, bonds, angles, impropers, dihedrals.
  !
  subroutine write_psf_file(at,psf_file,run_type_string,intrares_impropers,imp_filename,add_silica_23body,nneighb_only,alt_connect)

    character(len=*),           intent(in) :: psf_file
    type(atoms),                intent(in) :: at
    character(len=*), optional, intent(in) :: run_type_string
    type(Table),      optional, intent(in) :: intrares_impropers
    character(80),    optional, intent(in) :: imp_filename
    logical,          optional, intent(in) :: add_silica_23body
    logical, intent(in), optional :: nneighb_only
    type(Connection), intent(in), optional, target :: alt_connect

    type(Inoutput)          :: psf
    character(103)          :: sor
    character(*), parameter :: psf_format = '(I8,1X,A4,I5,1X,A4,1X,A4,1X,A4,1X,2G14.6,I8)'
    character(*), parameter :: title_format = '(I8,1X,A)'
    character(*), parameter :: int_format = 'I8'
    integer                 :: mm, i
    character(4)            :: QM_prefix_atom_mol_name
    integer                 :: atom_type_index, &
                               atom_type_PDB_index, &
                               atom_res_name_index, &
                               atom_mol_name_index, &
                               atom_res_number_index, &
                               atom_charge_index
    type(Table)             :: bonds
    type(Table)             :: angles
    type(Table)             :: dihedrals, impropers
    character(len=1024)     :: my_run_type_string
    !integer                 :: run_type
    logical                 :: do_add_silica_23body

    call system_timer('write_psf_file_pos')

    !intraresidual impropers: table or read in from file
    if (.not.present(intrares_impropers).and..not.present(imp_filename)) call print('WARNING!!! NO INTRARESIDUAL IMPROPERS USED!',verbosity=ERROR)
    if (present(imp_filename)) then
       call system_abort('Not yet implemented.')
    endif

    do_add_silica_23body = optional_default(.false.,add_silica_23body)
    if (do_add_silica_23body) then !!xxx there must be a SIO2 residue?
       if (at%cutoff.lt.SILICA_2body_CUTOFF) call system_abort('The connect cutoff '//at%cutoff//' is smaller than the required cutoff for silica. Cannot build connectivity to silica.')
    endif

    my_run_type_string = optional_default('',run_type_string)
!    if ((trim(my_run_type_string).eq.'QMMM_CORE') .or. &
!        (trim(my_run_type_string).eq.'QMMM_EXTENDED')) then
!       qm_flag_index = get_property(at,'cluster_mark')
!       if (trim(my_run_type_string).eq.'QMMM_EXTENDED') run_type=QMMM_RUN_EXTENDED
!       if (trim(my_run_type_string).eq.'QMMM_CORE') run_type=QMMM_RUN_CORE
!    endif
    atom_type_index = get_property(at,'atom_type')
    atom_type_PDB_index = get_property(at,'atom_type_PDB')
    atom_res_name_index = get_property(at,'atom_res_name')
    atom_mol_name_index = get_property(at,'atom_mol_name')
    atom_res_number_index = get_property(at,'atom_res_number')
    atom_charge_index = get_property(at,'atom_charge')
    call initialise(psf,trim(psf_file),action=OUTPUT)
    call print('   PSF file: '//trim(psf%filename))

    call print('PSF',file=psf)
    call print('',file=psf)

    write(sor,title_format) 1,'!NTITLE'
    call print(sor,file=psf)
    write(sor,'(A)') '  PSF file generated by libAtoms -- http://www.libatoms.org'
    call print(sor,file=psf)
    call print('',file=psf)

   ! ATOM section
    write(sor,title_format) at%N, '!NATOM'
    call print(sor,file=psf)
    do mm=1,at%N
       QM_prefix_atom_mol_name = ''
       QM_prefix_atom_mol_name = trim(at%data%str(atom_mol_name_index,mm))
!does not work for contiguous molecules, e.g. silica
!!!       if ((trim(my_run_type_string).eq.'QMMM_CORE') .or. &
!!!           (trim(my_run_type_string).eq.'QMMM_EXTENDED')) then
!!!          if (at%data%int(qm_flag_index,mm).ge.QMMM_RUN_CORE .and. at%data%int(qm_flag_index,mm).le.run_type) then
!!!             QM_prefix_atom_mol_name = 'QM'//trim(at%data%str(atom_mol_name_index,mm))
!!!!             call print('QM molecule '//QM_prefix_atom_mol_name)
!!!          endif
!!!       endif
!       call print('molecule '//QM_prefix_atom_mol_name)
!       call print('writing PSF file: atom type '//at%data%str(atom_type_index,mm))
       write(sor,psf_format) mm, QM_prefix_atom_mol_name, at%data%int(atom_res_number_index,mm), &
                     at%data%str(atom_res_name_index,mm),at%data%str(atom_type_PDB_index,mm),at%data%str(atom_type_index,mm), &
                     at%data%real(atom_charge_index,mm),ElementMass(at%Z(mm))/MASSCONVERT,0
       call print(sor,file=psf)
    enddo
    call print('',file=psf)
call print('PSF| '//at%n//' atoms')
   ! BOND section
    call create_bond_list(at,bonds,do_add_silica_23body,nneighb_only,alt_connect=alt_connect)
    if (any(bonds%int(1:2,1:bonds%N).le.0) .or. any(bonds%int(1:2,1:bonds%N).gt.at%N)) &
       call system_abort('write_psf_file_pos: element(s) of bonds not within (0;at%N]')
    call write_psf_section(data_table=bonds,psf=psf,section='BOND',int_format=int_format,title_format=title_format)
call print('PSF| '//bonds%n//' bonds')

   ! ANGLE section
    call create_angle_list(at,bonds,angles,do_add_silica_23body,nneighb_only,alt_connect=alt_connect)
    if (any(angles%int(1:3,1:angles%N).le.0) .or. any(angles%int(1:3,1:angles%N).gt.at%N)) then
       do i = 1, angles%N
          if (any(angles%int(1:3,i).le.0) .or. any(angles%int(1:3,i).gt.at%N)) &
          call print('angle: '//angles%int(1,i)//' -- '//angles%int(2,i)//' -- '//angles%int(3,i),verbosity=ERROR)
       enddo
       call system_abort('write_psf_file_pos: element(s) of angles not within (0;at%N]')
    endif
    call write_psf_section(data_table=angles,psf=psf,section='THETA',int_format=int_format,title_format=title_format)
call print('PSF| '//angles%n//' angles')

   ! DIHEDRAL section
    call create_dihedral_list(at,angles,dihedrals,do_add_silica_23body,nneighb_only,alt_connect=alt_connect)
    if (any(dihedrals%int(1:4,1:dihedrals%N).le.0) .or. any(dihedrals%int(1:4,1:dihedrals%N).gt.at%N)) &
       call system_abort('write_psf_file_pos: element(s) of dihedrals not within (0;at%N]')
    call write_psf_section(data_table=dihedrals,psf=psf,section='PHI',int_format=int_format,title_format=title_format)
call print('PSF| '//dihedrals%n//' dihedrals')

   ! IMPROPER section
    call create_improper_list(at,angles,impropers,intrares_impropers=intrares_impropers)
    if (any(impropers%int(1:4,1:impropers%N).le.0) .or. any(impropers%int(1:4,1:impropers%N).gt.at%N)) &
       call system_abort('write_psf_file_pos: element(s) of impropers not within (0;at%N]')
    call write_psf_section(data_table=impropers,psf=psf,section='IMPHI',int_format=int_format,title_format=title_format)
call print('PSF| '//impropers%n//' impropers')

   !empty DON, ACC and NNB sections
    write(sor,title_format) 0,'!NDON'
    call print(sor,file=psf)
    call print('',file=psf)

    write(sor,title_format) 0,'!NACC'
    call print(sor,file=psf)
    call print('',file=psf)

    write(sor,title_format) 0,'!NNB'
    call print(sor,file=psf)
    call print('',file=psf)

    call print('END',file=psf)

    call finalise(bonds)
    call finalise(angles)
    call finalise(dihedrals)
    call finalise(impropers)
    call finalise(psf)

    call system_timer('write_psf_file_pos')
!call print('psf written. stop.')
!stop

  end subroutine write_psf_file

  !% Writes a section of the PSF topology file.
  !% Given a $section$ (name of the section) and a table $data_table$ with the data,
  !% it writes them into $psf$ file with $int_format$ format.
  !% Used by write_psf_file.
  !
  subroutine write_psf_section(data_table,psf,section,int_format,title_format)

    type(Table),      intent(in) :: data_table
    type(InOutput),   intent(in) :: psf
    character(len=*), intent(in) :: section
    character(*),     intent(in) :: int_format
    character(*),     intent(in) :: title_format

    character(len=103)           :: sor
    integer                      :: mm, i, num_per_line

    if (.not.any(size(data_table%int,1).eq.(/2,3,4/))) &
       call system_abort('data table to print into psf file has wrong number of integers '//size(data_table%int,1))
    if (size(data_table%int,1).eq.2) num_per_line = 4
    if (size(data_table%int,1).eq.3) num_per_line = 3
    if (size(data_table%int,1).eq.4) num_per_line = 2

    sor = ''
    write(sor,title_format) data_table%N, '!N'//trim(section)
    call print(sor,file=psf)

    mm = 1
    do while (mm.le.(data_table%N-num_per_line+1-mod(data_table%N,num_per_line)))
       select case(size(data_table%int,1))
         case(2)
           write(sor,'(8'//trim(int_format)//')') &
              data_table%int(1,mm),   data_table%int(2,mm), &
              data_table%int(1,mm+1), data_table%int(2,mm+1), &
              data_table%int(1,mm+2), data_table%int(2,mm+2), &
              data_table%int(1,mm+3), data_table%int(2,mm+3)
           mm = mm + 4
         case(3)
           write(sor,'(9'//trim(int_format)//')') &
              data_table%int(1,mm),   data_table%int(2,mm),   data_table%int(3,mm), &
              data_table%int(1,mm+1), data_table%int(2,mm+1), data_table%int(3,mm+1), &
              data_table%int(1,mm+2), data_table%int(2,mm+2), data_table%int(3,mm+2)
           mm = mm + 3
         case(4)
           write(sor,'(8'//trim(int_format)//')') &
              data_table%int(1,mm),   data_table%int(2,mm),   data_table%int(3,mm), data_table%int(4,mm), &
              data_table%int(1,mm+1), data_table%int(2,mm+1), data_table%int(3,mm+1), data_table%int(4,mm+1)
           mm = mm + 2
       end select
       call print(sor,file=psf)
    enddo

   ! mm = data_table%N - mod(data_table%N,num_per_line) + 1
    sor = ''
    do i=1, mod(data_table%N,num_per_line) !if 0 then it does nothing
       select case(size(data_table%int,1))
         case(2)
           write(sor((i-1)*16+1:i*16),'(2'//trim(int_format)//')') data_table%int(1,mm), data_table%int(2,mm)
         case(3)
           write(sor((i-1)*24+1:i*24),'(3'//trim(int_format)//')') data_table%int(1,mm), data_table%int(2,mm), data_table%int(3,mm)
         case(4)
           write(sor((i-1)*32+1:i*32),'(4'//trim(int_format)//')') data_table%int(1,mm), data_table%int(2,mm), data_table%int(3,mm), data_table%int(4,mm)
       end select
       mm = mm + 1
    enddo

    if (mm .ne. data_table%N+1) call system_abort('psf writing: written '//(mm-1)//' of '//data_table%N)
    if (mod(data_table%N,num_per_line).ne.0) call print(sor,file=psf)
    call print('',file=psf)

  end subroutine write_psf_section

  !% Subroutine to create a $bonds$ table with all the nearest neighbour bonds
  !% according to the connectivity. If $add_silica_23body$ is true, include
  !% atom pairs that have SIO2 molecule type up to the silica_cutoff distance.
  !% Optionally use hysteretic_connect, if it is passed as alt_connect.
  !
  subroutine create_bond_list(at,bonds,add_silica_23body,nneighb_only,alt_connect)

  type(Atoms), intent(in)  :: at
  type(Table), intent(out) :: bonds
  logical,     intent(in)  :: add_silica_23body
  logical, intent(in), optional :: nneighb_only
  type(Connection), intent(in), optional, target :: alt_connect

    character(*), parameter  :: me = 'create_bond_list: '

  type(Table) :: atom_a,atom_b
  integer     :: i,j
  integer     :: atom_j
!  logical              :: do_qmmm
!  integer,dimension(3) :: pos_indices
!  integer              :: qm_flag_index
  integer              :: atom_mol_name_index
  logical :: add_bond

    call system_timer('create_bond_list')


    if (add_silica_23body) then
       if (at%cutoff.lt.SILICA_2body_CUTOFF) call system_abort('The connect cutoff '//at%cutoff//' is smaller than the required cutoff for silica. Cannot build connectivity to silica.')
       atom_mol_name_index = get_property(at,'atom_mol_name')
    endif
    call initialise(bonds,2,0,0,0,0)

!    do_qmmm = .true.
!    if (get_value(at%properties,trim('QM_flag'),pos_indices)) then
!       qm_flag_index = pos_indices(2)
!    else
!       do_qmmm = .false.
!    end if
!
    do i=1,at%N
       call initialise(atom_a,4,0,0,0,0)
       call append(atom_a,(/i,0,0,0/))
       if (add_silica_23body) then
          call bfs_step(at,atom_a,atom_b,nneighb_only=.false.,min_images_only=.true.,alt_connect=alt_connect) ! SILICON_2BODY_CUTOFF is not within nneigh_tol
       else
          call bfs_step(at,atom_a,atom_b,nneighb_only=nneighb_only,min_images_only=.true.,alt_connect=alt_connect)
       endif
       do j = 1,atom_b%N
          atom_j = atom_b%int(1,j)
          if (atom_j.gt.i) then
!      ! QM core atoms should be isolated atoms
!!call print('atom '//i//' (QM flag '//at%data%int(qm_flag_index,i)//')')
!!call print('atom '//atom_j//' (QM flag '//at%data%int(qm_flag_index,atom_j)//')')
!
!             if (do_qmmm) then
!                if (any((/at%data%int(qm_flag_index,i),at%data%int(qm_flag_index,atom_j)/).eq.1)) then
!!                   call print('not added '//i//' (QM flag '//at%data%int(qm_flag_index,i)//') -- '//atom_j//' (QM flag '//at%data%int(qm_flag_index,atom_j)//')')
!                   cycle
!                else
!!                   call print('added '//i//' (QM flag '//at%data%int(qm_flag_index,i)//') -- '//atom_j//' (QM flag '//at%data%int(qm_flag_index,atom_j)//')')
!                endif
!             endif

             add_bond = .true.

             if (add_silica_23body) then ! do not include H2O -- SiO2 bonds
                if ( ((trim(at%data%str(atom_mol_name_index,atom_j)) .eq.'SIO2' .and. trim(at%data%str(atom_mol_name_index,i)).ne.'SIO2')) .or. &
                     ((trim(at%data%str(atom_mol_name_index,atom_j)) .ne.'SIO2' .and. trim(at%data%str(atom_mol_name_index,i)).eq.'SIO2')) ) then !silica -- something
!call system_abort('should have not got here')
                   !add only nearest neighbours
                   if (.not.(is_nearest_neighbour_abs_index(at,i,atom_j,alt_connect=alt_connect))) add_bond = .false.
                elseif  ((trim(at%data%str(atom_mol_name_index,atom_j)) .eq.'SIO2' .and. trim(at%data%str(atom_mol_name_index,i)).eq.'SIO2')) then !silica -- silica
                   !add atom pairs within SILICON_2BODY_CUTOFF
!call system_abort('what not?')
                   if (.not.are_silica_2body_neighbours(at,i,atom_j,alt_connect=alt_connect)) add_bond = .false. !SILICON_2BODY_CUTOFF for Si-Si and Si-O, nearest neighbours otherwise(Si-H,O-O,O-H,H-H)
                else
!call system_abort('should have not got here')
                   !add only nearest neighbours
                   if (.not.(is_nearest_neighbour_abs_index(at,i,atom_j,alt_connect=alt_connect))) add_bond = .false.
                endif
             endif

             if (add_bond) then
                call append(bonds,(/i,atom_j/))
             endif

          else
!             call print('not added '//i//' -- '//atom_j)
          endif
       enddo
       call finalise(atom_a)
       call finalise(atom_b)
    enddo

    if (any(bonds%int(1:2,1:bonds%N).le.0) .or. any(bonds%int(1:2,1:bonds%N).gt.at%N)) &
       call system_abort('create_bond_list: element(s) of bonds not within (0;at%N]')

    call system_timer('create_bond_list')

  end subroutine create_bond_list


  !% Test if atom $j$ is a neighbour of atom $i$ for the silica 2 body terms.
  !% Optionally use hysteretic_connect, if it is passed as alt_connect.
  !%
  function are_silica_2body_neighbours(this,i,j,alt_connect)

    type(Atoms), intent(in), target :: this
    integer,     intent(in) :: i,j
    type(Connection), intent(in), optional, target :: alt_connect
    logical                 :: are_silica_2body_neighbours

    real(dp)                :: d
    integer :: neigh, Z_i, Z_j, j2

    are_silica_2body_neighbours = .false.

    do neigh = 1, atoms_n_neighbours(this,i)
       j2 = atoms_neighbour(this,i,neigh,distance=d,alt_connect=alt_connect)
       if (j2.eq.j) then
          Z_i = this%Z(i)
          Z_j = this%Z(j)
          if ( ((Z_i.eq.14).and.(Z_j.eq. 8)) .or. & !Si--O
               ((Z_i.eq. 8).and.(Z_j.eq.14)) .or. & !O--Si
               ((Z_i.eq. 8).and.(Z_j.eq. 8)) ) then !O--O
             if (d < SILICA_2BODY_CUTOFF) &
                are_silica_2body_neighbours = .true.
          elseif ((Z_i.eq.14).and.(Z_j.eq.14))  then !Si--Si
             if (d < SILICON_2BODY_CUTOFF) &
                are_silica_2body_neighbours = .true.
!call print(i//'--'//j//': '//d//' < 2.8? '//are_silica_2body_neighbours)
          else !Si--H, O--H, O--O, H--H
             if (d < (bond_length(Z_i,Z_j)*this%nneightol)) &
                are_silica_2body_neighbours = .true.
!call print(i//'--'//j//': '//d//' < nneigh_tol? '//are_silica_2body_neighbours)
          endif
          return
       endif
    enddo

    call system_abort('are_silica_2body_neighbours: atom '//j//'is not a neighbour of atom '//i)

  end function are_silica_2body_neighbours

  !% Test if atom $j$ is a nearest neighbour of atom $i$ for the silica 3 body terms.
  !% Optionally use hysteretic_connect, if it is passed as alt_connect.
  !%
  function are_silica_nearest_neighbours(this,i,j,alt_connect)

    type(Atoms), intent(in), target :: this
    integer,     intent(in) :: i,j
    type(Connection), intent(in), optional, target :: alt_connect
    logical                 :: are_silica_nearest_neighbours

    real(dp)                :: d
    integer :: neigh, Z_i, Z_j, j2

    are_silica_nearest_neighbours = .false.

    do neigh = 1, atoms_n_neighbours(this,i)
       j2 = atoms_neighbour(this,i,neigh,distance=d,alt_connect=alt_connect)
       if (j2.eq.j) then
          Z_i = this%Z(i)
          Z_j = this%Z(j)
          if ( ((Z_i.eq.14).and.(Z_j.eq. 8)) .or. & !Si--O
               ((Z_i.eq. 8).and.(Z_j.eq.14)) .or. & !O--Si
               ((Z_i.eq.14).and.(Z_j.eq.14)) ) then !Si--Si
             if (d < SILICON_2BODY_CUTOFF) &
                are_silica_nearest_neighbours = .true.
!call print(i//'--'//j//': '//d//' < 2.8? '//are_silica_nearest_neighbours)
          else !Si--H, O--H, O--O, H--H
             if (d < (bond_length(Z_i,Z_j)*this%nneightol)) &
                are_silica_nearest_neighbours = .true.
!call print(i//'--'//j//': '//d//' < nneigh_tol? '//are_silica_nearest_neighbours)
          endif
          return
       endif
    enddo

    call system_abort('are_silica_nearest_neighbours: atom '//j//'is not a neighbour of atom '//i)

  end function are_silica_nearest_neighbours

  !% Subroutine to create an $angles$ table with all the angles,
  !% according to the connectivity nearest neighbours. If $add_silica_23body$ is true, include
  !% atom triplets that have SIO2 molecule type and the distances are smaller than the cutoff
  !% for the corresponding 3 body cutoff.
  !% Optionally use hysteretic_connect if passed as alt_connect.
  !
  subroutine create_angle_list(at,bonds,angles,add_silica_23body,nneighb_only,alt_connect)

  type(Atoms), intent(in)  :: at
  type(Table), intent(in)  :: bonds
  type(Table), intent(out) :: angles
  logical,     intent(in)  :: add_silica_23body
  logical,     intent(in), optional  :: nneighb_only
  type(Connection), intent(in), optional, target :: alt_connect

    character(*), parameter  :: me = 'create_angle_list: '

  integer     :: i,j
  type(Table) :: atom_a, atom_b
  integer     :: atom_j, atom_1, atom_2
  logical     :: add_angle
  integer     :: atom_type_index
  integer     :: atom_mol_name_index

    call system_timer('create_angle_list')


    if (add_silica_23body) then
       atom_type_index = get_property(at,'atom_type') !different cutoff to define angles for different atom types
       atom_mol_name_index = get_property(at,'atom_mol_name')
    endif

    call initialise(angles,3,0,0,0,0)

    do i=1,bonds%N

       atom_1 = bonds%int(1,i)
       atom_2 = bonds%int(2,i)

      ! look for one more to the beginning: ??--1--2 where ??<2
       call initialise(atom_a,4,0,0,0,0)
       call append(atom_a,(/atom_1,0,0,0/))

       if (add_silica_23body) then
!call system_abort('stop!')
          call bfs_step(at,atom_a,atom_b,nneighb_only=.false.,min_images_only=.true.,alt_connect=alt_connect) ! SILICON_2BODY_CUTOFF is not within nneigh_tol
       else
          call bfs_step(at,atom_a,atom_b,nneighb_only=nneighb_only,min_images_only=.true.,alt_connect=alt_connect)
       endif

       do j = 1,atom_b%N
!call print('atom '//j//' out of '//atom_b%N//' which is atom '//atom_j)
          atom_j = atom_b%int(1,j)
          if (atom_j.ge.atom_2) cycle

          add_angle = .true.

          if (add_silica_23body) then ! do not include H2O -- SiO2 bonds
             if ( ((trim(at%data%str(atom_mol_name_index,atom_j)) .eq.'SIO2' .and. trim(at%data%str(atom_mol_name_index,atom_1)).ne.'SIO2')) .or. &
                  ((trim(at%data%str(atom_mol_name_index,atom_j)) .ne.'SIO2' .and. trim(at%data%str(atom_mol_name_index,atom_1)).eq.'SIO2')) ) then !silica -- something
                !add only nearest neighbours
                if (.not.(is_nearest_neighbour_abs_index(at,atom_1,atom_j,alt_connect=alt_connect))) add_angle = .false.
             elseif  ((trim(at%data%str(atom_mol_name_index,atom_j)) .eq.'SIO2' .and. trim(at%data%str(atom_mol_name_index,atom_1)).eq.'SIO2')) then !silica -- silica
                !add atom pairs within SILICON_2BODY_CUTOFF
                if (.not.are_silica_nearest_neighbours(at,atom_1,atom_j,alt_connect=alt_connect)) add_angle = .false.

                ! we kept 2.8 A cutoff for Si-Si and 5.5 A for Si-O and O-O: reduce it to 2.6 for Si-O-Si, O-Si-O, Si-O-H
                if (add_angle) then
                    if (any(at%Z(atom_j)*1000000+at%Z(atom_1)*1000+at%Z(atom_2).eq.(/14008014,8014008,14008001,1008014/))) then
                       add_angle =  (distance_min_image(at,atom_j,atom_1).le.SILICA_3BODY_CUTOFF).and. &
                                  (distance_min_image(at,atom_1,atom_2).le.SILICA_3BODY_CUTOFF)
                    elseif (any(at%Z(atom_j)*1000000+at%Z(atom_1)*1000+at%Z(atom_2).eq.(/14008008,8008014,1008008,8008001,8008008/))) then
                       add_angle = .false.
                    endif
                endif
             else !something -- something, neither are silica
                !add only nearest neighbours
                if (.not.(is_nearest_neighbour_abs_index(at,atom_1,atom_j,alt_connect=alt_connect))) add_angle = .false.
             endif
          endif


!          if (add_angle.and.(atom_j.lt.atom_2)) &
          if (add_angle) & !.and.(atom_j.lt.atom_2)) &
                call append(angles,(/atom_j,atom_1,atom_2/))

       enddo

       call finalise(atom_a)
       call finalise(atom_b)

      ! look for one more to the end: 1--2--?? where 1<??
       call initialise(atom_a,4,0,0,0,0)
       call append(atom_a,(/atom_2,0,0,0/))

       if (add_silica_23body) then
          call bfs_step(at,atom_a,atom_b,nneighb_only=.false.,min_images_only=.true.,alt_connect=alt_connect) ! SILICON_2BODY_CUTOFF is not within nneigh_tol
       else
          call bfs_step(at,atom_a,atom_b,nneighb_only=nneighb_only,min_images_only=.true.,alt_connect=alt_connect)
       endif

       do j = 1,atom_b%N

          atom_j = atom_b%int(1,j)
          if (atom_j.ge.atom_1) cycle

          add_angle = .true.

          if (add_silica_23body) then ! do not include H2O -- SiO2 bonds
             if ( ((trim(at%data%str(atom_mol_name_index,atom_j)) .eq.'SIO2' .and. trim(at%data%str(atom_mol_name_index,atom_2)).ne.'SIO2')) .or. &
                  ((trim(at%data%str(atom_mol_name_index,atom_j)) .ne.'SIO2' .and. trim(at%data%str(atom_mol_name_index,atom_2)).eq.'SIO2')) ) then !silica -- something
                !add only nearest neighbours
                if (.not.(is_nearest_neighbour_abs_index(at,atom_2,atom_j,alt_connect=alt_connect))) add_angle = .false.
             elseif  ((trim(at%data%str(atom_mol_name_index,atom_j)) .eq.'SIO2' .and. trim(at%data%str(atom_mol_name_index,atom_2)).eq.'SIO2')) then !silica -- silica
                !add atom pairs within SILICON_2BODY_CUTOFF
                if (.not.are_silica_nearest_neighbours(at,atom_2,atom_j,alt_connect=alt_connect)) add_angle = .false.

                ! we kept 2.8 A cutoff for Si-Si and 5.5 A for Si-O and O-O: reduce it to 2.6 for Si-O-Si, O-Si-O, Si-O-H
                if (add_angle) then
                    if (any(at%Z(atom_j)*1000000+at%Z(atom_2)*1000+at%Z(atom_1).eq.(/14008014,8014008,14008001,1008014/))) then
                       add_angle =  (distance_min_image(at,atom_j,atom_2).le.SILICA_3BODY_CUTOFF).and. &
                                  (distance_min_image(at,atom_1,atom_2).le.SILICA_3BODY_CUTOFF)
                    elseif (any(at%Z(atom_j)*1000000+at%Z(atom_2)*1000+at%Z(atom_1).eq.(/14008008,8008014,1008008,8008001,8008008/))) then
                       add_angle = .false.
                    endif
                endif
             else !something -- something, neither are silica
                !add only nearest neighbours
                if (.not.(is_nearest_neighbour_abs_index(at,atom_2,atom_j,alt_connect=alt_connect))) add_angle = .false.
             endif
          endif

!          if (add_angle.and.(atom_j.lt.atom_1)) &
          if (add_angle) & !.and.(atom_j.lt.atom_1)) &
             call append(angles,(/atom_1,atom_2,atom_j/))

       enddo
       call finalise(atom_a)
       call finalise(atom_b)

    enddo

    if (any(angles%int(1:3,1:angles%N).le.0) .or. any(angles%int(1:3,1:angles%N).gt.at%N)) then
       do i = 1, angles%N
          if (any(angles%int(1:3,i).le.0) .or. any(angles%int(1:3,i).gt.at%N)) &
          call print('angle: '//angles%int(1,i)//' -- '//angles%int(2,i)//' -- '//angles%int(3,i),verbosity=ERROR)
       enddo
       call system_abort('create_angle_list: element(s) of angles not within (0;at%N]')
    endif

    call system_timer('create_angle_list')

  end subroutine create_angle_list

  !% Create $dihedrals$ table simply use the $angles$ table and
  !% the connectivity nearest neighbours. If $add_silica_23body$ is true
  !% skip any dihedrals with Si atoms, as there are only 2 and 3 body terms.
  !% Optionally use hysteretic_connect if passed as alt_connect.
  !
  subroutine create_dihedral_list(at,angles,dihedrals,add_silica_23body,nneighb_only, alt_connect)

  type(Atoms), intent(in)  :: at
  type(Table), intent(in)  :: angles
  type(Table), intent(out) :: dihedrals
  logical,     intent(in)  :: add_silica_23body
  logical,     intent(in), optional  :: nneighb_only
  type(Connection), intent(in), optional, target :: alt_connect

    character(*), parameter  :: me = 'create_angle_list: '

  integer     :: i,j
  type(Table) :: atom_a, atom_b
  integer     :: atom_j

    call system_timer('create_dihedral_list')


    call initialise(dihedrals,4,0,0,0,0)

    do i=1,angles%N

       if (add_silica_23body) then !only add 2 and 3 body for silica, skip dihedrals
          if (any(at%Z(angles%int(1:3,i)).eq.14)) cycle
       endif

      ! look for one more to the beginning: ??--1--2--3
       call initialise(atom_a,4,0,0,0,0)
       call append(atom_a,(/angles%int(1,i),0,0,0/))
       call bfs_step(at,atom_a,atom_b,nneighb_only=nneighb_only,min_images_only=.true.,alt_connect=alt_connect)
       do j = 1,atom_b%N
          atom_j = atom_b%int(1,j)
          if (atom_j.ne.angles%int(2,i)) then
            ! make sure it's not included twice -- no need to O(N^2) check at the end
             if (atom_j.lt.angles%int(3,i)) then
                if (add_silica_23body) then !only add 2 and 3 body for silica, skip dihedrals
                   if (at%Z(atom_j).eq.14) cycle
                endif
                call append(dihedrals,(/atom_j,angles%int(1,i),angles%int(2,i),angles%int(3,i)/))
             endif
          endif
       enddo
       call finalise(atom_a)
       call finalise(atom_b)

      ! look for one more to the end: 1--2--3--??
       call initialise(atom_a,4,0,0,0,0)
       call append(atom_a,(/angles%int(3,i),0,0,0/))
       call bfs_step(at,atom_a,atom_b,nneighb_only=nneighb_only,min_images_only=.true.,alt_connect=alt_connect)
       do j = 1,atom_b%N
          atom_j = atom_b%int(1,j)
          if (atom_j.ne.angles%int(2,i)) then
            ! make sure it's not included twice -- no need to O(N^2) check at the end
             if (atom_j.lt.angles%int(1,i)) then
                if (add_silica_23body) then !only add 2 and 3 body for silica, skip dihedrals
                   if (at%Z(atom_j).eq.14) cycle
                endif
                call append(dihedrals,(/angles%int(1,i),angles%int(2,i),angles%int(3,i),atom_j/))
             endif
          endif
       enddo
       call finalise(atom_a)
       call finalise(atom_b)

    enddo

    if (any(dihedrals%int(1:4,1:dihedrals%N).le.0) .or. any(dihedrals%int(1:4,1:dihedrals%N).gt.at%N)) &
       call system_abort('create_dihedral_list: element(s) of dihedrals not within (0;at%N]')

    call system_timer('create_dihedral_list')

  end subroutine create_dihedral_list
  
  !% Create $impropers$ table simply use the $angles$ table.
  !% For intraresidual impropers, use the input table,
  !% the backbone residues can be calculated.
  !
  subroutine create_improper_list(at,angles,impropers,intrares_impropers)

  type(Atoms),           intent(in)  :: at
  type(Table),           intent(in)  :: angles
  type(Table),           intent(out) :: impropers
  type(Table), optional, intent(in)  :: intrares_impropers

  integer, dimension(4) :: imp_atoms
  integer               :: nn,mm
  logical               :: cont
  integer, allocatable, dimension(:) :: count_array ! to count number of bonds
  integer               :: i,j, i_impr
  integer               :: last, tmp
  integer               :: atom_res_name_index
  integer               :: atom_type_index
  integer               :: i_pro, tmp_atoms(3)
  logical               :: reordered

    call system_timer('create_improper_list')

    if (.not. at%connect%initialised) &
       call system_abort('create_bond_list: connectivity not initialised, call calc_connect first')

    call initialise(impropers,4,0,0,0,0)

    allocate (count_array(angles%N))

    do i = 1,at%N
      if (.not.any(trim(at%species(i)).eq.(/'C','N'/))) cycle
      count_array = 0
      where (angles%int(2,1:angles%N).eq.i) count_array = 1
      if (sum(count_array(1:size(count_array))).ne.3) cycle
     ! only N with a neighbour that has 3 neighbors can stay
      if (trim(at%species(i)).eq.'N') then
         cont = .false.
         !for the first X1-N-X2
         nn = find_in_array(angles%int(2,1:angles%N),i)
           !check X1
            count_array = 0
            where (angles%int(2,1:angles%N).eq.angles%int(1,nn)) count_array = 1
            if (sum(count_array(1:size(count_array))).ne.3) cont = .true.
           !check X2
            count_array = 0
            where (angles%int(2,1:angles%N).eq.angles%int(3,nn)) count_array = 1
            if (sum(count_array(1:size(count_array))).ne.3) cont = .true.
         !for the second X1-N-X3
         mm = find_in_array(angles%int(2,nn+1:angles%N),i)
           !check X1
            count_array = 0
            where (angles%int(2,1:angles%N).eq.angles%int(1,nn+mm)) count_array = 1
            if (sum(count_array(1:size(count_array))).ne.3) cont = .true.
           !check X3
            count_array = 0
            where (angles%int(2,1:angles%N).eq.angles%int(3,nn+mm)) count_array = 1
            if (sum(count_array(1:size(count_array))).ne.3) cont = .true.
         !no need to check X2-N-X3
         if (.not.cont) cycle
      endif

     ! add to impropers i and its neighbours
      imp_atoms = 0
      imp_atoms(1) = i
      !neighbours from first angle
      nn = find_in_array(angles%int(2,1:angles%N),i)
        j = nn
        imp_atoms(2) = angles%int(1,j)
        imp_atoms(3) = angles%int(3,j)
      !3rd neighbour from second angle
      mm = find_in_array(angles%int(2,nn+1:angles%N),i)
        j = nn+mm
        if (.not.any(angles%int(1,j).eq.imp_atoms(1:3))) then
           imp_atoms(4) = angles%int(1,j)
        else
           imp_atoms(4) = angles%int(3,j)
        endif

!VVV ORDER is done according to the topology file! - and is read in when finding motifs
!if you don't do this, you won't have only the backbone impropers!
      atom_res_name_index = get_property(at,'atom_res_name')
      if (all(at%data%str(atom_res_name_index,imp_atoms(2:4)).eq.at%data%str(atom_res_name_index,imp_atoms(1)))) &
         cycle ! these should be added when identifying the residues
!ORDER: check charmm.pot file - start with $i, end with  H or O or N, in this order -- for intraresidual residues this can be needed later on...
      reordered = .true.
      tmp = 0
      ! if there is H
      last = find_in_array(at%Z(imp_atoms(2:4)),1)
      if (last.gt.0) then
        last = last + 1
        tmp = imp_atoms(4)
        imp_atoms(4) = imp_atoms(last)
        imp_atoms(last) = tmp
!        call print('reordered H to the end in '// &
!                    trim(ElementName(at%Z(imp_atoms(1))))//imp_atoms(1)//'--'// &
!                    trim(ElementName(at%Z(imp_atoms(2))))//imp_atoms(2)//'--'// &
!                    trim(ElementName(at%Z(imp_atoms(3))))//imp_atoms(3)//'--'// &
!                    trim(ElementName(at%Z(imp_atoms(4))))//imp_atoms(4))
      else
        last = find_in_array(at%Z(imp_atoms(2:4)),8) ! at the C-terminal there should be one "CC X X OC", with the double bonding the last one
        if (last.gt.0) then
          last = last + 1
          tmp = imp_atoms(4)
          imp_atoms(4) = imp_atoms(last)
          imp_atoms(last) = tmp
!          call print('reordered O to the end in '// &
!                      trim(ElementName(at%Z(imp_atoms(1))))//imp_atoms(1)//'--'// &
!                      trim(ElementName(at%Z(imp_atoms(2))))//imp_atoms(2)//'--'// &
!                      trim(ElementName(at%Z(imp_atoms(3))))//imp_atoms(3)//'--'// &
!                      trim(ElementName(at%Z(imp_atoms(4))))//imp_atoms(4))
        else
          last = find_in_array(at%Z(imp_atoms(2:4)),7)
          if (last.gt.0) then
            last = last + 1
            tmp = imp_atoms(4)
            imp_atoms(4) = imp_atoms(last)
            imp_atoms(last) = tmp
!            call print('reordered N to the end in '// &
!                        trim(ElementName(at%Z(imp_atoms(1))))//imp_atoms(1)//'--'// &
!                        trim(ElementName(at%Z(imp_atoms(2))))//imp_atoms(2)//'--'// &
!                        trim(ElementName(at%Z(imp_atoms(3))))//imp_atoms(3)//'--'// &
!                        trim(ElementName(at%Z(imp_atoms(4))))//imp_atoms(4))
          else
            reordered = .false.
!            call print('not reordered improper '// &
!                        trim(ElementName(at%Z(imp_atoms(1))))//imp_atoms(1)//'--'// &
!                        trim(ElementName(at%Z(imp_atoms(2))))//imp_atoms(2)//'--'// &
!                        trim(ElementName(at%Z(imp_atoms(3))))//imp_atoms(3)//'--'// &
!                        trim(ElementName(at%Z(imp_atoms(4))))//imp_atoms(4))
          endif
        endif
      endif

      !checking and adding only backbone (i.e. not intraresidual impropers) where the order of the 2nd and 3rd atoms doesn't matter
      atom_res_name_index = get_property(at,'atom_res_name')
      if (all(at%data%str(atom_res_name_index,imp_atoms(2:4)).eq.at%data%str(atom_res_name_index,imp_atoms(1)))) &
         cycle ! these should be added when identifying the residues

      if (.not.reordered) then
        ! Found N-C-CP1-CP3 Pro backbone, reordering according to atomic types (could be also according to the H neighbours)
!         call print('|PRO Found Pro backbone')
         atom_type_index = get_property(at,'atom_type')
         if (trim(at%data%str(atom_type_index,imp_atoms(1))).ne.'N') call system_abort('something has gone wrong. what is this if not proline? '// &
                     trim(at%data%str(atom_type_index,imp_atoms(1)))//imp_atoms(1)//'--'// &
                     trim(at%data%str(atom_type_index,imp_atoms(2)))//imp_atoms(2)//'--'// &
                     trim(at%data%str(atom_type_index,imp_atoms(3)))//imp_atoms(3)//'--'// &
                     trim(at%data%str(atom_type_index,imp_atoms(4)))//imp_atoms(4))
         tmp_atoms = 0
         do i_pro = 2,4
            if (trim(at%data%str(atom_type_index,imp_atoms(i_pro))).eq.'C')   tmp_atoms(1) = imp_atoms(i_pro)
            if (trim(at%data%str(atom_type_index,imp_atoms(i_pro))).eq.'CP1') tmp_atoms(2) = imp_atoms(i_pro)
            if (trim(at%data%str(atom_type_index,imp_atoms(i_pro))).eq.'CP3') tmp_atoms(3) = imp_atoms(i_pro)
         enddo
         if (any(tmp_atoms(1:3).eq.0)) call system_abort('something has gone wrong. what is this if not proline?'// &
                     trim(at%data%str(atom_type_index,imp_atoms(1)))//imp_atoms(1)//'--'// &
                     trim(at%data%str(atom_type_index,imp_atoms(2)))//imp_atoms(2)//'--'// &
                     trim(at%data%str(atom_type_index,imp_atoms(3)))//imp_atoms(3)//'--'// &
                     trim(at%data%str(atom_type_index,imp_atoms(4)))//imp_atoms(4))
         imp_atoms(2:4) = tmp_atoms(1:3)
!         call print('Reordered Pro improper '// &
!                     trim(at%data%str(atom_type_index,imp_atoms(1)))//imp_atoms(1)//'--'// &
!                     trim(at%data%str(atom_type_index,imp_atoms(2)))//imp_atoms(2)//'--'// &
!                     trim(at%data%str(atom_type_index,imp_atoms(3)))//imp_atoms(3)//'--'// &
!                     trim(at%data%str(atom_type_index,imp_atoms(4)))//imp_atoms(4))
      endif
      call append(impropers,imp_atoms(1:4))
!      call print('Added backbone improper '// &
!                  trim(ElementName(at%Z(imp_atoms(1))))//imp_atoms(1)//'--'// &
!                  trim(ElementName(at%Z(imp_atoms(2))))//imp_atoms(2)//'--'// &
!                  trim(ElementName(at%Z(imp_atoms(3))))//imp_atoms(3)//'--'// &
!                  trim(ElementName(at%Z(imp_atoms(4))))//imp_atoms(4))
    enddo

   ! add intraresidual impropers from the given Table
    if (present(intrares_impropers)) then
       do i_impr = 1, intrares_impropers%N
          call append(impropers,intrares_impropers%int(1:4,i_impr))
!          call print('Added intraresidual improper '// &
!                      trim(ElementName(at%Z(intrares_impropers%int(1,i_impr))))//intrares_impropers%int(1,i_impr)//'--'// &
!                      trim(ElementName(at%Z(intrares_impropers%int(2,i_impr))))//intrares_impropers%int(2,i_impr)//'--'// &
!                      trim(ElementName(at%Z(intrares_impropers%int(3,i_impr))))//intrares_impropers%int(3,i_impr)//'--'// &
!                      trim(ElementName(at%Z(intrares_impropers%int(4,i_impr))))//intrares_impropers%int(4,i_impr))
       enddo
    else
       call print('WARNING!!! NO INTRARESIDUAL IMPROPERS USED!!!',verbosity=ERROR)
    endif

   ! final check
    if (any(impropers%int(1:4,1:impropers%N).le.0) .or. any(impropers%int(1:4,1:impropers%N).gt.at%N)) &
       call system_abort('create_improper_list: element(s) of impropers not within (0;at%N]')

    call system_timer('create_improper_list')

  end subroutine create_improper_list

  !% Returns the index of the first column of a property $prop$.
  !% Aborts if the property cannot be found in the atoms object.
  !
  function get_property(at,prop) result(prop_index)

    type(Atoms),      intent(in) :: at
    character(len=*), intent(in) :: prop
  
    integer,dimension(3) :: pos_indices
    integer              :: prop_index

    if (get_value(at%properties,trim(prop),pos_indices)) then
       prop_index = pos_indices(2)
    else
       call system_abort('get_property: No '//trim(prop)//' property assigned to the Atoms object!')
    end if

  end function get_property

  !% Calculates the position dependent charges for a silica molecule.
  !% See: Cole et al., J. Chem. Phys. 127 204704--204712 (2007)
  !
  subroutine create_pos_dep_charges(at,SiOH_list,charge) !,residue_names)

    type(Atoms), intent(in) :: at
    type(Table), intent(in) :: SiOH_list
    real(dp), allocatable, dimension(:), intent(out) :: charge
!    character(4), dimension(:), intent(in) :: residue_names

    real(dp) :: rij, fcq
    integer :: iatom, jatom
    integer :: atom_i, atom_j, jj
    logical, allocatable :: silica_mask(:)

    call system_timer('create_pos_dep_charges')

!    if (size(residue_names).ne.at%N) call system_abort('Residue names have a different size '//size(residue_names)//'then the number of atoms '//at%N)

    if (at%cutoff.lt.SILICON_2BODY_CUTOFF) call print('WARNING! The connection cutoff value '//at%cutoff// &
       ' is less than the silica_cutoff. The charges can be totally wrong. Check your connection object.',ERROR)

    allocate(charge(at%N))
    charge = 0._dp

    allocate(silica_mask(at%N)) !whether the atom is in the SiOH list
    silica_mask = .false.
    silica_mask(SiOH_list%int(1,1:SiOH_list%N)) = .true.

    do iatom = 1, SiOH_list%N
       atom_i = SiOH_list%int(1,iatom)
       do jatom = 1, atoms_n_neighbours(at,atom_i)
          atom_j = atoms_neighbour(at, atom_i, jatom)
         ! check if atom_j is in SiOH_list
          if (.not.silica_mask(atom_j)) then
!          if (.not.any(atom_j.eq.SiOH_list%int(1,1:SiOH_list%N))) then
!             ! if it's a water O/H, ignore it!
!             if (('TIP3'.eq.trim(residue_names(atom_j))) .or. &
!                 ('DOPA'.eq.trim(residue_names(atom_j))) .or. &
!                 ('TIP' .eq.trim(residue_names(atom_j)))) then
                 cycle
!             else
!                call print('Not found atom '//ElementName(at%Z(atom_j))//' '//atom_j//' in '//SiOH_list%int(1,1:SiOH_list%N),verbosity=ERROR)
!                call print('Unknown residue '//trim(residue_names(atom_j)))
!                call system_abort('Found a neighbour that is not part of the SiO2 residue')
!             endif
          endif
          if (atom_j.lt.atom_i) cycle

         ! get the distance between atom $i$ and $j$, apply (R_q + Delta_q) cutoff
          rij = distance_min_image(at, atom_i, atom_j)
          if (rij.gt.(Danny_R_q+Danny_Delta_q)) cycle

         !silanol endings (O--H and H--O)
          if (at%Z(atom_i).eq.1) then !H-O
             if (rij.lt.1.2_dp*(ElementCovRad(1)+ElementCovRad(8))) then
                charge(atom_i) = 0.2_dp
                charge(atom_j) = charge(atom_j) - 0.2_dp
             endif
          else
             if (at%Z(atom_j).eq.1) then !O-H
                if (rij.lt.1.2_dp*(ElementCovRad(1)+ElementCovRad(8))) then
                   charge(atom_i) = charge(atom_i) - 0.2_dp
                   charge(atom_j) = 0.2_dp
                endif
             else
                !Si--Si and O--O do not affect each other
                if (at%Z(atom_i).eq.at%Z(atom_j)) cycle
                !Si--O and O--Si
                fcq = calc_fc(rij)
                if (at%Z(atom_i).eq.14 .and. at%Z(atom_j).eq.8) then !Si-O
                   charge(atom_i) = charge(atom_i) + 0.4_dp * fcq
                   charge(atom_j) = charge(atom_j) - 0.4_dp * fcq
                else
                   if (at%Z(atom_i).eq.8 .and. at%Z(atom_j).eq.14) then !O-Si
                      charge(atom_i) = charge(atom_i) - 0.4_dp * fcq
                      charge(atom_j) = charge(atom_j) + 0.4_dp * fcq
                   else
                      call print(atom_i//trim(ElementName(at%Z(atom_i)))//'--'//atom_j//trim(ElementName(at%Z(atom_j)))//' of silica is none of H,O,Si!',verbosity=ERROR)
                      call system_abort('create_pos_dep_charge: unknown neighbour pair')
                   endif
                endif
             endif
          endif
       enddo
    enddo

    deallocate(silica_mask)

    call print('Calculated charge on atoms:',verbosity=ANAL)
    do jj = 1, at%N
       call print ('   atom '//jj//': '//charge(jj),verbosity=ANAL)
    enddo

    call system_timer('create_pos_dep_charges')

  end subroutine create_pos_dep_charges

  !% Cutoff function. Used by create_pos_dep_charges
  !
  function calc_fc(rij) result(fc)

    real(dp), intent(in) :: rij
    real(dp) :: fc

    real(dp) :: Danny_R, Danny_Delta

      Danny_R = Danny_R_q
      Danny_Delta = Danny_Delta_q

    if (rij.lt.(Danny_R-Danny_Delta)) then
       fc = 1._dp
    else
       if (rij.ge.(Danny_R+Danny_Delta)) then
          fc = 0._dp
       else
          fc = 1._dp - (rij - Danny_R + Danny_Delta) / (2._dp * Danny_Delta) + &
               sin(PI*(rij - Danny_R + Danny_Delta) / Danny_Delta) / (2._dp*PI)
       endif
    endif

  end function calc_fc

  !% Calculate silica cutoff for 2 body terms.
  !% The input $code$ = 1000 * Z(atom_i) + Z(atom_j)
  !
  function danny_cutoff(code) result(cutoff)

    integer, intent(in) :: code
    real(dp)            :: cutoff

    cutoff = 0._dp

    select case(code)
      case (14014)       !Si-Si
        cutoff = Danny_Si_Si_cutoff
      case (14008,8014) !Si-O
        cutoff = Danny_Si_Si_O_cutoff
      case (14001,1014)  !Si-H
        cutoff = Danny_Si_H_cutoff
      case (8008)       !O-O should never ever be less than 2.6
        cutoff = Danny_O_O_cutoff
      case (8001,1008)  !O-H
        cutoff = Danny_O_H_cutoff
      case (1001)        !H-H
        cutoff = Danny_H_H_cutoff
      case default
        call system_abort('danny_cutoff: unknown code '//code)
    end select

  end function danny_cutoff

  !% Remove bonds for metal ions - everything but H, C, N, O, Si, P, S, Cl
  !% Uses remove_bond in atoms_module.
  !
  subroutine delete_metal_connects(at)

    type(Atoms), intent(inout) :: at
    integer :: i, j

    do i = 1, at%N
       if (any(at%Z(i).eq.(/1,6,7,8,14,15,16,17/))) cycle
!       call print('found metal atom '//i)
       do j = 1, at%N
          if (i.ne.j) &
             call remove_bond(at%connect, i, j)
       enddo
    enddo

  end subroutine delete_metal_connects

end module topology_module
