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

module topology_module

  use system_module          
  use units_module           
  use extendable_str_module  
  use linearalgebra_module   
  use dictionary_module      
  use table_module           
  use periodictable_module   
  use connection_module           
  use atoms_types_module           
  use atoms_module           
  use clusters_module        
  use structures_module      


  implicit none
  private

public :: next_motif, find_motif_backbone
  !! private :: next_motif, write_psf_section, create_bond_list, write_pdb_file
  !! private :: create_angle_list, create_dihedral_list
  !! private :: create_improper_list
  !! private :: create_pos_dep_charges, calc_fc

  public  :: delete_metal_connects, &
             write_brookhaven_pdb_file, &
             write_cp2k_pdb_file, &
             write_psf_file, &
             write_psf_file_arb_pos, &
             create_residue_labels_arb_pos, &
             NONE_RUN, &
             QS_RUN, &
             MM_RUN, &
             QMMM_RUN_CORE, &
             QMMM_RUN_EXTENDED, &
             find_water_monomer, find_A2_monomer, find_AB_monomer, &
	     find_molecule_ids, &
             find_general_monomer, &
             find_monomer_pairs, &
             find_monomer_triplets, &
             calc_mean_pos


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
  real(dp), parameter, public :: TITANIA_2BODY_CUTOFF = 7.5_dp !Si-O, O-O 2body interaction
  real(dp), parameter, public :: SILICON_3BODY_CUTOFF = 3.8_dp !Si-Si-Si, Si-Si-O 3body cutoff
  real(dp), parameter, public :: SILICA_3BODY_CUTOFF = 3.6_dp !Si-O-Si, O-Si-O, Si-O-H 3body cutoff
  real(dp), parameter, public :: TITANIA_3BODY_CUTOFF = 7.5_dp !Si-O, O-O 2body interaction
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
  !% if EVB is being used, must pass exact connectivity object and set evb_nneighb_only to false
  !%
  subroutine create_residue_labels_arb_pos(at,do_CHARMM,intrares_impropers,find_silica_residue,pos_field_for_connectivity, &
       form_bond,break_bond, silica_pos_dep_charges, silica_charge_transfer, have_titania_potential, &
       find_molecules, evb_nneighb_only, remove_Si_H_silica_bonds, remove_Ti_H_titania_bonds, error)

    type(Atoms),           intent(inout),target :: at
    logical,     optional, intent(in)    :: do_CHARMM
    type(Table), optional, intent(out)   :: intrares_impropers
    logical,     optional, intent(in)    :: find_silica_residue
    character(len=*), optional, intent(in) :: pos_field_for_connectivity
    integer, optional, intent(in) :: form_bond(2), break_bond(2)
    logical,     optional, intent(in)    :: silica_pos_dep_charges
    real(dp), intent(in), optional :: silica_charge_transfer
    logical, intent(in), optional :: have_titania_potential
    logical, intent(in), optional :: find_molecules
    logical, intent(in), optional :: evb_nneighb_only
    logical, optional, intent(in) :: remove_Si_H_silica_bonds, remove_Ti_H_titania_bonds
    integer, optional, intent(out) :: error

    real(dp), pointer :: use_pos(:,:)
    type(Connection) :: t_connect
    type(Atoms) :: at_copy
    logical :: do_find_silica_residue, do_have_titania_potential
    logical :: use_pos_is_pos

    logical :: bond_exists, do_nneighb_only
    integer :: shift(3)
    real(dp) :: form_bond_dist

    integer :: ji, j

    INIT_ERROR(error)

    do_nneighb_only = optional_default(.true., evb_nneighb_only)

    ! save a copy
    at_copy = at

    use_pos_is_pos = .false.
    ! find desired position field (pos, avgpos, whatever)
    if (present(pos_field_for_connectivity)) then
      if (.not. assign_pointer(at, trim(pos_field_for_connectivity), use_pos)) then
	RAISE_ERROR("calc_topology can't find pos field '"//trim(pos_field_for_connectivity)//"'", error)
      endif
      if (trim(pos_field_for_connectivity) == 'pos') use_pos_is_pos = .true.
    else if (.not. assign_pointer(at, 'avgpos', use_pos)) then
      call print("WARNING: calc_topology can't find default pos field 'avgpos', trying to use pos instead")
      if (.not. assign_pointer(at, 'pos', use_pos)) then
	RAISE_ERROR("calc_topology can't find avgpos or pos fields",error)
      endif
      use_pos_is_pos = .true.
    endif

    do_find_silica_residue = optional_default(.false.,find_silica_residue)
    do_have_titania_potential = optional_default(.false.,have_titania_potential)
  
    ! copy desired pos to pos, and new connectivity
    !NB don't do if use_pos => pos
    if (.not. use_pos_is_pos) at_copy%pos = use_pos
    if (do_find_silica_residue) then
       call set_cutoff(at_copy,SILICA_2BODY_CUTOFF)
       call calc_connect(at_copy, alt_connect=t_connect)       
    elseif (do_have_titania_potential) then
       call set_cutoff(at_copy,1.3*bond_length(22, 8))
       call calc_connect(at_copy, alt_connect=t_connect)
    else
       ! use hysteretic connect to get nearest neighbour cutoff
       ! will use default cutoff, which is the same as heuristics_nneighb_only=.true.
       call calc_connect_hysteretic(at, DEFAULT_NNEIGHTOL, DEFAULT_NNEIGHTOL, alt_connect=t_connect)
    endif


    call break_form_bonds(at, t_connect, form_bond, break_bond, error=error)
    PASS_ERROR(error)

    ! now create labels using this connectivity object
    call create_residue_labels_internal(at,do_CHARMM,intrares_impropers,heuristics_nneighb_only=do_nneighb_only,alt_connect=t_connect,&
	 find_silica_residue=do_find_silica_residue, silica_pos_dep_charges=silica_pos_dep_charges, &
	 silica_charge_transfer=silica_charge_transfer, have_titania_potential=have_titania_potential, &
	 find_molecules=find_molecules, remove_Si_H_silica_bonds=remove_Si_H_silica_bonds, &
         remove_Ti_H_titania_bonds=remove_Ti_H_titania_bonds, error=error)
    PASS_ERROR(error)
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
  subroutine create_residue_labels_internal(at,do_CHARMM,intrares_impropers, heuristics_nneighb_only,alt_connect, &
       find_silica_residue, silica_pos_dep_charges, silica_charge_transfer, have_titania_potential, &
       find_molecules, remove_Si_H_silica_bonds, remove_Ti_H_titania_bonds, error) !, hysteretic_neighbours)

    type(Atoms),           intent(inout),target :: at
    logical,     optional, intent(in)    :: do_CHARMM
    type(Table), optional, intent(out)   :: intrares_impropers
    logical, intent(in), optional :: heuristics_nneighb_only
    type(Connection), intent(in), optional, target :: alt_connect
    logical,     optional, intent(in)    :: find_silica_residue, silica_pos_dep_charges
    real(dp), optional, intent(in) :: silica_charge_transfer
    logical,     optional, intent(in)    :: have_titania_potential
    logical, optional, intent(in) :: find_molecules
    integer, optional, intent(out) :: error
    logical, optional, intent(in) :: remove_Si_H_silica_bonds, remove_Ti_H_titania_bonds
!    logical, optional, intent(in) :: hysteretic_neighbours


    character(*), parameter  :: me = 'create_residue_labels_pos_internal: '

    type(Inoutput)                       :: lib
    character(4)                         :: cha_res_name(MAX_KNOWN_RESIDUES), Cres_name
    character(3)                         :: pdb_res_name(MAX_KNOWN_RESIDUES), pres_name
    integer                              :: residue_number(at%N)
    character(4)                         :: atom_name(at%N), atom_name_PDB(at%N)
    real(dp)                             :: atom_charge(at%N)
    integer                              :: atom_subgroup(at%N)
    type(Table)                          :: residue_type, list
    logical                              :: unidentified(at%N)
    integer, allocatable, dimension(:,:) :: motif
    integer                              :: i, m, n, nres
    character(4), allocatable            :: at_names(:), at_names_PDB(:)
    real(dp),     allocatable            :: at_charges(:)
    integer, allocatable                 :: at_subgroups(:)
    logical                              :: my_do_charmm
    character, dimension(:,:), pointer   :: atom_type, atom_type_PDB, atom_res_name, atom_mol_name
    integer, dimension(:), pointer       :: atom_res_number, atom_res_type, motif_atom_num, atom_subgroup_number
    real(dp), dimension(:), pointer      :: atom_charge_ptr
    logical                              :: ex
    type(extendable_str)                 :: residue_library
    integer                              :: i_impr, n_impr
    integer, allocatable                 :: imp_atoms(:,:)
    real(dp)                             :: mol_charge_sum
    logical                              :: found_residues
    type(Table)                          :: atom_Si, atom_SiO, SiOH_list, atom_Ti, atom_TiO, TiOH_list
    real(dp), dimension(:), allocatable  :: charge
integer :: j,atom_i, ji
!logical :: silanol
    type(Table) :: bondH,bondSi, bondTi
    integer :: bond_H,bond_Si, bond_Ti
!    integer                             :: qm_flag_index, pos_indices(3)
!    logical                             :: do_qmmm
    type(Connection), pointer :: use_connect
!    logical :: use_hysteretic_neighbours
    type(Table) :: O_atom, O_neighb
    integer :: hydrogen
    logical :: find_silica, titania_potential, do_silica_pos_dep_charges, do_find_molecules
    logical :: do_remove_Si_H_silica_bonds, do_remove_Ti_H_titania_bonds
    real(dp) :: do_silica_charge_transfer

    integer, pointer :: mol_id(:)
    type(allocatable_array_pointers), allocatable :: molecules(:)
    integer :: i_motif_at

    INIT_ERROR(error)

    call system_timer('create_residue_labels_pos_internal')

    find_silica = optional_default(.false.,find_silica_residue)
    do_silica_pos_dep_charges = optional_default(.true., silica_pos_dep_charges)
    do_silica_charge_transfer = optional_default(2.4_dp, silica_charge_transfer)
    titania_potential = optional_default(.false.,have_titania_potential)
    do_find_molecules = optional_default(.true., find_molecules)
    do_remove_Si_H_silica_bonds = optional_default(.true., remove_Si_H_silica_bonds)
    do_remove_Ti_H_titania_bonds = optional_default(.true., remove_Ti_H_titania_bonds)

    if (present(alt_connect)) then
      use_connect => alt_connect
    else
      use_connect => at%connect
    endif
    if (.not.use_connect%initialised) then
       RAISE_ERROR(me//'No connectivity data present in atoms structure', error)
    endif
!    use_hysteretic_neighbours = optional_default(.false.,hysteretic_neighbours)

    my_do_charmm = optional_default(.true.,do_CHARMM)

    call initialise(residue_library)
    call print_title('Creating CHARMM format')
    call get_param_value(at, 'Library', residue_library, error=error)
    PASS_ERROR_WITH_INFO('create_residue_labels_pos_internal: no residue library specified, but topology generation requested', error)
    call print('Library: '//string(residue_library))

    !Open the residue library
    !call print('Opening library...')
    call initialise(lib,string(residue_library),action=INPUT)

    !Set all atoms as initially unidentified
    unidentified = .true.

    ! add property for find_motif() to save the index of each atom in the motif
    call add_property(at,'motif_atom_num',0, ptr=motif_atom_num)

    !Read each of the residue motifs from the library
    n = 0
    nres = 0
    mol_charge_sum = 0._dp
    found_residues = .false.
    if (present(intrares_impropers)) call initialise(intrares_impropers,4,0,0,0,0)
    call allocate(residue_type,1,0,0,0,1000)
    call print('Identifying atoms...')

!!!!!!!!!!!  TITANIA POTENTIAL  !!!!!!!!!!
    if (titania_potential) then
  ! TIO residue if Ti atom is present in the atoms structure 
       if (any(at%Z(1:at%N).eq.22)) then
          call print('|-Looking for TIO residue, not from the library...')
          call print('| |-Found... will be treated as 1 molecule, 1 residue...')
          !all this bulk will be 1 residue
          n = n + 1
          cha_res_name(n) = 'TIO2'
          pdb_res_name(n) = '' !not used
          call append(residue_type,(/n/))
          nres = nres + 1

          !Add Ti atoms
          call initialise(atom_Ti,4,0,0,0,0)
          do i = 1,at%N
             if (at%Z(i).eq.22) then
                call append(atom_Ti,(/i,0,0,0/))
             endif
          enddo
          call print(atom_Ti%N//' Ti atoms found in total')
          !Add O atoms
          call bfs_step(at,atom_Ti,atom_TiO,nneighb_only=.true.,min_images_only=.true.,alt_connect=use_connect)
          call print(atom_TiO%N//' O atoms found in total')
!          if (any(at%Z(atom_TiO%int(1,1:atom_TiO%N)).eq.1)) call system_abort('Ti-H bond')
          if (do_remove_Ti_H_titania_bonds) then
             do i=1,atom_TiO%N !check Hs bonded to Ti. There shouldn't be any,removing the bond.
                 if (at%Z(atom_TiO%int(1,i)).eq.1) then
                    call print('WARNING! Ti and H are very close',verbosity=PRINT_ALWAYS)
                    bond_H = atom_TiO%int(1,i)
                    call initialise(bondH,4,0,0,0,0)
                    call append(bondH,(/bond_H,0,0,0/))
                    call bfs_step(at,bondH,bondTi,nneighb_only=.true.,min_images_only=.true.,alt_connect=use_connect)
                    do j = 1,bondTi%N
                       if (at%Z(bondTi%int(1,i)).eq.22) then
                          bond_Ti = bondTi%int(1,i)
                          call print('WARNING! Remove Ti '//bond_Ti//' and H '//bond_H//' bond ('//distance_min_image(at,bond_H,bond_Ti)//')',verbosity=PRINT_ALWAYS)
                          call remove_bond(use_connect,bond_H,bond_Ti)
                       endif
                    enddo
                 endif
             enddo
          endif

      !Add H atoms -- if .not.remove_Ti_H_bonds, we might include whole water molecules at this stage, adding the remaining -OH.
          call add_cut_hydrogens(at,atom_TiO,heuristics_nneighb_only=.true.,alt_connect=use_connect)
          call print(atom_TiO%N//' O/H atoms found in total')

!check if none of these atoms are identified yet
          if (any(.not.unidentified(atom_Ti%int(1,1:atom_Ti%N)))) then
             RAISE_ERROR('already identified atoms found again.', error)
          endif
          if (any(.not.unidentified(atom_TiO%int(1,1:atom_TiO%N)))) then
!             call system_abort('already identified atoms found again.')
             do i = 1,atom_TiO%N,-1
                if (.not.unidentified(atom_TiO%int(1,i))) then
                   call print('delete from TiO2 list already identified atom '//atom_TiO%int(1,1:atom_TiO%N))
                   call delete(atom_TiO,i)
                endif
             enddo
          endif
          unidentified(atom_Ti%int(1,1:atom_Ti%N)) = .false.
          unidentified(atom_TiO%int(1,1:atom_TiO%N)) = .false.

          !add atom, residue and molecule names
          !TIO
          do i = 1, atom_Ti%N                              !atom_Ti  only has Ti atoms
             atom_i = atom_Ti%int(1,i)
             atom_name(atom_i) = 'TIO'
             atom_name_PDB(atom_i) = 'TIO'
          enddo

          !OTB, OTI & HTI
          do i = 1, atom_TiO%N                             !atom_TiO only has O,H atoms
             atom_i = atom_TiO%int(1,i)
             !OTB & OTI
             if (at%Z(atom_i).eq.8) then
                call initialise(O_neighb,4,0,0,0)
                call initialise(O_atom,4,0,0,0)
                call append(O_atom,(/atom_i,0,0,0/))
                call bfs_step(at,O_atom,O_neighb,nneighb_only=.true.,min_images_only=.true.,alt_connect=use_connect)
                !check nearest neighbour number = 3
                if (O_neighb%N.ne.3) then
                   call print('WARNING! titania O '//atom_i//'has '//O_neighb%N//'/=3 nearest neighbours',PRINT_ALWAYS)
                   call print('neighbours: '//O_neighb%int(1,1:O_neighb%N))
                endif
                !check if it has a H nearest neighbour
                hydrogen = find_in_array(at%Z(O_neighb%int(1,1:O_neighb%N)),1)
                if (hydrogen.ne.0) then
                   atom_name(atom_i) = 'OTI' !silanol O
                   atom_name_PDB(atom_i) = 'OTI' !silanol O
!                   call print('Found OH silanol oxygen.'//atom_TiO%int(1,i)//' hydrogen: '//O_neighb%int(1,hydrogen))
                   !check if it has only 1 H nearest neighbour
                   if (hydrogen.lt.O_neighb%N) then
                      if(find_in_array(at%Z(O_neighb%int(1,hydrogen+1:O_neighb%N)),1).gt.0) then
                         RAISE_ERROR('More than 1 H neighbours of O '//atom_i, error)
                      endif
                   endif
                else
                   atom_name(atom_i) = 'OTB' !bridging O
                   atom_name_PDB(atom_i) = 'OTB' !bridging O
!                   call print('Found OB bridging oxygen.'//atom_TiO%int(1,i))
                endif
                call finalise(O_atom)
                call finalise(O_neighb)
             !HSI
             elseif (at%Z(atom_TiO%int(1,i)).eq.1) then
                atom_name(atom_TiO%int(1,i)) = 'HTI'
                atom_name_PDB(atom_TiO%int(1,i)) = 'HTI'
             else
                RAISE_ERROR('Non O/H atom '//atom_i//'!?', error)
             endif
          enddo
          
          !Add all the titanium atoms together
          call initialise(TiOH_list,4,0,0,0,0)
          call append (TiOH_list,atom_Ti)
          call append (TiOH_list,atom_TiO)

          !Residue numbers
          residue_number(TiOH_list%int(1,1:TiOH_list%N)) = nres

          !Charges
          !call create_pos_dep_charges(at,TiOH_list,charge) !,residue_names=cha_res_name(residue_type%int(1,residue_number(1:at%N))))
          atom_charge = 1.2_dp
          do i=1,at%n
             if(at%Z(i).eq.8) atom_charge(i) = -0.6_dp 
             call print('   '//i//'   '//atom_charge(i),verbosity=PRINT_ANAL)
          enddo
          call print("overall titania charge: "//sum(atom_charge(TiOH_list%int(1,1:TiOH_list%N))))

          call finalise(atom_Ti)
          call finalise(atom_TiO)
          call finalise(TiOH_list)

       else
          call print('WARNING! have_titania_potential is true, but found no Ti atoms in the atoms object!',PRINT_ALWAYS)
       endif

    endif
!!!!!!!!!!!!!! END TITANIA POTENTIAL !!!!!!!!!!!!!!!!
    

!!!!!!!!!!!!!!! DANNY POTENTIAL !!!!!!!!!!!!!!!!
    if (find_silica) then

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
          if (do_remove_Si_H_silica_bonds) then
             call initialise(bondH,4,0,0,0)
             call initialise(bondSi,4,0,0,0)

             do i=1,atom_SiO%N !check Hs bonded to Si. There shouldn't be any,removing the bond.
                 if (at%Z(atom_SiO%int(1,i)).eq.1) then
                    call print('WARNING! Si and H are very close',verbosity=PRINT_ALWAYS)
                    bond_H = atom_SiO%int(1,i)
                    call wipe(bondH)
                    call wipe(bondSi)
                    call append(bondH,(/bond_H,0,0,0/))
                    call bfs_step(at,bondH,bondSi,nneighb_only=.true.,min_images_only=.true.,alt_connect=use_connect)
                    do j = 1,bondSi%N
                       if (at%Z(bondSi%int(1,j)).eq.14) then
                          bond_Si = bondSi%int(1,j)
                          call print('WARNING! Remove Si '//bond_Si//' and H '//bond_H//' bond ('//distance_min_image(at,bond_H,bond_Si)//')',verbosity=PRINT_ALWAYS)
                          call remove_bond(use_connect,bond_H,bond_Si)
                       endif
                    enddo
                 endif
!                call print('atom_SiO has '//at%Z(atom_SiO%int(1,i)))
             enddo
          endif
          !Add H atoms -- if .not.remove_Si_H_bonds, we might include whole water molecules at this stage, adding the remaining -OH.
      !    call bfs_step(at,atom_SiO,atom_SIOH,nneighb_only=.true.,min_images_only=.true.)
          call add_cut_hydrogens(at,atom_SiO,heuristics_nneighb_only=.true.,alt_connect=use_connect)
          call print(atom_SiO%N//' O/H atoms found in total')

          !check if none of these atoms are identified yet
          if (any(.not.unidentified(atom_Si%int(1,1:atom_Si%N)))) then
             RAISE_ERROR('already identified atoms found again.', error)
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
                   call print('WARNING! silica O '//atom_i//'has '//O_neighb%N//'/=2 nearest neighbours',PRINT_ALWAYS)
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
                      if(find_in_array(at%Z(O_neighb%int(1,hydrogen+1:O_neighb%N)),1).gt.0) then
                         RAISE_ERROR('More than 1 H neighbours of O '//atom_i, error)
		      endif
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
                RAISE_ERROR('Non O/H atom '//atom_i//'!?', error)
             endif
          enddo

          !Add all the silica atoms together
          call initialise(SiOH_list,4,0,0,0,0)
          call append (SiOH_list,atom_Si)
          call append (SiOH_list,atom_SiO)

          !Residue numbers
          residue_number(SiOH_list%int(1,1:SiOH_list%N)) = nres

          !Charges
          if (do_silica_pos_dep_charges) then
             call create_pos_dep_charges(at,SiOH_list,charge) !,residue_names=cha_res_name(residue_type%int(1,residue_number(1:at%N))))
             atom_charge(SiOH_list%int(1,1:SiOH_list%N)) = 0._dp
             atom_charge(SiOH_list%int(1,1:SiOH_list%N)) = charge(SiOH_list%int(1,1:SiOH_list%N))
             call print("overall silica charge: "//sum(atom_charge(SiOH_list%int(1,1:SiOH_list%N))))
             call print('Atomic charges: ',PRINT_ANAL)
             call print('   ATOM     CHARGE',PRINT_ANAL)
             do i=1,at%N
                call print('   '//i//'   '//atom_charge(i),verbosity=PRINT_ANAL)
             enddo
          else
             where (at%z == 1) atom_charge = do_silica_charge_transfer/4.0_dp
             where (at%z == 8) atom_charge = -do_silica_charge_transfer/2.0_dp
             where (at%z ==14) atom_charge = do_silica_charge_transfer
          end if

          call finalise(atom_Si)
          call finalise(atom_SiO)
          call finalise(SiOH_list)

       else
          call print('WARNING! find_silica_residue is true, but found no silicon atoms in the atoms object!',PRINT_ALWAYS)
       endif

    endif
!!!!!!!!!!!!!!! END DANNY POTENTIAL !!!!!!!!!!!!!!!!

    do 

       ! Pull the next residue template from the library
       if (my_do_charmm) then
          call next_motif(lib,cres_name,pres_name,motif,atom_names=at_names,atom_charges=at_charges,atom_names_PDB=at_names_PDB, n_impr=n_impr,imp_atoms=imp_atoms,do_CHARMM=.true.,at_subgroups=at_subgroups)
       else
          call next_motif(lib,cres_name,pres_name,motif,atom_names=at_names,atom_names_PDB=at_names_PDB, do_CHARMM=.false.,at_subgroups=at_subgroups)
       endif

       if (cres_name=='NONE') then
          found_residues = .true.
          exit
       endif

! call print("motif")
!   call print("cres_name " // trim(cres_name)// " pres_name "//trim(pres_name) // " n_impr " // n_impr)
! do m=1, size(motif,1)
!   call print("motif  " // motif(m,1:7) // " name " // trim(at_names(m)) // " charge " // at_charges(m) // " PDB name " // trim(at_names_PDB(m)) // &
!    " subgroup " // at_subgroups(m))
! end do

       ! Store its CHARMM (3 char) and PDB (3 char) names
       ! e.g. cha/pdb_res_name(5) corresponds to the 5th residue found in the library
       n = n + 1
       cha_res_name(n) = cres_name
       pdb_res_name(n) = pres_name

       ! Search the atom structure for this residue
       call print('|-Looking for '//cres_name//'...',verbosity=PRINT_ANAL)
       call find_motif(at,motif,list,mask=unidentified,nneighb_only=heuristics_nneighb_only,alt_connect=use_connect) !,hysteretic_neighbours=use_hysteretic_neighbours)

       if (list%N > 0) then
          
          call print('| |-Found '//list%N//' occurrences of '//cres_name//' with charge '//(sum(at_charges(1:size(at_charges)))))
          mol_charge_sum = mol_charge_sum + list%N * sum(at_charges(1:size(at_charges)))

          ! Loop over all found instances of the residue

          do m = 1, list%N

             !Mark the newly identified atoms
!1 row is 1 residue
             unidentified(list%int(:,m)) = .false.

	     do i_motif_at=1, size(list%int,1)
	       motif_atom_num(list%int(i_motif_at,m)) = i_motif_at
	     end do

             !Store the residue info
             nres = nres + 1
             call append(residue_type,(/n/))
             ! if residue_type%int(1,2) == 3 then residue no. 2 matches the 3rd residue in the library
             residue_number(list%int(:,m)) = nres
                  !e.g. if residue_number(i) = j        j-th residue in the atoms object (1 to 8000 in case of 8000 H2O)
                  !        residue_type(1,j)   = k        k-th residue in the library file, in order
                  !        cha_res_name(k)   = 'ALA'    name of the k-th residue in the library file
                  !then atom 'i' is in a residue 'ALA'
	     atom_subgroup(list%int(:,m)) = at_subgroups
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

       call print('|',verbosity=PRINT_ANAL)

    end do

   ! check if residue library is empty
    if (.not.found_residues) then
      RAISE_ERROR('Residue library '//trim(lib%filename)//' does not contain any residues!', error)
    endif

    call print('Finished.')
    call print(nres//' residues found in total')

    if (any(unidentified)) then
       call print(count(unidentified)//' unidentified atoms',verbosity=PRINT_ALWAYS)
       call print(find(unidentified))
       do i=1,at%N
	  if (unidentified(i)) then
	    call print(ElementName(at%Z(i))//' atom '//i//' has avgpos: '//round(at%pos(1,i),5)//&
	      ' '//round(at%pos(2,i),5)//' '//round(at%pos(3,i),5),verbosity=PRINT_ALWAYS)
	    call print(ElementName(at%Z(i))//' atom '//i//' has number of neighbours: '//n_neighbours(at,i,alt_connect=alt_connect),verbosity=PRINT_ALWAYS)
	    do ji=1, n_neighbours(at, i, alt_connect=alt_connect)
	      j = neighbour(at, i, ji, alt_connect=alt_connect)
	      call print("  neighbour " // j // " is of type " // ElementName(at%Z(j)), verbosity=PRINT_ALWAYS)
	    end do
	  endif
       enddo

      call print(at%connect)
       ! THIS IS WHERE THE CALCULATION OF NEW PARAMETERS SHOULD GO
      RAISE_ERROR('create_residue_labels_pos_internal: Unidentified atoms', error)

    else
       call print('All atoms identified')
       call print('Total charge of the molecule: '//round(mol_charge_sum,5))
    end if

   ! add data to store CHARMM topology
    call add_property(at,'atom_type',repeat(' ',TABLE_STRING_LENGTH), ptr=atom_type)
    call add_property(at,'atom_type_PDB',repeat(' ',TABLE_STRING_LENGTH), ptr=atom_type_PDB)
    call add_property(at,'atom_res_name',repeat(' ',TABLE_STRING_LENGTH), ptr=atom_res_name)
    call add_property(at,'atom_mol_name',repeat(' ',TABLE_STRING_LENGTH), ptr=atom_mol_name)
    call add_property(at,'atom_res_number',0, ptr=atom_res_number)
    call add_property(at,'atom_subgroup_number',0, ptr=atom_subgroup_number)
    call add_property(at,'atom_res_type',0, ptr=atom_res_type)
    call add_property(at,'atom_charge',0._dp, ptr=atom_charge_ptr)

    atom_type(1,1:at%N) = 'X'
    atom_type_PDB(1,1:at%N) = 'X'
    atom_res_name(1,1:at%N) = 'X'
    atom_mol_name(1,1:at%N) = 'X'
    atom_res_number(1:at%N) = 0
    atom_subgroup_number(1:at%N) = 0
    atom_res_type(1:at%N) = 0
    atom_charge_ptr(1:at%N) = 0._dp

    do i=1, at%N
       atom_res_name(:,i) = pad(cha_res_name(residue_type%int(1,residue_number(i))), TABLE_STRING_LENGTH)
       atom_res_number(i) = residue_number(i)
       atom_subgroup_number(i) = atom_subgroup(i)
       atom_res_type(i) = residue_type%int(1,residue_number(i))
       atom_type(:,i) = pad(adjustl(atom_name(i)), TABLE_STRING_LENGTH)
       atom_type_PDB(:,i) = pad(adjustl(atom_name_PDB(i)), TABLE_STRING_LENGTH)
    end do

    if (my_do_charmm) then
       if (do_find_molecules) then
          allocate(molecules(at%N))
          call find_molecule_ids(at,molecules,heuristics_nneighb_only=heuristics_nneighb_only,alt_connect=alt_connect)
          do i=1, size(molecules)
             if (allocated(molecules(i)%i_a)) then
                ! special case for silica molecule
                if (find_silica .and. count(at%Z(molecules(i)%i_a) == 14) /=0) then
                   atom_mol_name(:,molecules(i)%i_a) = atom_res_name(:,molecules(i)%i_a)
                   ! special case for single atoms
                elseif (size(molecules(i)%i_a) == 1) then
                   atom_mol_name(:,molecules(i)%i_a) = atom_res_name(:,molecules(i)%i_a)
                   ! special case for H2O
                else if (size(molecules(i)%i_a) == 3) then
                   if (count(at%Z(molecules(i)%i_a) == 8) == 1 .and. &
                        count(at%Z(molecules(i)%i_a) == 1) == 2) then
                      atom_mol_name(:,molecules(i)%i_a) = atom_res_name(:,molecules(i)%i_a)
                   else ! default
                      call print("Found molecule containing "//size(molecules(i)%i_a)//" atoms and not water, single atom or silica")
                      do j=1,size(molecules(i)%i_a)
                         atom_mol_name(:,molecules(i)%i_a(j)) = pad("M"//i, TABLE_STRING_LENGTH)
                      end do
                   endif
                else ! default
                   call print("Found molecule containing "//size(molecules(i)%i_a)//" atoms and not water, single atom or silica")
                   do j=1,size(molecules(i)%i_a)
                      atom_mol_name(:,molecules(i)%i_a(j)) = pad("M"//i, TABLE_STRING_LENGTH)
                   end do
                endif
             end if ! allocated(molecules)
          end do ! i=1,size(molecules)
          deallocate(molecules)
       else 
          ! find_molecules is F, assume we have just one molecule
          ! we set molecule name to residue name of first atom

          call add_property(at,'mol_id',0,overwrite=.true.)
          if (.not. assign_pointer(at,'mol_id',mol_id)) &
               call system_abort("create_residue_labels_internal can't assign mol_id")
          mol_id(:) = 1
          do i=1, at%N
             atom_mol_name(:,i) = atom_res_name(:,1)
          end do
       end if
       atom_charge_ptr(1:at%N) = atom_charge(1:at%N)
    endif

    if (any(atom_res_number(1:at%N).le.0)) then
       RAISE_ERROR('create_residue_labels_pos_internal: atom_res_number is not >0 for every atom', error)
    endif
    if (any(atom_type(1,1:at%N).eq.'X')) then
       RAISE_ERROR('create_residue_labels_pos_internal: atom_type is not saved for at least one atom', error)
    endif
    if (any(atom_res_name(1,1:at%N).eq. 'X')) then
       RAISE_ERROR('create_residue_labels_pos_internal: atom_res_name is not saved for at least one atom', error)
    endif

    !Free up allocations
    call finalise(residue_type)
    call finalise(list)
    if (allocated(motif)) deallocate(motif)

    !Close the library
    call finalise(lib)

    call system_timer('create_residue_labels_pos_internal')

  end subroutine create_residue_labels_internal

  subroutine find_molecule_ids(at,molecules,heuristics_nneighb_only,alt_connect)
    type(Atoms), intent(inout) :: at
    type(allocatable_array_pointers), optional :: molecules(:)
    logical, intent(in), optional :: heuristics_nneighb_only
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
	   call bfs_step(at,cur_molec,next_atoms,nneighb_only = heuristics_nneighb_only, min_images_only = .true., alt_connect=alt_connect)
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
  subroutine next_motif(library,res_name,pdb_name,motif,atom_names,atom_charges,atom_names_PDB, n_impr,imp_atoms,do_CHARMM,at_subgroups,header_line)
    
    type(Inoutput),                   intent(in)  :: library
    character(4),                     intent(out) :: res_name
    character(3),                     intent(out) :: pdb_name
    integer,             allocatable, intent(out) :: motif(:,:)
    character(4),        allocatable, intent(out) :: atom_names(:), atom_names_PDB(:)
    real(dp), optional,  allocatable, intent(out) :: atom_charges(:)
    logical,  optional,               intent(in)  :: do_CHARMM
    integer,  optional,  allocatable, intent(out) :: imp_atoms(:,:)
    integer,  optional,               intent(out) :: n_impr
    integer, optional, allocatable, intent(out)   :: at_subgroups(:)
    character(len=*), optional, intent(out)       :: header_line

    character(20), dimension(11) :: fields
    integer                      :: status, num_fields, data(7), i, n_at, max_num_fields
    type(Table)                  :: motif_table
    integer :: tmp_at_subgroups(MAX_ATOMS_PER_RES)
    character(4)                 :: tmp_at_names(MAX_ATOMS_PER_RES), tmp_at_names_PDB(MAX_ATOMS_PER_RES)
    real(dp)                     :: tmp_at_charges(MAX_ATOMS_PER_RES),check_charge
    logical                      :: my_do_charmm
    character(len=STRING_LENGTH) :: line
   ! for improper generation
    integer                      :: imp_fields, tmp_imp_atoms(4,MAX_IMPROPERS_PER_RES)

    my_do_charmm = .true.
    if (present(do_CHARMM)) my_do_charmm = do_CHARMM

    if (my_do_charmm) then
       max_num_fields = 11
       imp_fields = 5
    else !do AMBER
       max_num_fields = 9
    endif

    status = 0

    if (present(header_line)) header_line =''
    do while(status==0)
       line = read_line(library,status)
       if (present(header_line)) then
	 if (len_trim(header_line) > 0) header_line=trim(header_line)//quip_new_line
	 header_line=trim(header_line)//trim(line)
       endif
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
          data(i) = string_to_int(fields(i+1))
       end do
       call append(motif_table,data)
       n_at = n_at + 1
       tmp_at_subgroups(n_at) = string_to_int(fields(1))
       tmp_at_names(n_at) = fields(9)
       if (my_do_charmm) then
          tmp_at_charges(n_at) = string_to_real(fields(10))
          check_charge=check_charge+tmp_at_charges(n_at)
          if(num_fields == 11) then
             tmp_at_names_PDB(n_at) = fields(11)
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
!       call print('WARNING next_motif: Charge of '//res_name//' residue is :'//round(check_charge,4),verbosity=PRINT_ALWAYS)
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

    if (present(at_subgroups)) then
      if (allocated(at_subgroups)) then
	if (size(at_subgroups) /= n_at) deallocate(at_subgroups)
      endif
      if (.not. allocated(at_subgroups)) allocate(at_subgroups(n_at))
      at_subgroups = tmp_at_subgroups(1:n_at)
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
    character, dimension(:,:), pointer :: atom_type, atom_res_name, atom_mol_name
    integer,   dimension(:),   pointer :: atom_res_number
    real(dp),  dimension(:),   pointer :: atom_charge
    character(len=STRING_LENGTH)       :: my_run_type_string
    real(dp)                 :: cell_lengths(3), cell_angles(3)

    external :: lattice_xyz_to_abc

    call system_timer('write_pdb_file')
    my_run_type_string = optional_default('',run_type_string)

    call initialise(pdb,trim(pdb_file),action=OUTPUT)
    call print('   PDB file: '//trim(pdb%filename))
!    call print('REMARK'//at%N,file=pdb)
!lattice information could be added in a line like this:
!CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1          
    call get_lattice_params(at%lattice, cell_lengths(1), cell_lengths(2), cell_lengths(3), &
					cell_angles(1), cell_angles(2), cell_angles(3))
    sor=''
    write(sor, '(a6,3f9.3,3f7.2,a16)') 'CRYST1', cell_lengths(:), DEGREES_PER_RADIAN*cell_angles(:), ' P 1           1'
    call print(sor, file=pdb)
    
!    if ((trim(my_run_type_string).eq.'QMMM_CORE') .or. &
!        (trim(my_run_type_string).eq.'QMMM_EXTENDED')) then
!       qm_flag_index = get_property(at,'cluster_mark')
!       if (trim(my_run_type_string).eq.'QMMM_EXTENDED') run_type=QMMM_RUN_EXTENDED
!       if (trim(my_run_type_string).eq.'QMMM_CORE') run_type=QMMM_RUN_CORE
!    endif

    if (.not. assign_pointer(at, 'atom_type_PDB', atom_type)) &
         call system_abort('Cannot assign pointer to "atom_type_PDB" property.')
    if (.not. assign_pointer(at, 'atom_res_name', atom_res_name)) &
         call system_abort('Cannot assign pointer to "atom_res_name" property.')
    if (.not. assign_pointer(at, 'atom_mol_name', atom_mol_name)) &
         call system_abort('Cannot assign pointer to "atom_mol_name" property.')
    if (.not. assign_pointer(at, 'atom_res_number', atom_res_number)) &
         call system_abort('Cannot assign pointer to "atom_res_numer" property.')
    if (.not. assign_pointer(at, 'atom_charge', atom_charge)) &
         call system_abort('Cannot assign pointer to "atom_charge" property.')

    do mm=1,at%N
      ! e.g. CP2K needs different name for QM molecules, if use isolated atoms
       sor = ''
       QM_prefix_atom_mol_name = ''
       QM_prefix_atom_mol_name = trim(a2s(atom_mol_name(:,mm)))
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
       write(sor,trim(pdb_format)) 'ATOM  ',mm,trim(a2s(atom_type(:,mm))),trim(a2s(atom_res_name(:,mm))),atom_res_number(mm), &
                             at%pos(1:3,mm),0._dp,0._dp,QM_prefix_atom_mol_name,ElementName(at%Z(mm)),atom_charge(mm)
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

  subroutine write_psf_file_arb_pos(at,psf_file,run_type_string,intrares_impropers,imp_filename,add_silica_23body,&
                                    pos_field_for_connectivity,form_bond,break_bond, remove_qmmm_link_bonds, &
                                    run_suffix, error)
    character(len=*),           intent(in) :: psf_file
    type(atoms),                intent(inout) :: at
    character(len=*), optional, intent(in) :: run_type_string
    type(Table),      optional, intent(in) :: intrares_impropers
    character(80),    optional, intent(in) :: imp_filename
    logical,          optional, intent(in) :: add_silica_23body
    character(len=*), optional, intent(in) :: pos_field_for_connectivity, run_suffix
    integer, optional, intent(in) :: form_bond(2), break_bond(2)
    logical, optional, intent(in) :: remove_qmmm_link_bonds
    integer, optional, intent(out) :: error

    real(dp), pointer :: use_pos(:,:)
    type(Connection) :: t_connect
    type(Atoms) :: at_copy
    !character(len=TABLE_STRING_LENGTH), pointer :: atom_res_name_p(:)
    logical :: use_pos_is_pos, do_add_silica_23body

    INIT_ERROR(error)

    ! save a copy
    at_copy = at

    do_add_silica_23body = optional_default(.false., add_silica_23body)

    ! find desired position field (pos, avgpos, whatever)
    use_pos_is_pos = .false.
    if (present(pos_field_for_connectivity)) then
      if (.not. assign_pointer(at, trim(pos_field_for_connectivity), use_pos)) then
	RAISE_ERROR("calc_topology can't find pos field '"//trim(pos_field_for_connectivity)//"'", error)
      endif
      if (trim(pos_field_for_connectivity) == 'pos') use_pos_is_pos = .true.
    else
      if (.not. assign_pointer(at, 'avgpos', use_pos)) then
	RAISE_ERROR("calc_topology can't find default pos field avgpos", error)
      endif
    endif

    ! copy desired pos to pos, and new connectivity
    !NB don't do if use_pos => pos
    if (.not. use_pos_is_pos) at_copy%pos = use_pos

    if (do_add_silica_23body) then
       call set_cutoff(at_copy,SILICA_2BODY_CUTOFF)
       call calc_connect(at_copy, alt_connect=t_connect)       
    else
       ! use hysteretic connect to get nearest neighbour cutoff
       ! will use default cutoff, which is the same as heuristics_nneighb_only=.true.
       call calc_connect_hysteretic(at, DEFAULT_NNEIGHTOL, DEFAULT_NNEIGHTOL, alt_connect=t_connect)
    endif
    

    call break_form_bonds(at, t_connect, form_bond, break_bond, error=error)
    PASS_ERROR(error)

    ! now create labels using this connectivity object
    ! if cutoff is set to 0, heuristics_nneighb_only doesn't matter
    ! if cutoff is large, then heuristics_nneighb_only must be true
    !    so might as well pass true
    if (do_add_silica_23body) then
       ! cutoff is large, must do heuristics_nneighb_only=.true., but EVB form bond won't work
       call write_psf_file (at,psf_file,run_type_string,intrares_impropers,imp_filename,add_silica_23body,heuristics_nneighb_only=.true.,alt_connect=t_connect,&
            remove_qmmm_link_bonds=remove_qmmm_link_bonds, run_suffix=run_suffix)
    else
       ! cutoff is set to 0, all bonds are already heuristics_nneighb_only except extra EVB form_bond bonds
       call write_psf_file (at,psf_file,run_type_string,intrares_impropers,imp_filename,add_silica_23body,heuristics_nneighb_only=.false.,alt_connect=t_connect,&
            remove_qmmm_link_bonds=remove_qmmm_link_bonds, run_suffix=run_suffix)
    endif
    PASS_ERROR(error)
    call finalise(t_connect)
  end subroutine write_psf_file_arb_pos

  !% Writes PSF topology file, to be used with the PDB coordinate file.
  !% PSF contains the list of atoms, bonds, angles, impropers, dihedrals.
  !
  subroutine write_psf_file(at,psf_file,run_type_string,intrares_impropers,imp_filename,add_silica_23body,heuristics_nneighb_only,alt_connect, &
       remove_qmmm_link_bonds, run_suffix, error)

    character(len=*),           intent(in) :: psf_file
    type(atoms),                intent(in) :: at
    character(len=*), optional, intent(in) :: run_type_string
    type(Table),      optional, intent(in) :: intrares_impropers
    character(80),    optional, intent(in) :: imp_filename
    logical,          optional, intent(in) :: add_silica_23body
    logical, intent(in), optional :: heuristics_nneighb_only
    type(Connection), intent(in), optional, target :: alt_connect
    character(len=*), optional, intent(in) :: run_suffix
    logical, optional, intent(in) :: remove_qmmm_link_bonds
    integer, intent(out), optional :: error

    type(Inoutput)          :: psf
    character(103)          :: sor
    character(*), parameter :: psf_format = '(I8,1X,A4,I5,1X,A4,1X,A4,1X,A4,1X,2G14.6,I8)'
    character(*), parameter :: title_format = '(I8,1X,A)'
    character(*), parameter :: int_format = 'I8'
    integer                 :: mm, i
    character(4)            :: QM_prefix_atom_mol_name
    character, dimension(:,:), pointer   :: atom_type, atom_type_PDB, atom_res_name, atom_mol_name
    integer, dimension(:), pointer       :: atom_res_number
    real(dp), dimension(:), pointer      :: atom_charge
    type(Table)             :: bonds
    type(Table)             :: angles
    type(Table)             :: dihedrals, impropers
    character(len=STRING_LENGTH) :: my_run_type_string
    !integer                 :: run_type
    logical                 :: do_add_silica_23body

    INIT_ERROR(error)

    call system_timer('write_psf_file_pos')

    !intraresidual impropers: table or read in from file
    if (.not.present(intrares_impropers).and..not.present(imp_filename)) call print('WARNING!!! NO INTRARESIDUAL IMPROPERS USED!',verbosity=PRINT_ALWAYS)
    if (present(imp_filename)) then
       call system_abort('Not yet implemented.')
    endif

    do_add_silica_23body = optional_default(.false.,add_silica_23body)
    if (do_add_silica_23body) then !!xxx there must be a SIO2 residue?
       if (at%cutoff.lt.SILICA_2BODY_CUTOFF) call system_abort('The connect cutoff '//at%cutoff//' is smaller than the required cutoff for silica. Cannot build connectivity to silica.')
    endif

    my_run_type_string = optional_default('',run_type_string)
!    if ((trim(my_run_type_string).eq.'QMMM_CORE') .or. &
!        (trim(my_run_type_string).eq.'QMMM_EXTENDED')) then
!       qm_flag_index = get_property(at,'cluster_mark')
!       if (trim(my_run_type_string).eq.'QMMM_EXTENDED') run_type=QMMM_RUN_EXTENDED
!       if (trim(my_run_type_string).eq.'QMMM_CORE') run_type=QMMM_RUN_CORE
!    endif
    if (.not. assign_pointer(at, 'atom_type', atom_type)) &
         call system_abort('Cannot assign pointer to "atom_type" property.')
    if (.not. assign_pointer(at, 'atom_type_PDB', atom_type_PDB)) &
         call system_abort('Cannot assign pointer to "atom_type_PDB" property.')
    if (.not. assign_pointer(at, 'atom_res_name', atom_res_name)) &
         call system_abort('Cannot assign pointer to "atom_res_name" property.')
    if (.not. assign_pointer(at, 'atom_mol_name', atom_mol_name)) &
         call system_abort('Cannot assign pointer to "atom_mol_name" property.')
    if (.not. assign_pointer(at, 'atom_res_number', atom_res_number)) &
         call system_abort('Cannot assign pointer to "atom_res_numer" property.')
    if (.not. assign_pointer(at, 'atom_charge', atom_charge)) &
         call system_abort('Cannot assign pointer to "atom_charge" property.')

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
       QM_prefix_atom_mol_name = trim(a2s(atom_mol_name(:,mm)))
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
       write(sor,psf_format) mm, QM_prefix_atom_mol_name, atom_res_number(mm), &
                     trim(a2s(atom_res_name(:,mm))),trim(a2s(atom_type_PDB(:,mm))),trim(a2s(atom_type(:,mm))), &
                     atom_charge(mm),ElementMass(at%Z(mm))/MASSCONVERT,0
       call print(sor,file=psf)
    enddo
    call print('',file=psf)
call print('PSF| '//at%n//' atoms')
   ! BOND section
    call create_bond_list(at,bonds,do_add_silica_23body,heuristics_nneighb_only,alt_connect=alt_connect, &
         remove_qmmm_link_bonds=remove_qmmm_link_bonds, run_suffix=run_suffix)
    if (any(bonds%int(1:2,1:bonds%N).le.0) .or. any(bonds%int(1:2,1:bonds%N).gt.at%N)) &
       call system_abort('write_psf_file_pos: element(s) of bonds not within (0;at%N]')
    call write_psf_section(data_table=bonds,psf=psf,section='BOND',int_format=int_format,title_format=title_format)
call print('PSF| '//bonds%n//' bonds')

   ! ANGLE section
    call create_angle_list(at,bonds,angles,do_add_silica_23body,heuristics_nneighb_only,alt_connect=alt_connect, &
         remove_qmmm_link_bonds=remove_qmmm_link_bonds, run_suffix=run_suffix)

    if (any(angles%int(1:3,1:angles%N).le.0) .or. any(angles%int(1:3,1:angles%N).gt.at%N)) then
       do i = 1, angles%N
          if (any(angles%int(1:3,i).le.0) .or. any(angles%int(1:3,i).gt.at%N)) &
          call print('angle: '//angles%int(1,i)//' -- '//angles%int(2,i)//' -- '//angles%int(3,i),verbosity=PRINT_ALWAYS)
       enddo
       call system_abort('write_psf_file_pos: element(s) of angles not within (0;at%N]')
    endif
    call write_psf_section(data_table=angles,psf=psf,section='THETA',int_format=int_format,title_format=title_format)
call print('PSF| '//angles%n//' angles')

   ! DIHEDRAL section
    call create_dihedral_list(at,angles,dihedrals,do_add_silica_23body,heuristics_nneighb_only,alt_connect=alt_connect, &
         remove_qmmm_link_bonds=remove_qmmm_link_bonds, run_suffix=run_suffix)
    if (any(dihedrals%int(1:4,1:dihedrals%N).le.0) .or. any(dihedrals%int(1:4,1:dihedrals%N).gt.at%N)) &
       call system_abort('write_psf_file_pos: element(s) of dihedrals not within (0;at%N]')
    call write_psf_section(data_table=dihedrals,psf=psf,section='PHI',int_format=int_format,title_format=title_format)
call print('PSF| '//dihedrals%n//' dihedrals')

   ! IMPROPER section
    call create_improper_list(at,angles,impropers,intrares_impropers=intrares_impropers, & !why not alt_connect?
         remove_qmmm_link_bonds=remove_qmmm_link_bonds, run_suffix=run_suffix)
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
  subroutine create_bond_list(at,bonds,add_silica_23body,heuristics_nneighb_only,alt_connect, &
       remove_qmmm_link_bonds, run_suffix)

  type(Atoms), intent(in)  :: at
  type(Table), intent(out) :: bonds
  logical,     intent(in)  :: add_silica_23body
  logical, intent(in), optional :: heuristics_nneighb_only
  type(Connection), intent(in), optional, target :: alt_connect
  logical, intent(in), optional :: remove_qmmm_link_bonds
  character(len=*), intent(in), optional :: run_suffix

    character(*), parameter  :: me = 'create_bond_list: '

  type(Table) :: atom_a,atom_b
  integer     :: i,j
  integer     :: atom_j
!  logical              :: do_qmmm
!  integer,dimension(3) :: pos_indices
!  integer              :: qm_flag_index
  character(STRING_LENGTH) :: my_run_suffix
  character, pointer, dimension(:,:) :: atom_mol_name
  integer, pointer, dimension(:) :: cluster_mark
  logical :: add_bond, do_remove_qmmm_link_bonds

    call system_timer('create_bond_list')

    do_remove_qmmm_link_bonds = optional_default(.false., remove_qmmm_link_bonds)
    my_run_suffix = optional_default('', run_suffix)

    if (do_remove_qmmm_link_bonds) then
       call assign_property_pointer(at, 'cluster_mark'//trim(my_run_suffix), cluster_mark)
    end if

    if (add_silica_23body) then
       if (at%cutoff.lt.SILICA_2BODY_CUTOFF) call system_abort('The connect cutoff '//at%cutoff//' is smaller than the required cutoff for silica. Cannot build connectivity to silica.')
       if (.not. assign_pointer(at, 'atom_mol_name', atom_mol_name)) &
            call system_abort('Cannot assign pointer to "atom_mol_name" property.')
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
          call bfs_step(at,atom_a,atom_b,nneighb_only=heuristics_nneighb_only,min_images_only=.true.,alt_connect=alt_connect)
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
                if ( ((trim(a2s(atom_mol_name(:,atom_j))) .eq.'SIO2' .and. trim(a2s(atom_mol_name(:,i))).ne.'SIO2')) .or. &
                     ((trim(a2s(atom_mol_name(:,atom_j))) .ne.'SIO2' .and. trim(a2s(atom_mol_name(:,i))).eq.'SIO2')) ) then !silica -- something
!call system_abort('should have not got here')
                   !add only nearest neighbours
                   if (.not.(is_nearest_neighbour_abs_index(at,i,atom_j,alt_connect=alt_connect))) add_bond = .false.
                elseif  ((trim(a2s(atom_mol_name(:,atom_j))) .eq.'SIO2' .and. trim(a2s(atom_mol_name(:,i))).eq.'SIO2')) then !silica -- silica
                   !add atom pairs within SILICON_2BODY_CUTOFF
!call system_abort('what not?')
                   if (.not.are_silica_2body_neighbours(at,i,atom_j,alt_connect=alt_connect)) add_bond = .false. !SILICON_2BODY_CUTOFF for Si-Si and Si-O, nearest neighbours otherwise(Si-H,O-O,O-H,H-H)
                else
!call system_abort('should have not got here')
                   !add only nearest neighbours
                   if (.not.(is_nearest_neighbour_abs_index(at,i,atom_j,alt_connect=alt_connect))) add_bond = .false.
                endif
             endif

             if (do_remove_qmmm_link_bonds) then
                ! skip bond if it is a QM-MM link
                if ((cluster_mark(i) == HYBRID_NO_MARK) .neqv. (cluster_mark(atom_j) == HYBRID_NO_MARK)) then
                   call print('create_bond_list: removing QM/MM link bond '//i//'--'//atom_j, PRINT_VERBOSE)
                   add_bond = .false.
                end if
             end if

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

    do neigh = 1, n_neighbours(this,i)
       j2 = neighbour(this,i,neigh,distance=d,alt_connect=alt_connect)
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

    do neigh = 1, n_neighbours(this,i)
       j2 = neighbour(this,i,neigh,distance=d,alt_connect=alt_connect)
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
  subroutine create_angle_list(at,bonds,angles,add_silica_23body,heuristics_nneighb_only,alt_connect, &
       remove_qmmm_link_bonds, run_suffix)

    type(Atoms), intent(in)  :: at
    type(Table), intent(in)  :: bonds
    type(Table), intent(out) :: angles
    logical,     intent(in)  :: add_silica_23body
    logical,     intent(in), optional  :: heuristics_nneighb_only
    type(Connection), intent(in), optional, target :: alt_connect
    logical, intent(in), optional :: remove_qmmm_link_bonds
    character(len=*), intent(in), optional :: run_suffix
    character(*), parameter  :: me = 'create_angle_list: '
    integer, pointer, dimension(:) :: cluster_mark


    integer     :: i,j
    type(Table) :: atom_a, atom_b
    integer     :: atom_j, atom_1, atom_2
    logical     :: add_angle, do_remove_qmmm_link_bonds
    character, pointer, dimension(:,:) :: atom_mol_name, atom_type
    character(STRING_LENGTH) :: my_run_suffix

    call system_timer('create_angle_list')

    do_remove_qmmm_link_bonds = optional_default(.false., remove_qmmm_link_bonds)
    my_run_suffix = optional_default('', run_suffix)

    if (do_remove_qmmm_link_bonds) then
       call assign_property_pointer(at, 'cluster_mark'//trim(my_run_suffix), cluster_mark)
    end if

    if (add_silica_23body) then
       if (.not. assign_pointer(at, 'atom_type', atom_type)) &
            call system_abort('Cannot assign pointer to "atom_type" property.')
       if (.not. assign_pointer(at, 'atom_mol_name', atom_mol_name)) &
            call system_abort('Cannot assign pointer to "atom_mol_name" property.')
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
          call bfs_step(at,atom_a,atom_b,nneighb_only=heuristics_nneighb_only,min_images_only=.true.,alt_connect=alt_connect)
       endif

       do j = 1,atom_b%N
!call print('atom '//j//' out of '//atom_b%N//' which is atom '//atom_j)
          atom_j = atom_b%int(1,j)
          if (atom_j.ge.atom_2) cycle

          add_angle = .true.

          if (add_silica_23body) then ! do not include H2O -- SiO2 bonds
             if ( ((trim(a2s(atom_mol_name(:,atom_j))) .eq.'SIO2' .and. trim(a2s(atom_mol_name(:,atom_1))).ne.'SIO2')) .or. &
                  ((trim(a2s(atom_mol_name(:,atom_j))) .ne.'SIO2' .and. trim(a2s(atom_mol_name(:,atom_1))).eq.'SIO2')) ) then !silica -- something
                !add only nearest neighbours
                if (.not.(is_nearest_neighbour_abs_index(at,atom_1,atom_j,alt_connect=alt_connect))) add_angle = .false.
             elseif  ((trim(a2s(atom_mol_name(:,atom_j))) .eq.'SIO2' .and. trim(a2s(atom_mol_name(:,atom_1))).eq.'SIO2')) then !silica -- silica
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

          if (do_remove_qmmm_link_bonds) then
             ! skip angle if it spans the QM-MM boundary
             if ( .not. (all(cluster_mark( (/atom_j, atom_1, atom_2/) ) == HYBRID_NO_MARK) .or. &
                         all(cluster_mark( (/atom_j, atom_1, atom_2/) ) /= HYBRID_NO_MARK))) then
                call print('create_angle_list: removing QM/MM spanning angle '//(/atom_j, atom_1, atom_2/), PRINT_VERBOSE)
                add_angle = .false.
             end if
          end if

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
          call bfs_step(at,atom_a,atom_b,nneighb_only=heuristics_nneighb_only,min_images_only=.true.,alt_connect=alt_connect)
       endif

       do j = 1,atom_b%N

          atom_j = atom_b%int(1,j)
          if (atom_j.ge.atom_1) cycle

          add_angle = .true.

          if (add_silica_23body) then ! do not include H2O -- SiO2 bonds
             if ( ((trim(a2s(atom_mol_name(:,atom_j))) .eq.'SIO2' .and. trim(a2s(atom_mol_name(:,atom_2))).ne.'SIO2')) .or. &
                  ((trim(a2s(atom_mol_name(:,atom_j))) .ne.'SIO2' .and. trim(a2s(atom_mol_name(:,atom_2))).eq.'SIO2')) ) then !silica -- something
                !add only nearest neighbours
                if (.not.(is_nearest_neighbour_abs_index(at,atom_2,atom_j,alt_connect=alt_connect))) add_angle = .false.
             elseif  ((trim(a2s(atom_mol_name(:,atom_j))) .eq.'SIO2' .and. trim(a2s(atom_mol_name(:,atom_2))).eq.'SIO2')) then !silica -- silica
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

          if (do_remove_qmmm_link_bonds) then
             ! skip angle if it spans the QM-MM boundary
             if ( .not. (all(cluster_mark( (/atom_j, atom_1, atom_2/) ) == HYBRID_NO_MARK) .or. &
                         all(cluster_mark( (/atom_j, atom_1, atom_2/) ) /= HYBRID_NO_MARK))) then
                call print('create_angle_list: removing QM/MM spanning angle '//(/atom_j, atom_1, atom_2/), PRINT_VERBOSE)
                add_angle = .false.
             end if
          end if

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
          call print('angle: '//angles%int(1,i)//' -- '//angles%int(2,i)//' -- '//angles%int(3,i),verbosity=PRINT_ALWAYS)
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
  subroutine create_dihedral_list(at,angles,dihedrals,add_silica_23body,heuristics_nneighb_only, alt_connect, &
       remove_qmmm_link_bonds, run_suffix)

  type(Atoms), intent(in)  :: at
  type(Table), intent(in)  :: angles
  type(Table), intent(out) :: dihedrals
  logical,     intent(in)  :: add_silica_23body
  logical,     intent(in), optional  :: heuristics_nneighb_only
  type(Connection), intent(in), optional, target :: alt_connect
  logical, intent(in), optional :: remove_qmmm_link_bonds
  character(len=*), intent(in), optional :: run_suffix
  character(*), parameter  :: me = 'create_angle_list: '

  logical     :: do_remove_qmmm_link_bonds
  character(STRING_LENGTH) :: my_run_suffix
  integer, pointer, dimension(:) :: cluster_mark

  integer     :: i,j
  type(Table) :: atom_a, atom_b
  integer     :: atom_j

    call system_timer('create_dihedral_list')

    do_remove_qmmm_link_bonds = optional_default(.false., remove_qmmm_link_bonds)
    my_run_suffix = optional_default('', run_suffix)

    if (do_remove_qmmm_link_bonds) then
       call assign_property_pointer(at, 'cluster_mark'//trim(my_run_suffix), cluster_mark)
    end if


    call initialise(dihedrals,4,0,0,0,0)

    do i=1,angles%N

       if (add_silica_23body) then !only add 2 and 3 body for silica, skip dihedrals
          if (any(at%Z(angles%int(1:3,i)).eq.14)) cycle
       endif

      ! look for one more to the beginning: ??--1--2--3
       call initialise(atom_a,4,0,0,0,0)
       call append(atom_a,(/angles%int(1,i),0,0,0/))
       call bfs_step(at,atom_a,atom_b,nneighb_only=heuristics_nneighb_only,min_images_only=.true.,alt_connect=alt_connect)
       do j = 1,atom_b%N
          atom_j = atom_b%int(1,j)
          if (atom_j.ne.angles%int(2,i)) then
            ! make sure it's not included twice -- no need to O(N^2) check at the end
             if (atom_j.lt.angles%int(3,i)) then
                if (add_silica_23body) then !only add 2 and 3 body for silica, skip dihedrals
                   if (at%Z(atom_j).eq.14) cycle
                endif
                
                if (do_remove_qmmm_link_bonds) then
                   ! skip angle if it spans the QM-MM boundary
                   if ( .not. (all(cluster_mark( (/atom_j, angles%int(1:3,i)/) ) == HYBRID_NO_MARK) .or. &
                               all(cluster_mark( (/atom_j, angles%int(1:3,i)/) ) /= HYBRID_NO_MARK))) then
                      call print('create_dihedral_list: removing QM/MM spanning dihedral '//(/atom_j, angles%int(1:3,i)/), PRINT_VERBOSE)
                      cycle
                   end if
                end if
          
                call append(dihedrals,(/atom_j,angles%int(1,i),angles%int(2,i),angles%int(3,i)/))
             endif
          endif
       enddo
       call finalise(atom_a)
       call finalise(atom_b)

      ! look for one more to the end: 1--2--3--??
       call initialise(atom_a,4,0,0,0,0)
       call append(atom_a,(/angles%int(3,i),0,0,0/))
       call bfs_step(at,atom_a,atom_b,nneighb_only=heuristics_nneighb_only,min_images_only=.true.,alt_connect=alt_connect)
       do j = 1,atom_b%N
          atom_j = atom_b%int(1,j)
          if (atom_j.ne.angles%int(2,i)) then
            ! make sure it's not included twice -- no need to O(N^2) check at the end
             if (atom_j.lt.angles%int(1,i)) then
                if (add_silica_23body) then !only add 2 and 3 body for silica, skip dihedrals
                   if (at%Z(atom_j).eq.14) cycle
                endif

                if (do_remove_qmmm_link_bonds) then
                   ! skip angle if it spans the QM-MM boundary
                   if ( .not. (all(cluster_mark( (/atom_j, angles%int(1:3,i)/) ) == HYBRID_NO_MARK) .or. &
                               all(cluster_mark( (/atom_j, angles%int(1:3,i)/) ) /= HYBRID_NO_MARK))) then
                      call print('create_dihedral_list: removing QM/MM spanning dihedral '//(/atom_j, angles%int(1:3,i)/), PRINT_VERBOSE)
                      cycle
                   end if
                end if

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
  subroutine create_improper_list(at,angles,impropers,intrares_impropers, &
       remove_qmmm_link_bonds, run_suffix)

  type(Atoms),           intent(in)  :: at
  type(Table),           intent(in)  :: angles
  type(Table),           intent(out) :: impropers
  type(Table), optional, intent(in)  :: intrares_impropers
  logical, intent(in), optional :: remove_qmmm_link_bonds  
  character(len=*), intent(in), optional :: run_suffix
  
  logical do_remove_qmmm_link_bonds  
  character(STRING_LENGTH) :: my_run_suffix
  integer, dimension(4) :: imp_atoms
  integer               :: nn,mm
  logical               :: cont
  integer, allocatable, dimension(:) :: count_array ! to count number of bonds
  integer               :: i,j, i_impr
  integer               :: last, tmp
  integer               :: i_pro, tmp_atoms(3)
  logical               :: reordered
  character, dimension(:,:), pointer   :: atom_res_name, atom_type
  integer, pointer :: atom_res_number(:), cluster_mark(:)

    call system_timer('create_improper_list')

    do_remove_qmmm_link_bonds = optional_default(.false., remove_qmmm_link_bonds)
    my_run_suffix = optional_default('', run_suffix)

    if (do_remove_qmmm_link_bonds) then
       call assign_property_pointer(at, 'cluster_mark'//trim(my_run_suffix), cluster_mark)
    end if

    if (.not. at%connect%initialised) &
       call system_abort('create_bond_list: connectivity not initialised, call calc_connect first')

    call initialise(impropers,4,0,0,0,0)

    allocate (count_array(angles%N))

call print("create_improper_list", PRINT_ANAL)
    do i = 1,at%N
call print("i " // i // " " // trim(a2s(at%species(:,i))), PRINT_ANAL)
      if (.not.any(trim(a2s(at%species(:,i))).eq.(/'C','N'/))) cycle
      count_array = 0
      where (angles%int(2,1:angles%N).eq.i) count_array = 1
      if (sum(count_array(1:size(count_array))).ne.3) cycle
     ! only N with a neighbour that has 3 neighbors can stay
      if (trim(a2s(at%species(:,i))).eq.'N') then
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
call print("i is N with 3 neighbours, continuing", PRINT_ANAL)

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

call print("original imp order " // imp_atoms // " Zs " // at%Z(imp_atoms), PRINT_ANAL)
!VVV ORDER is done according to the topology file! - and is read in when finding motifs
!if you don't do this, you won't have only the backbone impropers!
      if (.not. assign_pointer(at, 'atom_res_number', atom_res_number)) &
           call system_abort('Cannot assign point to "atom_res_number" property.')
      if (atom_res_number(imp_atoms(1)) == atom_res_number(imp_atoms(2)) .and. &
          atom_res_number(imp_atoms(1)) == atom_res_number(imp_atoms(3)) .and. &
          atom_res_number(imp_atoms(1)) == atom_res_number(imp_atoms(4))) cycle ! these should be added when identifying the residues
!      if (trim(a2s(atom_res_name(:,imp_atoms(1)))) .eq. trim(a2s(atom_res_name(:,imp_atoms(2)))) .and. &
!          trim(a2s(atom_res_name(:,imp_atoms(1)))) .eq. trim(a2s(atom_res_name(:,imp_atoms(3)))) .and. &
!          trim(a2s(atom_res_name(:,imp_atoms(1)))) .eq. trim(a2s(atom_res_name(:,imp_atoms(4))))) &
!          cycle  ! these should be added when identifying the residues

!ORDER: check charmm.pot file - start with $i, end with  H or O or N, in this order -- for intraresidual residues this can be needed later on...
      reordered = .true.
      tmp = 0
      ! if there is H
      last = find_in_array(at%Z(imp_atoms(2:4)),1)
      if (last.gt.0) then ! found a H
call print("reorder H to end", PRINT_ANAL)
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
        if (last.gt.0) then ! found an O
call print("reorder O to end", PRINT_ANAL)
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
          if (last.gt.0) then ! found an N
call print("reorder P to end", PRINT_ANAL)
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
call print("reordered " // reordered, PRINT_ANAL)
call print("final imp order " // imp_atoms // " Zs " // at%Z(imp_atoms), PRINT_ANAL)

      !checking and adding only backbone (i.e. not intraresidual impropers) where the order of the 2nd and 3rd atoms doesn't matter
      if (atom_res_number(imp_atoms(1)) == atom_res_number(imp_atoms(2)) .and. &
          atom_res_number(imp_atoms(1)) == atom_res_number(imp_atoms(3)) .and. &
          atom_res_number(imp_atoms(1)) == atom_res_number(imp_atoms(4))) cycle ! these should be added when identifying the residues
!      if (trim(a2s(atom_res_name(:,imp_atoms(1)))) .eq. trim(a2s(atom_res_name(:,imp_atoms(2)))) .and. &
!          trim(a2s(atom_res_name(:,imp_atoms(1)))) .eq. trim(a2s(atom_res_name(:,imp_atoms(3)))) .and. &
!          trim(a2s(atom_res_name(:,imp_atoms(1)))) .eq. trim(a2s(atom_res_name(:,imp_atoms(4))))) &
!          cycle  ! these should be added when identifying the residues

      if (.not.reordered) then
        ! Found N-C-CP1-CP3 Pro backbone, reordering according to atomic types (could be also according to the H neighbours)
!         call print('|PRO Found Pro backbone')
         if (.not. assign_pointer(at, 'atom_type', atom_type)) &
              call system_abort('Cannot assign pointer to "atom_type" property.')
call print("atom type " // trim(a2s(atom_type(:,imp_atoms(1)))), PRINT_ANAL)
call print("atom type " // trim(a2s(atom_type(:,imp_atoms(2)))), PRINT_ANAL)
call print("atom type " // trim(a2s(atom_type(:,imp_atoms(3)))), PRINT_ANAL)
call print("atom type " // trim(a2s(atom_type(:,imp_atoms(4)))), PRINT_ANAL)
         if (trim(a2s(atom_type(:,imp_atoms(1)))).ne.'N') call system_abort('something has gone wrong. what is this if not proline? '// &
                     trim(a2s(atom_type(:,imp_atoms(1))))//imp_atoms(1)//'--'// &
                     trim(a2s(atom_type(:,imp_atoms(2))))//imp_atoms(2)//'--'// &
                     trim(a2s(atom_type(:,imp_atoms(3))))//imp_atoms(3)//'--'// &
                     trim(a2s(atom_type(:,imp_atoms(4))))//imp_atoms(4))
         tmp_atoms = 0
         do i_pro = 2,4
            if (trim(a2s(atom_type(:,imp_atoms(i_pro)))).eq.'C')   tmp_atoms(1) = imp_atoms(i_pro)
            if (trim(a2s(atom_type(:,imp_atoms(i_pro)))).eq.'CP1') tmp_atoms(2) = imp_atoms(i_pro)
            if (trim(a2s(atom_type(:,imp_atoms(i_pro)))).eq.'CP3') tmp_atoms(3) = imp_atoms(i_pro)
         enddo
         if (any(tmp_atoms(1:3).eq.0)) call system_abort('something has gone wrong. what is this if not proline?'// &
                     trim(a2s(atom_type(:,imp_atoms(1))))//imp_atoms(1)//'--'// &
                     trim(a2s(atom_type(:,imp_atoms(2))))//imp_atoms(2)//'--'// &
                     trim(a2s(atom_type(:,imp_atoms(3))))//imp_atoms(3)//'--'// &
                     trim(a2s(atom_type(:,imp_atoms(4))))//imp_atoms(4))
         imp_atoms(2:4) = tmp_atoms(1:3)
!         call print('Reordered Pro improper '// &
!                     trim(a2s(atom_type(:,imp_atoms(1))))//imp_atoms(1)//'--'// &
!                     trim(a2s(atom_type(:,imp_atoms(2))))//imp_atoms(2)//'--'// &
!                     trim(a2s(atom_type(:,imp_atoms(3))))//imp_atoms(3)//'--'// &
!                     trim(a2s(atom_type(:,imp_atoms(4))))/imp_atoms(4))
      endif

      if (do_remove_qmmm_link_bonds) then
         ! skip improper if it spans the QM-MM boundary
         if ( .not. (all(cluster_mark(imp_atoms) == HYBRID_NO_MARK) .or. &
                     all(cluster_mark(imp_atoms) /= HYBRID_NO_MARK))) then
            call print('create_improper_list: removing QM/MM spanning improper '//imp_atoms)
            cycle
         end if
      end if


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
       call print('WARNING!!! NO INTRARESIDUAL IMPROPERS USED!!!',verbosity=PRINT_ALWAYS)
    endif

   ! final check
    if (any(impropers%int(1:4,1:impropers%N).le.0) .or. any(impropers%int(1:4,1:impropers%N).gt.at%N)) &
       call system_abort('create_improper_list: element(s) of impropers not within (0;at%N]')

    call system_timer('create_improper_list')

  end subroutine create_improper_list

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
       ' is less than the silica_cutoff. The charges can be totally wrong. Check your connection object.',PRINT_ALWAYS)

    allocate(charge(at%N))
    charge = 0._dp

    allocate(silica_mask(at%N)) !whether the atom is in the SiOH list
    silica_mask = .false.
    silica_mask(SiOH_list%int(1,1:SiOH_list%N)) = .true.

    do iatom = 1, SiOH_list%N
       atom_i = SiOH_list%int(1,iatom)
       do jatom = 1, n_neighbours(at,atom_i)
          atom_j = neighbour(at, atom_i, jatom)
         ! check if atom_j is in SiOH_list
          if (.not.silica_mask(atom_j)) then
!          if (.not.any(atom_j.eq.SiOH_list%int(1,1:SiOH_list%N))) then
!             ! if it's a water O/H, ignore it!
!             if (('TIP3'.eq.trim(residue_names(atom_j))) .or. &
!                 ('DOPA'.eq.trim(residue_names(atom_j))) .or. &
!                 ('TIP' .eq.trim(residue_names(atom_j)))) then
                 cycle
!             else
!                call print('Not found atom '//ElementName(at%Z(atom_j))//' '//atom_j//' in '//SiOH_list%int(1,1:SiOH_list%N),verbosity=PRINT_ALWAYS)
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
                      call print(atom_i//trim(ElementName(at%Z(atom_i)))//'--'//atom_j//trim(ElementName(at%Z(atom_j)))//' of silica is none of H,O,Si!',verbosity=PRINT_ALWAYS)
                      call system_abort('create_pos_dep_charge: unknown neighbour pair')
                   endif
                endif
             endif
          endif
       enddo
    enddo

    deallocate(silica_mask)

    call print('Calculated charge on atoms:',verbosity=PRINT_ANAL)
    do jj = 1, at%N
       call print ('   atom '//jj//': '//charge(jj),verbosity=PRINT_ANAL)
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

  subroutine find_water_monomer(at,water_index,OHH_ordercheck,monomer_cutoff,error)

     type(atoms), intent(in) :: at
     integer, dimension(3,at%N/3), intent(out) :: water_index
     logical, intent(in), optional :: OHH_ordercheck
     real(dp), intent(in), optional :: monomer_cutoff
     integer, intent(out), optional :: error

     real(dp) :: r, roh1, roh2, my_monomer_cutoff
     integer :: i, j, n, h1, h2, oindex, iO, iH1, iH2
     logical, dimension(at%N) :: H_associated
     logical :: do_OHH_ordercheck


     ! loop through atoms, find oxygens and their closest hydrogens
     oindex = 0
     H_associated = .false.

     do_OHH_ordercheck = optional_default(.true.,OHH_ordercheck)
     my_monomer_cutoff = optional_default(at%cutoff, monomer_cutoff)

     if(do_OHH_ordercheck) then
        do i = 1, at%N
           if(at%Z(i) == 8) then
              oindex = oindex + 1


              roh1 = my_monomer_cutoff ! initialise smaller oh distance
              roh2 = my_monomer_cutoff ! initialise larger oh distance
              h1 = 0
              h2 = 0
              do n = 1, n_neighbours(at, i)
                 j = neighbour(at, i, n, distance=r)

                 if( (at%Z(j) /= 1) .or. H_associated(j) ) cycle

                 if(r < roh1) then
                    if(h1 /= 0) then ! H1 was already found, so move it to H2
                       roh2 = roh1
                       h2 = h1
                    end if
                    roh1 = r
                    h1 = j
                 else if(r < roh2) then
                    roh2 = r
                    h2 = j
                 end if
              end do
              if(h1 == 0 .or. h2 == 0) then
                 RAISE_ERROR("Cannot find hydrogens for oxygen index "//i//". h1="//h1//", h2="//h2//" with cutoff of "//my_monomer_cutoff//" A", error)
!                 call write(at, "stdout")
              end if
              H_associated(h1) = .true.
              H_associated(h2) = .true.

              water_index(:,oindex) = (/i, h1, h2/)
           end if
        end do

        ! sanity check
        if(oindex /= at%N/3) then
           RAISE_ERROR("Number of oxygens is not equal to at%N/3", error)
        end if
     else
        do i = 1, at%N/3
           iO = (i-1)*3+1
           iH1 = iO+1
           iH2 = iO+2
           if( (at%Z(iO)/=8) .or. (at%Z(iH1)/=1) .or. (at%Z(iH2)/=1) ) then
              RAISE_ERROR("Atom order incorrect, does not match OHH order.", error)
           endif
           water_index(:,i) = (/iO, iH1, iH2/)
        enddo
     endif

  endsubroutine find_water_monomer

  subroutine find_A2_monomer(at,atomic_number,monomer_cutoff,A2_index,error)

     type(atoms), intent(in) :: at
     integer, intent(in) :: atomic_number
     real(dp), intent(in) :: monomer_cutoff
     integer, dimension(:), intent(out) :: A2_index
     integer, intent(out), optional :: error

     real(dp) :: r, r12
     integer :: i, j, n, jn, atom_2, atom_n1, atom_n2
     logical, dimension(at%N) :: atom_associated


     ! loop through atoms, find oxygens and their closest hydrogens
     atom_associated = .false.
     A2_index = 0

     do i = 1, at%N
        if(at%Z(i) == atomic_number .and. .not. atom_associated(i) ) then

           r12 = monomer_cutoff ! initialise distance
           atom_2 = 0
           do n = 1, n_neighbours(at, i)
              j = neighbour(at, i, n, distance=r, jn = jn, max_dist=monomer_cutoff)
              if( j == 0 ) cycle
              if( (at%Z(j) /= atomic_number) .or. atom_associated(j) ) cycle

              if(r < r12) then
                 r12 = r
                 atom_2 = j
                 atom_n2 = n
                 atom_n1 = jn
              endif
           end do
           if(atom_2 == 0) then
              RAISE_ERROR("Cannot find pair for atom index "//i, error)
           end if
           atom_associated(i) = .true.
           atom_associated(atom_2) = .true.

           A2_index(i) = atom_n2
           A2_index(atom_2) = atom_n1
        endif
     enddo

  endsubroutine find_A2_monomer

  subroutine find_AB_monomer(at,atomic_number,monomer_cutoff,AB_index,error)

     type(atoms), intent(in) :: at
     integer, dimension(2), intent(in) :: atomic_number 
     real(dp), intent(in) :: monomer_cutoff
     integer, dimension(:,:), intent(out) :: AB_index
     integer, intent(out), optional :: error

     real(dp) :: r, r12
     integer :: i, j, n, atom_2, atom_index
     logical, dimension(at%N) :: atom_associated


     ! loop through atoms, find oxygens and their closest hydrogens
     atom_index = 0
     atom_associated = .false.

     do i = 1, at%N
        if(at%Z(i) == atomic_number(1) .and. .not. atom_associated(i) ) then
           atom_index = atom_index + 1

           r12 = monomer_cutoff ! initialise distance
           atom_2 = 0
           do n = 1, n_neighbours(at, i)
              j = neighbour(at, i, n, distance=r, max_dist=monomer_cutoff)
              if( j == 0 ) cycle

              if( (at%Z(j) /= atomic_number(2)) .or. atom_associated(j) ) cycle

              if(r < r12) then
                 r12 = r
                 atom_2 = j
              endif
           end do
           if(atom_2 == 0) then
              RAISE_ERROR("Cannot find pair for atom index "//i, error)
           end if
           atom_associated(i) = .true.
           atom_associated(atom_2) = .true.

           if(atom_index <= size(AB_index,2)) then
              AB_index(:,atom_index) = (/i, atom_2/)
           else
              RAISE_ERROR("A2_index is too small: "//size(AB_index,2)//". Size required is at least "//atom_index,error)
           endif
        endif
     enddo

     ! sanity check
     if(atom_index /= count(at%Z == atomic_number(1)) .or. atom_index /= count(at%Z == atomic_number(2))) then
        RAISE_ERROR("Number of monomers is not equal to the number of monomer atoms.", error)
     end if

  endsubroutine find_AB_monomer

       
   subroutine break_form_bonds(at, conn, form_bond, break_bond, error)
      type(Atoms), intent(in) :: at
      type(Connection), intent(inout) :: conn
      integer, intent(in), optional :: form_bond(2), break_bond(2)
      integer, intent(out), optional :: error

      logical :: bond_exists
      real(dp) :: form_bond_dist
      integer :: j, ji
      integer :: shift(3)

      INIT_ERROR(error)

       if (present(form_bond)) then
	 if (all(form_bond >= 1) .and. all(form_bond <= at%N)) then
	    bond_exists = .false.
	    do ji=1, n_neighbours(at, form_bond(1), alt_connect=conn)
	       j = neighbour(at, form_bond(1), ji, shift=shift, alt_connect=conn)
	       if (j == form_bond(2)) then
		  bond_exists = .true.
		  exit
	       endif
	    end do
	    if (.not. bond_exists) then
	       form_bond_dist = distance_min_image(at, form_bond(1), form_bond(2), shift=shift)
	       call add_bond(conn, at%pos, at%lattice, form_bond(1), form_bond(2), shift, form_bond_dist, error=error)
	       PASS_ERROR(error)
	    endif
	 endif ! valid atom #s
       endif ! present(form_bond)
       if (present(break_bond)) then
	 if (all(break_bond >= 1) .and. all(break_bond <= at%N)) then
	    bond_exists = .false.
	    do ji=1, n_neighbours(at, break_bond(1), alt_connect=conn)
	       j = neighbour(at, break_bond(1), ji, shift=shift, alt_connect=conn)
	       if (j == break_bond(2)) then
		  bond_exists = .true.
		  exit
	       endif
	    end do
	    if (bond_exists) then
	       call remove_bond(conn, break_bond(1), break_bond(2), shift, error=error)
	       PASS_ERROR(error)
	    endif
	 endif ! valid atom #s
       endif ! present(break_bond)
   end subroutine break_form_bonds

   function find_motif_backbone(motif, is_proline) result(backbone)
   integer :: motif(:,:)
   logical :: is_proline
   integer :: backbone(5,3)

   integer :: i, ji, j, ki, k, li, l, mi, m, ii, ni, n, found_N, found_aC, found_cC

   is_proline = .false.

     backbone = 0
     do i=1, size(motif,1) ! look for N
!call print("check atom i " //i // " Z=" // motif(i,1))
       if (motif(i,1) /= 7) cycle
       ! i is an N
       found_N = i
!call print("atom i is a N, continuing")
       do ji=2, size(motif,2) ! look for alpha-C
	 j=motif(i,ji); if (j == 0) exit
!call print("check atom j "//j//" Z=" // motif(j,1))
	 if (motif(j,1) == 6) then ! j is a possible alpha-C
!call print("atom j is a C, continuing")
	   found_aC = j
	   do ki=2, size(motif,2) ! look for possible carboxyl-C
	     k = motif(j,ki); if (k == 0) exit
!call print("check atom k " //k//" Z=" // motif(k,1))
	     if (motif(k,1) == 6) then ! k is a possible carboxyl-C
!call print("atom k  is a C, continuing")
	       found_cC = k
	       do li=2, size(motif,2) ! loook for possible carboxyl-O
		 l = motif(k, li); if (l == 0) exit
!call print("check atom l "//l//" Z=" // motif(l,1))
		 if (motif(l,1) == 8 .and. count(motif(l,2:) /= 0) == 1) then ! found carboxyl O
!call print("atom l is an O with coord 1, continuing")

		     ii = 1
		     backbone(ii,1) = found_N
		     do mi=2, size(motif,2) ! find neighbours of N for backbone
		       m = motif(found_N, mi); if (m == 0) exit
!call print("neighbor m " // m // " Z="//motif(m,1)// " neighbor of N")
		       if (m == found_aC) cycle
		       if (motif(m,1) == 1) then
!call print("neighbor m is a H, adding to backbone")
			 ii = ii + 1
			 backbone(ii,1) = m
		       else if (motif(m,1) == 6) then 
			 is_proline = .true.
		       else
			 call print("find_motif_backbone confused by neighbor of N", PRINT_VERBOSE)
			 backbone = 0
			 return
		       endif
		     end do ! mi

		     ii = 1
		     backbone(ii,2) = found_aC
		     do mi=2, size(motif,2) ! find neighbours of alpha-C for backbone
		       m = motif(found_aC, mi); if (m == 0) exit
!call print("neighbor m " // m // " Z="//motif(m,1)// " neighbor of alpha-C")
		       if (m == found_N .or. m == found_cC) cycle
		       if (motif(m,1) == 1) then
!call print("neighbor m is a H, adding to backbone")
			 ii = ii + 1
			 backbone(ii,2) = m
		       else if (motif(m,1) /= 6) then 
			  call print("find_motif_backbone confused by neighbor of alpha-C, neither H or C", PRINT_VERBOSE)
			  backbone = 0
			  return
		       endif
		     end do ! mi

		     ii = 1
		     backbone(ii,3) = found_cC
		     do mi=2, size(motif,2) ! find neighbours of carboxyl-C for backbone
		       m = motif(found_cC, mi); if (m == 0) exit
!call print("neighbor m " // m // " Z="//motif(m,1)// " neighbor of carboxyl-C")
		       if (m == found_aC) cycle
		       if (motif(m,1) == 8) then
!call print("neighbor m is an O, adding to backbone")
			 ii = ii + 1
			 backbone(ii,3) = m
			 do ni=2, size(motif,2) ! look for OH
			   n = motif(m, ni); if (n == 0) exit
			   if (n == found_cC) cycle
!call print("neighbor n " // n // " Z="//motif(n,1)// " neighbor of carboxyl-C-O")
			   if (motif(n,1) == 1) then
!call print("neighbor n is an H, adding to backbone")
			     ii = ii + 1
			     backbone(ii,3) = n
			   else
			     call print("find_motif_backbone confused by non-H neighbor of carboxyl-O", PRINT_VERBOSE)
			     backbone = 0
			     return
			   endif
			 end do
		       else
			 call print("find_motif_backbone confused by neighbor of carboxyl-C", PRINT_VERBOSE)
			 backbone = 0
			 return
		       endif
		     end do ! mi
		 endif ! found carboxyl O
	       end do ! li
	     endif ! carboxyl-C
	   end do ! ki
	 endif ! found alpha C
       end do ! ji
     end  do ! i
   end function find_motif_backbone

   subroutine find_general_monomer(at,monomer_index,signature,is_associated,cutoff,general_ordercheck,use_smooth_cutoff,error)
     type(atoms), intent(in) :: at
     integer, intent(in), dimension(:) :: signature
     real(dp), dimension(:), allocatable :: r
     real(dp) :: r_ij, cutoff
     integer, dimension(1) :: temp
     integer, dimension(:), allocatable :: indices
     integer, dimension(:,:), allocatable :: monomer_index, monomer_index_working
     integer, intent(out), optional :: error
     logical, optional :: use_smooth_cutoff, general_ordercheck
     integer :: i, j, k, n, Z_uniq_index, Z_uniq, Z_uniq_pos, monomers_found
     logical, dimension(:), intent(inout) :: is_associated
     logical :: do_general_ordercheck, my_use_smooth_cutoff

     my_use_smooth_cutoff = optional_default(.false.,use_smooth_cutoff) 
     do_general_ordercheck = optional_default(.true.,general_ordercheck)

     allocate(monomer_index(size(signature),0))
     allocate(monomer_index_working(size(signature),0))
     allocate(r(size(signature)))
     allocate(indices(size(signature)))

     monomers_found = 0
     if(do_general_ordercheck) then
        do i = 1, at%N
           if(.not. is_associated(i) .and. any(signature .eq. at%Z(i))) then
              if (at%Z(i) == 1) cycle ! for now, don't let H be at the centre of a monomer

              r = cutoff ! initialise distances
              indices = 0 ! initialise indices
              temp = minloc(signature, signature .eq. at%Z(i)) ! will return location of first element of signature with correct Z
              indices(temp(1))=i
              r(temp(1))=0.0

              do n = 1, n_neighbours(at, i)
                 j = neighbour(at, i, n, distance=r_ij)
                 if(is_associated(j) .or. .not. any(signature .eq. at%Z(j)) ) cycle

                 do k=1,maxval(signature)
                   if (at%Z(j) .eq. k) then
                     if(r_ij .lt. maxval(r, mask = signature .eq. at%Z(j)) ) then
                       temp =  maxloc(r, mask = signature .eq. at%Z(j))
                       r(temp(1)) = r_ij
                       indices(temp(1)) = j
                       exit
                     end if
                   end if
                 end do
              end do

              if(any(indices .eq. 0)) cycle ! couldn't make a monomer with i as central atom

              monomers_found = monomers_found + 1 ! number of monomers found so far
              do k=1,size(signature)
                is_associated(indices(k))=.true.
              end do

              deallocate(monomer_index_working)
              allocate(monomer_index_working(size(monomer_index,1),size(monomer_index,2)))
              monomer_index_working =monomer_index
              deallocate(monomer_index)
              allocate(monomer_index(size(monomer_index_working,1),size(monomer_index_working,2)+1))
              monomer_index(:,:monomers_found-1) = monomer_index_working
              monomer_index(:,monomers_found) = indices
           end if
        end do
      else
        RAISE_ERROR("Finding general monomers without order checking not supported", error)
      end if

      deallocate(monomer_index_working)
      deallocate(r)
      deallocate(indices)

   end subroutine find_general_monomer

   subroutine find_monomer_pairs(at_in,monomer_pairs,mean_pos_diffs,pairs_diffs_map,monomer_one_index,monomer_two_index,monomers_identical,double_count,cutoff,pairs_shifts,error)
   ! finds pairs of monomers combined into dimer if *any pair* of atoms are within cutoff. Returns 2 by n array of pairs, 
   ! the elements of which refer to the second index of the monomer_one_index and monomer_two_index matrices respectively.
   ! also returns a set of diffs, which are vectors between the mean positions of the shifted versions of these pairs.
     type(atoms), intent(in) :: at_in
     integer, dimension(:,:), allocatable, intent(out) :: monomer_pairs
     real(dp),dimension(:,:), allocatable, intent(out) :: mean_pos_diffs
     integer, dimension(:), allocatable, intent(out) :: pairs_diffs_map
     integer, intent(in), dimension(:,:) :: monomer_one_index, monomer_two_index
     logical, intent(in) :: monomers_identical, double_count
     real(dp), intent(in) :: cutoff
     integer, dimension(:,:), allocatable, intent(out), optional :: pairs_shifts
     integer, intent(out), optional :: error

     type(atoms) :: at
     integer, dimension(:,:), allocatable ::  mean_pos_shifts
     integer, dimension(:), allocatable :: atomic_index_one, atomic_index_two
     integer :: i, j, i_atomic, j_atomic, n, i_neighbour, i_desc, k, m, n_pairs,monomer_one_size ,monomer_two_size,rep,n_repeats
     real(dp) :: r_one_two, mass_one, mass_two, temp_dist,min_dist

     real(dp), dimension(3) :: diff_one_two, mean_pos_one, mean_pos_two ! unweighted mean position of atoms in monomer, if monomer straddles cell boundary this will be in the middle of the cell
     integer, dimension(3) :: shift_one_two, shift_one, shift_two,min_image_shift
     integer, dimension(2) :: temp2d
     integer, dimension(1) :: temp1d
     logical, dimension(:), allocatable :: shifts_mask

     ! make a copy of the atoms object
     at = at_in
     call set_cutoff(at,cutoff+0.5_dp)
     call calc_connect(at)

     allocate(monomer_pairs(2,0))
     allocate(mean_pos_shifts(3,0))
     allocate(mean_pos_diffs(3,0))
     allocate(pairs_diffs_map(0))
     allocate(shifts_mask(0))

     monomer_one_size = size(monomer_one_index,1)
     monomer_two_size = size(monomer_two_index,1)

     allocate(atomic_index_one(monomer_one_size))
     allocate(atomic_index_two(monomer_two_size))

     i_desc=0
     n_pairs=0

     ! loop over all occurences of first monomer
     do i=1,size(monomer_one_index,2)
       atomic_index_one = monomer_one_index(:,i)
       mean_pos_one = calc_mean_pos(at,atomic_index_one)  

       !loop over neighbours of each atom, check at least one pair of heavy atoms within cutoff and belong to an occurrence of monomer_two
       do i_atomic =1,size(monomer_one_index,1)
         if (at%Z(atomic_index_one(i_atomic)) == 1) cycle
         do n = 1, n_neighbours(at,atomic_index_one(i_atomic))

           i_neighbour = neighbour(at,atomic_index_one(i_atomic),n,distance=r_one_two,shift=shift_one_two)
           if (at%Z(i_neighbour) == 1 ) cycle
           if( r_one_two >= cutoff ) cycle
           temp2d = maxloc(monomer_two_index, monomer_two_index .eq. i_neighbour)
           if (any(temp2d .eq. 0)) cycle                                     ! atom i_neighbour does not belong to the type of monomer we're looking for
           j = temp2d(2)                                                     ! atom i_neighbour belongs to monomer j
           atomic_index_two = monomer_two_index(:,j)
           mean_pos_two = calc_mean_pos(at,atomic_index_two)
  
           temp_dist = distance_min_image(at,i_neighbour,mean_pos_two,shift=shift_two)
           shift_one_two = shift_one_two + shift_two

           ! get the diff vector between the shifted mean positions                                                                                                                                                                                                             
           min_dist = distance_min_image(at,mean_pos_one,mean_pos_two,shift=min_image_shift)
           diff_one_two = diff_min_image(at,mean_pos_one,mean_pos_two) + matmul(at%lattice,shift_one_two-min_image_shift)

           ! this makes sure we don't double count. All monomer pairs in the unit cell are found, 
           ! as are half of the monomer pairs in which one of the monomers is from a neighbour cell
           if (monomers_identical) then
             if (all(shift_one_two .eq. 0)) then
               if (i .ge. j) cycle
             else
               if (i .gt. j) cycle
               if (i == j) then
                 temp1d = maxloc(monomer_pairs(1,:), monomer_pairs(1,:) .eq. i .and. monomer_pairs(2,:) .eq. i)                                                      
                 k=temp1d(1)
                 if (k > 0) then       
                   shifts_mask =(/ mean_pos_shifts(1,:) .eq. -shift_one_two(1)/) .and. &
                                (/ mean_pos_shifts(2,:) .eq. -shift_one_two(2)/) .and. &
                                (/ mean_pos_shifts(3,:) .eq. -shift_one_two(3)/)
                   if (count( (/pairs_diffs_map .eq. k /) .and. shifts_mask) > 0)   cycle !  cycle if negative shift already present  
                 end if
               end if
             end if
           end if

           ! check if this combination of monomers was already found
           temp1d = maxloc(monomer_pairs(1,:), monomer_pairs(2,:) .eq. j .and. monomer_pairs(1,:) .eq. i) 
           k=temp1d(1)                                                                                    ! if so, pair is at (:,k) in the monomer_pairs array , else equals 0

           if (k > 0 ) then   ! have found this pair before    
             shifts_mask =(/ mean_pos_shifts(1,:) .eq. shift_one_two(1)/) .and. &
                          (/ mean_pos_shifts(2,:) .eq. shift_one_two(2)/) .and. &
                          (/ mean_pos_shifts(3,:) .eq. shift_one_two(3)/)
             if (count( (/pairs_diffs_map .eq. k/) .and. shifts_mask) > 0) cycle ! cycle if we've already found this pair with this shift 

           else               ! new pair
             n_pairs = n_pairs + 1
             call reallocate(monomer_pairs,2,n_pairs,copy=.true.)
             monomer_pairs(:,n_pairs) = (/ i,j /)
             k=n_pairs
           end if

           n_repeats = 1
           if (double_count) n_repeats = 2
           do rep=1,n_repeats       ! save this shift, possibly twice if we're double-counting
             i_desc = i_desc + 1
             call reallocate(mean_pos_shifts,3,i_desc,copy=.true.)
             call reallocate(mean_pos_diffs,3,i_desc,copy=.true.)
             call reallocate(pairs_diffs_map,i_desc,copy=.true.)
             call reallocate(shifts_mask,i_desc)
             mean_pos_shifts(:,i_desc) = shift_one_two
             mean_pos_diffs(:,i_desc) = diff_one_two 
             pairs_diffs_map(i_desc) = k

           end do
         end do
       end do
     end do

     if (present(pairs_shifts)) then
       allocate(pairs_shifts(3,size(mean_pos_shifts,2)))
       pairs_shifts = mean_pos_shifts
     end if
     call finalise(at)
     deallocate(mean_pos_shifts)
     deallocate(atomic_index_one)
     deallocate(atomic_index_two)
     deallocate(shifts_mask)

   end subroutine find_monomer_pairs  

   subroutine find_monomer_triplets(at,monomer_triplets,triplets_diffs,triplets_diffs_map,monomer_one_index,monomer_two_index,monomer_three_index,one_two_identical,one_three_identical,two_three_identical,cutoff,error)
   ! loops through pairs of monomers made by find_monomer_pairs above and makes trimers if a third monomer is within cutoff. Returns 3 by n array, 
   ! the elements of which refer to the second index of the monomer_{one,two,three}_index matrices. 
   ! the columns of triplets shifts are (k,shift_one_two,shift_one_three), where k refers to a column of monomer_triplets and the shifts correspond to a set of periodic images which are within cutoff.
     type(atoms), intent(in) :: at
     integer, dimension(:,:), intent(out), allocatable :: monomer_triplets
     real(dp), dimension(:,:), intent(out), allocatable :: triplets_diffs
     integer, dimension(:), intent(out), allocatable ::  triplets_diffs_map
     integer, intent(in), dimension(:,:) :: monomer_one_index, monomer_two_index, monomer_three_index
     logical, intent(in) :: one_two_identical,one_three_identical, two_three_identical
     real(dp), intent(in) :: cutoff
     integer, intent(out), optional :: error

     integer, dimension(:,:), allocatable :: pairs_one_two, pairs_one_three,shifts_one_two,shifts_one_three,triplets_shifts
     real(dp), dimension(:,:), allocatable ::  diffs_one_two,diffs_one_three
     integer, dimension(:), allocatable :: map_one_two, map_one_three
     integer :: i, j, k,pos_ij,pos_ik,pos_jk,i_map,i_desc
     real(dp) :: dist, pairwise_cutoff,min_dist

     real(dp), dimension(3) :: diff_ij,diff_ik,diff_jk, min_distances
     integer, dimension(3) :: shift_ij, shift_ik, shift_jk
     integer, dimension(1) :: temp1d
     logical :: double_count

     double_count=.false.
     pairwise_cutoff = 2.0_dp*cutoff+ 0.1_dp ! plenty big to find all relevant dimers

     call find_monomer_pairs(at,pairs_one_two,  diffs_one_two,  map_one_two,  monomer_one_index,monomer_two_index  ,one_two_identical , double_count,pairwise_cutoff,shifts_one_two  ,error)
     call find_monomer_pairs(at,pairs_one_three,diffs_one_three,map_one_three,monomer_one_index,monomer_three_index,one_three_identical,double_count,pairwise_cutoff,shifts_one_three,error)
     ! for any equivalent monomer types we have to append the negative self-pairs since these are excluded by find_monomer_pairs

      if (.not. allocated(monomer_triplets)) allocate(monomer_triplets(3,0))
      if (.not. allocated(triplets_diffs)) allocate(triplets_diffs(6,0))
      if (.not. allocated(triplets_diffs_map)) allocate(triplets_diffs_map(0))
      if (.not. allocated(triplets_shifts)) allocate(triplets_shifts(6,0))

     !! make triplets
     do pos_ij=1,size(map_one_two)     ! loop over pairs (i,j) of 1 and 2
       i=pairs_one_two( 1, map_one_two(pos_ij) )
       j=pairs_one_two( 2, map_one_two(pos_ij) )
       diff_ij = diffs_one_two(:,pos_ij)
       shift_ij = shifts_one_two(:,pos_ij)
       
       do pos_ik=1,size(map_one_three)                   ! loop over monomers k of type 3 also paired with i 
         if (i /= pairs_one_three( 1, map_one_three(pos_ik) ) ) cycle

         k =pairs_one_three( 2 , map_one_three(pos_ik) )
         diff_ik = diffs_one_three(:,pos_ik)
         shift_ik = shifts_one_three(:,pos_ik)

         diff_jk = diff_ik - diff_ij
         shift_jk = shift_ik - shift_ij

         if (two_three_identical .and. all(shift_jk .eq. 0) .and. j==k ) cycle ! skip if  j and k are same monomer

!         if (two_three_identical .and. k .lt. j) cycle
!         if (two_three_identical .and. all(shift_jk .eq. 0) .and. j .ge. k ) cycle 
!         if (two_three_identical .and. j .ge. k ) cycle 

         ! check that is a valid triplet - i.e. two sides of the triangle are within cutoff length
         ! find minimum  intermolecular distance (excluding Hydrogen) between each monomer pair
         min_distances = (/ min_intermolecular_dist(at,monomer_one_index(:,i),monomer_two_index(:,j)   , diff_ij, pairwise_cutoff) , &
                            min_intermolecular_dist(at,monomer_one_index(:,i),monomer_three_index(:,k) , diff_ik, pairwise_cutoff) , &
                            min_intermolecular_dist(at,monomer_two_index(:,j),monomer_three_index(:,k) , diff_jk, pairwise_cutoff)   /)
         if ( count(min_distances .gt. cutoff) .gt. 1 ) cycle
         call add_triplet_if_new(monomer_triplets,triplets_diffs,triplets_diffs_map,triplets_shifts,i,j,k,diff_ij,diff_ik,shift_ij,shift_ik,two_three_identical)
       end do
     end do
     deallocate(pairs_one_two)
     deallocate(shifts_one_two)
     deallocate(diffs_one_two)
     deallocate(pairs_one_three)
     deallocate(shifts_one_three)
     deallocate(diffs_one_three)
     if(allocated(triplets_shifts)) deallocate(triplets_shifts)
!call print("triplets diffs")
!call print(triplets_diffs)
end subroutine find_monomer_triplets  


subroutine add_triplet_if_new(monomer_triplets,triplets_diffs,triplets_diffs_map,triplets_shifts,i,j,k,diff_ij,diff_ik,shift_ij,shift_ik,two_three_identical)
     ! this also excludes the case where the and three identical and this pair was already found in the opposite order
     integer, dimension(:,:), intent(inout), allocatable :: monomer_triplets
     real(dp), dimension(:,:), intent(inout), allocatable :: triplets_diffs
     integer, dimension(:), intent(inout), allocatable :: triplets_diffs_map
     integer, dimension(:,:), intent(inout), allocatable :: triplets_shifts
     integer, intent(in) :: i,j,k
     real(dp), dimension(3),intent(in) :: diff_ij,diff_ik
     integer, dimension(3),intent(in) :: shift_ij,shift_ik
     logical, intent(in) :: two_three_identical

     integer, dimension(3) ::  triplet
     logical, dimension(:), allocatable :: shift_mask
     integer, dimension(1):: temp1d
     integer :: pos_ijk, pos_ikj,n_triplets,n_shifts
     logical :: do_append

      if (.not. allocated(monomer_triplets)) allocate(monomer_triplets(3,0))
      if (.not. allocated(triplets_diffs)) allocate(triplets_diffs(6,0))
      if (.not. allocated(triplets_diffs_map)) allocate(triplets_diffs_map(0))
      if (.not. allocated(triplets_shifts)) allocate(triplets_shifts(6,0))

      n_triplets=size(monomer_triplets,2)
      n_shifts=size(triplets_diffs_map)
      triplet=(/i,j,k/)
      allocate(shift_mask(n_shifts))
      !call print("triplet has i, j,k = "//(/i,j,k/)//" and shifts "//(/shift_ij,shift_ik/))

      ! in case monomers two and three are identical, check this pair hasn't already been found in opposite order
      if (two_three_identical) then
        shift_mask =  (/  triplets_shifts(1,:) .eq. shift_ik(1)  /)   .and. &
                      (/  triplets_shifts(2,:) .eq. shift_ik(2)  /)   .and. &
                      (/  triplets_shifts(3,:) .eq. shift_ik(3)  /)   .and. &         

                      (/  triplets_shifts(4,:) .eq. shift_ij(1)  /)   .and. &
                      (/  triplets_shifts(5,:) .eq. shift_ij(2)  /)   .and. &
                      (/  triplets_shifts(6,:) .eq. shift_ij(3)  /) 

        temp1d = maxloc(monomer_triplets(1,:), monomer_triplets(1,:) .eq. i .and. monomer_triplets(2,:) .eq. k .and. monomer_triplets(3,:) .eq. j)
        pos_ikj = temp1d(1)
        temp1d=maxloc(triplets_diffs_map,triplets_diffs_map .eq. pos_ikj .and. shift_mask)
        if (temp1d(1) /= 0) then 
          deallocate(shift_mask)
          return
        end if
      end if

      ! find where this triplet is if we've found it before
      temp1d = maxloc(monomer_triplets(1,:), monomer_triplets(1,:) .eq. i .and. monomer_triplets(2,:) .eq. j .and. monomer_triplets(3,:) .eq. k)
      pos_ijk = temp1d(1)

      if (pos_ijk == 0) then ! append triplet if it's new
        n_triplets = n_triplets + 1
        call reallocate(monomer_triplets,3,n_triplets,copy=.true.)
        monomer_triplets(:,n_triplets) = triplet
        pos_ijk=n_triplets      
      end if

      shift_mask =  (/  triplets_shifts(1,:) .eq. shift_ij(1)  /)   .and. &       ! check if this shift already present for this triplet
                    (/  triplets_shifts(2,:) .eq. shift_ij(2)  /)   .and. &
                    (/  triplets_shifts(3,:) .eq. shift_ij(3)  /)   .and. &         

                    (/  triplets_shifts(4,:) .eq. shift_ik(1)  /)   .and. &
                    (/  triplets_shifts(5,:) .eq. shift_ik(2)  /)   .and. &
                    (/  triplets_shifts(6,:) .eq. shift_ik(3)  /) 

      temp1d=maxloc(triplets_diffs_map,triplets_diffs_map .eq. pos_ijk .and. shift_mask)
      deallocate(shift_mask)
      do_append = (temp1d(1) == 0)

      if (do_append) then  ! append shift if it's new
        n_shifts = n_shifts + 1
        call reallocate(triplets_shifts,6,n_shifts,copy=.true.)
        call reallocate(triplets_diffs,6,n_shifts,copy=.true.)
        call reallocate(triplets_diffs_map,n_shifts,copy=.true.)
        triplets_shifts(:,n_shifts) = (/shift_ij,shift_ik/)
        triplets_diffs(:,n_shifts) = (/diff_ij,diff_ik/)
        triplets_diffs_map(n_shifts) = pos_ijk
      end if      
end subroutine add_triplet_if_new

function calc_mean_pos(at,indices) result (pos)
! generates an absolute 'mean' position of a set of atoms - NB not the centre of mass
  type(Atoms), intent(in) :: at
  integer,dimension(:), intent(in) :: indices
  real(dp),dimension(3) :: pos
  integer :: i,n

  pos = 0.0_dp
  n=size(indices)

  do i =2,n
    pos = pos + diff_min_image(at,indices(1),indices(i))
  end do
  pos = pos / n
  pos = pos + at%pos(:,indices(1))

end function calc_mean_pos

function min_intermolecular_dist(at,mono1,mono2,mean_pos_diff,r_init) result (r_min)
  ! returns the minimum intermolecular distance between two monomers, excluding hydrogens
  type(Atoms), intent(in) :: at
  integer, dimension(:), intent(in) :: mono1,mono2
  real(dp), dimension(3), intent(in) :: mean_pos_diff
  real(dp), intent(in) :: r_init

  real(dp) :: r_min, temp_dist
  integer :: i_atomic, i_glob,j_atomic, j_glob
  real(dp), dimension(3) :: mean_pos_one, mean_pos_two

  mean_pos_one = calc_mean_pos(at,mono1)
  mean_pos_two = calc_mean_pos(at,mono2)

  r_min = r_init

  do i_atomic=1,size(mono1)
    i_glob=mono1(i_atomic)
    if (at%Z(i_glob) == 1) cycle

    do j_atomic=1,size(mono2)
      j_glob = mono2(j_atomic)
      if (at%Z(j_glob) == 1) cycle
      temp_dist = norm( diff_min_image(at,at%pos(:,i_glob),mean_pos_one) + mean_pos_diff + diff_min_image(at,mean_pos_two,at%pos(:,j_glob)) )
      r_min = min(r_min,temp_dist)
    end do
  end do

end function min_intermolecular_dist

end module topology_module
