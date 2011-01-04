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

! last modified -- 2009-02-19 -- Csilla
! takes an xyz file, identifies the residues and outputs a PSF and a PDB file.
#include "error.inc"

program xyz2pdb

!  use libatoms_module

  use atoms_module,            only: atoms, finalise, &
                                     set_cutoff, &
                                     calc_connect, print, &
                                     distance_min_image, bond_length, &
                                     map_into_cell, assign_pointer, add_property, &
				     has_property, atoms_sort
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
                                     PRINT_SILENT, PRINT_NORMAL, PRINT_ANAL, PRINT_NERD
  use table_module,            only: table, finalise, int_part, delete
  use topology_module,         only: create_residue_labels_arb_pos, delete_metal_connects, &
                                     write_brookhaven_pdb_file, &
                                     write_psf_file, &
                                     MM_RUN, silica_2body_cutoff
  use cinoutput_module,        only : read, write
  use error_module

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
    integer                     :: center_atom, i
    real(dp)                    :: shift(3)
    logical                     :: ex
    logical                     :: have_silica_potential
    logical                     :: use_avgpos
    integer :: at_i
    integer, pointer :: sort_index_p(:)
    logical                     :: do_sort
    integer :: error = ERROR_NONE
    character(80) :: sorted_name

    call system_initialise(verbosity=PRINT_silent,enable_timing=.true.)
    call verbosity_push(PRINT_NORMAL)
    call system_timer('program')

   ! reading in run parameters
    call initialise(params_in)
    call param_register(params_in, 'File', PARAM_MANDATORY, xyz_file, help_string="Input XYZ file containing species:pos or species:avgpos")
    call param_register(params_in, 'Residue_Library', PARAM_MANDATORY, Library, help_string="Residue library for pattern matching, e.g. $QUIP_DIR/libAtoms/protein_res.CHARMM.lib")
    call param_register(params_in, 'Neighbour_Tolerance', '1.2', Neighbour_Tolerance, help_string="Bond criterion tolerance factor.  Two atoms are bonded if d(a1,a2)<[CovRad(a1)+CovRad(a2)]*Neighbour_Tolerance.  This should normally not need changing.")
    call param_register(params_in, 'Delete_Metal_Connections', 'T', Delete_Metal_Connections, help_string="Whether not to take into account bonds of non(C,H,O,N,P) elements.")
    call param_register(params_in, 'Print_XSC', 'F', print_xsc, help_string="Whether to output XSC (NAMD cell file).")
    call param_register(params_in, 'Center_atom', '0', center_atom, help_string="Which atom should the system be centered around.  0 means no centering.  This is important for programs that can only treat nonperiodic systems not to cut the molecule in half.")
    call param_register(params_in, 'have_silica_potential', 'F', have_silica_potential, help_string="Whether there is a silica unit in the system to be treated with the Danny potential in CP2K.")
    call param_register(params_in, 'use_avgpos', 'T', use_avgpos, help_string="Whether to use the average positions (avgpos) for the pattern matching rather than the instantaneous positions (pos).")
    call param_register(params_in, 'sort', 'F', do_sort, help_string="Whether to sort atoms by molecule and/or residue.")
    if (.not. param_read_args(params_in)) then
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
    sorted_name = trim(Root)//'.sorted.xyz'

   ! print run parameters
    call print('Input:')
    call print('  XYZ file: '//trim(xyz_file))
    call print('  Residue Library: '//trim(Library))
    call print('  Print_XSC (NAMD cell file): '//print_xsc)
    call print('  Delete_Metal_Connections'//Delete_Metal_Connections)
    call print('  Neighbour Tolerance: '//Neighbour_Tolerance)
    call print('  have_silica_potential: '//have_silica_potential)
    call print('  use_avgpos to build connectivity: '//use_avgpos)
    call print('  sort atoms by molecule/residue : '//use_avgpos)
    if(center_atom > 0) &
         call print('  Center atom: '//center_atom)
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
    call read(my_atoms,trim(xyz_file))


   ! calculating connectivities
    call print('Calculating connectivities...')
   ! use nonuniform connect_cutoff, include only nearest neighbours, otherwise the adjacency matrix will contain disjoint regions
    if (have_silica_potential) then
       call set_cutoff(my_atoms,silica_2body_cutoff)
    else
       call set_cutoff(my_atoms,0.0_dp)
    endif
    my_atoms%nneightol = Neighbour_Tolerance
    call calc_connect(my_atoms)
   ! remove bonds for metal ions - everything but H, C, N, O, Si, P, S
    if (Delete_Metal_Connections) call delete_metal_connects(my_atoms)
!    call print(my_atoms%connect)

    ! center atoms if needed
    if(center_atom > 0) then
       shift = my_atoms%pos(:,center_atom)
       do i=1,my_atoms%N
          my_atoms%pos(:,i) = my_atoms%pos(:,i)-shift
       end do
    end if
    call map_into_cell(my_atoms)

   ! identify residues
    call print('Identifying residues...')
    call set_value(my_atoms%params,'Library',trim(Library))
    if (use_avgpos) then
       call create_residue_labels_arb_pos(my_atoms,do_CHARMM=.true.,intrares_impropers=intrares_impropers, &
	have_silica_potential=have_silica_potential,pos_field_for_connectivity="avgpos")
    else !use actual positions
       call create_residue_labels_arb_pos(my_atoms,do_CHARMM=.true.,intrares_impropers=intrares_impropers, &
	have_silica_potential=have_silica_potential,pos_field_for_connectivity="pos")
    endif

    if (do_sort) then
       ! sort by molecule, residue ID
       nullify(sort_index_p)
       if (.not. assign_pointer(my_atoms, 'sort_index', sort_index_p)) then
	 call add_property(my_atoms, 'sort_index', 0)
	 if (.not. assign_pointer(my_atoms, 'sort_index', sort_index_p)) &
	   call system_abort("WARNING: do_cp2k_calc failed to assign pointer for sort_index, not sorting")
       endif
       if (associated(sort_index_p)) then
	 do at_i=1, my_atoms%N
	   sort_index_p(at_i) = at_i
	 end do
       endif
       if (.not.(has_property(my_atoms,'mol_id')) .or. .not. has_property(my_atoms,'atom_res_number')) then
	 call system_abort("WARNING: can't do sort_by_molecule - need mol_id and atom_res_number")
       else
	 call atoms_sort(my_atoms, 'mol_id', 'atom_res_number', error=error)
	 HANDLE_ERROR(error)
	 if (associated(sort_index_p)) then
	   do at_i=1, my_atoms%N
	     if (sort_index_p(at_i) /= at_i) then
	       call print("sort() of my_atoms%data by mol_id, atom_res_number reordered some atoms")
	       exit
	     endif
	   end do
	 endif
       end if
       call calc_connect(my_atoms)
    end if ! do_sort

   ! print output PDB and PSF files
    call print('Writing files with CHARMM format...')
    call write_psf_file(my_atoms,psf_file=trim(psf_name),intrares_impropers=intrares_impropers,add_silica_23body=have_silica_potential)
    call write_brookhaven_pdb_file(my_atoms,trim(pdb_name))
    if (Print_XSC) call write_xsc_file(my_atoms,xsc_file=trim(xsc_name))
    if (do_sort) then
      call write(my_atoms, trim(sorted_name))
    endif

    call finalise(intrares_impropers)
    call finalise(my_atoms)
    call print ('Finished.')

    call system_timer('program')
    call system_finalise

contains

  subroutine print_usage

    call print('Usage: xyz2pdb File=filename.xyz [Residue_Library=library] [Print_XSC] [Delete_Metal_Connections=T] [Neighbour_tolerance=1.2] [Center_atom=num]')
    call print('')
    call print('  File=filename,        where your input file has extension .xyz')
    call print('  [Residue_Library=library],  optional, default is protein_res.CHARMM.lib')
    call print('  [Print_XSC],          optional, whether to print NAMD cell file, default is false')
    call print('  [Delete_Metal_Connections=T],  optional, default is true, only calculates connection for H,C,N,O,Si,P,S,Cl')
    call print('                        should work fine - only modify if you want bonds with other elements in your structure')
    call print('  [Neighbour_Tolerance=1.2],  optional, default is 1.2, should work fine - do not poke it ')
    call print('                        unless you know you have unusally long bonds in your structure')
    call print('  [Center_atom=num], optionally shifts atoms so that the atom with index num is at')
    call print('                        the origin, before dropping periodic boundary condition info')
    call print('  [use_avgpos=T], optionally use instantaneous positions instead of average positions, if this is set to F')
    call print('  [sort=F], optionally sort atoms by residue/molecule')
    call print('  [have_silica_potential=F], optionally switches on the use of silica potential described by Cole et al., JChemPhys, 2007.')

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
