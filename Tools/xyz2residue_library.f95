program xyz2residue_library
use libatoms_module
implicit none

   type(Atoms) :: at, m_at
   integer i_at, j_at, nn_ii, nn_i
   integer, pointer :: mol_id(:)
   character(len=10) :: specie

   call system_initialise(verbosity=PRINT_SILENT)
   call verbosity_push(PRINT_NORMAL)

   call read(at, "stdin")

   call calc_connect(at)
   call find_molecule_ids(at)

   if (.not. assign_pointer(at, "mol_id", mol_id)) &
      call system_abort("Failed to find mol_id property")

   do i_at=1, at%N
      if (mol_id(i_at) /= 0) then
	 call select(m_at, at, mask=(mol_id(:) == mol_id(i_at)))
	 call calc_connect(m_at)

	 call print("%residue X"//count(mol_id == mol_id(i_at))//" R"//mol_id(i_at)//" R"//mol_id(i_at))
	 do j_at=1, m_at%N
	    line = "0 "//m_at%Z(j_at)
	    if (atoms_n_neighbours(m_at, j_at) > 6) &
	       call system_abort("Too many neighbours "//atoms_n_neighbours(m_at, j_at)//" for atom "//j_at)
	    do nn_ii=1, atoms_n_neighbours(m_at, j_at)
	       nn_i = atoms_neighbour(m_at, j_at, nn_ii)
	       line = trim(line)//" "//nn_i
	    end do
	    do nn_ii=atoms_n_neighbours(m_at, j_at)+1, 6
	       line = trim(line)//" "//0
	    end do
	    specie = string_cat_string_array(" ",m_at%species(:,j_at))
	    line = trim(line)//" "//trim(specie)//"  0.0 "//trim(specie)
	    call print(trim(line))
	 end do
	 call print("")

      endif
      where (mol_id == mol_id(i_at))
	 mol_id = 0
      end where
   end do

   call system_finalise()
end program

