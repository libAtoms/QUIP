program test_surf_cell
use libatoms_module
implicit none

type(Atoms) :: at_prim, at_surf, at_slab, at_slab_cut

real(dp) :: surf_v(3), tol
logical :: third_vec_normal
integer :: max_n

integer :: i
integer :: surf_i(3,3)
type(Dictionary) :: cli_params
real(dp) :: F(3,3), old_axes(3,3), old_axes_inv(3,3), new_axes(3,3), z_min, z_max, z_vacuum
integer :: n3_min, n3_max
character(len=STRING_LENGTH) :: comment, Z_order
character(len=STRING_LENGTH) :: verbosity

call system_initialise(verbosity=PRINT_SILENT)

call initialise(cli_params)
call param_register(cli_params, "surf_v", PARAM_MANDATORY, surf_v, help_string="vector normal to surface in Cartesian coords")
call param_register(cli_params, "z_min", PARAM_MANDATORY, z_min, help_string="minimum Z position to include in slab")
call param_register(cli_params, "z_max", PARAM_MANDATORY, z_max, help_string="maximum Z position to include in slab")
call param_register(cli_params, "z_vacuum", PARAM_MANDATORY, z_vacuum, help_string="amount of vacuum to create")
call param_register(cli_params, "third_vec_normal", "T", third_vec_normal, help_string="if true force the third vector to be normal to the surface")
call param_register(cli_params, "tol", "1.0e-6", tol, help_string="tolerance for various normal/parallel tests")
call param_register(cli_params, "max_n", "4", max_n, help_string="maximum supercell component to search for matching vectors")
call param_register(cli_params, "verbosity", "NORMAL", verbosity, help_string="verbosity level.  VERBOSE for more output about surface cell")
if (.not. param_read_args(cli_params)) then
   call param_print_help(cli_params)
   call system_abort("Error reading params from command line")
endif
call finalise(cli_params)

call verbosity_push(verbosity_of_str(trim(verbosity)))

call read_vasp(at_prim, "stdin")

call print("at_prim", verbosity=PRINT_VERBOSE)
if (current_verbosity() >= PRINT_VERBOSE) call write(at_prim, "stdout")

call print("surface unit cell for lattice", verbosity=PRINT_VERBOSE)
call print(at_prim%lattice, verbosity=PRINT_VERBOSE)
call print ("normal to "//surf_v, verbosity=PRINT_VERBOSE)

! create the surface unit cell lattice
call surface_unit_cell(surf_i, surf_v, at_prim%lattice, third_vec_normal=third_vec_normal, tol=tol, max_n=max_n)

call print("", verbosity=PRINT_VERBOSE)
call print("in lattice coordinates", verbosity=PRINT_VERBOSE)
call print(surf_i, verbosity=PRINT_VERBOSE)
call print("in cartesian coordinates", verbosity=PRINT_VERBOSE)
call print(matmul(at_prim%lattice,surf_i), verbosity=PRINT_VERBOSE)

! create the actual surface unit cell with atoms
call arbitrary_supercell(at_surf, at_prim, surf_i)

! rotate to desired axes (surface normal to z)
old_axes(:,1) = at_surf%lattice(:,1) ! map to x
old_axes(:,2) = at_surf%lattice(:,1) .cross. at_surf%lattice(:,2) ! surf_v, map to z
old_axes(:,3) = at_surf%lattice(:,1) .cross. (at_surf%lattice(:,1) .cross. at_surf%lattice(:,2)) ! map to y
old_axes(:,1) = old_axes(:,1) / norm(old_axes(:,1))
old_axes(:,2) = old_axes(:,2) / norm(old_axes(:,2))
old_axes(:,3) = old_axes(:,3) / norm(old_axes(:,3))
new_axes(:,1) = (/ 1, 0, 0 /)
new_axes(:,2) = (/ 0, 0, 1 /)
new_axes(:,3) = (/ 0, 1, 0 /)

call print("old_axes", verbosity=PRINT_VERBOSE)
call print(old_axes, verbosity=PRINT_VERBOSE)
call print("new_axes", verbosity=PRINT_VERBOSE)
call print(new_axes, verbosity=PRINT_VERBOSE)

call matrix3x3_inverse(old_axes, old_axes_inv)
F = matmul(new_axes, old_axes_inv)
new_axes = matmul(F, at_surf%lattice)
if (new_axes(3,3) < 0.0_dp) then
   F(3,:) = -F(3,:)
endif

call print("F", verbosity=PRINT_VERBOSE)
call print(F, verbosity=PRINT_VERBOSE)

call print("F.old_axes", verbosity=PRINT_VERBOSE)
call print(matmul(F, old_axes), verbosity=PRINT_VERBOSE)

call set_lattice(at_surf, matmul(F, at_surf%lattice), .true.)

n3_min = floor((z_min-minval(at_surf%pos(3,:)))/at_surf%lattice(3,3))-1
n3_max = floor((z_max-maxval(at_surf%pos(3,:)))/at_surf%lattice(3,3))+2
call print("n3 min max "//n3_min//" "//n3_max, verbosity=PRINT_VERBOSE)
call supercell(at_slab, at_surf, 1, 1, n3_max-n3_min+1)
do i=1, at_slab%N
   at_slab%pos(:,i) = at_slab%pos(:,i) + n3_min*at_surf%lattice(:,3)
end do

call print("at_slab", verbosity=PRINT_VERBOSE)
if (current_verbosity() >= PRINT_VERBOSE) call write(at_slab, "stdout")

call select(at_slab_cut, at_slab, mask=(at_slab%pos(3,:) >= z_min .and. at_slab%pos(3,:) <= z_max))

old_axes = at_slab_cut%lattice
old_axes(:,3) = old_axes(:,3)/norm(old_axes(:,3)) * (z_max-z_min + z_vacuum)
call set_lattice(at_slab_cut, old_axes, .false.)

! print out final slab

call get_param_value(at_prim, "VASP_Comment", comment)
comment = trim(comment)//" surf_v "//surf_v
call set_param_value(at_slab_cut, "VASP_Comment", comment)

call get_param_value(at_prim, "VASP_Z_order", Z_order)
call set_param_value(at_slab_cut, "VASP_Z_order", Z_order)

call print("at_slab_cut", verbosity=PRINT_VERBOSE)
if (current_verbosity() >= PRINT_VERBOSE) call write(at_slab_cut, "stdout")

call write_vasp(at_slab_cut, "stdout", fix_order=.true., cartesian=.false.)

call verbosity_pop()
call system_finalise()

contains

subroutine write_vasp(at, filename, fix_order, cartesian)
   type(Atoms), intent(in) :: at
   character(*), intent(in) :: filename
   logical, intent(in) :: fix_order, cartesian

   type(inoutput) :: io
   integer :: sorted_Zs(at%N)
   integer, allocatable :: uniq_Zs(:)
   integer :: i, j
   character(len=STRING_LENGTH) :: comment, Z_order
   integer :: n_species, i_species, Ns(100)

   call initialise(io, filename)
   call get_param_value(at, "VASP_Comment", comment)

   if (.not. get_value(at%params, 'VASP_Z_order', Z_order)) then
      Z_order = ''
   endif

   call print(trim(comment), file=io, verbosity=PRINT_ALWAYS)
   call print(1.0_dp, file=io, verbosity=PRINT_ALWAYS)
   call print(at%lattice(:,1), file=io, verbosity=PRINT_ALWAYS)
   call print(at%lattice(:,2), file=io, verbosity=PRINT_ALWAYS)
   call print(at%lattice(:,3), file=io, verbosity=PRINT_ALWAYS)
   if (fix_order) then
      sorted_Zs = at%Z
      call sort_array(sorted_Zs)
      call uniq(sorted_Zs, uniq_Zs)
      if (len_trim(Z_order) > 0) then
	 read(unit=Z_order, fmt=*) line, uniq_Zs
      endif

      line = ""
      do i=1, size(uniq_Zs)
	 line=trim(line)//" "//ElementName(uniq_Zs(i))
      end do
      call print(trim(line), file=io, verbosity=PRINT_ALWAYS)

      line = ""
      do i=1, size(uniq_Zs)
	 line=trim(line)//" "//count(at%Z == uniq_Zs(i))
      end do

      call print(trim(line), file=io, verbosity=PRINT_ALWAYS)

      if (cartesian) then
	 call print("Cartesian")
	 do i=1, size(uniq_Zs)
	    do j=1, at%N
	       if (at%Z(j) == uniq_Zs(i)) call print(at%pos(:,j), file=io, verbosity=PRINT_ALWAYS)
	    end do
	 end do
      else
	 call print("Direct")
	 do i=1, size(uniq_Zs)
	    do j=1, at%N
	       if (at%Z(j) == uniq_Zs(i)) call print(matmul(at%g,at%pos(:,j)), file=io, verbosity=PRINT_ALWAYS)
	    end do
	 end do
      endif

      deallocate(uniq_Zs)
   else
      line=ElementName(at%Z(1))
      n_species = 1
      do i=2, at%N
	 if (at%Z(i) /= at%Z(i-1)) then
	    line=trim(line)//" "//ElementName(at%Z(i))
	    n_species = n_species + 1
	 endif
      end do
      call print(trim(line), file=io, verbosity=PRINT_ALWAYS)

      line=""//at%Z(1)
      i_species = 1
      Ns(i_species) = 1
      do i=2, at%N
	 if (at%Z(i) /= at%Z(i-1)) then
	    i_species = i_species + 1
	 endif
	 Ns(i_species) = Ns(i_species) + 1
      end do
      call print(Ns(1:n_species), file=io, verbosity=PRINT_ALWAYS)

      if (cartesian) then
	 call print("Cartesian")
	 do i=1, at%N
	    call print(at%pos(:,i), file=io, verbosity=PRINT_ALWAYS)
	 end do
      else
	 call print("Direct")
	 do i=1, at%N
	    call print(matmul(at%g,at%pos(:,i)), file=io, verbosity=PRINT_ALWAYS)
	 end do
      endif
   end if


   call finalise(io)

end subroutine write_vasp


subroutine read_vasp(at, filename)
   type(Atoms), intent(inout) :: at
   character(*), intent(in) :: filename

   type(Inoutput) :: io
   character(len=STRING_LENGTH) :: comment, Z_order
   character(len=10) :: species(100)
   integer :: n_species, i, j, N_tot
   integer :: Ns(100), species_Zs(100)
   real(dp) :: lattice(3,3)
   real(dp) :: lattice_scale
   logical :: cartesian
   integer :: status

   call initialise(io, filename)
   line=read_line(io); comment=trim(line)
   line=read_line(io); read (unit=line, fmt=*) lattice_scale
   line=read_line(io); read (unit=line, fmt=*) lattice(:,1)
   line=read_line(io); read (unit=line, fmt=*) lattice(:,2)
   line=read_line(io); read (unit=line, fmt=*) lattice(:,3)
   lattice = lattice * lattice_scale
   ! species line is optional
   line=read_line(io)
   call split_string_simple(line, species, n_species, " 	")
   read(unit=species(1), fmt=*, iostat=status) Ns(1)
   if (status == 0) then ! it's a number
      read(unit=line, fmt=*) Ns(1:n_species)
      do i=1, n_species
	 species_Zs(i) = i
      end do
   else ! actually species
      line=read_line(io)
      read(unit=line, fmt=*) Ns(1:n_species)
      do i=1, n_species
	 species_Zs(i) = atomic_number(species(i))
      end do
   endif
   Z_order='S '//species_Zs(1:n_species)
   line=read_line(io)
   if (line(1:1) == 'S') line=read_line(io)
   if (line(1:1) == 'c' .or. line(1:1) == 'C') then 
      cartesian = .true.
   else if (line(1:1) == 'd' .or. line(1:1) == 'D') then 
      cartesian = .false.
   else
      call system_abort ("confused by cartesian/direct line of '"//trim(line)//"'")
   endif

   call initialise(at, sum(Ns(1:n_species)), lattice)
   N_tot = 0
   do i=1, n_species
   do j=1, Ns(i)
      line=read_line(io)
      N_tot = N_tot + 1
      read (unit=line, fmt=*) at%pos(:,N_tot)
      at%Z(N_tot) = species_Zs(i)
   end do
   end do

   call set_atoms(at, at%Z)

   if (.not. cartesian) then
      at%pos = matmul(at%lattice, at%pos)
   endif

   if (len_trim(comment) == 0) then
      call set_param_value(at, "VASP_Comment", "none")
   else
      call set_param_value(at, "VASP_Comment", trim(comment))
   endif

   call set_param_value(at, "VASP_Z_order", trim(Z_order))

   call finalise(io)

end subroutine read_vasp

end program
