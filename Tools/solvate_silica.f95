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

! solvate a molecule. Needs a solute and a solvent XYZ input file.
! shifts and rotates the solvent file so that
!   atom1_point has position atom1_pos
!   atom1_point -> atom2_vector is` aligned with atom12_vector
!   atom3 is in the atom1_pos - atom1_pos+atom12_vector - atom123_plane plane.
! don't include water molecules closer to the solute than exclusion
! rotation_sign to play around with

program solvate_silica

  use libatoms_module

  implicit none

  type(Atoms)                 :: at, wat, at2
  type(dictionary)            :: cli_params
  character(len=STRING_LENGTH) :: filename, waterfilename
  integer                     :: i, j, stat

  !solvation of solvent
  real(dp)                    :: exclusion
  integer                     :: silica_atoms
  logical                     :: add

  !orientation of solvent
  integer                     :: atom1_point, atom2_vector, atom3_plane
  real(dp)                    :: atom1_pos(3), atom12_vector(3), atom123_plane_point(3)
  integer                     :: rotation_sign1
  integer                     :: rotation_sign2
  real(dp)                    :: shift(3)
  real(dp)                    :: vector_1(3), vector_2(3)
  real(dp)                    :: theta, axis(3), origin(3)
  real(dp)                    :: unit_vector12(3), projection_of_atom3_on_atom12_vector(3)
  real(dp)                    :: surface_normal(3)
  logical                     :: rotation2_pluspi

  call system_initialise(verbosity=PRINT_SILENT)
!  call verbosity_push(PRINT_NORMAL)

  call initialise(cli_params)
  call param_register(cli_params,"file","stdin", filename, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params,"waterfile","stdin", waterfilename, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params,"silica_atoms",'0', silica_atoms, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params,"atom1_point",PARAM_MANDATORY, atom1_point, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params,"atom2_vector",PARAM_MANDATORY, atom2_vector, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params,"atom3_plane",PARAM_MANDATORY, atom3_plane, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params, 'atom1_pos', '(/0.0 0.0 0.0/)', atom1_pos, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params, 'atom12_vector', PARAM_MANDATORY, atom12_vector, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params, 'atom123_plane_point', PARAM_MANDATORY, atom123_plane_point, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params,"exclusion","0.0", exclusion, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params,"rotation1_sign",'1', rotation_sign1, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params,"rotation2_sign",'1', rotation_sign2, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params,"rotation2_pluspi",'F', rotation2_pluspi, help_string="No help yet.  This source file was $LastChangedBy$")
  if (.not. param_read_args(cli_params)) then
    !call system_abort("Usage: decimate [file=(stdin)] [waterfile=(stdin)] silica_atoms atom1_point atom2_vector atom3_plane atom1_pos atom12_vector atom123_plane_point exclusion rotation1_sign rotation2_sign2 rotation2_pluspi")
    call verbosity_push(PRINT_NORMAL)
    call print_usage
    call system_abort('could not parse argument line')
  endif
  call finalise(cli_params)

  if (abs(rotation_sign1) /= 1) call system_abort("rotation1_sign must be +/-1")
  if (abs(rotation_sign2) /= 1) call system_abort("rotation2_sign must be +/-1")

  !Read in the 2 files
  call read(at, filename, error=stat)
  if (stat /= 0) call system_abort('Could not read file to solvate.')
  call map_into_cell(at)
  call read(wat, waterfilename, error=stat)
  if (stat /= 0) call system_abort('Could not read water file.')
  call map_into_cell(wat)

  !Initialise the solvated atoms object
  call initialise(at2,0,reshape(reshape(wat%lattice, (/9/)), (/3,3/)))

!!!position and orientation of solvent

    !shift atom1_point to the origin (centre of rotation)

       shift = at%pos(1:3,atom1_point)
       do i=1,at%N
          at%pos(1:3,i) = at%pos(1:3,i) - shift(1:3)
       enddo
       call calc_connect(at)
       call map_into_cell(at)
       call calc_dists(at)

    !rotate solvent to align atom1->atom2 vector with atom12_vector

       !vectors we want to lay on top of each other
       vector_1(1:3) = diff_min_image(at,atom1_point,atom2_vector) !pointing from atom1 -> atom2
       vector_2(1:3) = atom12_vector(1:3) !pointing from atom1 -> atom2, with arbitrary length
     
       call print('vector_1('//atom1_point//'->'//atom2_vector//') = '//vector_1(1:3))
       call print('vector_2('//atom1_point//'->'//atom2_vector//') = '//vector_2(1:3))
     
       !calc. angle between the two vectors
       theta = acos( dot_product(vector_1(1:3),vector_2(1:3)) / &
                     sqrt(dot_product(vector_1(1:3),vector_1(1:3))) / &
                     sqrt(dot_product(vector_2(1:3),vector_2(1:3))) )
     
       call print('angle = '//theta)
     
       !calc. axis that is normal to both vectors
       axis(1:3) = vector_1(1:3) .cross. vector_2(1:3)
     
       call print('axis = '//axis)
     
       origin(1:3) = 0._dp
     
       !rotate only the 2nd config
       call rotate(at%pos, real(rotation_sign1,dp)*axis, theta, origin)

call print('atom1--atom2 vector should be all right:')
call write(at, 'stdout')

    !rotate solvent around atom12_vector axis to move atom3_plane into the plane that is defined by
    !    - origin(1:3) = (/0.,0.,0./),
    !    - atom12_vector(1:3) and
    !    - atom123_planepoint(1:3)-atom1_point(1:3)

       !projection of atom3 onto atom12_vector
       !  (at1->at3 .dot. unit_vector) * unit_vector  : the unit vector in atom12_vector direction
       unit_vector12 = atom12_vector(1:3) / norm(atom12_vector(1:3))
       projection_of_atom3_on_atom12_vector(1:3) = &
            ( (at%pos(1:3,atom3_plane)-at%pos(1:3,atom1_point)) .dot. unit_vector12(1:3) ) * unit_vector12(1:3)

call print("Projection of atom3 to atom12 vector: "//projection_of_atom3_on_atom12_vector(1:3))

       !we want to rotate this
       vector_1(1:3) = at%pos(1:3,atom3_plane) - projection_of_atom3_on_atom12_vector(1:3)
       !onto vector_2, i.e. vector from projection-point, _|_ to vector12, on the plane
         surface_normal(1:3) = ( atom123_plane_point(1:3) - at%pos(1:3,atom1_point) ) .cross. atom12_vector(1:3)
       vector_2(1:3) = surface_normal(1:3) .cross. atom12_vector(1:3)
!OR (-1)* the same __/\
       vector_2(1:3) = vector_2(1:3) / norm(vector_2(1:3)) * norm(vector_1(1:3))
     
       call print('vector_1() = '//vector_1(1:3))
       call print('vector_2() = '//vector_2(1:3))
     
       !calc. angle between the two vectors
       theta = acos( dot_product(vector_1(1:3),vector_2(1:3)) / &
                     sqrt(dot_product(vector_1(1:3),vector_1(1:3))) / &
                     sqrt(dot_product(vector_2(1:3),vector_2(1:3))) )
     
       call print('angle = '//theta)
     
       !calc. axis that is normal to both vectors
       axis(1:3) = vector_1(1:3) .cross. vector_2(1:3)
     
       call print('axis = '//axis)
       call print('atom12_vector should be the same = '//atom12_vector)
     
       origin(1:3) = 0._dp
     
       !rotate only the 2nd config
       if (rotation2_pluspi) theta = theta + PI
       call rotate(at%pos, (real(rotation_sign2,dp)) * axis, theta, origin)

call print('atom3 should be in plane with vector12 and atom3_plane_point:')
call print('v13 x v12 . v1planepoint = '//(diff_min_image(at, atom1_point, atom3_plane) .cross. atom12_vector(1:3) .dot. atom123_plane_point(1:3)))
call write(at, 'stdout')

  !shifting it to the required position
     do i=1,at%N
        at%pos(1:3,i) = at%pos(1:3,i) + atom1_pos(1:3)
     enddo

  !Add atoms to solvate to at2
  do i=1,at%N
     call add_atoms(at2,at%pos(1:3,i),at%Z(i))
  enddo

  !add silica
  do i=1,silica_atoms
     call add_atoms(at2,wat%pos(1:3,i  ),wat%Z(i  ))
  enddo
  !Add water molecules between the given limits to at2
  do i=silica_atoms+1,wat%N,3
     add = .true.
     do j = 1,at%N
        if (distance_min_image(at2,j,wat%pos(1:3,i  )).lt.exclusion .or. &
            distance_min_image(at2,j,wat%pos(1:3,i+1)).lt.exclusion .or. &
            distance_min_image(at2,j,wat%pos(1:3,i+2)).lt.exclusion ) then
           add = .false.
        endif
     enddo
     if (add) then
        call add_atoms(at2,wat%pos(1:3,i  ),wat%Z(i  ))
        call add_atoms(at2,wat%pos(1:3,i+1),wat%Z(i+1))
        call add_atoms(at2,wat%pos(1:3,i+2),wat%Z(i+2))
     endif
  enddo

  call map_into_cell(at2)
  call verbosity_push(PRINT_NORMAL)
  if (stat == 0) call write(at2, 'stdout')

  call verbosity_pop()
  call system_finalise()

contains

    !call system_abort("Usage: decimate [file=(stdin)] [waterfile=(stdin)] silica_atoms atom1_point atom2_vector atom3_plane atom1_pos atom12_vector atom123_plane_point exclusion rotation1_sign rotation2_sign2 rotation2_pluspi")
  subroutine print_usage

    call print("Usage: decimate [file=(stdin)] [waterfile=(stdin)] silica_atoms atom1_point atom2_vector atom3_plane atom1_pos atom12_vector atom123_plane_point exclusion rotation1_sign rotation2_sign2 rotation2_pluspi")
    call print('')
    call print('  file=filename,           the solute file')
    call print('  waterfile=filename,      the solvent (containing silica and water) file, with cell size at least the size of the solute file the size of the solvated molecule will be the same as the cell of the waterfile')
    call print('  silica_atoms,            the number of silica atoms in the solvent file')
    call print('  atom1_point,             the solute atom with this index to set (atom1_point) coordinate in the solvent file')
    call print('  atom2_vector,            the solute atom with this index to set (atom1_point,atom2_vector) direction in the solvent file')
    call print('  atom3_plane,             the solute atom with this index to set (atom1_point,atom2_vector,atom3_plane) plane in the solvent file')
    call print('  atom1_pos,               the position of atom1_point in the solvated file')
    call print('  atom12_vector,           the direction of atom1_point -> atom2_vector vector in the solvated file')
    call print('  atom123_plane point,     a point on the plane that contains atom1_point, atom2_vector and atom3_plane in the solvated file')
    call print('  [exclusion=0.],          optionally keeps a layer around the solvent where there will be no solvent molecule (valid for water, not silica)')
    call print('  [rotation1_sign=1],      optionally rotate in the negative direction, when rotating the solute molecule to align atom1&2 with atom12_vector')
    call print('  [rotation2_sign=1],      optionally rotate in the negative direction, when rotating the solute molecule to rotate atom3_plane into the (atom1_pos,atom12_vector,atom123_plane_point) plane')
    call print('  [rotation2_pluspi=F],    optionally rotate pi more, to have atom3_plane on the other side of (atom1_point,atom2_vector) vector')

    call print('')

  end subroutine

end program solvate_silica
