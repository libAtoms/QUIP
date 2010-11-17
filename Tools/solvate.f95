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

! solvate a molecule. Needs a solute and a solvent XYZ input file (only water).
! [x/y/z][min/max] specify the new lattice, which must not be smaller than the solute's lattice and larger than the solvate lattice.
! could be improved: e.g. shift the water file to cover *min -- *max, at the moment it works only if centred around the origin,
!                         multiply water file to cover any box size,
!                         check whether water molecules touch any other added water molecules.
! center_around_atom : to center the solvent around atom and shift that atom to the origin before mapping into cell

program solvate

  use libatoms_module

  implicit none

  type(Atoms)                 :: at, wat, at2
  integer                     :: i, j, stat
  type(dictionary)            :: cli_params
  character(len=FIELD_LENGTH) :: filename, waterfilename, solvated_filename
!  type(inoutput)              :: xyzfile, waterfile
  real(dp)                    :: xmin, xmax, ymin, ymax, zmin, zmax, exclusion, security_zone
  logical                     :: add
  integer                     :: center_around_atom
  real(dp)                    :: shift(3)
  real(dp), pointer           :: pos_p(:,:), avgpos_p(:,:)

  call system_initialise(verbosity=PRINT_SILENT)
  call verbosity_push(PRINT_NORMAL)

  call initialise(cli_params)
  call param_register(cli_params,"file","stdin", filename, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params,"waterfile","stdin", waterfilename, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params,"outfile","stdout", solvated_filename, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params,"xmin",PARAM_MANDATORY, xmin, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params,"ymin",PARAM_MANDATORY, ymin, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params,"zmin",PARAM_MANDATORY, zmin, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params,"xmax",PARAM_MANDATORY, xmax, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params,"ymax",PARAM_MANDATORY, ymax, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params,"zmax",PARAM_MANDATORY, zmax, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params,"exclusion","2.4_dp", exclusion, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params,"center_around_atom","0", center_around_atom, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(cli_params,"security_zone","1._dp", security_zone, help_string="No help yet.  This source file was $LastChangedBy$")
  if (.not. param_read_args(cli_params)) then
    !call system_abort("Usage: decimate [file=(stdin)] [waterfile=(stdin)] xmin xmax ymin ymax zmin zmax [exclusion=(2.4)] [security_zone=(1.0)] [center_around_atom=0]")
    call print_usage
    call system_abort('could not parse argument line')
  endif
  call finalise(cli_params)

  !Read in the 2 files
  call print("Reading file to solvate: "//trim(filename))
  call read(at, filename, error=stat)
  if (stat /= 0) call system_abort('Could not read file to solvate.')
  !center around given atom, if required
  if (center_around_atom > 0) then
     if (center_around_atom > at%N) call system_abort("center_around_atom must be less than the number of solvent atoms "//at%N)
     shift = at%pos(1:3,center_around_atom)
     do i=1,at%N
        at%pos(1:3,i) = at%pos(1:3,i) - shift(1:3)
     enddo
  endif
  call map_into_cell(at)

  call print("Reading waterfile: "//trim(waterfilename))
  call read(wat, waterfilename, error=stat)
  if (stat /= 0) call system_abort('Could not read water file.')
  call map_into_cell(wat)

  !check for cell size -- if the water file is at least as big as the solute.
  do i=1,3
     do j=1,3
        if (at%lattice(i,j).gt.wat%lattice(i,j)) call system_abort('Too small water file')
     enddo
  enddo

  call print("Starting solvation.")
  !Initialise the solvated atoms object
  call initialise(at2,0,reshape((/xmax-xmin,0._dp,0._dp,0._dp,ymax-ymin,0._dp,0._dp,0._dp,zmax-zmin/),(/3,3/)))

  !Add solute atoms to at2
  do i=1,at%N
     call add_atoms(at2,at%pos(1:3,i),at%Z(i))
  enddo

  !Add water molecules between the given limits to at2
  do i=1,wat%N,3
     add = .false.
     if (minval(wat%pos(1,i:i+2)).gt.xmin+security_zone .and. &
         maxval(wat%pos(1,i:i+2)).lt.xmax-security_zone .and. &
         minval(wat%pos(2,i:i+2)).gt.ymin+security_zone .and. &
         maxval(wat%pos(2,i:i+2)).lt.ymax-security_zone .and. &
         minval(wat%pos(3,i:i+2)).gt.zmin+security_zone .and. &
         maxval(wat%pos(3,i:i+2)).lt.zmax-security_zone ) then
        add = .true.
     endif
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

  !add avgpos property that can be used for the topology calc.
  call add_property(at2,'avgpos',0._dp, n_cols=3)
  if (.not.(assign_pointer(at2, "avgpos", avgpos_p))) call system_abort('??')
  if (.not.(assign_pointer(at2, "pos", pos_p))) call system_abort('??')
  avgpos_p(1:3,1:at2%N) = pos_p(1:3,1:at2%N)
  
  if (stat == 0) then
     call print("Printing solvated file into "//trim(solvated_filename))
     call write(at2, solvated_filename)
  endif

  call verbosity_pop()
  call system_finalise()

contains

    !call system_abort("Usage: decimate [file=(stdin)] [waterfile=(stdin)] xmin xmax ymin ymax zmin zmax [exclusion=(2.4)] [security_zone=(1.0)] [center_around_atom=0]")
  subroutine print_usage

    call print("Usage: decimate [file=(stdin)] [waterfile=(stdin)] xmin xmax ymin ymax zmin zmax [exclusion=(2.4)] [security_zone=(1.0)] [center_around_atom=0]")
    call print('')
    call print('  file=filename,           the solute file')
    call print('  waterfile=filename,      the solvent (only water) file, with cell size at least the size of the solute file')
    call print('  xmin,                    the lower limit of the solvated file in the x direction (should be >= 0.5*solvent file cell edge in x direction)')
    call print('  xmax,                    the lower limit of the solvated file in the y direction (should be <= 0.5*solvent file cell edge in x direction)')
    call print('  ymin,                    the lower limit of the solvated file in the z direction (should be >= 0.5*solvent file cell edge in y direction)')
    call print('  ymax,                    the lower limit of the solvated file in the x direction (should be <= 0.5*solvent file cell edge in y direction)')
    call print('  zmin,                    the lower limit of the solvated file in the y direction (should be >= 0.5*solvent file cell edge in z direction)')
    call print('  zmax,                    the lower limit of the solvated file in the z direction (should be <= 0.5*solvent file cell edge in z direction)')
    call print('  [exclusion=0.],          optionally keeps a layer around the solvent where there will be no solvent molecule')
    call print('  [security_zone=1.],      optionally keeps a layer on the edge of the solvated box where there will be no water molecule placed to avoid clashes in PBC')
    call print('  [center_around_atom=1],  optionally shifts the solute molecule so that the atom with this index will be at the origin')

    call print('')

  end subroutine

end program solvate
