! solvate a molecule. Needs a solute and a solvate XYZ input file.
! [x/y/z][min/max] specify the new lattice, which must not be smaller than the solute's lattice and larger than the solvate lattice.
! could be improved: e.g. shift the water file to cover *min -- *max, at the moment it works only if centred around the origin,
!                         multiply water file to cover any box size,
!                         check whether water molecules touch any other added water molecules.

program solvate

  use libatoms_module

  implicit none

  type(Atoms)                 :: at, wat, at2
  integer                     :: i, j, stat
  type(dictionary)            :: cli_params
  character(len=FIELD_LENGTH) :: filename, waterfilename
  type(inoutput)              :: xyzfile, waterfile
  real(dp)                    :: xmin, xmax, ymin, ymax, zmin, zmax, exclusion, security_zone
  logical                     :: add

  call system_initialise(verbosity=SILENT)
  call verbosity_push(NORMAL)

  call initialise(cli_params)
  call param_register(cli_params,"file","stdin", filename)
  call param_register(cli_params,"waterfile","stdin", waterfilename)
  call param_register(cli_params,"xmin",PARAM_MANDATORY, xmin)
  call param_register(cli_params,"ymin",PARAM_MANDATORY, ymin)
  call param_register(cli_params,"zmin",PARAM_MANDATORY, zmin)
  call param_register(cli_params,"xmax",PARAM_MANDATORY, xmax)
  call param_register(cli_params,"ymax",PARAM_MANDATORY, ymax)
  call param_register(cli_params,"zmax",PARAM_MANDATORY, zmax)
  call param_register(cli_params,"exclusion","2.4_dp", exclusion)
  call param_register(cli_params,"security_zone","1._dp", security_zone) !in all directions, so it comes to 2.0 Angstroms
  if (.not. param_read_args(cli_params, do_check = .true.)) then
    call system_abort("Usage: decimate [file=(stdin)] [waterfile=(stdin)] xmin xmax ymin ymax zmin zmax [exclusion=(2.4)] [security_zone=(1.0)]")
  endif
  call finalise(cli_params)

  call initialise(xyzfile, filename, INPUT)
  call initialise(waterfile, waterfilename, INPUT)

  !Read in the 2 files
  call read_xyz(at, xyzfile, status=stat)
  if (stat /= 0) call system_abort('Could not read file to solvate.')
  call map_into_cell(at)
  call read_xyz(wat, waterfile, status=stat)
  if (stat /= 0) call system_abort('Could not read water file.')
  call map_into_cell(wat)
  do i=1,3
     do j=1,3
        if (at%lattice(i,j).gt.wat%lattice(i,j)) call system_abort('Too small water file')
     enddo
  enddo

  !Initialise the solvated atoms object
  call initialise(at2,0,reshape((/xmax-xmin,0._dp,0._dp,0._dp,ymax-ymin,0._dp,0._dp,0._dp,zmax-zmin/),(/3,3/)))

  !Add atoms of at to at2
  do i=1,at%N
     call add_atom_single(at2,at%pos(1:3,i),at%Z(i))
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
        call add_atom_single(at2,wat%pos(1:3,i  ),wat%Z(i  ))
        call add_atom_single(at2,wat%pos(1:3,i+1),wat%Z(i+1))
        call add_atom_single(at2,wat%pos(1:3,i+2),wat%Z(i+2))
     endif
  enddo

  call map_into_cell(at2)
  if (stat == 0) call print_xyz(at2, mainlog, all_properties=.true.)

  call verbosity_pop()
  call system_finalise()

end program solvate
