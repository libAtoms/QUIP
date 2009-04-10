module netcdf_module

  !% This module reads and writes Atoms objects to and from disk in NetCDF format, using a 
  !% superset of the AMBER conventions described at {\tt http://amber.scripps.edu/netcdf/nctraj.html}
  !% In addition to the coordinates, lattice and velocites which are included in the
  !% AMBER NetCDF convention, this module also reads and writes all other atomic parameters
  !% (per frame data) and properties (per atom data). The file format is defined by the
  !% an Atoms object used to initialise a NetCDFFile for the first time and can't be changed
  !% thereafter. This means that the number of atoms, set of parameters and set of properties
  !% is fixed for the duration of a simulations. This shouldn't be a problem in practice,
  !% but it's important to remember to add all properties and parameters before the 
  !% call to 'initialise(netcdf,action=OUTPUT)'.

  use System_module
  use Atoms_module
  use netcdf

  implicit none

  private

  public :: NetCDFFile

  type NetCDFFile
     integer :: ncid  !% NetCDF ID of this file

     character(255) :: filename !% Filename with which this NetCDFFile is associated

     integer :: frame_dim_id, spatial_dim_id, atom_dim_id, &
          cell_spatial_dim_id, cell_angular_dim_id, label_dim_id, spatial_var_id, &
          cell_spatial_var_id, cell_angular_var_id, string_dim_id

     type(Dictionary) :: params   !% Mapping from parameter names to types and variable IDs
     type(Dictionary) :: properties !% Mapping from property names to types and variable IDs
     
     integer :: n_frame   !% Number of frames in file - unlimited
     integer :: n_spatial !% Number of spatial dimensions (i.e. 3)
     integer :: n_atom    !% Number of atoms - fixed once file is initialised 
     integer :: n_cell_spatial  !% Number of cell lengths (i.e. 3)
     integer :: n_cell_angular  !% Number of cell angles (i.e. 3)
     integer :: n_label !% Length of string data. Equal to TABLE_STRING_LENGTH

     integer :: frame !% Current frame for reading and writing data

     logical :: initialised = .false.
     
  end type NetCDFFile


  public :: initialise
  interface initialise
     module procedure netcdffile_initialise
  end interface

  public :: finalise
  interface finalise
     module procedure netcdffile_finalise
  end interface

  public :: read_netcdf
  interface read_netcdf
     module procedure  netcdffile_read_atoms
  end interface

  public :: write_netcdf
  interface write_netcdf
     module procedure  netcdffile_write_atoms
  end interface

  public :: write_prmtop

contains

  !% OMIT
  subroutine netcdf_check(status)
    integer, intent(IN) :: status
    
    if (status /= nf90_noerr) &
         call system_abort('netcdf error '//nf90_strerror(status))

  end subroutine netcdf_check

  !% Write an AMBER 'prmtop' file for the Atoms object 'at' to
  !% the file 'filename'. This is a header file containing atom
  !% species and massesm, and enables NetCDF trajectories produced 
  !% by this module to be directly visualised in VMD.
  subroutine write_prmtop(at, filename)
    type(Atoms), intent(in) :: at
    character(*), intent(in) :: filename

    integer :: i
    type(Inoutput) :: file

    call initialise(file, filename, action=OUTPUT)

    write (file%unit, '(A)')          "%VERSION  libAtoms"
    write (file%unit, '(A)')          "%FLAG TITLE"
    write (file%unit, '(A)')          "%FORMAT(20a4)"
    write (file%unit, '(A4,75X,A1)')  "NASN", " "

    write (file%unit, '(A)')          "%FLAG POINTERS"
    write (file%unit, '(A)')          "%FORMAT(10I8)"
    write (file%unit, '(10I8)')       at%N, (0, i = 1, 11)
    write (file%unit, '(10I8)')       (0, i = 1, 12)
    write (file%unit, '(6I8)')        (0, i = 1, 6)

    write (file%unit, '(A)')          "%FLAG ATOM_NAME"
    write (file%unit, '(A)')          "%FORMAT(20a4)"
    if (has_property(at, 'species')) then
       write (file%unit, '(20A4)')       (at%species(i), i=1,at%N)
    else
       write (file%unit, '(20A4)')       (ElementMass(at%Z(i)), i=1,at%N)
    end if

    if (has_property(at, 'mass')) then
       write (file%unit, '(A)')          "%FLAG MASS"
       write (file%unit, '(A)')          "%FORMAT(5E16.5)"
       write (file%unit, '(5E16.5)')     (at%mass(i), i = 1,at%N)
    end if

    call finalise(file)

  end subroutine write_prmtop

  !% Initialise a NetCDFFile using the file 'filename'. 'action'
  !% should be one of 'INPUT', 'OUTPUT' or 'INOUT'. If 'filename'
  !% already exists then the default is 'INOUT', otherwise it is
  !% 'OUTPUT'. In 'INPUT' mode the file is opened for reading only.
  !% 'OUTPUT' mode will overwrite an existing file, while 'INOUT'
  !% will open it in read/write mode if it already exists, or in write
  !% only mode if it does not. The 'append' argument only makes
  !% sense when used with 'INOUT' mode: if it is set to true then
  !% subsequent writes will go to the end of the file; if false
  !% (the default) then writes will overwrite existing data by default
  !% (the 'frame' optional argument to the 'read_atoms' and 'write_atoms'
  !% methods can be used to override this behaviour, or 'this%frame'
  !% can be set explicity).
  subroutine netcdffile_initialise(this, filename, at, action, append)
    type(NetCDFFile), intent(inout) :: this
    character(*), intent(in) :: filename
    type(Atoms), intent(inout), optional :: at
    integer, intent(in), optional :: action
    logical, intent(in), optional :: append
    
    logical :: exist

    this%filename = filename

    if (present(action)) then
       if (action == INPUT) then
          call netcdffile_open_read(this, filename, read_only=.true.,append=append)
       else if (action == OUTPUT) then
          if (.not. present(at)) &
               call system_abort('netcdffile_initialise: at must be present when writing new file')
          call netcdffile_open_write(this, filename, at)
       else if (action == INOUT) then

          inquire (file=filename, exist=exist)

          if (exist) then
             call netcdffile_open_read(this, filename, read_only=.false.,append=append)
          else
             if (.not. present(at)) &
                  call system_abort('netcdffile_initialise: at must be present when writing new file')
             call netcdffile_open_write(this, filename, at)
          end if
          
       end if
    else
       ! Guess whether to open the file in INOUT or OUTPUT mode
       inquire (file=filename, exist=exist)

       if (exist) then
          call netcdffile_open_read(this, filename, read_only=.false.,append=append)
       else
          if (.not. present(at)) &
               call system_abort('netcdffile_initialise: at must be present when writing new file')
          call netcdffile_open_write(this, filename, at)
       end if
    end if

  end subroutine netcdffile_initialise


  !% Close the file associated with this NetCDF object and finalise
  !% internal data. 
  subroutine netcdffile_finalise(this)
    type(NetCDFFile), intent(inout) :: this

    if (this%initialised) then
       call netcdf_check(nf90_close(this%ncid))
       this%initialised = .false.
       this%n_frame = 0
       this%n_atom = 0
       this%n_spatial = 0
       this%n_cell_spatial = 0
       this%n_cell_angular = 0
       this%n_label = 0
       
       call finalise(this%params)
       call finalise(this%properties)

    end if

  end subroutine netcdffile_finalise
    
  !% OMIT
  subroutine netcdffile_open_write(this, filename, at)
    type(NetCDFFile), intent(inout) :: this
    character(*), intent(in) :: filename
    type(Atoms), intent(inout) :: at

    integer :: var_id, ndims, lookup(3), i
    character(key_len) :: name
    real(dp) :: lattice_lengths(3), lattice_angles(3)

    if (this%initialised) &
         call system_abort('netcdffile_open_write - this NetCDFFile is already open')

    this%n_frame = 0
    this%n_atom = at%N
    this%n_spatial = 3
    this%n_cell_spatial = 3
    this%n_cell_angular = 3
    this%n_label = TABLE_STRING_LENGTH
   
    call netcdf_check(nf90_create(filename, ior(NF90_64BIT_OFFSET,NF90_CLOBBER), this%ncid))
    
    ! Global attributes
    call netcdf_check(nf90_put_att(this%ncid, NF90_GLOBAL, 'Conventions', 'AMBER'))
    call netcdf_check(nf90_put_att(this%ncid, NF90_GLOBAL, 'ConventionVersion', '1.0'))
    call netcdf_check(nf90_put_att(this%ncid, NF90_GLOBAL, 'application', 'libAtoms'))
    call netcdf_check(nf90_put_att(this%ncid, NF90_GLOBAL, 'program', 'libAtoms'))
    call netcdf_check(nf90_put_att(this%ncid, NF90_GLOBAL, 'programVersion', 'svn: '//SVN_VERSION)) 
    call netcdf_check(nf90_put_att(this%ncid, NF90_GLOBAL, 'title', 'Atoms Object'))

    ! Dimensions
    call netcdf_check(nf90_def_dim(this%ncid, 'frame', NF90_UNLIMITED, this%frame_dim_id))
    call netcdf_check(nf90_def_dim(this%ncid, 'spatial', this%n_spatial, this%spatial_dim_id))
    call netcdf_check(nf90_def_dim(this%ncid, 'atom', at%N, this%atom_dim_id))
    call netcdf_check(nf90_def_dim(this%ncid, 'cell_spatial', this%n_cell_spatial, this%cell_spatial_dim_id))
    call netcdf_check(nf90_def_dim(this%ncid, 'cell_angular', this%n_cell_angular, this%cell_angular_dim_id))
    call netcdf_check(nf90_def_dim(this%ncid, 'label', TABLE_STRING_LENGTH, this%label_dim_id))
    call netcdf_check(nf90_def_dim(this%ncid, 'string', VALUE_LEN, this%string_dim_id))


    ! Label variables
    call netcdf_check(nf90_def_var(this%ncid, 'spatial', NF90_CHAR, (/ this%spatial_dim_id /), &
         this%spatial_var_id))
    call netcdf_check(nf90_def_var(this%ncid, 'cell_spatial', NF90_CHAR, (/ this%cell_spatial_dim_id /), &
         this%cell_spatial_var_id))
    call netcdf_check(nf90_def_var(this%ncid, 'cell_angular', NF90_CHAR, &
         (/ this%label_dim_id, this%cell_angular_dim_id /), this%cell_angular_var_id))

    ! add entries for lattice to at%params

    call get_lattice_params(at%lattice, lattice_lengths(1), &
         lattice_lengths(2), lattice_lengths(3), lattice_angles(1), &
         lattice_angles(2), lattice_angles(3))

    lattice_angles = lattice_angles * DEGREES_PER_RADIAN
    
    call set_value(at%params, 'cell_lengths', lattice_lengths)
    call set_value(at%params, 'cell_angles',  lattice_angles)

    call initialise(this%params)

    do i=1,at%params%N

       if (trim(at%params%keys(i)) == 'Lattice' .or. &
           trim(at%params%keys(i)) == 'Properties') cycle
       
       select case(at%params%entries(i)%type)

          case(T_INTEGER)
             call netcdf_check(nf90_def_var(this%ncid, at%params%keys(i), NF90_INT, &
                  (/ this%frame_dim_id /), var_id))

             call set_value(this%params, at%params%keys(i), (/ 1, T_INTEGER, var_id /))

          case(T_REAL)
             call netcdf_check(nf90_def_var(this%ncid, at%params%keys(i), NF90_DOUBLE, &
                  (/ this%frame_dim_id /), var_id))

             call set_value(this%params, at%params%keys(i), (/ 1, T_REAL, var_id /))

          case (T_CHAR)
             call netcdf_check(nf90_def_var(this%ncid, at%params%keys(i), NF90_CHAR, &
                  (/ this%string_dim_id, this%frame_dim_id /), var_id))

             call set_value(this%params, at%params%keys(i), (/ 1, T_CHAR, var_id /))
          
          case(T_INTEGER_A)
             if (size(at%params%entries(i)%i_a) /= 3) cycle

             call netcdf_check(nf90_def_var(this%ncid, at%params%keys(i), NF90_INT, &
                  (/ this%spatial_dim_id, this%frame_dim_id /), var_id))

             call set_value(this%params, at%params%keys(i), &
                  (/ 3, T_INTEGER_A, var_id /))

          case(T_REAL_A)
             if (size(at%params%entries(i)%r_a) /= 3) cycle

             call netcdf_check(nf90_def_var(this%ncid, at%params%keys(i), NF90_DOUBLE, &
                  (/ this%spatial_dim_id, this%frame_dim_id /), var_id))

             call set_value(this%params, at%params%keys(i), &
                  (/ 3, T_REAL_A, var_id /))

          case (T_CHAR_A)

             if (size(at%params%entries(i)%s_a) /= 3) cycle

             call netcdf_check(nf90_def_var(this%ncid, at%params%keys(i), NF90_CHAR, &
                  (/ this%string_dim_id, this%spatial_dim_id, this%frame_dim_id /), var_id))

             call set_value(this%params, at%params%keys(i), &
                  (/ 3, T_CHAR_A, var_id /))

          case default

             cycle
             
       end select
       call netcdf_check(nf90_put_att(this%ncid, var_id, 'type', at%params%entries(i)%type))

    end do

    ! Copy per-atom properties
    call initialise(this%properties)
    do i=1,at%properties%N
       if (.not. get_value(at%properties, at%properties%keys(i), lookup)) &
            call system_abort('netcdf_open_write: missing property '//at%properties%keys(i))

       ndims = (lookup(3)-lookup(2)+1)
       if (ndims /= 1 .and. ndims /= 3) cycle
       
       select case(lookup(1))
          
          case(PROPERTY_INT)

             if (ndims == 1) then
                call netcdf_check(nf90_def_var(this%ncid, at%properties%keys(i), NF90_INT, &
                     (/ this%atom_dim_id, this%frame_dim_id /), var_id))

                call set_value(this%properties, at%properties%keys(i), &
                     (/ 1, PROPERTY_INT, var_id /))
             else
                call netcdf_check(nf90_def_var(this%ncid, at%properties%keys(i), NF90_INT, &
                     (/ this%spatial_dim_id, this%atom_dim_id, this%frame_dim_id /), var_id))

                call set_value(this%properties, at%properties%keys(i), &
                     (/ 3, PROPERTY_INT, var_id /))
             end if
             
          
          case(PROPERTY_REAL)
             ! Do some name mangling to maintain AMBER and VMD compatibility
             name = at%properties%keys(i)
             if (trim(name) == 'pos')  name = 'coordinates'
             if (trim(name) == 'velo') name = 'velocities'

             if (ndims == 1) then
                call netcdf_check(nf90_def_var(this%ncid, name, NF90_DOUBLE, &
                     (/ this%atom_dim_id, this%frame_dim_id /), var_id))

                call set_value(this%properties, at%properties%keys(i), &
                     (/ 1, PROPERTY_REAL, var_id /))
             else
                call netcdf_check(nf90_def_var(this%ncid, name, NF90_DOUBLE, &
                     (/ this%spatial_dim_id, this%atom_dim_id, this%frame_dim_id /), var_id))

                call set_value(this%properties, at%properties%keys(i), &
                     (/ 3, PROPERTY_REAL, var_id /))
             end if
                       
          case(PROPERTY_STR)

             if (ndims == 1) then
                call netcdf_check(nf90_def_var(this%ncid, at%properties%keys(i), NF90_CHAR, &
                     (/ this%label_dim_id, this%atom_dim_id, this%frame_dim_id /), var_id))

                call set_value(this%properties, at%properties%keys(i), &
                     (/ 1, PROPERTY_STR, var_id /))
             else
                call netcdf_check(nf90_def_var(this%ncid, at%properties%keys(i), NF90_CHAR, &
                     (/ this%label_dim_id, this%spatial_dim_id, this%atom_dim_id, this%frame_dim_id /), var_id))

                call set_value(this%properties, at%properties%keys(i), &
                     (/ 3, PROPERTY_STR, var_id /))
             end if

          case(PROPERTY_LOGICAL)

             if (ndims == 1) then
                call netcdf_check(nf90_def_var(this%ncid, at%properties%keys(i), NF90_INT, &
                     (/ this%atom_dim_id, this%frame_dim_id /), var_id))

                call set_value(this%properties, at%properties%keys(i), &
                     (/ 1, PROPERTY_LOGICAL, var_id /))

             else
                call netcdf_check(nf90_def_var(this%ncid, at%properties%keys(i), NF90_INT, &
                     (/ this%spatial_dim_id, this%atom_dim_id, this%frame_dim_id /), var_id))

                call set_value(this%properties, at%properties%keys(i), &
                     (/ 3, PROPERTY_LOGICAL, var_id /))
             end if

       end select
       call netcdf_check(nf90_put_att(this%ncid, var_id, 'type', lookup(1)))

    end do

    call netcdf_check(nf90_enddef(this%ncid))

    ! Put label variables
    call netcdf_check(nf90_put_var(this%ncid, this%spatial_var_id, 'xyz'))
    call netcdf_check(nf90_put_var(this%ncid, this%cell_spatial_var_id, 'abc'))
    call netcdf_check(nf90_put_var(this%ncid, this%cell_angular_var_id, (/'alpha', 'beta ', 'gamma'/)))

    !**TODO** set units correctly

    this%n_frame = 0
    this%frame = 1
    this%initialised = .true.

    call netcdf_check(nf90_sync(this%ncid))

  end subroutine netcdffile_open_write

  !% OMIT
  subroutine netcdffile_open_read(this, filename, read_only, append)
    type(NetCDFFile), intent(inout) :: this
    character(*), intent(in) :: filename
    logical, optional :: read_only, append

    logical :: do_read_only, do_append
    character(TABLE_STRING_LENGTH) :: spatial_var, cell_spatial_var, cell_angular_var(3)
    character(key_len) :: name
    integer :: ndims, nvars, type, mytype, size, varid
    integer, allocatable, dimension(:) :: dimids

    if (this%initialised) call system_abort('netcdffile_open_read: this NetCDFFile is already open')

    do_read_only = optional_default(.true., read_only)
    do_append = optional_default(.false., append)

    if (do_read_only) then
       call netcdf_check(nf90_open(filename, ior(NF90_64BIT_OFFSET,NF90_NOWRITE), this%ncid))
    else
       call netcdf_check(nf90_open(filename, ior(NF90_64BIT_OFFSET,NF90_WRITE), this%ncid))
    end if

    ! Dimensions
    call netcdf_check(nf90_inq_dimid(this%ncid, 'frame', this%frame_dim_id))
    call netcdf_check(nf90_inq_dimid(this%ncid, 'spatial', this%spatial_dim_id))
    call netcdf_check(nf90_inq_dimid(this%ncid, 'atom', this%atom_dim_id))
    call netcdf_check(nf90_inq_dimid(this%ncid, 'cell_spatial', this%cell_spatial_dim_id))
    call netcdf_check(nf90_inq_dimid(this%ncid, 'cell_angular', this%cell_angular_dim_id))
    call netcdf_check(nf90_inq_dimid(this%ncid, 'label', this%label_dim_id))
    call netcdf_check(nf90_inq_dimid(this%ncid, 'string', this%string_dim_id))

    call netcdf_check(nf90_inquire_dimension(this%ncid, this%frame_dim_id, len=this%n_frame))

    call netcdf_check(nf90_inquire_dimension(this%ncid, this%atom_dim_id, len=this%n_atom))

    call netcdf_check(nf90_inquire_dimension(this%ncid, this%spatial_dim_id, len=this%n_spatial))
    if (this%n_spatial /= 3) call system_abort('netcdffile_open_read: number of spatial dimensions /= 3')

    call netcdf_check(nf90_inquire_dimension(this%ncid, this%cell_spatial_dim_id, len=this%n_cell_spatial))
    if (this%n_cell_spatial /= 3) call system_abort('netcdffile_open_read: number of cell spatial dimensions /= 3')

    call netcdf_check(nf90_inquire_dimension(this%ncid, this%cell_angular_dim_id, len=this%n_cell_angular))
    if (this%n_cell_angular /= 3) call system_abort('netcdffile_open_read: number of cell angular dimensions /= 3')

    call netcdf_check(nf90_inquire_dimension(this%ncid, this%label_dim_id, len=this%n_label))
    if (this%n_label /= TABLE_STRING_LENGTH) &
         call system_abort('netcdffile_open_read: string length /= '//TABLE_STRING_LENGTH)

    ! Label variables
    call netcdf_check(nf90_inq_varid(this%ncid, 'spatial', this%spatial_var_id))
    call netcdf_check(nf90_inq_varid(this%ncid, 'cell_spatial', this%cell_spatial_var_id))
    call netcdf_check(nf90_inq_varid(this%ncid, 'cell_angular', this%cell_angular_var_id))

    ! Check the label variables
    call netcdf_check(nf90_get_var(this%ncid, this%spatial_var_id, spatial_var, &
         start=(/1/), count=(/this%n_spatial/)))

    if (spatial_var(1:3) /= 'xyz') &
         call system_abort('netcdffile_open_read: spatial('//trim(spatial_var)//') /= xyz')

    call netcdf_check(nf90_get_var(this%ncid, this%cell_spatial_var_id, cell_spatial_var, &
         start=(/1/), count=(/this%n_cell_spatial/)))
    if (cell_spatial_var(1:3) /= 'abc') &
         call system_abort('netcdffile_open_read: cell_spatial('//trim(cell_spatial_var)//') /= abc')

    call netcdf_check(nf90_get_var(this%ncid, this%cell_angular_var_id, cell_angular_var, &
         start=(/1,1/), count=(/this%n_label,this%n_cell_angular/)))

    if (cell_angular_var(1)(1:5) /= 'alpha' .or. &
        cell_angular_var(2)(1:4) /= 'beta' .or. &
        cell_angular_var(3)(1:5) /= 'gamma') &
          call system_abort('netcdffile_open_read: cell_angular /= (/alpha,beta,gamma/)')


    call netcdf_check(nf90_inquire(this%ncid, nvariables=nvars))

    call initialise(this%params)
    call initialise(this%properties)

    do varid=1,nvars
       call netcdf_check(nf90_inquire_variable(this%ncid, varid, ndims=ndims))
       allocate(dimids(ndims))

       call netcdf_check(nf90_inquire_variable(this%ncid, varid, dimids=dimids, &
            name=name, xtype=type))

       if (dimids(ndims) == this%frame_dim_id) then
          if (ndims == 0) then
             ! one-off label variable, do nothing
             
          else if (ndims == 1) then
             if (dimids(1) == this%frame_dim_id) then
                ! it's a scalar per-frame property, i.e. a param
                call netcdf_check(nf90_inquire_dimension(this%ncid, dimids(1), len=size))
                call netcdf_check(nf90_get_att(this%ncid, varid, 'type', mytype))
                call set_value(this%params, name, (/size, mytype, varid/))
             else
                call system_abort('netcdffile_open_read: Unknown one dimensional variable '//trim(name))
             end if
          else if (ndims == 2) then
             if (dimids(ndims) == this%frame_dim_id .and. dimids(ndims-1) == this%atom_dim_id) then
                ! it's a scalar per-atom property
                call netcdf_check(nf90_get_att(this%ncid, varid, 'type', mytype))
                call set_value(this%properties, name, (/1, mytype, varid/))
               
             else if (dimids(ndims) == this%frame_dim_id) then
                ! it's an array per-frame parameter
                call netcdf_check(nf90_inquire_dimension(this%ncid, dimids(1), len=size))
                call netcdf_check(nf90_get_att(this%ncid, varid, 'type', mytype))
                call set_value(this%params, name, (/size, mytype, varid/))
             else
                call system_abort('netcdffile_open_read: Unknown two dimensional variable '//trim(name))
             end if
          else if (ndims == 3) then
             if (dimids(ndims) == this%frame_dim_id .and. dimids(ndims-1) == this%atom_dim_id) then
                if (dimids(1) == this%label_dim_id) then
                   ! it's a per-atom scalar string
                   call netcdf_check(nf90_get_att(this%ncid, varid, 'type', mytype))
                   call set_value(this%properties, name, (/1, mytype, varid/))
                else
                   ! it's an array per-atom property
                   call netcdf_check(nf90_inquire_dimension(this%ncid, dimids(1), len=size))
                   call netcdf_check(nf90_get_att(this%ncid, varid, 'type', mytype))
                
                   ! Reverse name mangling required for AMBER/VMD compatibility
                   if (trim(name) == 'coordinates') name = 'pos'
                   if (trim(name) == 'velocities') name = 'velo'
                   call set_value(this%properties, name, (/size, mytype, varid/))
                end if
             else
                call system_abort('netcdffile_open_read: Unknown three dimensional variable '//trim(name))
             end if
          else if (ndims == 4) then
             if (dimids(ndims) == this%frame_dim_id .and. dimids(ndims-1) == this%atom_dim_id &
                  .and. dimids(1) == this%label_dim_id) then
                ! it's a per-atom string array property
                call netcdf_check(nf90_inquire_dimension(this%ncid, dimids(2), len=size))
                call netcdf_check(nf90_get_att(this%ncid, varid, 'type', mytype))
                call set_value(this%properties, name, (/size, mytype, varid/))
             else
                call system_abort('netcdffile_open_read: Unknown three dimensional variable '//trim(name))
             end if
          else
             call system_abort('netcdffile_open_read: Unknown '//ndims//'-dimensional variable '//trim(name))
             
          end if
       end if
       deallocate(dimids)
    end do

    call print('Properties read from file:', VERBOSE)
    call print(this%properties, VERBOSE)
    call print('', VERBOSE)
    call print('Parameters read from file:', VERBOSE)
    call print(this%params, VERBOSE)
    call print('', VERBOSE)

    this%initialised = .true.
    if (do_append) then
       this%frame = this%n_frame + 1
    else
       this%frame = 1
    end if

  end subroutine netcdffile_open_read


  !% Write an Atoms object 'at' to an open NetCDFFile 'this'.
  !% The Atoms object must have the same number of atoms, parameters and properties
  !% as the one used to initialise the NetCDFFile. The optional argument 'frame'
  !% can be used to specify where in the file this Atoms object should be written;
  !% if it is not present then we use 'this%frame' and then advance it by one.
  !% If 'prmtop' is true then we also write an AMBER style prmtop file.
  subroutine netcdffile_write_atoms(this, at, frame, prmtop)
    type(NetCDFFile), intent(inout) :: this
    type(Atoms), intent(inout) :: at
    integer, optional :: frame
    logical, optional :: prmtop

    type(Dictionary) :: tmp_params
    real(dp) :: lattice_lengths(3), lattice_angles(3)
    integer :: use_frame, i, j, lookup(3), alookup(3)
    integer, allocatable, dimension(:,:) :: tmp_integer
    integer :: var_id, ndims
    character(key_len) :: name
    logical :: mismatch

    if (present(frame)) then
       use_frame = frame
    else
       use_frame = this%frame
       this%frame = this%frame + 1
       this%n_frame = this%n_frame + 1
    end if

    if (.not. this%initialised) &
         call system_abort('netcdffile_write_atoms: NetCDFFile not initialised')

    if (at%N /= this%n_atom) call system_abort('netcdffile_write_atoms: at%N('//(at%N)//') /= ('&
         //this%n_atom//')')

    call get_lattice_params(at%lattice, lattice_lengths(1), &
         lattice_lengths(2), lattice_lengths(3), lattice_angles(1), &
         lattice_angles(2), lattice_angles(3))
    lattice_angles = lattice_angles * DEGREES_PER_RADIAN

    tmp_params = at%params
    if (has_key(tmp_params, 'Lattice')) call remove_value(tmp_params, 'Lattice')
    if (has_key(tmp_params, 'Properties')) call remove_value(tmp_params, 'Properties')

    call set_value(tmp_params, 'cell_lengths', lattice_lengths)
    call set_value(tmp_params, 'cell_angles',  lattice_angles)

    mismatch = this%params%N /= tmp_params%N
    if (.not. mismatch) mismatch = any(this%params%keys(1:this%params%N) /= tmp_params%keys(1:this%params%N))
    if (mismatch) then
       call Print('NetCDF: Entering REDEF mode')
       call print('this%params')
       call print(this%params)
       call print('')
       call print('tmp_params')
       call print(tmp_params)

       call netcdf_check(nf90_redef(this%ncid))

       do i=1,tmp_params%N
          
          if (has_key(this%params, trim(tmp_params%keys(i)))) cycle

          if (trim(tmp_params%keys(i)) == 'Lattice' .or. &
               trim(tmp_params%keys(i)) == 'Properties') cycle
       
          select case(tmp_params%entries(i)%type)
             
          case(T_INTEGER)
             call netcdf_check(nf90_def_var(this%ncid, tmp_params%keys(i), NF90_INT, &
                  (/ this%frame_dim_id /), var_id))
             
             call set_value(this%params, tmp_params%keys(i), (/ 1, T_INTEGER, var_id /))

          case(T_REAL)
             call netcdf_check(nf90_def_var(this%ncid, tmp_params%keys(i), NF90_DOUBLE, &
                  (/ this%frame_dim_id /), var_id))

             call set_value(this%params, tmp_params%keys(i), (/ 1, T_REAL, var_id /))

          case (T_CHAR)
             call netcdf_check(nf90_def_var(this%ncid, tmp_params%keys(i), NF90_CHAR, &
                  (/ this%string_dim_id, this%frame_dim_id /), var_id))

             call set_value(this%params, tmp_params%keys(i), (/ 1, T_CHAR, var_id /))
          
          case(T_INTEGER_A)
             if (size(tmp_params%entries(i)%i_a) /= 3) cycle

             call netcdf_check(nf90_def_var(this%ncid, tmp_params%keys(i), NF90_INT, &
                  (/ this%spatial_dim_id, this%frame_dim_id /), var_id))

             call set_value(this%params, tmp_params%keys(i), &
                  (/ 3, T_INTEGER_A, var_id /))

          case(T_REAL_A)
             if (size(tmp_params%entries(i)%r_a) /= 3) cycle

             call netcdf_check(nf90_def_var(this%ncid, tmp_params%keys(i), NF90_DOUBLE, &
                  (/ this%spatial_dim_id, this%frame_dim_id /), var_id))

             call set_value(this%params, tmp_params%keys(i), &
                  (/ 3, T_REAL_A, var_id /))

          case (T_CHAR_A)

             if (size(tmp_params%entries(i)%s_a) /= 3) cycle

             call netcdf_check(nf90_def_var(this%ncid, tmp_params%keys(i), NF90_CHAR, &
                  (/ this%string_dim_id, this%spatial_dim_id, this%frame_dim_id /), var_id))

             call set_value(this%params, tmp_params%keys(i), &
                  (/ 3, T_CHAR_A, var_id /))

          case default

             cycle
             
          end select
          call netcdf_check(nf90_put_att(this%ncid, var_id, 'type', tmp_params%entries(i)%type))
       end do

       call netcdf_check(nf90_enddef(this%ncid))

    end if

    do i=1,tmp_params%N

       if (trim(tmp_params%keys(i)) == 'Lattice' .or. &
           trim(tmp_params%keys(i)) == 'Properties') cycle

       if (.not. get_value(this%params, tmp_params%keys(i),lookup)) &
            call system_abort('netcdffile_write_atoms: unknown param '//trim(tmp_params%keys(i))//' in Atoms object')
       
       select case(lookup(2))

          case(T_INTEGER)
             call netcdf_check(nf90_put_var(this%ncid, lookup(3), tmp_params%entries(i)%i, &
                  start=(/use_frame/)))

          case(T_REAL)
             call netcdf_check(nf90_put_var(this%ncid, lookup(3), tmp_params%entries(i)%r, &
                  start=(/use_frame/)))

          case(T_CHAR)
             call netcdf_check(nf90_put_var(this%ncid, lookup(3), tmp_params%entries(i)%s, &
                  start=(/1,use_frame/), count=(/len(trim(tmp_params%entries(i)%s)),1/)))

          case(T_INTEGER_A)
             call netcdf_check(nf90_put_var(this%ncid, lookup(3), tmp_params%entries(i)%i_a, &
                  start=(/1,use_frame/), count=(/lookup(1),1/)))

          case(T_REAL_A)
             call netcdf_check(nf90_put_var(this%ncid, lookup(3), tmp_params%entries(i)%r_a, &
                  start=(/1,use_frame/), count=(/lookup(1),1/)))

          case(T_CHAR_A)
             do j=1,lookup(1) 
                call netcdf_check(nf90_put_var(this%ncid, lookup(3), tmp_params%entries(i)%s_a(j), &
                     start=(/1,j,use_frame/), count=(/len(trim(tmp_params%entries(i)%s_a(j))),1,1/)))
             end do

       end select
    end do

    mismatch = this%properties%N /= at%properties%N
    if (.not. mismatch) mismatch = any(this%properties%keys(1:this%properties%N) /= at%properties%keys(1:this%properties%N))
    if (mismatch) then
       call Print('NetCDF: entering REDEF mode')
       call print('this%properties')
       call print(this%properties)
       call print('')
       call print('at%properties')
       call print(at%properties)

       call netcdf_check(nf90_redef(this%ncid))

       do i=1,at%properties%N

          if (has_key(this%properties, at%properties%keys(i))) cycle

          if (.not. get_value(at%properties, at%properties%keys(i), lookup)) &
               call system_abort('netcdffile_open_write: missing property '//at%properties%keys(i))

          ndims = (lookup(3)-lookup(2)+1)
          if (ndims /= 1 .and. ndims /= 3) cycle
       
          select case(lookup(1))
          
          case(PROPERTY_INT)

             if (ndims == 1) then
                call netcdf_check(nf90_def_var(this%ncid, at%properties%keys(i), NF90_INT, &
                     (/ this%atom_dim_id, this%frame_dim_id /), var_id))

                call set_value(this%properties, at%properties%keys(i), &
                     (/ 1, PROPERTY_INT, var_id /))
             else
                call netcdf_check(nf90_def_var(this%ncid, at%properties%keys(i), NF90_INT, &
                     (/ this%spatial_dim_id, this%atom_dim_id, this%frame_dim_id /), var_id))

                call set_value(this%properties, at%properties%keys(i), &
                     (/ 3, PROPERTY_INT, var_id /))
             end if
             
          
          case(PROPERTY_REAL)
             ! Do some name mangling to maintain AMBER and VMD compatibility
             name = at%properties%keys(i)
             if (trim(name) == 'pos')  name = 'coordinates'
             if (trim(name) == 'velo') name = 'velocities'

             if (ndims == 1) then
                call netcdf_check(nf90_def_var(this%ncid, name, NF90_DOUBLE, &
                     (/ this%atom_dim_id, this%frame_dim_id /), var_id))

                call set_value(this%properties, at%properties%keys(i), &
                     (/ 1, PROPERTY_REAL, var_id /))
             else
                call netcdf_check(nf90_def_var(this%ncid, name, NF90_DOUBLE, &
                     (/ this%spatial_dim_id, this%atom_dim_id, this%frame_dim_id /), var_id))

                call set_value(this%properties, at%properties%keys(i), &
                     (/ 3, PROPERTY_REAL, var_id /))
             end if
                       
          case(PROPERTY_STR)

             if (ndims == 1) then
                call netcdf_check(nf90_def_var(this%ncid, at%properties%keys(i), NF90_CHAR, &
                     (/ this%label_dim_id, this%atom_dim_id, this%frame_dim_id /), var_id))

                call set_value(this%properties, at%properties%keys(i), &
                     (/ 1, PROPERTY_STR, var_id /))
             else
                call netcdf_check(nf90_def_var(this%ncid, at%properties%keys(i), NF90_CHAR, &
                     (/ this%label_dim_id, this%spatial_dim_id, this%atom_dim_id, this%frame_dim_id /), var_id))

                call set_value(this%properties, at%properties%keys(i), &
                     (/ 3, PROPERTY_STR, var_id /))
             end if

          case(PROPERTY_LOGICAL)

             if (ndims == 1) then
                call netcdf_check(nf90_def_var(this%ncid, at%properties%keys(i), NF90_INT, &
                     (/ this%atom_dim_id, this%frame_dim_id /), var_id))

                call set_value(this%properties, at%properties%keys(i), &
                     (/ 1, PROPERTY_LOGICAL, var_id /))

             else
                call netcdf_check(nf90_def_var(this%ncid, at%properties%keys(i), NF90_INT, &
                     (/ this%spatial_dim_id, this%atom_dim_id, this%frame_dim_id /), var_id))

                call set_value(this%properties, at%properties%keys(i), &
                     (/ 3, PROPERTY_LOGICAL, var_id /))
             end if

          end select
          call netcdf_check(nf90_put_att(this%ncid, var_id, 'type', lookup(1)))
       end do
       
       call netcdf_check(nf90_enddef(this%ncid))

    end if

    do i=1,at%properties%N

       if (.not. get_value(this%properties, at%properties%keys(i), lookup)) &
            call system_abort('netcdffile_write_atoms: property '//&
            trim(at%properties%keys(i))//' not found in this%properties')

       if (.not. get_value(at%properties, at%properties%keys(i), alookup)) &
            call system_abort('netcdffile_write_atoms: unknown property '//&
            trim(at%properties%keys(i))//' in Atoms object')

       select case(lookup(2))
          
          case(PROPERTY_INT)
             if (lookup(1) == 1) then
                call netcdf_check(nf90_put_var(this%ncid, lookup(3), &
                     at%data%int(alookup(2):alookup(3),1:at%N), &
                     start=(/1,use_frame/), count=(/at%N,1/)))
             else
                call netcdf_check(nf90_put_var(this%ncid, lookup(3), &
                     at%data%int(alookup(2):alookup(3),1:at%N), &
                     start=(/1,1,use_frame/), count=(/lookup(1),at%N,1/)))
             end if

          case(PROPERTY_REAL)
             if (lookup(1) == 1) then
                call netcdf_check(nf90_put_var(this%ncid, lookup(3), &
                     at%data%real(alookup(2):alookup(3),1:at%N), &
                     start=(/1,use_frame/), count=(/at%N,1/)))
             else
                call netcdf_check(nf90_put_var(this%ncid, lookup(3), &
                     at%data%real(alookup(2):alookup(3),1:at%N), &
                     start=(/1,1,use_frame/), count=(/lookup(1),at%N,1/)))
             end if

          case(PROPERTY_STR)
             if (lookup(1) == 1) then
                call netcdf_check(nf90_put_var(this%ncid, lookup(3), &
                     at%data%str(alookup(2):alookup(3),1:at%N), &
                     start=(/1,1,use_frame/), count=(/TABLE_STRING_LENGTH,at%N,1/)))
             else
                call netcdf_check(nf90_put_var(this%ncid, lookup(3), &
                     at%data%str(alookup(2):alookup(3),1:at%N), &
                     start=(/1,1,1,use_frame/), count=(/TABLE_STRING_LENGTH,lookup(1),at%N,1/)))
             end if
                
          case(PROPERTY_LOGICAL)
             if (lookup(1) == 1) then
                allocate(tmp_integer(alookup(2):alookup(3),1:at%N))
                tmp_integer = at%data%logical(alookup(2):alookup(3),1:at%N)
                call netcdf_check(nf90_put_var(this%ncid, lookup(3), &
                     tmp_integer, start=(/1,use_frame/), count=(/at%N,1/)))
                deallocate(tmp_integer)
             else
                allocate(tmp_integer(alookup(2):alookup(3),1:at%N))
                tmp_integer = at%data%logical(alookup(2):alookup(3),1:at%N)
                call netcdf_check(nf90_put_var(this%ncid, lookup(3), &
                     tmp_integer, start=(/1,1,use_frame/), count=(/lookup(1),at%N,1/)))
                deallocate(tmp_integer)
             end if

       end select

    end do
    
    call netcdf_check(nf90_sync(this%ncid))

    if (present(prmtop)) then
       if (prmtop) then
          ! Change file extension to .prmtop
          call write_prmtop(at, this%filename(1:index(this%filename,'.',.true.))//'prmtop')
       end if
    end if

  end subroutine netcdffile_write_atoms


  !% Read an Atoms object 'at' from an open NetCDFFile 'this'. The optional
  !% argument 'frame' can be used to seek to an arbitary position in the file.
  !% If it is not present we default to sequential access, with the frame
  !% pointer 'this%frame' advancing by one on each call to 'read_atoms'.
  subroutine netcdffile_read_atoms(this, at, frame)
    type(NetCDFFile), intent(inout) :: this
    type(Atoms), intent(out) :: at
    integer, optional :: frame

    integer :: i, use_frame, lookup(3), alookup(3), tmp_i
    real(dp) :: lattice(3,3), lattice_lengths(3), lattice_angles(3), tmp_r
    character(VALUE_LEN) :: tmp_s
    integer, allocatable, dimension(:) :: tmp_i_a
    real(dp), allocatable, dimension(:) :: tmp_r_a
    character(VALUE_LEN), allocatable, dimension(:) :: tmp_s_a
    integer, allocatable, dimension(:,:) :: tmp_integer
    type(Dictionary) :: tmp_params

    if (present(frame)) then
       use_frame = frame
    else
       use_frame = this%frame
       this%frame = this%frame + 1
    end if

    call initialise(tmp_params)
    do i=1,this%params%N
       if (.not. get_value(this%params, this%params%keys(i),lookup)) &
            call system_abort('netcdffile_read_atoms: unknown param '//this%params%keys(i))
       
       select case(lookup(2))

          case(T_INTEGER)
             call netcdf_check(nf90_get_var(this%ncid, lookup(3), tmp_i, &
                  start=(/use_frame/)))
             call set_value(tmp_params, this%params%keys(i), tmp_i)

          case(T_REAL)
             call netcdf_check(nf90_get_var(this%ncid, lookup(3), tmp_r, &
                  start=(/use_frame/)))
             call set_value(tmp_params, this%params%keys(i), tmp_r)

          case(T_CHAR)
             call netcdf_check(nf90_get_var(this%ncid, lookup(3), tmp_s, &
                  start=(/1,use_frame/), count=(/VALUE_LEN,1/)))
             call set_value(tmp_params, this%params%keys(i), tmp_s)

          case(T_INTEGER_A)
             allocate(tmp_i_a(lookup(1)))
             call netcdf_check(nf90_get_var(this%ncid, lookup(3), tmp_i_a, &
                  start=(/1,use_frame/), count=(/lookup(1),1/)))
             call set_value(tmp_params, this%params%keys(i), tmp_i_a)
             deallocate(tmp_i_a)
             
          case(T_REAL_A)
             allocate(tmp_r_a(lookup(1)))
             call netcdf_check(nf90_get_var(this%ncid, lookup(3), tmp_r_a, &
                  start=(/1,use_frame/), count=(/lookup(1),1/)))
             call set_value(tmp_params, this%params%keys(i), tmp_r_a)
             deallocate(tmp_r_a)

          case(T_CHAR_A)
             allocate(tmp_s_a(lookup(1)))
             call netcdf_check(nf90_get_var(this%ncid, lookup(3), tmp_s_a, &
                  start=(/1,1,use_frame/), count=(/VALUE_LEN,lookup(1),1/)))
             call set_value(tmp_params, this%params%keys(i), tmp_s_a)
             deallocate(tmp_s_a)

       end select
    end do

    if (.not. get_value(tmp_params, 'cell_lengths', lattice_lengths)) &
         call system_abort('netcdffile_read_atoms: missing cell_lengths variable')

    if (.not. get_value(tmp_params, 'cell_angles', lattice_angles)) &
         call system_abort('netcdffile_read_atoms: missing cell_angles variable')

    lattice_angles = lattice_angles * RADIANS_PER_DEGREE

    lattice = make_lattice(lattice_lengths(1), lattice_lengths(2), lattice_lengths(3), &
         lattice_angles(1), lattice_angles(2), lattice_angles(3))

    call initialise(at, this%n_atom, lattice)

    at%params = tmp_params

    ! Now read the properties, creating them in the Atoms object as we go
    do i=1,this%properties%N
       
       if (.not. get_value(this%properties, this%properties%keys(i), lookup)) &
            call system_abort('netcdffile_read_atoms: property '//&
            trim(at%properties%keys(i))//' not found in this%properties')

       select case(lookup(2))
          
          case(PROPERTY_INT)
             call add_property(at, this%properties%keys(i), 0, n_cols=lookup(1), lookup=alookup)
             if (lookup(1) == 1) then
                call netcdf_check(nf90_get_var(this%ncid, lookup(3), &
                     at%data%int(alookup(2):alookup(3),1:at%N), &
                     start=(/1,use_frame/), count=(/at%N,1/)))
             else
                call netcdf_check(nf90_get_var(this%ncid, lookup(3), &
                     at%data%int(alookup(2):alookup(3),1:at%N), &
                     start=(/1,1,use_frame/), count=(/lookup(1),at%N,1/)))
             end if

          case(PROPERTY_REAL)
             call add_property(at, this%properties%keys(i), 0.0_dp, n_cols=lookup(1), lookup=alookup)
             if (lookup(1) == 1) then
                call netcdf_check(nf90_get_var(this%ncid, lookup(3), &
                     at%data%real(alookup(2):alookup(3),1:at%N), &
                     start=(/1,use_frame/), count=(/at%N,1/)))
             else
                call netcdf_check(nf90_get_var(this%ncid, lookup(3), &
                     at%data%real(alookup(2):alookup(3),1:at%N), &
                     start=(/1,1,use_frame/), count=(/lookup(1),at%N,1/)))
             end if

          case(PROPERTY_STR)
             call add_property(at, this%properties%keys(i), repeat(" ",TABLE_STRING_LENGTH), &
                  n_cols=lookup(1), lookup=alookup)
             if (lookup(1) == 1) then
                call netcdf_check(nf90_get_var(this%ncid, lookup(3), &
                     at%data%str(alookup(2):alookup(3),1:at%N), &
                     start=(/1,1,use_frame/), count=(/TABLE_STRING_LENGTH,at%N,1/)))
             else
                call netcdf_check(nf90_get_var(this%ncid, lookup(3), &
                     at%data%str(alookup(2):alookup(3),1:at%N), &
                     start=(/1,1,1,use_frame/), count=(/TABLE_STRING_LENGTH,lookup(1),at%N,1/)))
             end if
                
          case(PROPERTY_LOGICAL)
             call add_property(at, this%properties%keys(i), .false., n_cols=lookup(1), lookup=alookup)
             allocate(tmp_integer(alookup(2):alookup(3),1:at%N))
             if (lookup(1) == 1) then
                call netcdf_check(nf90_get_var(this%ncid, lookup(3), &
                     tmp_integer, start=(/1,use_frame/), count=(/at%N,1/)))
             else
                allocate(tmp_integer(alookup(2):alookup(3),1:at%N))
                call netcdf_check(nf90_get_var(this%ncid, lookup(3), &
                     tmp_integer, start=(/1,1,use_frame/), count=(/lookup(1),at%N,1/)))
             end if
             at%data%logical(alookup(2):alookup(3),1:at%N) = tmp_integer
             deallocate(tmp_integer)

       end select
    end do

  end subroutine netcdffile_read_atoms

  
end module netcdf_module
