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

module CInOutput_module

  !% Interface to C routines for reading and writing Atoms objects to and from XYZ and NetCDF files.

  use iso_c_binding
  use error_module
  use System_module, only: dp, optional_default, s2a, a2s, parse_string, print, verbosity_push, PRINT_ANAL, PRINT_VERBOSE, PRINT_ALWAYS, INPUT, OUTPUT, INOUT
  use Atoms_module, only: Atoms, initialise, finalise, add_property, bcast, has_property, set_lattice, atoms_repoint
  use PeriodicTable_module, only: atomic_number_from_symbol, ElementName
  use Extendable_str_module, only: Extendable_str, operator(//), string
  use Dictionary_module, only: Dictionary, has_key, get_value, set_value, print, print_keys, lookup_entry_i, lower_case, subset, &
       T_INTEGER, T_CHAR, T_REAL, T_LOGICAL, T_INTEGER_A, T_REAL_A, T_INTEGER_A2, T_REAL_A2, T_LOGICAL_A, T_CHAR_A
  use Table_module, only: Table, allocate, append, TABLE_STRING_LENGTH
  use MPI_Context_module, only: MPI_context

  implicit none

  private

  interface

     subroutine ciofree(at) bind(c)
       use iso_c_binding, only: C_PTR
       type(C_PTR), intent(in), value :: at
     end subroutine ciofree

     function cioinit(at, filename, action, append, netcdf4, no_compute_index, &
          n_frame, n_atom, n_int, n_real, n_str, n_logical, n_param, n_property, &
          property_name, property_type, property_ncols, property_start, property_filter, &
          param_name, param_type, param_size, param_value, param_int, param_real, param_logical, param_int_a, &
          param_real_a, param_logical_a, param_int_a2, param_real_a2, param_filter, lattice, got_index, pnetcdf4) bind(c)
       use iso_c_binding, only: C_INT, C_CHAR, C_PTR, C_LONG
       type(C_PTR), intent(out) :: at
       character(kind=C_CHAR,len=1), dimension(*), intent(in) :: filename
       integer(kind=C_INT), intent(in) :: action, append, netcdf4, no_compute_index
       type(C_PTR) :: n_frame, n_atom, n_int, n_real, n_str, n_logical, n_param, n_property, &
            property_name, property_type, property_ncols, property_start, property_filter, &
            param_name, param_type, param_size, param_value, &
            param_int, param_real, param_logical, param_int_a, param_real_a, param_logical_a, &
            param_int_a2, param_real_a2, &
            param_filter, lattice, got_index, pnetcdf4
       integer(kind=C_INT) :: cioinit
     end function cioinit

     function cioskip(at, n_skip) bind(c)
       use iso_c_binding, only: C_PTR, C_INT, C_LONG
       type(C_PTR), intent(in), value :: at
       integer(C_INT), intent(in) :: n_skip
       integer(C_INT) :: cioskip
     end function cioskip

     function cioquery(at, frame) bind(c)
       use iso_c_binding, only: C_PTR, C_INT, C_LONG
       type(C_PTR), intent(in), value :: at
       integer(C_INT), intent(in) :: frame
       integer(C_INT) :: cioquery
     end function cioquery

     function cioread(at, frame, int_data, real_data, str_data, logical_data, zero) bind(c)
       use iso_c_binding, only: C_INT, C_PTR, C_DOUBLE, C_CHAR, C_LONG, C_INT
       type(C_PTR), intent(in), value :: at
       integer(C_INT), intent(in) :: frame
       type(C_PTR), intent(in), value :: int_data, real_data, str_data, logical_data
       integer(C_INT), intent(in) :: zero
       integer(C_INT) :: cioread
     end function cioread

     function ciowrite(at, int_data, real_data, str_data, logical_data, prefix, intformat, realformat, frame, &
       shuffle, deflate, deflate_level, swap) bind(c)
       use iso_c_binding, only: C_INT, C_PTR, C_DOUBLE, C_CHAR, C_LONG, C_INT
       type(C_PTR), intent(in), value :: at
       type(C_PTR), intent(in), value :: int_data, real_data, str_data, logical_data
       character(kind=C_CHAR,len=1), dimension(*), intent(in) :: prefix, intformat, realformat
       integer(C_INT) :: frame
       integer(C_INT) :: shuffle, deflate, deflate_level, swap
       integer(C_INT) :: ciowrite
     end function ciowrite

  end interface

  type CInOutput
     type(C_PTR) :: c_at
     logical :: initialised = .false.
     integer :: action
     
     type(C_PTR) :: c_n_frame, c_n_atom, c_n_int, c_n_real, c_n_str, c_n_logical, c_n_param, c_n_property, &
          c_property_name, c_property_type, c_property_start, c_property_ncols, c_property_filter, &
          c_param_name, c_param_type, c_param_size, c_param_value, c_pint, c_preal, c_plogical, c_pint_a, c_preal_a, c_plogical_a, &
          c_pint_a2, c_preal_a2, c_param_filter, c_lattice, c_got_index, c_netcdf4

     integer, pointer :: n_atom, n_frame
     integer, pointer :: n_int, n_real, n_str, n_logical, n_param, n_property, got_index, netcdf4
     integer :: current_frame

     integer, pointer, dimension(:) :: property_type, property_ncols, property_start, &
          param_type, param_size, pint, property_filter, param_filter
     logical, pointer, dimension(:) :: plogical
     integer, pointer, dimension(:,:) :: pint_a
     logical, pointer, dimension(:,:) :: plogical_a
     integer, pointer, dimension(:,:) :: pint_a2
     real(dp), pointer, dimension(:) :: preal
     real(dp), pointer, dimension(:,:) :: preal_a, lattice
     real(dp), pointer, dimension(:,:) :: preal_a2
     character(1), dimension(:,:), pointer :: property_name, param_name, param_value
     type(MPI_Context) :: mpi

  end type CInOutput

  interface initialise
     !% Open a file for reading or writing. File type (Extended XYZ or NetCDF) is guessed from
     !% filename extension. Filename 'stdout' with action 'OUTPUT' to write XYZ to stdout.
     module procedure CInOutput_initialise
  end interface

  interface finalise
     !% Close file and free memory. Calls close() interface.
     module procedure CInOutput_finalise
  end interface

  interface close
     !% Close file. After a call to close(), you can can call initialise() again
     !% to reopen a new file.
     module procedure CInOutput_close
  end interface

  interface read
     !% Read an Atoms object from this CInOutput stream.
     !%
     !% Important properties which may be present (non-exhaustive list)
     !% \begin{itemize}
     !% \item {\bf species}, str, 1 col -- atomic species, e.g. Si or H
     !% \item {\bf pos}, real, 3 cols -- cartesian positions, in A
     !% \item {\bf Z}, int, 1 col -- atomic numbers
     !% \item {\bf mass}, real, 1 col -- atomic masses, in A,eV,fs units system
     !% \item {\bf velo}, real, 3 cols -- velocities, in A/fs
     !% \item {\bf acc}, real, 3 cols -- accelerations, in A/fs$^2$
     !% \item {\bf hybrid}, int, 1 col -- one for QM atoms and zero for hybrid atoms
     !% \item {\bf frac_pos}, real, 3 cols -- fractional positions of atoms
     !% \end{itemize}

     module procedure CInOutput_read
     module procedure atoms_read
     module procedure atoms_read_cinoutput
  end interface

  interface write
     !% Write an Atoms object to this CInOutput stream.
     module procedure CInOutput_write
     module procedure atoms_write
     module procedure atoms_write_cinoutput
  end interface

  interface query
     !% Query the CInOutput stream to fill in values of n_frame, n_atom, n_property etc.
     !% Called once automatically by read(), but if the structure of the XYZ file changes
     !% between frames you need to call it again each time.
     module procedure CInOutput_query
  end interface

  public :: CInOutput, initialise, finalise, close, read, write, query

  ! temporary hack...
  integer, parameter :: PROPERTY_INT     = 1 !% Property types, used by Atoms and Table
  integer, parameter :: PROPERTY_REAL    = 2
  integer, parameter :: PROPERTY_STR     = 3
  integer, parameter :: PROPERTY_LOGICAL = 4
  integer, parameter :: KEY_LEN = 256
  integer, parameter :: VALUE_LEN = 1024

contains

  subroutine cinoutput_initialise(this, filename, action, append, netcdf4, no_compute_index, mpi, error)
    type(CInOutput), intent(inout)  :: this
    character(*), intent(in) :: filename
    integer, intent(in), optional :: action
    logical, intent(in), optional :: append
    logical, optional, intent(in) :: netcdf4
    logical, optional, intent(in) :: no_compute_index
    type(MPI_context), optional, intent(in) :: mpi
    integer, intent(out), optional :: error

    integer :: do_append, do_netcdf4, do_no_compute_index

    INIT_ERROR(error)

    this%action = optional_default(INPUT, action)
    do_netcdf4 = 1
    if (present(netcdf4)) then
       if (.not. netcdf4) do_netcdf4 = 0
    end if
    do_append = 0
    if (present(append)) then
       if (append) do_append = 1
    end if
    do_no_compute_index = 0
    if (present(no_compute_index)) then
       if (no_compute_index) do_no_compute_index = 1
    end if

    if (present(mpi)) then
       this%mpi = mpi
    end if

    if (.not. this%mpi%active .or. (this%mpi%active .and. (this%action /= INPUT .or. (this%action == INPUT .and. this%mpi%my_proc == 0)))) then
       if (cioinit(this%c_at, trim(filename)//C_NULL_CHAR, this%action, do_append, do_netcdf4, do_no_compute_index, &
            this%c_n_frame, this%c_n_atom, this%c_n_int, this%c_n_real, this%c_n_str, this%c_n_logical, this%c_n_param, this%c_n_property, &
            this%c_property_name, this%c_property_type, this%c_property_ncols, this%c_property_start, this%c_property_filter, &
            this%c_param_name, this%c_param_type, this%c_param_size, this%c_param_value, &
            this%c_pint, this%c_preal, this%c_plogical, this%c_pint_a, this%c_preal_a, this%c_plogical_a, this%c_pint_a2, this%c_preal_a2, &
            this%c_param_filter, this%c_lattice, this%c_got_index, this%c_netcdf4) == 0) then
          RAISE_ERROR_WITH_KIND(ERROR_IO,"Error opening file "//filename, error)
       endif

       call c_f_pointer(this%c_n_frame, this%n_frame)
       call c_f_pointer(this%c_n_atom, this%n_atom)
       call c_f_pointer(this%c_n_int, this%n_int)
       call c_f_pointer(this%c_n_real, this%n_real)
       call c_f_pointer(this%c_n_str, this%n_str)
       call c_f_pointer(this%c_n_logical, this%n_logical)
       call c_f_pointer(this%c_n_param, this%n_param)
       call c_f_pointer(this%c_n_property, this%n_property)
       call c_f_pointer(this%c_got_index, this%got_index)
       call c_f_pointer(this%c_netcdf4, this%netcdf4)

       call c_f_pointer(this%c_lattice, this%lattice, (/3,3/))

       call c_f_pointer(this%c_property_name, this%property_name, (/KEY_LEN,this%n_property/))
       call c_f_pointer(this%c_property_type, this%property_type, (/this%n_property/))
       call c_f_pointer(this%c_property_ncols, this%property_ncols, (/this%n_property/))
       call c_f_pointer(this%c_property_start, this%property_start, (/this%n_property/))
       call c_f_pointer(this%c_property_filter, this%property_filter, (/this%n_property/))

       call c_f_pointer(this%c_param_name, this%param_name, (/KEY_LEN,this%n_param/))
       call c_f_pointer(this%c_param_type, this%param_type, (/this%n_param/))
       call c_f_pointer(this%c_param_size, this%param_size, (/this%n_param/))
       call c_f_pointer(this%c_param_value, this%param_value, (/VALUE_LEN,this%n_param/))
       call c_f_pointer(this%c_pint, this%pint, (/this%n_param/))
       call c_f_pointer(this%c_preal, this%preal, (/this%n_param/))
       call c_f_pointer(this%c_plogical, this%plogical, (/this%n_param/))
       call c_f_pointer(this%c_pint_a, this%pint_a, (/3,this%n_param/))
       call c_f_pointer(this%c_preal_a, this%preal_a, (/3,this%n_param/))
       call c_f_pointer(this%c_plogical_a, this%plogical_a, (/3,this%n_param/))
       call c_f_pointer(this%c_pint_a, this%pint_a, (/3,this%n_param/))
       call c_f_pointer(this%c_pint_a2, this%pint_a2, (/9,this%n_param/))
       call c_f_pointer(this%c_preal_a2, this%preal_a2, (/9,this%n_param/))
       call c_f_pointer(this%c_param_filter, this%param_filter, (/this%n_param/))
       
       if (this%action /= INPUT .and. do_append /= 0) then
          this%current_frame = this%n_frame
       else
          this%current_frame = 0
       end if
    end if
    this%initialised = .true.

  end subroutine cinoutput_initialise

  subroutine cinoutput_close(this)
    type(CInOutput), intent(inout) :: this

    if (this%initialised) then
       call ciofree(this%c_at) ! closes files and frees C struct
       this%initialised = .false.
    end if

  end subroutine cinoutput_close

  subroutine cinoutput_finalise(this)
    type(CInOutput), intent(inout) :: this
    call cinoutput_close(this)
  end subroutine cinoutput_finalise

  subroutine cinoutput_query(this, frame, error)
    type(CInOutput), intent(inout) :: this
    integer, optional, intent(in) :: frame
    integer, intent(out), optional :: error
    integer(C_INT) :: do_frame

    integer :: cioquery_status, cioskip_status

    INIT_ERROR(error)

    if (.not. this%initialised) then
       RAISE_ERROR("This CInOutput object is not initialised", error)
    endif
    if (present(frame) .and. this%got_index == 0) then
      if (frame /= this%current_frame) then
	if (frame < this%current_frame) then
	  RAISE_ERROR("cinoutput_query: CInOutput object not seekable and frame argument passed " // frame // " < this%current_frame " // this%current_frame, error)
	endif
	cioskip_status = cioskip(this%c_at, frame-this%current_frame);
	if (cioskip_status == 0) then
          RAISE_ERROR_WITH_KIND(ERROR_IO,"Error querying CInOutput file while skipping to desired frame", error)
	endif
      endif
    endif

    do_frame = optional_default(this%current_frame, frame)

    cioquery_status = cioquery(this%c_at, do_frame)
    if (cioquery_status == 0) then
       RAISE_ERROR_WITH_KIND(ERROR_IO,"Error querying CInOutput file", error)
    endif

    call c_f_pointer(this%c_property_name, this%property_name, (/KEY_LEN,this%n_property/))
    call c_f_pointer(this%c_property_type, this%property_type, (/this%n_property/))
    call c_f_pointer(this%c_property_ncols, this%property_ncols, (/this%n_property/))
    call c_f_pointer(this%c_property_start, this%property_start, (/this%n_property/))
    call c_f_pointer(this%c_property_filter, this%property_filter, (/this%n_property/))

    call c_f_pointer(this%c_param_name, this%param_name, (/KEY_LEN,this%n_param/))
    call c_f_pointer(this%c_param_type, this%param_type, (/this%n_param/))
    call c_f_pointer(this%c_param_size, this%param_size, (/this%n_param/))
    call c_f_pointer(this%c_param_value, this%param_value, (/VALUE_LEN,this%n_param/))
    call c_f_pointer(this%c_pint, this%pint, (/this%n_param/))
    call c_f_pointer(this%c_preal, this%preal, (/this%n_param/))
    call c_f_pointer(this%c_plogical, this%plogical, (/this%n_param/))
    call c_f_pointer(this%c_pint_a, this%pint_a, (/3,this%n_param/))
    call c_f_pointer(this%c_preal_a, this%preal_a, (/3,this%n_param/))
    call c_f_pointer(this%c_plogical_a, this%plogical_a, (/3,this%n_param/))
    call c_f_pointer(this%c_pint_a2, this%pint_a2, (/9,this%n_param/))
    call c_f_pointer(this%c_preal_a2, this%preal_a2, (/9,this%n_param/))
    call c_f_pointer(this%c_param_filter, this%param_filter, (/this%n_param/))

  end subroutine cinoutput_query

  subroutine cinoutput_read(this, at, frame, zero, error)
    use iso_c_binding, only: C_INT
    type(CInOutput), intent(inout) :: this
    type(Atoms), target, intent(out) :: at
    integer, optional, intent(in) :: frame
    logical, optional, intent(in) :: zero
    integer, intent(out), optional :: error

    type(Table), target :: tmp_data
    type(Dictionary) :: empty_dictionary
    type(C_PTR) :: int_ptr, real_ptr, str_ptr, log_ptr
    integer :: i
    character(len=KEY_LEN) :: namestr
    integer :: do_zero
    integer(C_INT) :: do_frame
    integer :: tmp_do_frame
    integer :: n_skip
    integer :: cioskip_status

    INIT_ERROR(error)
call verbosity_push(PRINT_ANAL)

    if (.not. this%initialised) then
       RAISE_ERROR_WITH_KIND(ERROR_IO,"This CInOutput object is not initialised", error)
    endif
    if (this%action /= INPUT .and. this%action /= INOUT) then
       RAISE_ERROR_WITH_KIND(ERROR_IO,"Cannot read from action=OUTPUT CInOutput object", error)
    endif

    if (.not. this%mpi%active .or. (this%mpi%active .and. this%mpi%my_proc == 0)) then

       call print('cinoutput_read: doing read on proc '//this%mpi%my_proc, PRINT_VERBOSE)

       if (present(frame) .and. this%got_index == 0) then
          if (frame /= this%current_frame) then
             if (frame < this%current_frame) then
                  BCAST_RAISE_ERROR_WITH_KIND(ERROR_IO,"cinoutput_read: CInOutput object not seekable and frame argument passed " // frame // " < this%current_frame " // this%current_frame, error, this%mpi)
	     endif
             n_skip = frame-this%current_frame
             cioskip_status = cioskip(this%c_at, n_skip)
             if (cioskip_status == 0) then
                BCAST_RAISE_ERROR_WITH_KIND(ERROR_IO,"Error querying CInOutput file while skipping to desired frame", error, this%mpi)
             endif
             this%current_frame = this%current_frame + n_skip
          endif
       endif

       do_frame = optional_default(this%current_frame, frame)

       if (this%got_index == 1) then
          if (do_frame < 0) do_frame = this%n_frame + do_frame ! negative frames count backwards from end
          if (do_frame < 0 .or. do_frame >= this%n_frame) then
             call finalise(tmp_data)
             call finalise(at)
             BCAST_RAISE_ERROR_WITH_KIND(ERROR_IO_EOF,"cinoutput_read: frame "//int(do_frame)//" out of range 0 <= frame < "//int(this%n_frame), error, this%mpi)
          end if
       end if

       do_zero = 0
       if (present(zero)) then
          if (zero) do_zero = 1
       end if

       if (this%got_index == 1) then
          tmp_do_frame = do_frame
          call cinoutput_query(this, tmp_do_frame, error=error)
       else
          call cinoutput_query(this, error=error)
       end if
       BCAST_PASS_ERROR(error, this%mpi)


       ! Make a blank data table, then initialise atoms object from it. We will
       ! get C routine to read directly into final Atoms structure to save copying.

       call allocate(tmp_data, this%n_int, this%n_real, this%n_str, this%n_logical, int(this%n_atom))
       call append(tmp_data, blank_rows=int(this%n_atom))
       this%lattice = 0.0_dp
       this%lattice(1,1) = 1.0_dp
       this%lattice(2,2) = 1.0_dp
       this%lattice(3,3) = 1.0_dp

       int_ptr = C_NULL_PTR
       if (this%n_int /= 0) int_ptr = c_loc(tmp_data%int(1,1))

       real_ptr = C_NULL_PTR
       if (this%n_real /= 0) real_ptr = c_loc(tmp_data%real(1,1))

       str_ptr = C_NULL_PTR
       if (this%n_str /= 0) str_ptr = c_loc(tmp_data%str(1,1))

       log_ptr = C_NULL_PTR
       if (this%n_logical /= 0) log_ptr = c_loc(tmp_data%logical(1,1))

       if (cioread(this%c_at, do_frame, int_ptr, real_ptr, str_ptr, log_ptr, do_zero) == 0) then
          call finalise(tmp_data)
          BCAST_RAISE_ERROR_WITH_KIND(ERROR_IO,"Error reading from file", error, this%mpi)
       end if

       call initialise(empty_dictionary) ! pass an empty dictionary to prevent pos, species, z properties being created
       call initialise(at, int(this%n_atom), transpose(this%lattice), properties=empty_dictionary)
       call finalise(empty_dictionary)

       ! temporary hack - we copy from data Table into new properties Dictionary
       do i=1,this%n_property
          select case(this%property_type(i))
          case(PROPERTY_INT)
             if (this%property_ncols(i) == 1) then
                call add_property(at, c_array_to_f_string(this%property_name(:,i)), &
                     tmp_data%int(this%property_start(i)+1,1:at%n), overwrite=.true.,error=error)
                BCAST_PASS_ERROR(error, this%mpi)
             else
                call add_property(at, c_array_to_f_string(this%property_name(:,i)), &
                     tmp_data%int(this%property_start(i)+1:this%property_start(i)+this%property_ncols(i),1:at%n), &
                     overwrite=.true.,error=error)
                BCAST_PASS_ERROR(error, this%mpi)
             end if
             
          case(PROPERTY_REAL)
             if (this%property_ncols(i) == 1) then
                call add_property(at, c_array_to_f_string(this%property_name(:,i)), &
                     tmp_data%real(this%property_start(i)+1,1:at%n), overwrite=.true.,error=error)
                BCAST_PASS_ERROR(error, this%mpi)
             else
                call add_property(at, c_array_to_f_string(this%property_name(:,i)), &
                     tmp_data%real(this%property_start(i)+1:this%property_start(i)+this%property_ncols(i),1:at%n), &
                     overwrite=.true., error=error)
                BCAST_PASS_ERROR(error, this%mpi)
             end if
             
          case(PROPERTY_LOGICAL)
             if (this%property_ncols(i) == 1) then
                call add_property(at, c_array_to_f_string(this%property_name(:,i)), &
                     tmp_data%logical(this%property_start(i)+1,1:at%n), overwrite=.true.,error=error)
                BCAST_PASS_ERROR(error, this%mpi)
             else
                BCAST_RAISE_ERROR("CInOutput_read: logical properties with n_cols != 1 no longer supported", error, this%mpi)
             end if
             
          case(PROPERTY_STR)
             if (this%property_ncols(i) == 1) then
                call add_property(at, c_array_to_f_string(this%property_name(:,i)), &
                     tmp_data%str(this%property_start(i)+1,1:at%n), overwrite=.true., error=error)
                BCAST_PASS_ERROR(error, this%mpi)
             else
                BCAST_RAISE_ERROR("CInOutput_read: string properties with n_cols != 1 no longer supported", error, this%mpi)
             end if
             
          case default
             BCAST_RAISE_ERROR("CInOutput_read: unknown property type "//this%property_type(i), error, this%mpi)
          end select
       end do
       
       do i=1,this%n_param
          namestr = c_array_to_f_string(this%param_name(:,i))
          if (trim(namestr) == "Lattice" .or. trim(namestr) == "Properties") cycle
          select case(this%param_type(i))
          case(T_INTEGER)
             call set_value(at%params, namestr, this%pint(i))
          case(T_REAL)
             call set_value(at%params, namestr, this%preal(i))
          case(T_LOGICAL)
             call set_value(at%params, namestr, this%plogical(i))
          case(T_INTEGER_A)
             call set_value(at%params, namestr, this%pint_a(:,i))
          case (T_REAL_A)
             call set_value(at%params, namestr, this%preal_a(:,i))
          case (T_LOGICAL_A)
             call set_value(at%params, namestr, this%plogical_a(:,i))
          case(T_CHAR)
             call set_value(at%params, namestr, c_array_to_f_string(this%param_value(:,i)))
          case(T_INTEGER_A2)
             call set_value(at%params, namestr, reshape(this%pint_a2(:,i), (/3,3/)))
          case(T_REAL_A2)
             call set_value(at%params, namestr, reshape(this%preal_a2(:,i), (/3,3/)))
          case default
             BCAST_RAISE_ERROR('cinoutput_read: unsupported parameter i='//i//' key='//trim(namestr)//' type='//this%param_type(i), error, this%mpi)
          end select
       end do

       call atoms_repoint(at)
       call set_lattice(at, transpose(this%lattice), scale_positions=.false.)

       if (.not. has_property(at,"Z") .and. .not. has_property(at, "species")) then
          call print ("at%properties keys", PRINT_ALWAYS)
          call print_keys(at%properties, PRINT_ALWAYS)
          BCAST_RAISE_ERROR('cinoutput_read: atoms object read from file has neither Z nor species', error, this%mpi)
       else if (.not. has_property(at,"species") .and. has_property(at,"Z")) then
          call add_property(at, "species", repeat(" ",TABLE_STRING_LENGTH))
          call atoms_repoint(at)
          do i=1,at%N
             at%species(:,i) = s2a(ElementName(at%Z(i)))
          end do
       else if (.not. has_property(at,"Z") .and. has_property(at, "species")) then
          call add_property(at, "Z", 0)
          call atoms_repoint(at)
          do i=1,at%n
             at%Z(i) = atomic_number_from_symbol(a2s(at%species(:,i)))
          end do
       end if

       call finalise(tmp_data)
       this%current_frame = this%current_frame + 1

    end if

    if (this%mpi%active) then
       BCAST_CHECK_ERROR(error,this%mpi)
       call bcast(this%mpi, at)
    end if

  end subroutine cinoutput_read


  subroutine cinoutput_write(this, at, properties, properties_array, prefix, int_format, real_format, frame, &
       shuffle, deflate, deflate_level, error)
    type(CInOutput), intent(inout) :: this
    type(Atoms), target, intent(in) :: at
    character(*), intent(in), optional :: properties
    character(*), intent(in), optional :: properties_array(:)
    character(*), intent(in), optional :: prefix
    character(*), intent(in), optional :: int_format, real_format
    integer, intent(in), optional :: frame
    logical, intent(in), optional :: shuffle, deflate
    integer, intent(in), optional :: deflate_level
    integer, intent(out), optional :: error

    type(C_PTR) :: int_ptr, real_ptr, str_ptr, log_ptr
    logical :: dum, got_properties
    character(len=VALUE_LEN) :: valuestr
    integer :: i, j, k, extras, n, tmp_int_a2(3,3)
    real(dp) :: tmp_real_a2(3,3)
    integer(C_INT) :: do_frame
    type(Dictionary) :: selected_properties
    character(len=KEY_LEN) :: do_prefix, do_int_format, do_real_format
    integer :: do_shuffle, do_deflate
    integer :: do_deflate_level, do_swap
    type(Table), target :: tmp_data
    character(len=KEY_LEN) :: tmp_properties_array(at%properties%n)
    integer n_properties

    INIT_ERROR(error)

    if (.not. this%initialised) then
      RAISE_ERROR("This CInOutput object is not initialised", error)
    endif
    if (this%action /= OUTPUT .and. this%action /= INOUT) then
      RAISE_ERROR("Cannot write to action=INPUT CInOutput object", error)
    endif

    if (this%mpi%active .and. this%mpi%my_proc /= 0) return

    do_prefix = optional_default('', prefix)
    do_int_format = optional_default('%8d', int_format)
    do_real_format = optional_default('%16.8f', real_format)

    do_frame = optional_default(this%current_frame, frame)

    do_shuffle = transfer(optional_default(.true., shuffle),do_shuffle)
    do_deflate = transfer(optional_default(.true., deflate),do_shuffle)
    do_deflate_level = optional_default(6, deflate_level)
    got_properties = present(properties) .or. present(properties_array)
    do_swap = transfer(.not. got_properties, do_swap)

    if (present(properties) .and. present(properties_array)) then
       RAISE_ERROR('CInOutput_write: "properties" and "properties_array" cannot both be present.', error)
    end if
    
    call initialise(selected_properties)
    if (.not. got_properties) then
       do i=1,at%properties%N
          call set_value(selected_properties, string(at%properties%keys(i)), 1)
       end do
    else
       if (present(properties_array)) then
          call subset(at%properties, properties_array, selected_properties, error=error)
          PASS_ERROR(error)
       else ! we've got a colon-separated string
          call parse_string(properties, ':', tmp_properties_array, n_properties, error=error)
          PASS_ERROR(error)
          call subset(at%properties, tmp_properties_array(1:n_properties), selected_properties, error=error)
          PASS_ERROR(error)
       end if
    end if

    this%n_atom = at%n
    this%lattice = transpose(at%lattice)

    extras = 0
    if (has_key(at%params, "Lattice"))    extras = extras + 1
    if (has_key(at%params, "Properties")) extras = extras + 1
    
    this%n_param = at%params%N - extras
    call c_f_pointer(this%c_param_name, this%param_name, (/KEY_LEN,this%n_param/))
    call c_f_pointer(this%c_param_type, this%param_type, (/this%n_param/))
    call c_f_pointer(this%c_param_size, this%param_size, (/this%n_param/))
    call c_f_pointer(this%c_param_value, this%param_value, (/VALUE_LEN,this%n_param/))
    call c_f_pointer(this%c_pint, this%pint, (/this%n_param/))
    call c_f_pointer(this%c_preal, this%preal, (/this%n_param/))
    call c_f_pointer(this%c_plogical, this%plogical, (/this%n_param/))
    call c_f_pointer(this%c_pint_a, this%pint_a, (/3,this%n_param/))
    call c_f_pointer(this%c_preal_a, this%preal_a, (/3,this%n_param/))
    call c_f_pointer(this%c_plogical_a, this%plogical_a, (/3,this%n_param/))
    call c_f_pointer(this%c_pint_a2, this%pint_a2, (/9,this%n_param/))
    call c_f_pointer(this%c_preal_a2, this%preal_a2, (/9,this%n_param/))
    call c_f_pointer(this%c_param_filter, this%param_filter, (/this%n_param/))
    n = 1
    do i=1, at%params%N
       if (lower_case(string(at%params%keys(i))) == lower_case('Lattice') .or. &
           lower_case(string(at%params%keys(i))) == lower_case('Properties')) cycle
       call f_string_to_c_array(string(at%params%keys(i)), this%param_name(:,n))
       this%param_filter(n) = 1
       select case(at%params%entries(i)%type)
       case(T_INTEGER)
          dum = get_value(at%params, string(at%params%keys(i)), this%pint(n))
          this%param_size(n) = 1
          this%param_type(n) = T_INTEGER
       case(T_REAL)
          dum = get_value(at%params, string(at%params%keys(i)), this%preal(n))
          this%param_size(n) = 1
          this%param_type(n) = T_REAL
       case(T_LOGICAL)
          dum = get_value(at%params, string(at%params%keys(i)), this%plogical(n))
          this%param_size(n) = 1
          this%param_type(n) = T_LOGICAL
       case(T_INTEGER_A)
          dum = get_value(at%params, string(at%params%keys(i)), this%pint_a(:,n))
          this%param_size(n) = 3
          this%param_type(n) = T_INTEGER_A
       case (T_REAL_A)
          dum = get_value(at%params, string(at%params%keys(i)), this%preal_a(:,n))
          this%param_size(n) = 3
          this%param_type(n) = T_REAL_A
       case(T_LOGICAL_A)
          dum = get_value(at%params, string(at%params%keys(i)), this%plogical_a(:,n))
          this%param_size(n) = 3
          this%param_type(n) = T_LOGICAL_A
       case(T_CHAR)
          dum = get_value(at%params, string(at%params%keys(i)), valuestr)
          call f_string_to_c_array(valuestr, this%param_value(:,n))
          this%param_size(n) = 1
          this%param_type(n) = T_CHAR
       case(T_INTEGER_A2)
          if (.not. get_value(at%params, string(at%params%keys(i)), tmp_int_a2)) then
               RAISE_ERROR('Bad atoms parameter '//string(at%params%keys(i))//' should be 3x3 int array',error)
	  endif
          this%pint_a2(:,n) = reshape(tmp_int_a2, (/9/))
          this%param_size(n) = 9
          this%param_type(n) = T_INTEGER_A2
       case (T_REAL_A2)
          if (.not. get_value(at%params, string(at%params%keys(i)), tmp_real_a2)) then
               RAISE_ERROR('Bad atoms parameter '//string(at%params%keys(i))//' should be 3x3 real array',error)
	  endif
          this%preal_a2(:,n) = reshape(tmp_real_a2, (/9/))
          this%param_size(n) = 9
          this%param_type(n) = T_REAL_A2
       case default
          RAISE_ERROR('cinoutput_write: unsupported parameter i='//i//' '//string(at%params%keys(i))//' type='//at%params%entries(i)%type, error)
       end select
       n = n + 1
    end do

    this%n_property = at%properties%N
    call c_f_pointer(this%c_property_name, this%property_name, (/KEY_LEN,this%n_property/))
    call c_f_pointer(this%c_property_type, this%property_type, (/this%n_property/))
    call c_f_pointer(this%c_property_ncols, this%property_ncols, (/this%n_property/))
    call c_f_pointer(this%c_property_start, this%property_start, (/this%n_property/))
    call c_f_pointer(this%c_property_filter, this%property_filter, (/this%n_property/))

    ! temporary hack
    this%n_int = 0
    this%n_real = 0
    this%n_logical = 0
    this%n_str = 0
    do i=1,selected_properties%n
       j = lookup_entry_i(at%properties, string(selected_properties%keys(i)))
       call f_string_to_c_array(string(at%properties%keys(j)), this%property_name(:,i))

       select case(at%properties%entries(j)%type)
       case(T_INTEGER_A)
          this%property_type(i) = PROPERTY_INT
          this%property_start(i) = this%n_int
          this%property_ncols(i) = 1
          this%n_int = this%n_int + 1
          this%property_filter(i) = 1
          
       case(T_REAL_A)
          this%property_type(i) = PROPERTY_REAL
          this%property_start(i) = this%n_real
          this%property_ncols(i) = 1
          this%n_real = this%n_real + 1
          this%property_filter(i) = 1

       case(T_LOGICAL_A)
          this%property_type(i) = PROPERTY_LOGICAL
          this%property_start(i) = this%n_logical
          this%property_ncols(i) = 1
          this%n_logical = this%n_logical + 1
          this%property_filter(i) = 1
          
       case(T_INTEGER_A2)
          this%property_type(i) = PROPERTY_INT
          this%property_start(i) = this%n_int
          this%property_ncols(i) = at%properties%entries(j)%len2(1)
          this%n_int = this%n_int + this%property_ncols(i)
          this%property_filter(i) = 1

       case(T_REAL_A2)
          this%property_type(i) = PROPERTY_REAL
          this%property_start(i) = this%n_real
          this%property_ncols(i) = at%properties%entries(j)%len2(1)
          this%n_real = this%n_real + this%property_ncols(i)
          this%property_filter(i) = 1

       case(T_CHAR_A)
          this%property_type(i) = PROPERTY_STR
          this%property_start(i) = this%n_str
          this%property_ncols(i) = 1
          this%n_str = this%n_str + 1
          this%property_filter(i) = 1          

       case default
          RAISE_ERROR("CInOutput_write: bad property type "//at%properties%entries(j)%type//" for property "//selected_properties%keys(i), error)

       end select
    end do
    j = selected_properties%n+1
    do i=1,at%properties%n
       if (has_key(selected_properties, string(at%properties%keys(i)))) cycle
       call f_string_to_c_array(string(at%properties%keys(i)), this%property_name(:,j))

       select case(at%properties%entries(i)%type)
       case(T_INTEGER_A)
          this%property_type(j) = PROPERTY_INT
          this%property_start(j) = this%n_int
          this%property_ncols(j) = 1
          this%n_int = this%n_int + 1
          this%property_filter(j) = 0
          
       case(T_REAL_A)
          this%property_type(j) = PROPERTY_REAL
          this%property_start(j) = this%n_real
          this%property_ncols(j) = 1
          this%n_real = this%n_real + 1
          this%property_filter(j) = 0

       case(T_LOGICAL_A)
          this%property_type(j) = PROPERTY_LOGICAL
          this%property_start(j) = this%n_logical
          this%property_ncols(j) = 1
          this%n_logical = this%n_logical + 1
          this%property_filter(j) = 0
          
       case(T_INTEGER_A2)
          this%property_type(j) = PROPERTY_INT
          this%property_start(j) = this%n_int
          this%property_ncols(j) = at%properties%entries(i)%len2(1)
          this%n_int = this%n_int + this%property_ncols(j)
          this%property_filter(j) = 0

       case(T_REAL_A2)
          this%property_type(j) = PROPERTY_REAL
          this%property_start(j) = this%n_real
          this%property_ncols(j) = at%properties%entries(i)%len2(1)
          this%n_real = this%n_real + this%property_ncols(j)
          this%property_filter(j) = 0

       case(T_CHAR_A)
          this%property_type(j) = PROPERTY_STR
          this%property_start(j) = this%n_str
          this%property_ncols(j) = 1
          this%n_str = this%n_str + 1
          this%property_filter(j) = 0

       case default
          RAISE_ERROR("CInOutput_write: bad property type "//at%properties%entries(i)%type//" for property "//at%properties%keys(i), error)

       end select
       j = j + 1
    end do

    call allocate(tmp_data, this%n_int, this%n_real, this%n_str, this%n_logical)
    call append(tmp_data, blank_rows=at%n)

    ! now copy data from at%properties Dictionary into temporaray data Table
    
    do i=1,this%n_property
       j = lookup_entry_i(at%properties, c_array_to_f_string(this%property_name(:,i)))

       select case(at%properties%entries(j)%type)
       case(T_INTEGER_A)
          tmp_data%int(this%property_start(i)+1,1:at%n) = at%properties%entries(j)%i_a(:)
          
       case(T_REAL_A)
          tmp_data%real(this%property_start(i)+1,1:at%n) = at%properties%entries(j)%r_a(:)

       case(T_LOGICAL_A)
          tmp_data%logical(this%property_start(i)+1,1:at%n) = at%properties%entries(j)%l_a(:)
          
       case(T_INTEGER_A2)
          tmp_data%int(this%property_start(i)+1:this%property_start(i)+this%property_ncols(i),1:at%n) = at%properties%entries(j)%i_a2(:,:)

       case(T_REAL_A2)
          tmp_data%real(this%property_start(i)+1:this%property_start(i)+this%property_ncols(i),1:at%n) = at%properties%entries(j)%r_a2(:,:)
          
       case(T_CHAR_A)
          do k=1,at%n
             tmp_data%str(this%property_start(i)+1:this%property_start(i)+this%property_ncols(i),k) = a2s(at%properties%entries(j)%s_a(:,k))
          end do

       case default
          RAISE_ERROR("CInOutput_write: bad property type "//at%properties%entries(j)%type//" for property "//at%properties%keys(j), error)

       end select
    end do

    int_ptr = C_NULL_PTR
    if (this%n_int /= 0) int_ptr = c_loc(tmp_data%int(1,1))

    real_ptr = C_NULL_PTR
    if (this%n_real /= 0) real_ptr = c_loc(tmp_data%real(1,1))
    
    str_ptr = C_NULL_PTR
    if (this%n_str /= 0) str_ptr = c_loc(tmp_data%str(1,1))

    log_ptr = C_NULL_PTR
    if (this%n_logical /= 0) log_ptr = c_loc(tmp_data%logical(1,1))

    if (ciowrite(this%c_at, int_ptr, real_ptr, str_ptr, log_ptr, &
         trim(do_prefix)//C_NULL_CHAR, &
         trim(do_int_format)//C_NULL_CHAR, trim(do_real_format)//C_NULL_CHAR, do_frame, &
         do_shuffle, do_deflate, do_deflate_level, do_swap) == 0) then
       RAISE_ERROR_WITH_KIND(ERROR_IO,"Error writing file.", error)
    end if

    call finalise(selected_properties)
    call finalise(tmp_data)
    this%current_frame = this%current_frame + 1

  end subroutine cinoutput_write

  function c_array_to_f_string(carray) result(fstring)
    !% Convert a null-terminated array of characters in a fixed length buffer
    !% to a blank-padded fortran string

    character(1), dimension(:), intent(in) :: carray
    character(len=size(carray)) :: fstring

    integer :: i

    fstring = repeat(' ',len(carray))
    do i=1,size(carray)
       if (carray(i) == C_NULL_CHAR) exit
       fstring(i:i) = carray(i)
    end do

  end function c_array_to_f_string

  function c_string_to_f_string(cstring) result(fstring)
    !% Convert a null-terminated string to a blank-padded fortran string

    character(*) :: cstring
    character(len=len(cstring)) :: fstring

    integer :: i

    fstring = repeat(' ',len(cstring))
    do i=1,len(cstring)
       if (cstring(i:i) == C_NULL_CHAR) exit
       fstring(i:i) = cstring(i:i)
    end do

  end function c_string_to_f_string

  subroutine f_string_to_c_array(fstring, carray)
    !% Convert a blank-padded fortran string to a null terminated char array of the same length
    character(*), intent(in) :: fstring
    character(1), dimension(len(fstring)+1), intent(out) :: carray

    integer :: i

    do i=1,len_trim(fstring)
       carray(i) = fstring(i:i)
    end do
    carray(len_trim(fstring)+1) = C_NULL_CHAR

  end subroutine f_string_to_c_array


  subroutine atoms_read(this, filename, frame, zero, mpi, error)
    !% Read Atoms object from XYZ or NetCDF file.
    type(Atoms), intent(inout) :: this
    character(len=*), intent(in) :: filename
    integer, optional, intent(in) :: frame
    logical, optional, intent(in) :: zero
    type(MPI_context), optional, intent(inout) :: mpi
    integer, intent(out), optional :: error

    type(CInOutput) :: cio

    INIT_ERROR(error)

    call initialise(cio, filename, INPUT, mpi=mpi, error=error)
    PASS_ERROR_WITH_INFO('While reading "' // filename // '".', error)
    call read(cio, this, frame, zero, error=error)
    PASS_ERROR_WITH_INFO('While reading "' // filename // '".', error)
    call finalise(cio)

  end subroutine atoms_read

  subroutine atoms_read_cinoutput(this, cio, frame, zero, error)
    type(Atoms), target, intent(inout) :: this
    type(CInOutput), intent(inout) :: cio
    integer, optional, intent(in) :: frame
    logical, optional, intent(in) :: zero
    integer, intent(out), optional :: error

    INIT_ERROR(error)
    call cinoutput_read(cio, this, frame, zero, error=error)
    PASS_ERROR(error)

  end subroutine atoms_read_cinoutput

  subroutine atoms_write(this, filename, append, properties, properties_array, prefix, int_format, real_format, error)
    !% Write Atoms object to XYZ or NetCDF file. Use filename "stdout" to write to terminal.
    use iso_fortran_env
    type(Atoms), intent(in) :: this
    character(len=*), intent(in) :: filename
    logical, optional, intent(in) :: append
    character(*), intent(in), optional :: properties
    character(*), intent(in), optional :: properties_array(:)
    character(*), intent(in), optional :: prefix
    character(*), intent(in), optional :: int_format, real_format
    integer, intent(out), optional :: error

    type(CInOutput) :: cio
    type(MPI_context) :: mpi

    if (trim(filename) == 'stdout') then
       call initialise(mpi)
       if (mpi%active .and. mpi%my_proc /= 0) return
       flush(output_unit)
    end if

    INIT_ERROR(error)
    call initialise(cio, filename, OUTPUT, append, error=error)
    PASS_ERROR_WITH_INFO('While writing "' // filename // '".', error)
    call write(cio, this, properties, properties_array, prefix, int_format, real_format, error=error)
    PASS_ERROR_WITH_INFO('While writing "' // filename // '".', error)
    call finalise(cio)

  end subroutine atoms_write

  subroutine atoms_write_cinoutput(this, cio, properties, properties_array, prefix, int_format, real_format, frame, shuffle, deflate, deflate_level, error)
    type(Atoms), target, intent(inout) :: this
    type(CInOutput), intent(inout) :: cio
    character(*), intent(in), optional :: properties
    character(*), intent(in), optional :: properties_array(:)    
    character(*), intent(in), optional :: prefix, int_format, real_format
    integer, intent(in), optional :: frame
    logical, intent(in), optional :: shuffle, deflate
    integer, intent(in), optional :: deflate_level
    integer, intent(out), optional :: error

    INIT_ERROR(error)
    call cinoutput_write(cio, this, properties, properties_array, prefix, int_format, real_format, frame, shuffle, deflate, deflate_level, error=error)
    PASS_ERROR(error)

  end subroutine atoms_write_cinoutput

end module CInOutput_module
