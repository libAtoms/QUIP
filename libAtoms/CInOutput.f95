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
  !% The XYZ routines are much faster (around a factor of a hundred!) than atoms_print_xyz and 
  !% atoms_read_xyz. 

  use iso_c_binding
  use error_module
  use Atoms_module
  use Dictionary_module
  use Table_module
  use MPI_Context_module

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

     function ciowrite(at, int_data, real_data, str_data, logical_data, intformat, realformat, frame, &
       shuffle, deflate, deflate_level, swap) bind(c)
       use iso_c_binding, only: C_INT, C_PTR, C_DOUBLE, C_CHAR, C_LONG, C_INT
       type(C_PTR), intent(in), value :: at
       type(C_PTR), intent(in), value :: int_data, real_data, str_data, logical_data
       character(kind=C_CHAR,len=1), dimension(*), intent(in) :: intformat, realformat
       integer(C_INT) :: frame
       integer(C_INT) :: shuffle, deflate, deflate_level, swap
       integer(C_INT) :: ciowrite
     end function ciowrite

     function cioupdate(at, int_data, real_data, str_data, logical_data) bind(c)
       use iso_c_binding, only: C_INT, C_PTR, C_DOUBLE, C_CHAR, C_LONG
       type(C_PTR), intent(in), value :: at
       type(C_PTR), intent(in), value :: int_data, real_data, str_data, logical_data
       integer(C_INT) :: cioupdate
     end function cioupdate

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

contains

  subroutine cinoutput_initialise(this, filename, action, append, netcdf4, no_compute_index, mpi, ierror)
    type(CInOutput), intent(inout)  :: this
    character(*), intent(in), optional :: filename
    integer, intent(in), optional :: action
    logical, intent(in), optional :: append
    logical, optional, intent(in) :: netcdf4
    logical, optional, intent(in) :: no_compute_index
    type(MPI_context), optional, intent(in) :: mpi
    integer, optional, intent(inout) :: ierror

    integer :: do_append, do_netcdf4, do_no_compute_index

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

    if (present(filename)) then
       if (cioinit(this%c_at, trim(filename)//C_NULL_CHAR, this%action, do_append, do_netcdf4, do_no_compute_index, &
            this%c_n_frame, this%c_n_atom, this%c_n_int, this%c_n_real, this%c_n_str, this%c_n_logical, this%c_n_param, this%c_n_property, &
            this%c_property_name, this%c_property_type, this%c_property_ncols, this%c_property_start, this%c_property_filter, &
            this%c_param_name, this%c_param_type, this%c_param_size, this%c_param_value, &
            this%c_pint, this%c_preal, this%c_plogical, this%c_pint_a, this%c_preal_a, this%c_plogical_a, this%c_pint_a2, this%c_preal_a2, &
            this%c_param_filter, this%c_lattice, this%c_got_index, this%c_netcdf4) == 0) &
            RAISE_ERROR("Error opening file "//filename, ierror)
    else
       if (cioinit(this%c_at, ""//C_NULL_CHAR, this%action, do_append, do_netcdf4, do_no_compute_index, &
            this%c_n_frame, this%c_n_atom, this%c_n_int, this%c_n_real, this%c_n_str, this%c_n_logical, this%c_n_param, this%c_n_property, &
            this%c_property_name, this%c_property_type, this%c_property_ncols, this%c_property_start, this%c_property_filter, &
            this%c_param_name, this%c_param_type, this%c_param_size, this%c_param_value, &
            this%c_pint, this%c_preal, this%c_plogical, this%c_pint_a, this%c_preal_a, this%c_plogical_a, this%c_pint_a2, this%c_preal_a2, &
            this%c_param_filter, this%c_lattice, this%c_got_index, this%c_netcdf4) == 0) &
            RAISE_ERROR("Error allocating C structure", ierror)
    end if

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

    this%initialised = .true.
    
    if (this%action /= INPUT .and. do_append /= 0) then
       this%current_frame = this%n_frame
    else
       this%current_frame = 0
    end if

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

  subroutine cinoutput_query(this, frame, ierror)
    type(CInOutput), intent(inout) :: this
    integer, optional, intent(in) :: frame
    integer, optional, intent(out) :: ierror
    integer(C_INT) :: do_frame

    integer :: cioquery_status, cioskip_status

    if (.not. this%initialised) call system_abort("This CInOutput object is not initialised")
    if (present(frame) .and. this%got_index == 0) then
      if (frame /= this%current_frame) then
	if (frame < this%current_frame) &
	  call system_abort("cinoutput_query: CInOutput object not seekable and frame argument passed " // frame // " < this%current_frame " // this%current_frame)
	cioskip_status = cioskip(this%c_at, frame-this%current_frame);
	if (cioskip_status == 0) then
! Lars
!	  if (present(status)) then
!	    status = 1
!	    return
!	  else
!	    call system_abort("Error querying CInOutput file while skipping to desired frame")
!	  endif
          RAISE_ERROR("Error querying CInOutput file while skipping to desired frame", ierror)
	endif
      endif
    endif

    do_frame = optional_default(this%current_frame, frame)

    cioquery_status = cioquery(this%c_at, do_frame)
    if (cioquery_status == 0) then
! Lars
!      if (present(status)) then
!	status = 1
!	return
!      else
!        call system_abort("Error querying CInOutput file")
!      endif
       RAISE_ERROR("Error querying CInOutput file", ierror)
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

  subroutine cinoutput_read(this, at, frame, zero, ierror)
    use iso_c_binding, only: C_INT
    type(CInOutput), intent(inout) :: this
    type(Atoms), target, intent(out) :: at
    integer, optional, intent(in) :: frame
    logical, optional, intent(in) :: zero
    integer, optional, intent(inout) :: ierror

    type(Dictionary) :: properties
    type(Table) :: data
    type(C_PTR) :: int_ptr, real_ptr, str_ptr, log_ptr
    integer :: i
    character(len=KEY_LEN) :: namestr
    integer :: do_zero
    integer(C_INT) :: do_frame
    integer :: tmp_do_frame
    integer :: n_skip
    integer :: cioskip_status

! Lars
!    if (present(status)) status = 0

    if (.not. this%initialised) call system_abort("This CInOutput object is not initialised")
    if (this%action /= INPUT .and. this%action /= INOUT) call system_abort("Cannot read from action=OUTPUT CInOutput object")

    if (.not. this%mpi%active .or. (this%mpi%active .and. this%mpi%my_proc == 0)) then

       call print('cinoutput_read: doing read on proc '//this%mpi%my_proc, VERBOSE)

       if (present(frame) .and. this%got_index == 0) then
          if (frame /= this%current_frame) then
             if (frame < this%current_frame) &
                  call system_abort("cinoutput_read: CInOutput object not seekable and frame argument passed " // frame // " < this%current_frame " // this%current_frame)
             n_skip = frame-this%current_frame
             cioskip_status = cioskip(this%c_at, n_skip)
             if (cioskip_status == 0) then
                RAISE_ERROR("Error querying CInOutput file while skipp!ing to desired frame", ierror)
! Lars
!                if (present(status)) then
!                   status = 1
!                   return
!                else
!                   call system_abort("Error querying CInOutput file while skipp!ing to desired frame")
!                endif
             endif
             this%current_frame = this%current_frame + n_skip
          endif
       endif

       do_frame = optional_default(this%current_frame, frame)

       if (this%got_index == 1) then
          if (do_frame < 0) do_frame = this%n_frame + do_frame ! negative frames count backwards from end
          if (do_frame < 0 .or. do_frame >= this%n_frame) then
             call finalise(properties)
             call finalise(data)
             call finalise(at)
             RAISE_ERROR("cinoutput_read: frame "//int(do_frame)//" out of range 0 <= frame < "//int(this%n_frame), ierror)


! Lars
!             if (present(status)) then
!                call finalise(properties)
!                call finalise(data)
!                call finalise(at)
!                status = 1
!                return
!             else
!                call system_abort("cinoutput_read: frame "//int(do_frame)//" out of range 0 <= frame < "//int(this%n_frame))
!             end if
          end if
       end if

       do_zero = 0
       if (present(zero)) then
          if (zero) do_zero = 1
       end if

       if (this%got_index == 1) then
          tmp_do_frame = do_frame
! Lars
!          call cinoutput_query(this, tmp_do_frame, status=ierror)
          call cinoutput_query(this, tmp_do_frame, ierror=ierror)
       else
! Lars
!          call cinoutput_query(this, status=ierror)
          call cinoutput_query(this, ierror=ierror)
       end if
! Lars
!       if (present(status)) then
!          if (status /= 0) return
!       endif
       PASS_ERROR(ierror)

       call initialise(properties)
       do i=1,this%n_property
          call set_value(properties, c_array_to_f_string(this%property_name(:,i)), &
               (/ this%property_type(i), this%property_start(i)+1, this%property_start(i)+this%property_ncols(i) /))
       end do

       ! Make a blank data table, then initialise atoms object from it. We will
       ! get C routine to read directly into final Atoms structure to save copying.
       call allocate(data, this%n_int, this%n_real, this%n_str, this%n_logical, int(this%n_atom))
       call append(data, blank_rows=int(this%n_atom))
       this%lattice = 0.0_dp
       this%lattice(1,1) = 1.0_dp
       this%lattice(2,2) = 1.0_dp
       this%lattice(3,3) = 1.0_dp
       call initialise(at, int(this%n_atom), transpose(this%lattice), data, properties)

       int_ptr = C_NULL_PTR
       if (this%n_int /= 0) int_ptr = c_loc(at%data%int(1,1))

       real_ptr = C_NULL_PTR
       if (this%n_real /= 0) real_ptr = c_loc(at%data%real(1,1))

       str_ptr = C_NULL_PTR
       if (this%n_str /= 0) str_ptr = c_loc(at%data%str(1,1))

       log_ptr = C_NULL_PTR
       if (this%n_logical /= 0) log_ptr = c_loc(at%data%logical(1,1))

       if (cioread(this%c_at, do_frame, int_ptr, real_ptr, str_ptr, log_ptr, do_zero) == 0) then
          call finalise(properties)
          call finalise(data)
          call finalise(at)
          RAISE_ERROR("Error reading from file", ierror)
! Lars
!          if (present(status)) then
!             call finalise(properties)
!             call finalise(data)
!             call finalise(at)
!             status = 1
!             return
!          else
!             call system_abort("Error reading from file")
!          end if
       end if

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
             call system_abort('cinoutput_read: unsupported parameter i='//i//' key='//trim(namestr)//' type='//this%param_type(i))
          end select
       end do

       call atoms_repoint(at)
       call set_lattice(at, transpose(this%lattice), scale_positions=.false.)

       if (.not. has_property(at,"Z") .and. .not. has_property(at, "species")) then
          call print ("at%properties", ERROR)
          call print(at%properties, ERROR)
          call system_abort('cinoutput_read: atoms object read from file has neither Z nor species')
       else if (.not. has_property(at,"species") .and. has_property(at,"Z")) then
          call add_property(at, "species", repeat(" ",TABLE_STRING_LENGTH))
          call atoms_repoint(at)
          at%species = ElementName(at%Z)
       else if (.not. has_property(at,"Z") .and. has_property(at, "species")) then
          call add_property(at, "Z", 0)
          call atoms_repoint(at)
          do i=1,at%n
             at%Z(i) = atomic_number_from_symbol(at%species(i))
          end do
       end if

       call finalise(properties)
       call finalise(data)
       this%current_frame = this%current_frame + 1
    end if

    if (this%mpi%active) call bcast(this%mpi, at)

  end subroutine cinoutput_read


  subroutine cinoutput_write(this, at, properties, int_format, real_format, frame, &
       shuffle, deflate, deflate_level, status)
    type(CInOutput), intent(inout) :: this
    type(Atoms), target, intent(in) :: at
    character(*), intent(in), optional :: properties(:)
    character(*), intent(in), optional :: int_format, real_format
    integer, intent(in), optional :: frame
    logical, intent(in), optional :: shuffle, deflate
    integer, intent(in), optional :: deflate_level
    integer, optional, intent(out) :: status

    type(C_PTR) :: int_ptr, real_ptr, str_ptr, log_ptr
    logical :: dum
    character(len=VALUE_LEN) :: valuestr
    integer :: i, j, lookup(3), extras, n, tmp_int_a2(3,3)
    real(dp) :: tmp_real_a2(3,3)
    integer(C_INT) :: do_frame
    type(Dictionary) :: selected_properties
    character(len=KEY_LEN) :: do_int_format, do_real_format
    integer :: do_shuffle, do_deflate
    integer :: do_deflate_level, do_swap

    if (present(status)) status = 0

    if (.not. this%initialised) call system_abort("This CInOutput object is not initialised")
    if (this%action /= OUTPUT .and. this%action /= INOUT) call system_abort("Cannot write to action=INPUT CInOutput object")

    if (this%mpi%active .and. this%mpi%my_proc /= 0) return

    do_int_format = optional_default('%8d', int_format)
    do_real_format = optional_default('%16.8f', real_format)

    do_frame = optional_default(this%current_frame, frame)

    do_shuffle = transfer(optional_default(.true., shuffle),do_shuffle)
    do_deflate = transfer(optional_default(.true., deflate),do_shuffle)
    do_deflate_level = optional_default(6, deflate_level)
    do_swap = transfer(.not. present(properties), do_swap)
    
    call initialise(selected_properties)
    if (.not. present(properties)) then
       do i=1,at%properties%N
          call set_value(selected_properties, at%properties%keys(i), 1)
       end do
    else
       do i=1,size(properties)
          if (.not. has_key(at%properties, properties(i))) then
             call print('cinoutput_write: skipping unknown property '//trim(properties(i)))
             cycle
          end if
          call set_value(selected_properties, properties(i), 1)
       end do
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
       if (lower_case(trim(at%params%keys(i))) == lower_case('Lattice') .or. &
           lower_case(trim(at%params%keys(i))) == lower_case('Properties')) cycle
       call f_string_to_c_array(at%params%keys(i), this%param_name(:,n))
       this%param_filter(n) = 1
       select case(at%params%entries(i)%type)
       case(T_INTEGER)
          dum = get_value(at%params, at%params%keys(i), this%pint(n))
          this%param_size(n) = 1
          this%param_type(n) = T_INTEGER
       case(T_REAL)
          dum = get_value(at%params, at%params%keys(i), this%preal(n))
          this%param_size(n) = 1
          this%param_type(n) = T_REAL
       case(T_LOGICAL)
          dum = get_value(at%params, at%params%keys(i), this%plogical(n))
          this%param_size(n) = 1
          this%param_type(n) = T_LOGICAL
       case(T_INTEGER_A)
          dum = get_value(at%params, at%params%keys(i), this%pint_a(:,n))
          this%param_size(n) = 3
          this%param_type(n) = T_INTEGER_A
       case (T_REAL_A)
          dum = get_value(at%params, at%params%keys(i), this%preal_a(:,n))
          this%param_size(n) = 3
          this%param_type(n) = T_REAL_A
       case(T_LOGICAL_A)
          dum = get_value(at%params, at%params%keys(i), this%plogical_a(:,n))
          this%param_size(n) = 3
          this%param_type(n) = T_LOGICAL_A
       case(T_CHAR)
          dum = get_value(at%params, at%params%keys(i), valuestr)
          call f_string_to_c_array(valuestr, this%param_value(:,n))
          this%param_size(n) = 1
          this%param_type(n) = T_CHAR
       case(T_INTEGER_A2)
          if (.not. get_value(at%params, at%params%keys(i), tmp_int_a2)) &
               call system_abort('Bad atoms parameter '//trim(at%params%keys(i))//' should be 3x3 int array')
          this%pint_a2(:,n) = reshape(tmp_int_a2, (/9/))
          this%param_size(n) = 9
          this%param_type(n) = T_INTEGER_A2
       case (T_REAL_A2)
          if (.not. get_value(at%params, at%params%keys(i), tmp_real_a2)) &
               call system_abort('Bad atoms parameter '//trim(at%params%keys(i))//' should be 3x3 real array')
          this%preal_a2(:,n) = reshape(tmp_real_a2, (/9/))
          this%param_size(n) = 9
          this%param_type(n) = T_REAL_A2
       case default
          call system_abort('cinoutput_write: unsupported parameter i='//i//' '//trim(at%params%keys(i))//' type='//at%params%entries(i)%type)
       end select
       n = n + 1
    end do

    this%n_property = at%properties%N
    call c_f_pointer(this%c_property_name, this%property_name, (/KEY_LEN,this%n_property/))
    call c_f_pointer(this%c_property_type, this%property_type, (/this%n_property/))
    call c_f_pointer(this%c_property_ncols, this%property_ncols, (/this%n_property/))
    call c_f_pointer(this%c_property_start, this%property_start, (/this%n_property/))
    call c_f_pointer(this%c_property_filter, this%property_filter, (/this%n_property/))

    ! First add the selected properties in the order they were specified...
    do i=1,selected_properties%n
       j = lookup_entry_i(at%properties, selected_properties%keys(i))
       call f_string_to_c_array(at%properties%keys(j), this%property_name(:,i))
       dum = get_value(at%properties, at%properties%keys(j), lookup)
       this%property_type(i) = lookup(1)
       this%property_start(i) = lookup(2)-1
       this%property_ncols(i) = lookup(3)-lookup(2)+1
       this%property_filter(i) = 1
    end do
    ! ... then the non-selected ones
    j = selected_properties%n+1
    do i=1,at%properties%n
       if (has_key(selected_properties, at%properties%keys(i))) cycle
       call f_string_to_c_array(at%properties%keys(i), this%property_name(:,j))
       dum = get_value(at%properties, at%properties%keys(i), lookup)
       this%property_type(j) = lookup(1)
       this%property_start(j) = lookup(2)-1
       this%property_ncols(j) = lookup(3)-lookup(2)+1
       this%property_filter(j) = 0
       j = j +1
    end do

    int_ptr = C_NULL_PTR
    this%n_int = at%data%intsize
    if (this%n_int /= 0) int_ptr = c_loc(at%data%int(1,1))

    real_ptr = C_NULL_PTR
    this%n_real = at%data%realsize
    if (this%n_real /= 0) real_ptr = c_loc(at%data%real(1,1))
    
    str_ptr = C_NULL_PTR
    this%n_str = at%data%strsize
    if (this%n_str /= 0) str_ptr = c_loc(at%data%str(1,1))

    log_ptr = C_NULL_PTR
    this%n_logical = at%data%logicalsize
    if (this%n_logical /= 0) log_ptr = c_loc(at%data%logical(1,1))

    if (ciowrite(this%c_at, int_ptr, real_ptr, str_ptr, log_ptr, &
         trim(do_int_format)//C_NULL_CHAR, trim(do_real_format)//C_NULL_CHAR, do_frame, &
         do_shuffle, do_deflate, do_deflate_level, do_swap) == 0) then
       if (present(status)) then
          status = 1
       else
          call system_abort("Error writing file")
       end if
    end if

    call finalise(selected_properties)
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
    character(1), dimension(len(fstring)), intent(out) :: carray

    integer :: i

    do i=1,len_trim(fstring)
       carray(i) = fstring(i:i)
    end do
    carray(len_trim(fstring)+1) = C_NULL_CHAR

  end subroutine f_string_to_c_array


  subroutine atoms_read(this, filename, frame, zero, ierror)
    !% Read Atoms object from XYZ or NetCDF file.
    type(Atoms), intent(inout) :: this
    character(len=*), intent(in) :: filename
    integer, optional, intent(in) :: frame
    logical, optional, intent(in) :: zero
    integer, optional, intent(inout) :: ierror

    type(CInOutput) :: cio

    call initialise(cio, filename, INPUT, ierror=ierror)
    PASS_ERROR_WITH_INFO('While reading "' // filename // '".', ierror)
    call read(cio, this, frame, zero, ierror=ierror)
    PASS_ERROR_WITH_INFO('While reading "' // filename // '".', ierror)
    call finalise(cio)

  end subroutine atoms_read

  subroutine atoms_read_cinoutput(this, cio, frame, zero, ierror)
    type(Atoms), target, intent(inout) :: this
    type(CInOutput), intent(inout) :: cio
    integer, optional, intent(in) :: frame
    logical, optional, intent(in) :: zero
    integer, optional, intent(inout) :: ierror

    call cinoutput_read(cio, this, frame, zero, ierror=ierror)
    PASS_ERROR(ierror)

  end subroutine atoms_read_cinoutput

  subroutine atoms_write(this, filename, append, properties, status)
    !% Write Atoms object to XYZ or NetCDF file. Use filename "stdout" to write to terminal.
    type(Atoms), intent(in) :: this
    character(len=*), intent(in) :: filename
    logical, optional, intent(in) :: append
    character(*), intent(in), optional :: properties(:)    
    integer, optional, intent(out) :: status

    type(CInOutput) :: cio

    call initialise(cio, filename, OUTPUT, append)
    call write(cio, this, properties, status=status)
    call finalise(cio)
    
  end subroutine atoms_write

  subroutine atoms_write_cinoutput(this, cio, properties, int_format, real_format, frame, shuffle, deflate, deflate_level, status)
    type(Atoms), target, intent(inout) :: this
    type(CInOutput), intent(inout) :: cio
    character(*), intent(in), optional :: properties(:)    
    character(*), intent(in), optional :: int_format, real_format
    integer, intent(in), optional :: frame
    logical, intent(in), optional :: shuffle, deflate
    integer, intent(in), optional :: deflate_level
    integer, optional, intent(out) :: status

    call cinoutput_write(cio, this, properties, int_format, real_format, frame, shuffle, deflate, deflate_level, status)

  end subroutine atoms_write_cinoutput



end module CInOutput_module
