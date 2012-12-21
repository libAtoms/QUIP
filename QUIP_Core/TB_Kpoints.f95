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

!X
!X TB_KPoints Module 
!X
!% This is a module for doing (weighted and unweighted) k-points sums, 
!% sorting, computing phase factors.  
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

module TB_KPoints_module

use system_module, only : dp, inoutput, initialise, INPUT, is_file_readable, NUMERICAL_ZERO, system_abort
use units_module
use extendable_str_module
use linearalgebra_module
use mpi_context_module
use quip_Common_module

implicit none

private

character(len=20) :: kpt_default_file = "k_pt_mesh.data.xml"

!% Structure for storing k-points and doing k-points sums (in parallel)
public :: KPoints
type KPoints
  integer :: N = 0, g_N = 0
  real(dp), allocatable  :: k_pts(:,:), weights(:), g_k_pts(:,:), g_weights(:)
  logical :: non_gamma = .false.
  logical :: no_sum_over_my_kpt = .false.
  type(MPI_context) :: mpi_global, mpi_my_kpt, mpi_across_kpts
end type

type(extendable_str), save :: cur_data
logical :: parse_in_kp
integer :: parse_cur_kp
type(KPoints), pointer :: parse_kp

!% Initialise by reading parameters from file/inoutput
public :: Initialise
public :: KPoints_Initialise_filename
interface Initialise
  module procedure KPoints_Initialise_inoutput, KPoints_Initialise_str, KPoints_Initialise_kp
  module procedure KPoints_Initialise_mesh, KPoints_Initialise_density
end interface Initialise

!% Prepare splitting according to mpi. 
public :: init_mpi
interface init_mpi
  module procedure KPoints_init_mpi
end interface init_mpi

!% Interface to deallocate everything.
public :: Finalise
interface Finalise
  module procedure KPoints_Finalise
end interface Finalise

!% Free mpi_contexts generated in init_mpi.
public :: end_mpi
interface end_mpi
  module procedure KPoints_end_mpi
end interface end_mpi

!% Interface to print everything.
public :: Print
interface Print
  module procedure KPoints_Print
end interface Print

!% Calculate a phase factor for a given k-point and cell shift.
public :: calc_phase
interface calc_phase
  module procedure KPoints_calc_phase
end interface

!% Read params in XML format.
private :: KPoints_read_points_xml
interface read_points_xml
  module procedure KPoints_read_points_xml
end interface read_points_xml

!% Minimum over k-points of a quantity that's duplicated among processors of each k-point.
public :: min
interface min
  module procedure KPoints_kmin_real, KPoints_kmin_real1 
end interface min

!% Maximum over k-points of a quantity that's duplicated among processors of each k-point.
public :: max
interface max
  module procedure KPoints_kmax_real, KPoints_kmax_real1
end interface max

!% Weighted sum over k-points of a quantity that's duplicated among processors of each k-point.
public :: ksum_dup
interface ksum_dup
  module procedure KPoints_ksum_dup_r1, KPoints_ksum_dup_c1
end interface ksum_dup

!% Local weighted sum over k-points of a quantity.
public :: local_ksum
interface local_ksum
  module procedure KPoints_local_ksum_real1, KPoints_local_ksum_real2
  module procedure KPoints_local_ksum_complex1, KPoints_local_ksum_complex2, KPoints_local_ksum_complex4
end interface local_ksum

!% UNweighted sum over k-points of an array that's distributed among processors of each k-point.
public :: ksum_distrib
interface ksum_distrib
  module procedure KPoints_ksum_distrib_real
end interface ksum_distrib

!% in-place UNweighted sum over k-points of an array that's distributed among processors of each k-point.
public :: ksum_distrib_inplace
interface ksum_distrib_inplace
  module procedure KPoints_ksum_distrib_inplace_real1, KPoints_ksum_distrib_inplace_real2
  module procedure KPoints_ksum_distrib_inplace_complex1
end interface ksum_distrib_inplace

public :: collect
interface collect
  module procedure Kpoints_collect_real2
end interface collect

contains

subroutine KPoints_Initialise_filename(this, filename, mpi_obj)
  type(KPoints), intent(inout) :: this
  character(len=*), intent(in) :: filename
  type(MPI_context), intent(in), optional :: mpi_obj

  type(inoutput) io

  call Initialise(io,filename,INPUT)
  call Initialise(this,io,mpi_obj)
  call Finalise(io)

end subroutine

subroutine KPoints_Initialise_inoutput(this, from_io, mpi_obj)
  type(KPoints), intent(inout) :: this
  type(inoutput), intent(inout), target, optional :: from_io
  type(MPI_context), intent(in), optional :: mpi_obj

  logical file_present
  type(inoutput), pointer :: in
  type(extendable_str) :: ss

  file_present = .true.

  if (present(from_io)) then
    in => from_io
  else
    if (is_file_readable(trim(kpt_default_file))) then
      allocate(in)
      call Initialise(in, trim(kpt_default_file), INPUT)
    else
      file_present = .false.
    endif
  endif

  call Initialise(ss)
  if (file_present) then
    if (present(mpi_obj)) then
      call read(ss, in%unit, convert_to_string=.true., mpi_comm=mpi_obj%communicator)
    else
      call read(ss, in%unit, convert_to_string=.true.)
    endif
  endif

  call Initialise(this, string(ss), mpi_obj)

  call Finalise(ss)

end subroutine

subroutine KPoints_Initialise_str(this, from_str, mpi_obj)
  type(KPoints), intent(inout) :: this
  character(len=*), intent(in) :: from_str
  type(MPI_context), intent(in), optional :: mpi_obj

  call Finalise(this)

  call KPoints_read_points_xml(this, from_str)

  call finish_initialise(this, mpi_obj)

end subroutine KPoints_Initialise_str

subroutine KPoints_Initialise_kp(this, from_kpoints, mpi_obj)
  type(KPoints), intent(inout) :: this
  type(KPoints), intent(in) :: from_kpoints
  type(MPI_context), intent(in), optional :: mpi_obj

  call Finalise(this)

  this%non_gamma = from_kpoints%non_gamma
  this%N = from_kpoints%g_N
  if (this%N > 0) then
    allocate(this%k_pts(3,this%N))
    allocate(this%weights(this%N))
    this%k_pts = from_kpoints%g_k_pts
    this%weights = from_kpoints%g_weights
  endif

  call finish_initialise(this, mpi_obj)

end subroutine KPoints_Initialise_kp

function reciprocal_lattice(lattice)
  real(dp), intent(in) :: lattice(3,3)
  real(dp) :: reciprocal_lattice(3,3)

  call matrix3x3_inverse(lattice,reciprocal_lattice)
  reciprocal_lattice = -2.0_dp*PI*transpose(reciprocal_lattice)

end function reciprocal_lattice

subroutine KPoints_Initialise_density(this, lattice, k_space_density, monkhorst_pack, mpi_obj)
  type(KPoints), intent(inout) :: this
  real(dp), intent(in) :: lattice(3,3)
  real(dp), intent(in) :: k_space_density
  logical, intent(in), optional :: monkhorst_pack
  type(MPI_context), intent(in), optional :: mpi_obj

  integer :: n(3)
  real(dp) :: recip_lattice(3,3)

  if (k_space_density > 0.0_dp) then
    recip_lattice = reciprocal_lattice(lattice)

    call matrix3x3_inverse(lattice,recip_lattice)
    recip_lattice = -2.0_dp*PI*transpose(recip_lattice)

    n(1) = norm(recip_lattice(:,1))*k_space_density+1.0_dp-NUMERICAL_ZERO
    n(2) = norm(recip_lattice(:,2))*k_space_density+1.0_dp-NUMERICAL_ZERO
    n(3) = norm(recip_lattice(:,3))*k_space_density+1.0_dp-NUMERICAL_ZERO
    where (n == 0)
      n = 1
    end where
    where (n > 1)
      n = n + mod(n,2)
    end where
  else
    n = 1
  endif
  call initialise(this, n, monkhorst_pack, mpi_obj)

end subroutine KPoints_initialise_density

subroutine KPoints_Initialise_mesh(this, n, monkhorst_pack, mpi_obj)
  type(KPoints), intent(inout) :: this
  integer, intent(in) :: n(3)
  logical, intent(in), optional :: monkhorst_pack
  type(MPI_context), intent(in), optional :: mpi_obj

  integer :: i, i1, i2, i3, ik
  integer :: min(3), max(3), step(3)
  real(dp) :: denom(3)
  logical :: shift(3)
  real(dp), allocatable :: mesh_wt(:,:,:)

  call Finalise(this)

  if (minval(n) < 1) call system_abort("Called KPoints_initialise_mesh with garbage mesh size " // n)

  if (maxval(n) > 1) then
    shift = .false.
    if (present(monkhorst_pack)) then
      if (monkhorst_pack .and. mod(n(1),2) == 0) then
	shift(1) = .true.
      end if
      if (monkhorst_pack .and. mod(n(2),2) == 0) then
	shift(2) = .true.
      end if
      if (monkhorst_pack .and. mod(n(3),2) == 0) then
	shift(3) = .true.
      end if
    endif

    this%N = 0
    do i=1, 3
      if (shift(i)) then
        if (mod(n(i),2) == 0) then
          min(i) = -(n(i)-1)
          max(i) = n(i)-1
          step(i) = 2
          denom(i) = 2*n(i)
        else
          call system_abort("kpoints_initialise_mesh tried to Monkhorst-Pack shift odd dimension " // i)
        endif
      else
        if (mod(n(i),2) == 0) then
          min(i) = -n(i)/2+1
          max(i) = n(i)/2
          step(i) = 1
          denom(i) = n(i)
        else
          if (n(i) == 1) then
            min(i) = 0
            max(i) = 0
            step(i) = 1
            denom(i) = 1
          else
            min(i) = -(n(i)-1)/2
            max(i) = (n(i)-1)/2
            step(i) = 1
            denom(i) = n(i)
          end if
        end if
      endif
    end do

    allocate(mesh_wt(min(1):max(1),min(2):max(2),min(3):max(3)))
    mesh_wt = 1.0_dp
    do i1=min(1), max(1), step(1)
    do i2=min(2), max(2), step(2)
    do i3=min(3), max(3), step(3)
      if (-i1 >= min(1) .and. -i1 <= max(1) .and. &
          -i2 >= min(2) .and. -i2 <= max(2) .and. &
          -i3 >= min(3) .and. -i3 <= max(3)) then ! reflected k-point is in range
        if ((i1 /= 0 .or. i2 /= 0 .or. i3 /= 0) .and. & 
            (mesh_wt(i1,i2,i3) > 0.0_dp .and. mesh_wt(-i1,-i2,-i3) > 0.0_dp)) then ! not gamma and both points have weight
          mesh_wt(-i1,-i2,-i3) = mesh_wt(-i1,-i2,-i3) + mesh_wt(i1,i2,i3)
          mesh_wt(i1,i2,i3) = 0.0_dp
        endif
      endif
    end do
    end do
    end do

    this%N = 0
    do i1=min(1), max(1), step(1)
    do i2=min(2), max(2), step(2)
    do i3=min(3), max(3), step(3)
      if (mesh_wt(i1,i2,i3) > 0.0_dp) this%N = this%N + 1
    end do
    end do
    end do
    allocate(this%k_pts(3,this%N))
    allocate(this%weights(this%N))

    ik = 1
    do i1=min(1), max(1), step(1)
    do i2=min(2), max(2), step(2)
    do i3=min(3), max(3), step(3)
      if (mesh_wt(i1,i2,i3) > 0.0_dp) then
        this%k_pts(1,ik) = real(i1,dp)/denom(1)
        this%k_pts(2,ik) = real(i2,dp)/denom(2)
        this%k_pts(3,ik) = real(i3,dp)/denom(3)
        this%weights(ik) = mesh_wt(i1,i2,i3)
        ik = ik + 1
      endif
    end do
    end do
    end do

    deallocate(mesh_wt)
  else ! maxval(n) <= 1
    this%non_gamma = .false.
    this%N = 1
    allocate(this%k_pts(3,this%N))
    allocate(this%weights(this%N))
    this%k_pts = 0.0_dp
    this%weights = 1.0_dp
  endif

  call finish_initialise(this, mpi_obj)

end subroutine KPoints_initialise_mesh

subroutine finish_initialise(this, mpi_obj)
  type(KPoints), intent(inout) :: this
  type(MPI_Context), intent(in), optional :: mpi_obj

  if (this%N > 0) then
    this%weights = this%weights / sum(this%weights)
    this%non_gamma = (maxval(abs(this%k_pts(1:3,1:this%N))) > 0.0_dp)
  else
    this%N = 1
    allocate(this%k_pts(3,1))
    allocate(this%weights(1))
    this%k_pts(1:3,1) = 0.0_dp
    this%weights(1) = 1.0_dp
    this%non_gamma = .false.
  endif

  this%g_N = this%N
  allocate(this%g_k_pts(3,this%g_N))
  allocate(this%g_weights(this%g_N))
  this%g_k_pts = this%k_pts
  this%g_weights = this%weights

  if (present(mpi_obj)) then
    call KPoints_init_mpi(this, mpi_obj)
  endif

end subroutine finish_initialise

subroutine KPoints_init_mpi(this, mpi_obj)
  type(KPoints), intent(inout) :: this
  type(MPI_context), intent(in) :: mpi_obj

  integer i
  integer t, local_N
  integer, allocatable :: local_index(:)
  integer my_split_ident
  integer n_procs_per_matrix
  real(dp), allocatable :: t_k_pts(:,:), t_weights(:)

  this%mpi_global = mpi_obj

  if (this%g_N >= this%mpi_global%n_procs) then
    n_procs_per_matrix = 1
  else
    n_procs_per_matrix = this%mpi_global%n_procs/this%g_N
  endif

  allocate(local_index(this%g_N))
  local_N = 0
  do i=0, this%g_N-1
    t = this%mpi_global%my_proc/n_procs_per_matrix
    if (mod(i,this%mpi_global%n_procs) == t) then
      local_N = local_N + 1
      local_index(local_N) = i+1
    endif
  end do

  if (this%mpi_global%active) then

    my_split_ident = this%mpi_global%my_proc/n_procs_per_matrix
    call split_context(this%mpi_global, my_split_ident, this%mpi_my_kpt)

    my_split_ident = this%mpi_my_kpt%my_proc
    call split_context(this%mpi_global, my_split_ident, this%mpi_across_kpts)
  endif

  allocate(t_k_pts(3,local_N))
  allocate(t_weights(local_N))
  t_k_pts(1:3,1:local_N) = this%g_k_pts(1:3,local_index(1:local_N))
  t_weights(1:local_N) = this%g_weights(local_index(1:local_N))

  deallocate(this%k_pts)
  deallocate(this%weights)
  this%N = local_N
  allocate(this%k_pts(3,this%N))
  allocate(this%weights(this%N))
  this%k_pts = t_k_pts
  this%weights = t_weights
  deallocate(t_k_pts)
  deallocate(t_weights)

end subroutine KPoints_init_mpi

subroutine KP_characters_handler(in)
  character(len=*), intent(in) :: in

  if (parse_in_kp) call concat(cur_data, in, keep_lf=.false.)
end subroutine

subroutine KP_startElement_handler(URI, localname, name, attributes)
  character(len=*), intent(in)   :: URI  
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name 
  type(dictionary_t), intent(in) :: attributes

  integer status
  character(len=1024) :: val

  if (name == "KPoints") then
    parse_in_kp = .true.
    parse_cur_kp=1
    call QUIP_FoX_get_value(attributes, "N", val, status)
    if (status /= 0) &
      call system_abort("parse kpoints xml: KPoints doesn't have N attribute")
    read(val, *) parse_kp%N
    if (parse_kp%N > 0) then
      allocate(parse_kp%k_pts(3,parse_kp%N))
      allocate(parse_kp%weights(parse_kp%N))
    endif
  elseif (parse_in_kp .and. name == "point") then
    if (parse_cur_kp > parse_kp%N) &
      call system_abort("parse kpoints xml found too many point elements " // parse_cur_kp // " > " // parse_kp%N)
    call QUIP_FoX_get_value(attributes, "weight", val, status)
    if (status /= 0) &
      call system_abort("parse kpoints xml: point " // parse_cur_kp // " doesn't have weight attribute")
    read(val, *) parse_kp%weights(parse_cur_kp)
  endif

  call zero(cur_data)

end subroutine

subroutine KP_endElement_handler(URI, localname, name)
  character(len=*), intent(in)   :: URI  
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name 

  character(len=1024) :: str

  if (name == "KPoints") then
    parse_in_kp = .false.
  else if (parse_in_kp .and. name == "point") then
    if (parse_cur_kp > parse_kp%N) &
      call system_abort("parse kpoints xml: too many points specified, at least " // parse_cur_kp)
    str = string(cur_data)
    read(str, *) parse_kp%k_pts(1:3,parse_cur_kp)
    parse_cur_kp = parse_cur_kp + 1
    call zero(cur_data)
  endif

end subroutine

subroutine KPoints_read_points_xml(this, str)
  type(KPoints), intent(inout), target :: this
  character(len=*), intent(in) :: str

  type(xml_t) :: fxml

  if (len(trim(str)) <= 0) then
    this%N = 0
    return
  endif

  parse_cur_kp = 1
  parse_kp => this

  call open_xml_string(fxml, str)

  call parse(fxml, &
    characters_handler=KP_characters_handler, &
    startElement_handler=KP_startElement_handler, &
    endElement_handler=KP_endElement_handler)

  call close_xml_t(fxml)

end subroutine

subroutine KPoints_Finalise(this)
  type(KPoints), intent(inout) :: this

  call end_mpi(this)

  if (allocated(this%k_pts)) deallocate(this%k_pts)
  if (allocated(this%weights)) deallocate(this%weights)
  if (allocated(this%g_k_pts)) deallocate(this%g_k_pts)
  if (allocated(this%g_weights)) deallocate(this%g_weights)

  this%N = 0
  this%g_N = 0
  this%non_gamma = .false.

end subroutine KPoints_Finalise

subroutine KPoints_end_mpi(this)
  type(KPoints), intent(inout) :: this

  call free_context(this%mpi_my_kpt)
  call free_context(this%mpi_across_kpts)
end subroutine KPoints_end_mpi

subroutine KPoints_Print(this,file)
  type(KPoints),    intent(in)           :: this
  type(Inoutput), intent(inout),optional:: file

  integer::i

  call Print('KPoints : ', file=file)

  call Print ('KPoints : non_gamma ' // this%non_gamma, file=file)
  call Print ('KPoints : g_N ' // this%g_N, file=file)
  do i=1, this%g_N
    call Print ('KPoints : g_k_pt ' // i // " " // this%g_k_pts(1:3,i) // " " // this%g_weights(i), file=file)
  end do

end subroutine KPoints_Print

function KPoints_calc_phase(this, ik, shift)
  type(KPoints), intent(in) :: this
  integer, intent(in) :: ik
  integer, intent(in) :: shift(3)
  complex(dp) :: KPoints_calc_phase

  real(dp) phase

  phase = 2.0_dp*PI*sum(shift(1:3)*this%k_pts(1:3,ik))
  KPoints_calc_phase = cmplx(cos(phase), sin(phase), KIND(phase))

end function KPoints_calc_phase


function KPoints_kmin_real(this, v)
  type(KPoints), intent(in) :: this
  real(dp), intent(in) :: v
  real(dp) :: KPoints_kmin_real

  if (this%mpi_my_kpt%my_proc == 0) then
    KPoints_kmin_real = min(this%mpi_across_kpts, v)
  endif

  call bcast(this%mpi_my_kpt, KPoints_kmin_real)
end function KPoints_kmin_real

function KPoints_kmin_real1(this, v)
  type(KPoints), intent(in) :: this
  real(dp), intent(in) :: v(:)
  real(dp) :: KPoints_kmin_real1

  if (this%mpi_my_kpt%my_proc == 0) then
    KPoints_kmin_real1 = min(this%mpi_across_kpts, minval(v))
  endif

  call bcast(this%mpi_my_kpt, KPoints_kmin_real1)
end function  KPoints_kmin_real1

function KPoints_kmax_real(this, v)
  type(KPoints), intent(in) :: this
  real(dp), intent(in) :: v
  real(dp) :: KPoints_kmax_real

  if (this%mpi_my_kpt%my_proc == 0) then
    KPoints_kmax_real = max(this%mpi_across_kpts, v)
  endif

  call bcast(this%mpi_my_kpt, KPoints_kmax_real)
end function KPoints_kmax_real

function KPoints_kmax_real1(this, v)
  type(KPoints), intent(in) :: this
  real(dp), intent(in) :: v(:)
  real(dp) :: KPoints_kmax_real1

  if (this%mpi_my_kpt%my_proc == 0) then
    KPoints_kmax_real1 = max(this%mpi_across_kpts, maxval(v))
  endif

  call bcast(this%mpi_my_kpt, KPoints_kmax_real1)
end function KPoints_kmax_real1

function KPoints_ksum_dup_r1(this, v)
  type(KPoints), intent(in) :: this
  real(dp), intent(in) :: v(:)
  real(dp) :: KPoints_ksum_dup_r1

  if (this%mpi_my_kpt%my_proc == 0) then
    if (size(v) /= size(this%weights)) call System_abort("size mismatch in KPoints_ksum_dup_r1")
    KPoints_ksum_dup_r1 = sum(this%mpi_across_kpts, dot_product(v,this%weights))
  endif

  call bcast(this%mpi_my_kpt, KPoints_ksum_dup_r1)
end function KPoints_ksum_dup_r1

function KPoints_ksum_dup_c1(this, v)
  type(KPoints), intent(in) :: this
  complex(dp), intent(in) :: v(:)
  complex(dp) :: KPoints_ksum_dup_c1

  if (this%mpi_my_kpt%my_proc == 0) then
    if (size(v) /= size(this%weights)) call System_abort("size mismatch in KPoints_ksum_dup_c1")
    KPoints_ksum_dup_c1 = sum(this%mpi_across_kpts, dot_product(v,this%weights))
  endif

  call bcast(this%mpi_my_kpt, KPoints_ksum_dup_c1)
end function KPoints_ksum_dup_c1

function KPoints_local_ksum_real1(this, v)
  type(KPoints), intent(in) :: this
  real(dp), intent(in) :: v(:)
  real(dp) :: KPoints_local_ksum_real1

  KPoints_local_ksum_real1 = dot_product(v,this%weights)
end function KPoints_local_ksum_real1

function KPoints_local_ksum_real2(this, v)
  type(KPoints), intent(in) :: this
  real(dp), intent(in) :: v(:,:)
  real(dp) :: KPoints_local_ksum_real2(size(v,1))

  integer ik

  KPoints_local_ksum_real2(:) = 0.0_dp
  do ik=1, this%N
    KPoints_local_ksum_real2(:) = KPoints_local_ksum_real2(:) + v(:,ik)*this%weights(ik)
  end do
end function KPoints_local_ksum_real2

function KPoints_local_ksum_complex1(this, v)
  type(KPoints), intent(in) :: this
  complex(dp), intent(in) :: v(:)
  complex(dp) :: KPoints_local_ksum_complex1

  KPoints_local_ksum_complex1 = dot_product(v,this%weights)
end function KPoints_local_ksum_complex1

function KPoints_local_ksum_complex2(this, v)
  type(KPoints), intent(in) :: this
  complex(dp), intent(in) :: v(:,:)
  complex(dp) :: KPoints_local_ksum_complex2(size(v,1))

  integer ik

  KPoints_local_ksum_complex2(:) = 0.0_dp
  do ik=1, this%N
    KPoints_local_ksum_complex2(:) = KPoints_local_ksum_complex2(:) + v(:,ik)*this%weights(ik)
  end do
end function KPoints_local_ksum_complex2

function KPoints_local_ksum_complex4(this, v)
  type(KPoints), intent(in) :: this
  complex(dp), intent(in) :: v(:,:,:,:)
  complex(dp) :: KPoints_local_ksum_complex4(size(v,1),size(v,2),size(v,3))

  integer ik

  KPoints_local_ksum_complex4(:,:,:) = 0.0_dp
  do ik=1, this%N
    KPoints_local_ksum_complex4(:,:,:) = KPoints_local_ksum_complex4(:,:,:) + v(:,:,:,ik)*this%weights(ik)
  end do
end function KPoints_local_ksum_complex4

function KPoints_ksum_distrib_real(this, v)
  type(KPoints), intent(in) :: this
  real(dp), intent(in) :: v
  real(dp) :: KPoints_ksum_distrib_real

  if (this%no_sum_over_my_kpt) then
    KPoints_ksum_distrib_real = sum(this%mpi_across_kpts, v)
  else
    KPoints_ksum_distrib_real = sum(this%mpi_across_kpts, sum(this%mpi_my_kpt, v))
  endif
end function KPoints_ksum_distrib_real

subroutine KPoints_ksum_distrib_inplace_real2(this, v)
  type(KPoints), intent(in) :: this
  real(dp), intent(inout) :: v(:,:)

  if (.not. this%no_sum_over_my_kpt) then
    call sum_in_place(this%mpi_my_kpt, v)
  endif
  call sum_in_place(this%mpi_across_kpts, v)
end subroutine KPoints_ksum_distrib_inplace_real2

subroutine KPoints_ksum_distrib_inplace_real1(this, v)
  type(KPoints), intent(in) :: this
  real(dp), intent(inout) :: v(:)

  if (.not. this%no_sum_over_my_kpt) then
    call sum_in_place(this%mpi_my_kpt, v)
  endif
  call sum_in_place(this%mpi_across_kpts, v)
end subroutine KPoints_ksum_distrib_inplace_real1

subroutine KPoints_ksum_distrib_inplace_complex1(this, v)
  type(KPoints), intent(in) :: this
  complex(dp), intent(inout) :: v(:)

  if (.not. this%no_sum_over_my_kpt) then
    call sum_in_place(this%mpi_my_kpt, v)
  endif
  call sum_in_place(this%mpi_across_kpts, v)
end subroutine KPoints_ksum_distrib_inplace_complex1

subroutine KPoints_collect_real2(this, v_in, v_out)
  type(KPoints), intent(in) :: this
  real(dp), intent(in) :: v_in(:,:)
  real(dp), intent(out) :: v_out(:,:)

  if (size(v_in,1) /= size(v_out,1)) then
    call system_abort("Called KPoints_collect_real2 with vector size mismatch v_in " // &
      size(v_in,1) // " v_out " // size(v_out,1))
  endif

  if (size(v_in,2) /= this%N)  then
    call system_abort("Called KPoints_collect_real2 with local N_k mismatch v_in " // &
      size(v_in,2) // " this%N " // this%N)
  end if

  if (size(v_out,2) /= this%g_N)  then
    call system_abort("Called KPoints_collect_real2 with local N_k mismatch v_in " // &
      size(v_in,2) // " this%N " // this%N)
  end if

  call collect(this%mpi_across_kpts, v_in, v_out)
end subroutine KPoints_collect_real2

end module TB_KPoints_module

