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
!X IPModel_MBX
!X
!% MBX interatomic potential available at https://github.com/paesanilab/MBX/, (branch: https://github.com/chemphys/MBX.git ; git checkout master-dev) 
!%
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#include "error.inc"

module IPModel_MBX_module

use error_module
use system_module, only : dp, inoutput, print, verbosity_push_decrement, verbosity_pop, operator(//), split_string,string_to_int
use dictionary_module
use paramreader_module
use linearalgebra_module
use atoms_types_module
use atoms_module

use mpi_context_module
use QUIP_Common_module
use units_module, only : KCAL_MOL

implicit none

! MBX functions
!external get_energy
!external get_energy_g
!external finalize_system

private

include 'IPModel_interface.h'

public :: IPModel_MBX
type IPModel_MBX
  integer :: n_types = 0
  integer, allocatable :: atomic_num(:), type_of_atomic_num(:)
  integer,allocatable,dimension(:) :: nats
  character(len=5),allocatable,dimension(:) :: at_name
  character(len=5),allocatable,dimension(:) :: monomers
  integer :: nmon
  character(len=STRING_LENGTH) :: json_file ! 20 char.s was short for me...


  real(dp) :: cutoff = 0.0_dp

  character(len=STRING_LENGTH) :: label

end type IPModel_MBX

logical, private :: parse_in_ip, parse_matched_label
type(IPModel_MBX), private, pointer :: parse_ip

interface Initialise
  module procedure IPModel_MBX_Initialise_str
end interface Initialise

interface Finalise
  module procedure IPModel_MBX_Finalise
end interface Finalise

interface Print
  module procedure IPModel_MBX_Print
end interface Print

interface Calc
  module procedure IPModel_MBX_Calc
end interface Calc

contains



subroutine IPModel_MBX_Initialise_str(this, args_str, param_str, error)
  type(IPModel_MBX), intent(inout) :: this
  character(len=*), intent(in) :: args_str, param_str
  type(Dictionary) :: params
  integer, optional, intent(out):: error
  character(len=STRING_LENGTH) :: nats_string, at_name_string,monomers_string,json_file_string,n_monomers_types_string,tmp_monomer_string
  !character(len=*) :: tmp_monomer_string
  character(len=STRING_LENGTH), allocatable :: nats_fields(:), at_name_fields(:), monomers_fields(:),n_monomers_types_fields(:)
  !character(len=STRING_LENGTH), dimension(1000) :: n_monomers_types_fields(:)

  !character(len=10),allocatable,dimension(:) :: n_monomers_types

  integer:: i,j,k,mon_ind,at_ind,tmp_num_monomer,n_group1,n_group2,n_group3,n_group4,sum_nats
  real(dp),allocatable  :: coord(:)


  external initialize_system
  external finalize_system


  INIT_ERROR(error)
  call Finalise(this)
  call finalize_system()

  call initialise(params)
  !write(*,*)  "DEBUGGING: initialised params reader"
  !call param_register(params, 'nats', PARAM_MANDATORY, nats_string, help_string="Number of atoms in the monomers, format {i1 i2 i3 ...}")
  !write(*,*)  "DEBUGGING: param_register nats done"
  !call param_register(params, 'monomers', PARAM_MANDATORY, monomers_string, help_string="Monomer names, format {ch4 h2o ...}")
  !write(*,*)  "DEBUGGING: param_register monomers done"
  !call param_register(params, 'at_name', PARAM_MANDATORY,at_name_string, help_string="Name of the atoms, format {C H H H H O H H ...}")
  !write(*,*)  "DEBUGGING: param_register at_name done"
  call param_register(params, 'n_monomers_types', PARAM_MANDATORY, n_monomers_types_string, &
          help_string="Numbers and types of monomers, format {2_CHHHH_3_OHH_1_CHHHH ...}")
  !write(*,*)  "DEBUGGING: param_register monomers done"
  call param_register(params, 'nmon', PARAM_MANDATORY, this%nmon, help_string="Number of monomers")
  !write(*,*)  "DEBUGGING: param_register nmon done"
  call param_register(params, 'json_file', PARAM_MANDATORY, json_file_string, help_string="Name of json file (to be already in the directory)")
  !write(*,*)  "DEBUGGING: param_register json_file done"

  if(.not. param_read_line(params, args_str, ignore_unknown=.true., task='IPModel_MBX_Initialise args_str')) then
     RAISE_ERROR("IPModel_MBX_Init failed to parse args_str='"//trim(args_str)//"'", error)
  end if

  call finalise(params)

  allocate(this%nats(this%nmon))
  allocate(this%monomers(this%nmon))
  allocate(n_monomers_types_fields((this%nmon)*2))
  call split_string(n_monomers_types_string,'_','{}',n_monomers_types_fields(:),n_group4,matching=.true.)



  !do i=1, n_group4
  !  !write(*,*) ("IPModel_MBX : nats " // this%nats(i) // " monomers " // this%monomers(i))
  !  write(*,*) ("IPModel_MBX :  " // trim(adjustl(n_monomers_types_fields(i))))
  !end do



  tmp_num_monomer=0
  mon_ind=0
  do i=1,n_group4
    if (mod(i,2) == 1) then
      tmp_num_monomer = string_to_int(n_monomers_types_fields(i)) 
      !write(*,*) ("DEBUGGING tmp_num_monomer"//tmp_num_monomer)
    else
      !tmp_monomer_string=""
      tmp_monomer_string = trim(adjustl(n_monomers_types_fields(i))) !//CHAR(0)
      !write(*,*) ("DEBUGGING tmp_monomer_string"//tmp_monomer_string)
      do j=1,tmp_num_monomer
        mon_ind=mon_ind+1
        if (tmp_monomer_string .eq. 'OHH') then
          !write(*,*) "DEBUGGING if statement comparing to OHH"
          this%monomers(mon_ind) = "h2o"//CHAR(0)
          this%nats(mon_ind) = 3
        else if (tmp_monomer_string .eq. 'CHHHH') then
          !write(*,*) "DEBUGGING if statement comparing to CHHHH"
          this%monomers(mon_ind) = "ch4"//CHAR(0)
          this%nats(mon_ind) = 5 
        else if (tmp_monomer_string .eq. 'COO') then
          this%monomers(mon_ind) = "co2"//CHAR(0)
          this%nats(mon_ind) = 3
        endif 
      end do
    end if
  end do
  
  sum_nats = int(sum(this%nats))
  allocate(this%at_name(sum_nats))
  !write(*,*) ("DEBUGGING allocated at_name")
  at_ind=0
  mon_ind=0
  do i=1,n_group4
    if (mod(i,2) == 1) then
      tmp_num_monomer = string_to_int(n_monomers_types_fields(i))
    else
      tmp_monomer_string = trim(adjustl(n_monomers_types_fields(i)))
      !write(*,*) ("this tmp_monomer_string: "//tmp_monomer_string)
      do j=1,tmp_num_monomer
        !write(*,*) ("do j = "//j)
        mon_ind=mon_ind+1
        do k=1,this%nats(mon_ind)  
          !write(*,*) ("do k = "//k)
          at_ind=at_ind+1
          this%at_name(at_ind)=tmp_monomer_string(k:(k+1))
        end do
      end do
    end if
  end do

  !call split_string(nats_string,' ','{}',nats_fields(:),n_group1,matching=.true.)
  !allocate(this%nats(n_group1))
  !do i=1,n_group1
  !  this%nats(i) = string_to_int(nats_fields(i))
  !end do

  !call split_string(monomers_string,' ','{}',monomers_fields(:),n_group2,matching=.true.)
  !allocate(this%monomers(n_group2))
  !do i=1,n_group2
  !  !this%monomers(i) = monomers_fields(i)
  !  write(*,*) trim(monomers_fields(i))
  !  this%monomers(i) = trim(adjustl(monomers_fields(i)))//CHAR(0)
  !end do

  !call split_string(at_name_string,' ','{}',at_name_fields(:),n_group3,matching=.true.)
  !allocate(this%at_name(n_group3))
  !do i=1,n_group3
  !  !this%at_name(i) = at_name_fields(i)
  !  this%at_name(i) = trim(adjustl(at_name_fields(i)))//CHAR(0)
  !end do

  this%json_file = trim(adjustl(json_file_string))//CHAR(0)

  !write(*,*) ("IPModel_MBX : MBX Potential")
  !write(*,*) ("IPModel_MBX : at_name = " // this%at_name)
  !write(*,*) ("IPModel_MBX : nmon = " // this%nmon)
  !write(*,*) ("IPModel_MBX : json_file = " // this%json_file)

  !do i=1, this%nmon
  !  write(*,*) ("IPModel_MBX : nats " // this%nats(i) // " monomers " // this%monomers(i))
  !  write(*,*) ("IPModel_MBX : at_name" // this%at_name(i))
  !end do



  allocate( coord((3*sum_nats)) ) !j
  coord=0.0_dp
  call randomise(coord, 1.0_dp)

  !do i=1,int(sum_nats*3)
  !  write(*,*) ("coords" // coord(i))
  !enddo


  call initialize_system(coord, this%nats, this%at_name, this%monomers, this%nmon, this%json_file)

  !write(*,*) "Initialisation done"


  if (allocated(nats_fields)) deallocate(nats_fields)
  if (allocated(at_name_fields)) deallocate(at_name_fields)
  if (allocated(monomers_fields)) deallocate(monomers_fields)
  if (allocated(n_monomers_types_fields)) deallocate(n_monomers_types_fields)
  n_monomers_types_string=""

  write(*,*) "Initialisation: deallocated variables "


end subroutine IPModel_MBX_Initialise_str


subroutine IPModel_MBX_Finalise(this)
  type(IPModel_MBX), intent(inout) :: this
  !external finalize_system

  ! Add finalisation code here

  if (allocated(this%atomic_num)) deallocate(this%atomic_num)
  if (allocated(this%type_of_atomic_num)) deallocate(this%type_of_atomic_num)

  if (allocated(this%nats)) deallocate(this%nats)
  if (allocated(this%at_name)) deallocate(this%at_name)
  if (allocated(this%monomers)) deallocate(this%monomers)

  this%nmon = 0
  this%json_file = ''

  !call finalize_system()

end subroutine IPModel_MBX_Finalise


subroutine IPModel_MBX_Calc(this, at, e, local_e, f, virial, local_virial, args_str, mpi, error)
   type(IPModel_MBX), intent(inout):: this
   type(Atoms), intent(inout)      :: at
   real(dp), intent(out), optional :: e, local_e(:)
   real(dp), intent(out), optional :: f(:,:), local_virial(:,:)   !% Forces, dimensioned as \texttt{f(3,at%N)}, local virials, dimensioned as \texttt{local_virial(9,at%N)} 
   real(dp) :: e_kcal_mol, e_eV
   real(dp),allocatable :: force_eV_A(:,:), grads_kcal_mol_A(:)
   real(dp), intent(out), optional :: virial(3,3)
   real(dp),allocatable  :: coord(:)
   !real(dp), dimension(3) :: lattice
   character(len=*), optional      :: args_str
   type(MPI_Context), intent(in), optional :: mpi
   integer, intent(out), optional :: error
   integer :: i, j
   integer :: sum_nats

   external get_energy
   external get_energy_g

   ! Add calc() code here

   write(*,*) ("number of atoms" // at%N)
   allocate( coord(3*at%N) ) !j
   do i = 1, at%N
      do j=1,3
        write(*,*) ("positions "//j//" "//i//" "//at%pos(j,i))
        coord((i-1)*3+j) = at%pos(j,i) 
      enddo
   enddo

   !if( is_diagonal(at%lattice) ) then
   !   do i = 1, 3
   !      lattice(i) = at%lattice(i,i)
   !   enddo
   !else
   !   RAISE_ERROR('IPModel_MBX_Calc - lattice must be orthorhombic for now',error)
   !endif
   


   !write(*,*) ("nats" // this%nats)
   sum_nats = (sum(this%nats))
   write(*,*) ("sum_nats" // sum_nats)
   write(*,*) ("coord" // coord)
   ! Energy call no gradients no pbc
   if (present(f)) then
     allocate(grads_kcal_mol_A(3*at%N))
     allocate(force_eV_A(3,at%N))
     call get_energy_g(coord, sum_nats, e_kcal_mol, grads_kcal_mol_A)
     do i=1,at%N
        do j=1,3
          force_eV_A(j,i) = -grads_kcal_mol_A(3*(i-1)+j)*KCAL_MOL
          write(*,*) ("force "//j//","//i)
          write(*,*) ("force value" //force_eV_A(j,i))
        enddo
     enddo
     write(*,*) ("E / kcal_mol"//e_kcal_mol)
     e_eV = e_kcal_mol*KCAL_MOL
     write(*,*) ("E / eV"//e_eV)
   else
     if (present(e)) then
       call get_energy(coord, sum_nats, e_kcal_mol)
       write(*,*) ("E / kcal_mol"//e_kcal_mol)
       e_eV = e_kcal_mol*KCAL_MOL
       write(*,*) ("E / eV"//e_eV)
     endif
   endif
   !!write(*,*) "Testing functions that use an array of the coordinates "
   !!write(*,*)
   !!write(*,*) "Energy (no gradients) = " , 
   !!write(*,*)

   INIT_ERROR(error)

   if (present(e)) e = e_eV ! if (present(e)) e = 0.0_dp
   if (present(local_e)) then
      RAISE_ERROR('IPModel_MBX_Calc - local energies not implemented',error)
      call check_size('Local_E',local_e,(/at%N/),'IPModel_MBX_Calc', error)
      local_e = 0.0_dp
   endif
   if (present(f)) then
      call check_size('Force',f,(/3,at%Nbuffer/),'IPModel_MBX_Calc', error)
      f = force_eV_A 
   end if
   if (present(virial)) virial = 0.0_dp
   if (present(local_virial)) then
      RAISE_ERROR('IPModel_MBX_Calc - local virials not implemented',error)
      call check_size('Local_virial',local_virial,(/9,at%Nbuffer/),'IPModel_MBX_Calc', error)
      local_virial = 0.0_dp
   endif

   !RAISE_ERROR('IPModel_Calc - not implemented',error)
   if (allocated(coord)) deallocate(coord)
   if (allocated(grads_kcal_mol_A)) deallocate(grads_kcal_mol_A)
   if (allocated(force_eV_A))  deallocate(force_eV_A)

end subroutine IPModel_MBX_Calc


subroutine IPModel_MBX_Print(this, file)
  type(IPModel_MBX), intent(in) :: this
  type(Inoutput), intent(inout),optional :: file

  integer :: ti

  call Print("IPModel_MBX : MBX Potential", file=file)
  call Print("IPModel_MBX : at_name = " // this%at_name, file=file)
  call Print("IPModel_MBX : nmon = " // this%nmon, file=file)
  call Print("IPModel_MBX : json_file = " // this%json_file, file=file)

  do ti=1, this%nmon
    call Print ("IPModel_MBX : nats " // this%nats(ti) // " monomers " // this%monomers(ti), file=file)
  end do

end subroutine IPModel_MBX_Print

subroutine IPModel_MBX_read_params_xml(this, param_str)
  type(IPModel_MBX), intent(inout), target :: this
  character(len=*), intent(in) :: param_str

  type(xml_t) :: fxml

  if (len(trim(param_str)) <= 0) return

  parse_in_ip = .false.
  parse_matched_label = .false.
  parse_ip => this

  call open_xml_string(fxml, param_str)
  call parse(fxml,  &
    startElement_handler = IPModel_startElement_handler, &
    endElement_handler = IPModel_endElement_handler)
  call close_xml_t(fxml)

  if (this%n_types == 0) then
    call system_abort("IPModel_MBX_read_params_xml parsed file, but n_types = 0")
  endif

end subroutine IPModel_MBX_read_params_xml

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X 
!% XML param reader functions
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
subroutine IPModel_startElement_handler(URI, localname, name, attributes)
  character(len=*), intent(in)   :: URI
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name
  type(dictionary_t), intent(in) :: attributes

  integer :: status
  character(len=STRING_LENGTH) :: value

  integer ti

  if (name == 'MBX_params') then ! new MBX stanza

    if (parse_matched_label) return ! we already found an exact match for this label

    call QUIP_FoX_get_value(attributes, 'label', value, status)
    if (status /= 0) value = ''

    if (len(trim(parse_ip%label)) > 0) then ! we were passed in a label
      if (value == parse_ip%label) then ! exact match
        parse_matched_label = .true.
        parse_in_ip = .true.
      else ! no match
        parse_in_ip = .false.
      endif
    else ! no label passed in
      parse_in_ip = .true.
    endif

    if (parse_in_ip) then
      if (parse_ip%n_types /= 0) then
        call finalise(parse_ip)
      endif

      call QUIP_FoX_get_value(attributes, 'n_types', value, status)
      if (status == 0) then
        read (value, *) parse_ip%n_types
      else
        call system_abort("Can't find n_types in MBX_params")
      endif

      allocate(parse_ip%atomic_num(parse_ip%n_types))
      parse_ip%atomic_num = 0

      call QUIP_FoX_get_value(attributes, "cutoff", value, status)
      if (status /= 0) call system_abort ("IPModel_MBX_read_params_xml cannot find cutoff")
      read (value, *) parse_ip%cutoff
    endif


  elseif (parse_in_ip .and. name == 'per_type_data') then

    call QUIP_FoX_get_value(attributes, "type", value, status)
    if (status /= 0) call system_abort ("IPModel_MBX_read_params_xml cannot find type")
    read (value, *) ti

    call QUIP_FoX_get_value(attributes, "atomic_num", value, status)
    if (status /= 0) call system_abort ("IPModel_MBX_read_params_xml cannot find atomic_num")
    read (value, *) parse_ip%atomic_num(ti)

    if (allocated(parse_ip%type_of_atomic_num)) deallocate(parse_ip%type_of_atomic_num)
    allocate(parse_ip%type_of_atomic_num(maxval(parse_ip%atomic_num)))
    parse_ip%type_of_atomic_num = 0
    do ti=1, parse_ip%n_types
      if (parse_ip%atomic_num(ti) > 0) &
        parse_ip%type_of_atomic_num(parse_ip%atomic_num(ti)) = ti
    end do


  endif

end subroutine IPModel_startElement_handler

subroutine IPModel_endElement_handler(URI, localname, name)
  character(len=*), intent(in)   :: URI
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name

  if (parse_in_ip) then
    if (name == 'MBX_params') then
      parse_in_ip = .false.
    end if
  endif

end subroutine IPModel_endElement_handler

end module IPModel_MBX_module
