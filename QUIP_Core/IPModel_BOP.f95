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
!X IPModel_BOP module  
!X
!% Module for Bond-Order potential  
!%
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#include "error.inc"

module IPModel_BOP_module

use libatoms_module

use mpi_context_module
use QUIP_Common_module

implicit none
private

include 'IPModel_interface.h'

public :: IPModel_BOP
type IPModel_BOP
  integer :: n_types = 0
  integer, allocatable :: atomic_num(:), type_of_atomic_num(:)

  real(dp) :: cutoff = 0.0_dp

  integer     ::  n
  type(Atoms) :: at
  integer, allocatable, dimension(:)  :: aptr
  integer, allocatable, dimension(:)  :: bptr
  integer     :: bptr_num

  character(len=FIELD_LENGTH) :: label

end type IPModel_BOP

logical, private :: parse_in_ip, parse_matched_label
type(IPModel_BOP), private, pointer :: parse_ip

interface Initialise
  module procedure IPModel_BOP_Initialise_str
end interface Initialise

interface Finalise
  module procedure IPModel_BOP_Finalise
end interface Finalise

interface Print
  module procedure IPModel_BOP_Print, IPModel_BOP_Print_ptr
end interface Print

interface Calc
  module procedure IPModel_BOP_Calc
end interface Calc

contains

subroutine IPModel_BOP_Initialise_str(this, args_str, param_str)
  type(IPModel_BOP), intent(inout) :: this
  character(len=*), intent(in) :: args_str, param_str

  type(Dictionary) :: params

  call Finalise(this)

  call initialise(params)
  this%label=''
  call param_register(params, 'label', '', this%label, help_string="No help yet.  This source file was $LastChangedBy$")
  if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='IPModel_BOP_Initialise_str args_str')) then
    call system_abort("IPModel_BOP_Initialise_str failed to parse label from args_str="//trim(args_str))
  endif
  call finalise(params)

  call IPModel_BOP_read_params_xml(this, param_str)

end subroutine IPModel_BOP_Initialise_str

subroutine IPModel_BOP_Finalise(this)
  type(IPModel_BOP), intent(inout) :: this

  if (allocated(this%atomic_num)) deallocate(this%atomic_num)
  if (allocated(this%type_of_atomic_num)) deallocate(this%type_of_atomic_num)

  this%n_types = 0
  this%label = ''
end subroutine IPModel_BOP_Finalise


! nat is the number of atoms whose energy and forces have to be computed with bop library
subroutine IPModel_BOP_Calc(this, at, e, local_e, f, virial, local_virial, args_str, mpi, error)
   type(IPModel_BOP), intent(inout):: this
   type(Atoms), intent(inout)      :: at   !Active + buffer atoms
   real(dp), intent(out), optional :: e, local_e(:)
   real(dp), intent(out), optional :: f(:,:), local_virial(:,:)   !% Forces, dimensioned as \texttt{f(3,at%N)}, local virials, dimensioned as \texttt{local_virial(9,at%N)} 
   real(dp), intent(out), optional :: virial(3,3)
   character(len=*), optional      :: args_str
   type(MPI_Context), intent(in), optional :: mpi
   integer, intent(out), optional :: error

   real(dp), allocatable, dimension(:,:) :: f_bop
!   real(dp), allocatable, dimension(:)   :: local_e_bop 
   integer, allocatable, dimension(:)    :: map_tobop
   type(Atoms)                     :: at_bop   !Active + buffer atoms
   type(Atoms)                     :: at_tmp   
   real(dp)                        :: e_bop
   integer                         :: nat  !Number of active atoms, the first nat of at
   type(Dictionary)                :: params
   integer                         :: i, nreplicate_x, nreplicate_y, nreplicate_z
   logical                         :: lpbc, lreplicate_cell, lprint_xyz, lcrack,latoms
   integer, pointer, dimension(:)  :: map,map_bop, active, original_cell, hybrid 
   real(dp), pointer, dimension(:,:):: forces_bop
!   real(dp), pointer, dimension(:) :: local_energy 
   integer :: nsafe, natoms_active

   logical :: has_atom_mask_name
   character(FIELD_LENGTH) :: atom_mask_name
   real(dp) :: r_scale, E_scale
   logical :: do_rescale_r, do_rescale_E

   INIT_ERROR(error)

   if (present(e)) e = 0.0_dp
   if (present(local_e)) then
      call check_size('Local_E',local_e,(/at%N/),'IPModel_BOP_Calc', error)
      local_e = 0.0_dp
   endif
   if (present(f)) then
      call check_size('Force',f,(/3,at%Nbuffer/),'IPModel_BOP_Calc', error)
      f = 0.0_dp
   end if
   if (present(virial)) virial = 0.0_dp
   if (present(local_virial)) then
      call check_size('Local_virial',local_virial,(/9,at%Nbuffer/),'IPModel_BOP_Calc', error)
      local_virial = 0.0_dp
      RAISE_ERROR("IPModel_BOP_Calc: local_virial calculation requested but not supported yet.", error)
   endif

   lpbc = .false.
   if (present(args_str)) then
     call initialise(params)  
     call param_register(params, 'nat', '0', nat, help_string="No help yet.  This source file was $LastChangedBy$")
     call param_register(params, 'pbc', '.false.', lpbc, help_string="No help yet.  This source file was $LastChangedBy$")
     call param_register(params, 'replicate_x', '1', nreplicate_x, help_string="No help yet.  This source file was $LastChangedBy$")
     call param_register(params, 'replicate_y', '1', nreplicate_y, help_string="No help yet.  This source file was $LastChangedBy$")
     call param_register(params, 'replicate_z', '1', nreplicate_z, help_string="No help yet.  This source file was $LastChangedBy$")
     call param_register(params, 'print_xyz', 'F', lprint_xyz, help_string="No help yet.  This source file was $LastChangedBy$")
     call param_register(params, 'crack', 'F', lcrack, help_string="No help yet.  This source file was $LastChangedBy$")
     call param_register(params, 'all_atoms', 'F', latoms, help_string="No help yet.  This source file was $LastChangedBy$")
     call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, help_string="No help yet.  This source file was $LastChangedBy$")
     call param_register(params, 'r_scale', '1.0',r_scale, has_value_target=do_rescale_r, help_string="Recaling factor for distances. Default 1.0.")
     call param_register(params, 'E_scale', '1.0',E_scale, has_value_target=do_rescale_E, help_string="Recaling factor for energy. Default 1.0.")

     if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='IPModel_BOP_Calc args_str')) then
       call system_abort("IPModel_BOP_Calc failed to parse args_str='"//trim(args_str)//"'")
     endif
     call finalise(params)

     if(has_atom_mask_name) then
        RAISE_ERROR('IPModel_BOP_Calc: atom_mask_name found, but not supported', error)
     endif
     if (do_rescale_r .or. do_rescale_E) then
        RAISE_ERROR("IPModel_BOP_Calc: rescaling of potential with r_scale and E_scale not yet implemented!", error)
     end if

   else
     nat = 1 
   endif


   if(latoms.and.lcrack) call system_abort('IPModel_BOP: all_atoms ' // latoms// " and crack " // lcrack) 

   if(nat.eq.0) then
     if(lcrack) then
       if (.not. assign_pointer(at, 'hybrid', hybrid)) &
         call system_abort('IPModel_BOP: atoms structure is missing hybrid property')
       natoms_active = count(hybrid/=0)
       this%n = at%N
     elseif(latoms) then
       this%n = at%N
     else
       this%n = 1
     endif
   else
      this%n = nat 
   endif

   lreplicate_cell = .false.
   if(nreplicate_x.ne.1.or.nreplicate_y.ne.1.or.nreplicate_z.ne.1) then
      lreplicate_cell = .true.
   endif

   if(lreplicate_cell.and.lpbc) then
        call system_abort("IPModel_BOP_Calc: or PBC or replicates the unit cell") 
   endif

   at_tmp = at
   call add_property(at_tmp, 'map', 0)
   if (.not. assign_pointer(at_tmp, 'map', map)) &
       call system_abort('active pointer assignment failed')
   do i=1,at_tmp%N
      map(i) = i
   enddo

   if(lreplicate_cell) then
      nsafe = 1
      call print("Replicate cell : "//nreplicate_x//" "// nreplicate_y //" "// nreplicate_z, PRINT_NERD)
!      call supercell_minus_plus(at_bop, at_tmp, nreplicate_x, nsafe, nsafe)
!      call supercell_minus_plus(at_tmp, at_bop, nsafe, nreplicate_y, nsafe)
!      call supercell_minus_plus(at_bop, at_tmp, nsafe, nsafe,nreplicate_z)
      call supercell(at_bop, at_tmp, nreplicate_x, nsafe, nsafe)
      call supercell(at_tmp, at_bop, nsafe, nreplicate_y, nsafe)
      call supercell(at_bop, at_tmp, nsafe, nsafe,nreplicate_z)
   elseif(lpbc) then
      call print ("Periodic boundary condition applied", PRINT_NERD)
      call IPModel_BOP_compute_buffer(this,at,at_bop)
      at_bop = at_tmp
   else   ! no pbc
      at_bop = at_tmp
   endif 
!  Assign map property for computing forces with BOP library
   if (.not. has_property(at_bop, 'map')) &
         call system_abort('IPModel_BOP: atoms structure has no "map" property')
   if (.not. assign_pointer(at_bop, 'map', map_bop)) &
         call system_abort('IPModel_BOP passed atoms structure with no map property')
   allocate(map_tobop(at_bop%N))
   map_tobop = map_bop

!  A vacuum region surrounding the cluster has to be applied before computing the connecting
   at_bop%lattice = at_bop%lattice * 3.0_dp
!  Cutoff set to the potential value in order to reduce the number of neighbours for the BOP library 
   call set_cutoff(at_bop, this%cutoff)
!  Recompute connectivity
   call print('IPModel_BOP : cutoff used for connectivity' // at_bop%cutoff, PRINT_NERD)
   call set_lattice(at_bop,at_bop%lattice,scale_positions=.false., remap=.true.,reconnect=.true.)

!  Atoms to be passed to the Bop library
   this%at = at_bop 

   if(this%n.gt.this%at%N) then
     call system_abort("IPModel_BOP_Calc: number of atoms for BOP library "// this%n // " > at%N = " // this%at%N) 
   endif

   if (present(e)) e = 0.0_dp
   if (present(local_e)) then
      local_e = 0.0_dp
      if(size(local_e).gt.this%at%N) then
        call system_abort("IPModel_BOP_Calc: reallocate local_e to the right value = (" // this%n // ") !")
      endif 
   endif
   if (present(f)) then 
     f = 0.0_dp
     if(size(f,2).gt.this%at%N) then
        call system_abort("IPModel_BOP_Calc: reallocate f to the right value = (3," // this%n // ") !")
     endif 
   endif
   if (present(virial)) virial = 0.0_dp
  
!  Compute the two pointer vectors required by the BOP library for the at_bop object
   call print('Computing pointers for BOP library', PRINT_NERD)
   call IPModel_BOP_Calc_ptr(this)

   if(allocated(f_bop)) deallocate(f_bop)
   allocate(f_bop(3,this%n))
!   if(allocated(local_e_bop)) deallocate(local_e_bop)
!   allocate(local_e_bop(this%n))
   f_bop = 0.0_dp
   e_bop = 0.0_dp

!  Compute forces and energy with BOP library 
#ifdef HAVE_BOP
   call print('Calling BOP library : active atoms='// this%n// " on "//this%at%N )
   call system_timer('BOP')
   call bop(this%n,this%at%N,this%at%pos,this%aptr,this%bptr_num,this%bptr,map_tobop,e_bop,f_bop)
   call system_timer('BOP')
   call print('BOP energy : ' // e_bop)
#else
   call print('No BOP library!!!')
#endif

!  Print properties if print_xyz=T is passed in args_str
   if(lprint_xyz) then
     call add_property(at_bop, 'active', 0)
     if (.not. assign_pointer(at_bop, 'active', active)) &
         call system_abort('active pointer assignment failed')
     active(1:this%n) = 1
     call add_property(at_bop, 'original_cell', 0)
     if (.not. assign_pointer(at_bop, 'original_cell', original_cell)) &
         call system_abort('original_cell pointer assignment failed')
     original_cell(1:at%n) = 1
     call add_property(at_bop, 'forces_bop', 0.0_dp, n_cols=3)
     if (.not. assign_pointer(at_bop, 'forces_bop', forces_bop)) &
         call system_abort('forces_bop: failed to assign pointer to forces_bop')
     forces_bop(:,1:this%n) = f_bop(:,1:this%n)
!     call add_property(at_bop, 'local_energy', 0.0_dp)
!     if (.not. assign_pointer(at_bop, 'local_energy', local_energy)) &
!         call system_abort('local_energy: failed to assign pointer to local_bop')
!     local_energy(1:this%n) = local_e_bop(1:this%n)
     call write(at_bop, "bop_cell.xyz", properties='active:original_cell:map:forces_bop')
   endif

   if(present(f)) then
!     if (present(mpi)) call sum_in_place(mpi, f_bop)
     if(lcrack) then
       f(:,1:natoms_active) = f_bop(:,1:natoms_active)  
     else
       f(:,1:this%n) = f_bop(:,1:this%n)  
     endif
   endif
   if(present(e)) e = e_bop 
!   if (present(mpi) .and. present(e)) e = sum(mpi, e_bop)
   call print('Energy computed with BOP library : ' // e_bop // " eV ",PRINT_NERD)
   if (present(mpi)) then
      if (present(local_e)) call sum_in_place(mpi, local_e)
   endif
   if(current_verbosity()  >= PRINT_NERD) then
     do i =1, this%n
       call print('Forces computed with BOP library on atom ' // i // " : " // f_bop(1,i) // " "// f_bop(2,i) // " " // f_bop(3,i) )
     enddo
   endif
    
   call finalise(at_bop)
   call finalise(at_tmp)

end subroutine IPModel_BOP_Calc

! compute aptr and bptr vectors required by the BOP library.
subroutine IPModel_BOP_Calc_ptr(this)
   type(IPModel_BOP), intent(inout) :: this
   integer                          :: nn, iat
   integer                          :: ik, k, itemp
   real(dp)                         :: rik 

   if(allocated(this%aptr)) deallocate(this%aptr)
   allocate(this%aptr(this%at%N))
   if(allocated(this%bptr)) deallocate(this%bptr)
   nn = 0
   do iat = 1, this%at%N 
     nn = nn + atoms_n_neighbours(this%at, iat) + 2
   enddo
   nn = nn + 1
   this%bptr_num = nn
   allocate(this%bptr(this%bptr_num))

   itemp =  -1000000
   nn = 0
   do iat = 1, this%at%N
      nn = nn + 1
      this%bptr(nn) = iat 
      this%aptr(iat) = nn
      if(nn.ge.this%bptr_num) then
         call system_abort("IPModel_BOP_Calc_ptr: nn > this%bptr_num")
      endif
      do ik = 1, atoms_n_neighbours(this%at, iat)
         nn = nn + 1
         k = atoms_neighbour(this%at, iat, ik, rik)
         this%bptr(nn) = k
         if(nn.ge.this%bptr_num) then
           call system_abort("IPModel_BOP_Calc_ptr: nn > this%bptr_num")
         endif
      enddo
      nn = nn + 1
      this%bptr(nn) = itemp 
      if(nn.ge.this%bptr_num) then
         call system_abort("IPModel_BOP_Calc_ptr: nn > this%bptr_num")
      endif
   enddo
   this%bptr(this%bptr_num) = 0
   if(current_verbosity()  >= PRINT_VERBOSE) then
     call print(this,this%aptr,this%bptr)
   endif

end subroutine IPModel_BOP_Calc_ptr

subroutine IPModel_BOP_Print(this, file)
  type(IPModel_BOP), intent(in) :: this
  type(Inoutput), intent(inout),optional :: file

  integer :: ti

  call Print("IPModel_BOP : Bond-Order-Potential", file=file)
  call Print("IPModel_BOP : n_types = " // this%n_types // " cutoff = " // this%cutoff, file=file)

  do ti=1, this%n_types
    call Print ("IPModel_BOP : type " // ti // " atomic_num " // this%atomic_num(ti), file=file)
    call verbosity_push_decrement()
    call Print ("IPModel_BOP : " // &
        "cutoff " // this%cutoff , &
         file=file)
    call verbosity_pop()
  end do

end subroutine IPModel_BOP_Print

subroutine IPModel_BOP_Print_ptr(this, aptr, bptr)
  type(IPModel_BOP), intent(in) :: this
  integer, intent(inout), dimension(:)  :: aptr
  integer, intent(inout), dimension(:)  :: bptr
  type(inoutput)   :: file_a
  type(inoutput)   :: file_b
  integer          :: i

  call print("Printing aptr and bptr")
  call initialise(file_a,'aptr.in', OUTPUT)
  call initialise(file_b,'bptr.in', OUTPUT)
  call print(this%at%N, file=file_a)
  do i = 1, this%at%N
    call print("  " // i // "    " // this%aptr(i), file=file_a)
  enddo
  call print(this%bptr_num, file=file_b)
  do i = 1, this%bptr_num
    call print("  " // i // "    " // this%bptr(i), file=file_b)
  enddo
  call finalise(file_a)
  call finalise(file_b)

end subroutine IPModel_BOP_Print_ptr

subroutine IPModel_BOP_read_params_xml(this, param_str)
  type(IPModel_BOP), intent(inout), target :: this
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
    call system_abort("IPModel_BOP_read_params_xml parsed file, but n_types = 0")
  endif

end subroutine IPModel_BOP_read_params_xml

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
  character(len=FIELD_LENGTH) :: value

  integer ti

  if (name == 'BOP_params') then ! new BOP stanza

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
        read (value, *), parse_ip%n_types
      else
        call system_abort("Can't find n_types in BOP_params")
      endif

      allocate(parse_ip%atomic_num(parse_ip%n_types))
      parse_ip%atomic_num = 0

      call QUIP_FoX_get_value(attributes, "cutoff", value, status)
      if (status /= 0) call system_abort ("IPModel_BOP_read_params_xml cannot find cutoff")
      read (value, *) parse_ip%cutoff
    endif


  elseif (parse_in_ip .and. name == 'per_type_data') then

    call QUIP_FoX_get_value(attributes, "type", value, status)
    if (status /= 0) call system_abort ("IPModel_BOP_read_params_xml cannot find type")
    read (value, *) ti

    call QUIP_FoX_get_value(attributes, "atomic_num", value, status)
    if (status /= 0) call system_abort ("IPModel_BOP_read_params_xml cannot find atomic_num")
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
    if (name == 'BOP_params') then
      parse_in_ip = .false.
    end if
  endif

end subroutine IPModel_endElement_handler

subroutine IPModel_BOP_compute_buffer(this,at,at_bop)
  type(IPModel_BOP), intent(inout):: this
  type(Atoms), intent(inout)  :: at
  type(Atoms), intent(inout) :: at_bop
  integer                  :: cellsNa,cellsNb,cellsNc, n_atom
  integer                  :: cellsNa2,cellsNb2,cellsNc2
  real(dp)                 :: cutoff
  integer                  :: i_atom, i, j, k
  type(table)              :: list

  cutoff = this%cutoff
  call print("calc_connect: cutoff calc_connect " // cutoff)
  cellsNa = at%connect%cellsNa
  cellsNb = at%connect%cellsNb
  cellsNc = at%connect%cellsNc
  cellsNa2 =  nint(real(cellsNa,dp) / 2.0_dp )
  cellsNb2 =  nint(real(cellsNb,dp) / 2.0_dp )
  cellsNc2 =  nint(real(cellsNc,dp) / 2.0_dp )
  call print("calc_connect: cells_N[abc] " // cellsNa // " " // cellsNb // " " // cellsNc)

  call calc_connect(at)

! pick an atom in the middle of the unit box 
  i_atom = -1
  loop_c : do k = cellsNc2, cellsNc
    loop_b : do j = cellsNb2, cellsNb
       loop_a : do i = cellsNa2, cellsNa
         n_atom =  at%connect%cell(i,j,k)%N
         if(n_atom.ne.0) then
            i_atom = at%connect%cell(i,j,k)%int(1,1)
            exit loop_c 
         endif
       enddo  loop_a
    enddo  loop_b
  enddo loop_c
  if(i_atom == -1) call system_abort("IPModel_BOP_compute_buffer : no atom in the halved simulation cell")

  call append(list, (/i_atom,0,0,0/)) 
  call BFS_grow(at, list, 4)
  call print('---------------cluster--------------------', PRINT_NERD)
  call print(list, PRINT_NERD)
  call print('--------------end cluster-----------------', PRINT_NERD)

  call build_cluster(at, list, at_bop)

end subroutine IPModel_BOP_compute_buffer
 
subroutine  build_cluster(at, list, at_bop)
  type(Atoms), intent(inout) :: at
  type(Atoms), intent(inout) :: at_bop
  type(table)                :: list
  real(dp), dimension(3)     :: dist, rmax, length
  integer, dimension(3)      :: n_cell_rep
  integer                    :: iatom, jatom, i

  call print ('Number of atoms belonging to the cluster : ' // list%N )
  call print ('Cluster center atom                      : ' // list%int(1,1) )   
  call print ('               shift                     : ' // list%int(2:4,1) )   

! Find how many cells have to be added 
  iatom = list%int(1,1)
  rmax = 0.0_dp 
  do i = 2, list%N
    jatom = list%int(1,i) 
    dist = abs( at%pos(:,iatom) -  pos_j(at, jatom, list%int(2:4,i)) ) 
    if(dist(1) .gt. rmax(1)) rmax(1) = dist(1) 
    if(dist(2) .gt. rmax(2)) rmax(2) = dist(2) 
    if(dist(3) .gt. rmax(3)) rmax(3) = dist(3) 
  enddo
  length(1) = at%lattice(1,1) / real(at%connect%cellsNa,dp)
  length(2) = at%lattice(2,2) / real(at%connect%cellsNb,dp)
  length(3) = at%lattice(3,3) / real(at%connect%cellsNc,dp)

  n_cell_rep = int(rmax/length) + 1
  call print("Number of cell repeated along the three directions : " // n_cell_rep(1:3) ) 

! Add to the simulation cell the just computed number of cells
  call add_cells(at, n_cell_rep, at_bop) 
  
end subroutine  build_cluster

subroutine add_cells(at, n_cell_rep, at_bop)
  type(Atoms), intent(inout) :: at
  type(Atoms), intent(inout) :: at_bop
  integer, dimension(3)      :: cells, n_cell_rep, periodicity
  integer                    :: i, j, k, ic, ic_right_p, ic_right_m, ii, iatom
  real(dp)                   :: pos(3), length

  cells(1) = at%connect%cellsNa
  cells(2) = at%connect%cellsNb
  cells(3) = at%connect%cellsNc
  
  at_bop = at
  do i = 1, 3    ! x, y, z
    if(n_cell_rep(i).eq.0) cycle  !skip this direction
    do ic = 1, n_cell_rep(i)  
       ic_right_p = cells(i) - ic + 1 + int((ic - 1) / cells(i)) * cells(i)
       ic_right_m = ic - 1 - int((ic - 1) / cells(i)) * cells(i)

       periodicity    = 0
       periodicity(i) = int( (ic-1)/cells(i) + 1 )
       call print ("ic_right_p   " // ic_right_p  // "  ; ic    " // ic, PRINT_NERD)
       call print ("ic_right_m   " // ic_right_m  // "  ; ic    " // ic, PRINT_NERD)
       call print ("periodicity  " // periodicity, PRINT_NERD)

       if(i.eq.1) then
          do j = 1, cells(2)
           do k = 1, cells(3)
            do ii = 1, at%connect%cell(ic_right_p,j,k)%N
              iatom = at%connect%cell(ic_right_p,j,k)%int(1,ii)
              call print( j //" "// k // " atom : " // iatom // " " // at%connect%cell(ic_right_p,j,k)%N ) 
              pos = at%pos(:,iatom) + real(periodicity,dp) * at%lattice(i,i) 
              call add_atoms(at_bop, pos, at%Z(iatom))  
            enddo

            do ii = 1, at%connect%cell(ic_right_m,j,k)%N  
              iatom = at%connect%cell(ic_right_m,j,k)%int(1,ii)
              call print( j //" "// k // " atom : " // iatom // " " // at%connect%cell(ic_right_m,j,k)%N ) 
              pos = at%pos(:,iatom) + real(periodicity,dp) * at%lattice(i,i)
              call add_atoms(at_bop, pos, at%Z(iatom))
             enddo
           enddo
          enddo
       elseif(i.eq.2) then 

       elseif(i.eq.3) then 

       endif

    enddo
    length = at%lattice(i,i) / cells(i) 
    at_bop%lattice(i,i) = at_bop%lattice(i,i) + n_cell_rep(i) * length * 2.0_dp 

  enddo

  stop 

end subroutine add_cells

function pos_j(at, j, shift)
  type(Atoms), intent(inout) :: at
  real(dp), dimension(3)     :: pos_j
  integer,  dimension(3)     :: shift
  integer                    :: j

  pos_j = (at%lattice .mult. shift) + at%pos(:,j)
  
end function pos_j

!!$  subroutine supercell_minus_plus(aa, a, n1, n2, n3)
!!$    type(Atoms), intent(out)::aa  !% Output (big) cell
!!$    type(Atoms), intent(in)::a    !% Input cell
!!$    integer, intent(inout)::n1, n2, n3
!!$    real(dp)::lattice(3,3), p(3)
!!$    integer::i,j,k,n
!!$    type(Table) :: big_data
!!$    integer     :: nn1, n1save, n2save, n3save
!!$    integer, pointer, dimension(:)  :: map_new, map
!!$
!!$    n1save = n1
!!$    n2save = n2
!!$    n3save = n3
!!$    n1 = 2 * n1 - 1
!!$    n2 = 2 * n2 - 1
!!$    n3 = 2 * n3 - 1
!!$    
!!$    call allocate(big_data,a%data%intsize,a%data%realsize,&
!!$         a%data%strsize, a%data%logicalsize, a%N*n1*n2*n3)
!!$
!!$    ! Replicate atomic data n1*n2*n3 times
!!$    do i=1,n1*n2*n3
!!$       call append(big_data,a%data)
!!$    end do
!!$
!!$    lattice(:,1) = a%lattice(:,1)*n1save
!!$    lattice(:,2) = a%lattice(:,2)*n2save
!!$    lattice(:,3) = a%lattice(:,3)*n3save
!!$    call atoms_initialise(aa, a%N*n1*n2*n3, lattice, data=big_data, properties=a%properties)
!!$    if (a%use_uniform_cutoff) then
!!$       call set_cutoff(aa, a%cutoff)
!!$    else
!!$       call set_cutoff_factor(aa, a%cutoff)
!!$    end if
!!$
!!$   if (.not. assign_pointer(a, 'map', map)) &
!!$       call system_abort('active pointer assignment failed')
!!$   if (.not. assign_pointer(aa, 'map', map_new)) &
!!$       call system_abort('active pointer assignment failed')
!!$
!!$    nn1 = 0
!!$    do i = 0,n1save-1
!!$       do j = 0,n2save-1
!!$          do k = 0,n3save-1
!!$             p = a%lattice .mult. (/i,j,k/)
!!$             do n = 1,a%N
!!$                ! overwrite position with shifted pos
!!$                nn1 = nn1 +  1
!!$                aa%pos(:,nn1) = a%pos(:,n)+p
!!$                map_new(nn1) = map(n) 
!!$             end do
!!$             if(i.eq.0.and.j.eq.0.and.k.eq.0) cycle
!!$             do n = 1,a%N
!!$                ! overwrite position with shifted pos
!!$                nn1 = nn1 +  1
!!$                aa%pos(:,nn1) = a%pos(:,n)-p
!!$                map_new(nn1) = map(n) 
!!$             end do
!!$          end do
!!$       end do
!!$    end do
!!$
!!$    call finalise(big_data)
!!$
!!$  end subroutine supercell_minus_plus

end module IPModel_BOP_module
