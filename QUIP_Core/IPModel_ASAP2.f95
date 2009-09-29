!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X     QUIP: quantum mechanical and interatomic potential simulation package
!X     
!X     Portions written by Noam Bernstein, while working at the
!X     Naval Research Laboratory, Washington DC. 
!X
!X     Portions written by Gabor Csanyi, Copyright 2006-2007.   
!X
!X     When using this software,  please cite the following reference:
!X
!X     reference
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X IPModel_ASAP2
!X
!% Interface to ASAP2 potential.
!% P. Tangney and S. Scandolo,
!% An ab initio parametrized interatomic force field for silica
!% J. Chem. Phys, 117, 8898 (2002). 
!%
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
module IPModel_ASAP2_module

use libatoms_module

use mpi_context_module
use QUIP_Common_module

implicit none
private

include 'IPModel_interface.h'

logical, private :: asap_initialised = .false.

public :: IPModel_ASAP2
type IPModel_ASAP2
  integer :: n_types = 0
  real(dp) :: betapol, tolpol, yukalpha, yuksmoothlength
  integer :: maxipol, pred_order
  integer, allocatable :: atomic_num(:), type_of_atomic_num(:)
  real(dp), allocatable, dimension(:) :: pol, z
  real(dp), allocatable, dimension(:,:) :: D_ms, gamma_ms, R_ms, B_pol, C_pol
  integer :: iesr(3)

  real(dp) :: cutoff_coulomb, cutoff_ms

  character(len=FIELD_LENGTH) :: label
  type(mpi_context) :: mpi
  logical :: initialised

end type IPModel_ASAP2

logical :: parse_in_ip, parse_matched_label
type(IPModel_ASAP2), pointer :: parse_ip

interface Initialise
  module procedure IPModel_ASAP2_Initialise_str
end interface Initialise

interface Finalise
  module procedure IPModel_ASAP2_Finalise
end interface Finalise

interface Print
  module procedure IPModel_ASAP2_Print
end interface Print

interface Calc
  module procedure IPModel_ASAP2_Calc
end interface Calc

contains


subroutine IPModel_ASAP2_Initialise_str(this, args_str, param_str, mpi)
  type(IPModel_ASAP2), intent(inout) :: this
  character(len=*), intent(in) :: args_str, param_str
  type(mpi_context), intent(in), optional :: mpi

  type(Dictionary) :: params

  this%initialised = .false.

  call Finalise(this)

  call initialise(params)
  this%label=''
  call param_register(params, 'label', '', this%label)
  if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='IPModel_ASAP2_Initialise_str args_str')) then
    call system_abort("IPModel_ASAP2_Initialise_str failed to parse label from args_str="//trim(args_str))
  endif
  call finalise(params)

  call IPModel_ASAP2_read_params_xml(this, param_str)
  this%initialised = .true.
  
  if (present(mpi)) this%mpi = mpi

end subroutine IPModel_ASAP2_Initialise_str

subroutine IPModel_ASAP2_Finalise(this)
  type(IPModel_ASAP2), intent(inout) :: this

  this%initialised = .false.

  if (allocated(this%atomic_num)) deallocate(this%atomic_num)
  if (allocated(this%type_of_atomic_num)) deallocate(this%type_of_atomic_num)
  if (allocated(this%pol)) deallocate(this%pol)
  if (allocated(this%z)) deallocate(this%z)
  if (allocated(this%D_ms)) deallocate(this%D_ms)
  if (allocated(this%gamma_ms)) deallocate(this%gamma_ms)
  if (allocated(this%R_ms)) deallocate(this%R_ms)
  if (allocated(this%B_pol)) deallocate(this%B_pol)
  if (allocated(this%C_pol)) deallocate(this%C_pol)
  this%n_types = 0
  this%label = ''
end subroutine IPModel_ASAP2_Finalise

!% Smooth cutoff function from D.J. Cole \emph{et al.}, J. Chem. Phys. {\bf 127}, 204704 (2007).
subroutine smooth_cutoff(x,R,D,fc,dfc_dx)

  real(dp) x,R,D,fc,dfc_dx
  real(dp), parameter :: pi = dacos(-1.0d0)

  if (x .lt. (R-D)) then
     fc = 1.0d0
     dfc_dx = 0.0d0
     return
  else if (x .gt. (R+D)) then
     fc = 0.0d0
     dfc_dx = 0.0d0
     return
  else
     fc = 1.0d0 - (x-R+D)/(2.0d0*D) + 1.0d0/(2.0d0*pi)*dsin((pi/D)*(x-R+D))
     dfc_dx = 1.0d0/(2.0d0*D)* (dcos((pi/D)*(x-R+D)) - 1.0d0)
     return
  end if
end subroutine smooth_cutoff


subroutine asap_rs_charges(this, at, e, local_e, f, virial, field, args_str)
   type(IPModel_ASAP2), intent(inout):: this
   type(Atoms), intent(inout)      :: at
   real(dp), intent(out), optional :: e, local_e(:)
   real(dp), intent(out), optional :: f(:,:)
   real(dp), intent(out), optional :: virial(3,3)
   real(dp), intent(out), optional :: field(:,:)
   character(len=*), optional, intent(in) :: args_str

   integer i, j, m, ti, tj
   real(dp) :: r_ij, u_ij(3), zv2, gamjir, gamjir3, gamjir2, fc, dfc_dr
   real(dp) :: de, dforce, expfactor

   do i=1, at%n
      ti = get_type(this%type_of_atomic_num, at%Z(i))
      do m = 1, atoms_n_neighbours(at, i, max_dist=this%cutoff_coulomb)
         
         j = atoms_neighbour(at, i, m, distance=r_ij, cosines=u_ij, max_dist=this%cutoff_coulomb)
         if (j <= 0) cycle
         if (r_ij .feq. 0.0_dp) cycle

         r_ij = r_ij/BOHR
         tj = get_type(this%type_of_atomic_num, at%Z(j))
         zv2 = this%z(ti)*this%z(tj)

         gamjir = zv2/r_ij
         gamjir3 = gamjir/(r_ij**2.0_dp)
         gamjir2 = zv2/(r_ij**2.0_dp)
         expfactor = exp(-this%yukalpha*r_ij)

         call smooth_cutoff(r_ij, this%cutoff_coulomb-this%yuksmoothlength, &
              this%yuksmoothlength, fc, dfc_dr)

         de = gamjir

         if (present(e) .or. present(local_e)) then
            if (present(e))       e = e + 0.5_dp*de*expfactor*fc
            if (present(local_e)) local_e(i) = local_e(i) + 0.5_dp*de*expfactor*fc
         end if

         if (present(f) .or. present(virial) .or. present(field)) then
            dforce = gamjir3*expfactor*fc*r_ij + de*(this%yukalpha*fc - dfc_dr)*expfactor

            if (present(f))      f(:,i) = f(:,i) - dforce*u_ij
            if (present(virial)) virial = virial + 0.5_dp*dforce*(u_ij .outer. u_ij)*r_ij
            if (present(field))  field(:,i) = field(:,i) + gamjir3*expfactor*fc*r_ij/this%z(ti)
         end if
      end do
   end do
end subroutine asap_rs_charges


subroutine asap_rs_dipoles(this, at, e, local_e, f, virial, args_str)
   type(IPModel_ASAP2), intent(inout):: this
   type(Atoms), intent(inout)      :: at
   real(dp), intent(out), optional :: e, local_e(:)
   real(dp), intent(out), optional :: f(:,:)
   real(dp), intent(out), optional :: virial(3,3)
   character(len=*), optional, intent(in) :: args_str

   
end subroutine asap_rs_dipoles


!% Morse-stretch potential, defined by
!% \begin{displaymath}
!% U_ij = D_ij \left[ e^{\gamma_{ij}\left( 1 - r_{ij}/r^0_{ij} \right)} 
!%                  - 2 e^{\left( 1 - r_{ij}/r^0_{ij} \right)} \right]
!% \end{displaymath}
subroutine asap_morse_stretch(this, at, e, local_e, f, virial, args_str)
   type(IPModel_ASAP2), intent(inout):: this
   type(Atoms), intent(inout)      :: at
   real(dp), intent(out), optional :: e, local_e(:)
   real(dp), intent(out), optional :: f(:,:)
   real(dp), intent(out), optional :: virial(3,3)
   character(len=*), optional, intent(in) :: args_str

   integer i, j, m, ti, tj
   real(dp) :: r_ij, u_ij(3), dms, gammams, rms
   real(dp) :: exponentms, factorms, phi, de
   real(dp) :: dforce
   real(dp) :: elimitij(this%n_types, this%n_types)

   ! Evaluate potential at cutoff. Will be subtracted from total energy.
   elimitij = 0.0_dp
   do ti=1,this%n_types
      do tj=1,this%n_types
         phi = exp(this%gamma_ms(ti,tj)*(1.0_dp - this%cutoff_ms/this%r_ms(ti,tj)))
         elimitij(ti,tj) = this%d_ms(ti,tj)*(phi - 2.0_dp*sqrt(phi))
      end do
   end do
 
   do i=1, at%n
      ti = get_type(this%type_of_atomic_num, at%Z(i))
      do m = 1, atoms_n_neighbours(at, i, max_dist=this%cutoff_ms)
         
         j = atoms_neighbour(at, i, m, distance=r_ij, cosines=u_ij, max_dist=this%cutoff_ms)
         if (j <= 0) cycle
         if (r_ij .feq. 0.0_dp) cycle

         r_ij = r_ij/BOHR
         tj = get_type(this%type_of_atomic_num, at%Z(j))

         dms = this%d_ms(ti,tj)
         gammams = this%gamma_ms(ti,tj)
         rms = this%r_ms(ti, tj)

         exponentms = gammams*(1.0_dp-r_ij/rms)
         factorms = gammams/rms
         phi = exp(exponentms)

         if (present(e) .or. present(local_e)) then
            de = dms*(phi-2.0_dp*sqrt(phi)) - elimitij(ti, tj)

            if (present(e))       e = e + 0.5_dp*de
            if (present(local_e)) local_e(i) = local_e(i) + 0.5_dp*de
         end if

         if (present(f) .or. present(virial)) then
            dforce  = -dms*(factorms*phi - factorms*dsqrt(phi))

            if (present(f))      f(:,i) = f(:,i) + dforce*u_ij
            if (present(virial)) virial = virial - 0.5_dp*dforce*(u_ij .outer. u_ij)*r_ij
         end if
      end do
   end do

end subroutine asap_morse_stretch


subroutine IPModel_ASAP2_Calc(this, at, e, local_e, f, virial, args_str)
   type(IPModel_ASAP2), intent(inout):: this
   type(Atoms), intent(inout)      :: at
   real(dp), intent(out), optional :: e, local_e(:)
   real(dp), intent(out), optional :: f(:,:)
   real(dp), intent(out), optional :: virial(3,3)
   character(len=*), optional, intent(in) :: args_str

   type(Dictionary) :: params
   logical :: restart, calc_dipoles, calc_field
   real(dp), allocatable, target :: thefield(:,:)
   real(dp), pointer :: field(:,:)
   

   call initialise(params)
   call param_register(params, 'restart', 'F', restart)
   call param_register(params, 'calc_dipoles', 'F', calc_dipoles)
   call param_register(params, 'calc_field', 'F', calc_field)
   if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='IPModel_ASAP2_Calc args_str')) then
      call system_abort("IPModel_ASAP2_Calc failed to parse args_str="//trim(args_str))
   endif
   call finalise(params)

   if (present(e)) e = 0.0_dp
   if (present(f)) f = 0.0_dp
   if (present(local_e)) local_e = 0.0_dp
   if (present(virial)) virial = 0.0_dp

   if (calc_field) then
      if (.not. has_property(at, 'field')) call add_property(at, 'field', 0.0_dp, n_cols=3)
      if (.not. assign_pointer(at, 'field', field)) &
           call system_abort('IPModel_ASAP2_calc failed to assign pointer to "field" property')
   else
      allocate(thefield(3,at%n))
      field => thefield
   end if

   field = 0.0_dp
   call asap_rs_charges(this, at, e, local_e, f, virial, field, args_str)
   call asap_rs_dipoles(this, at, e, local_e, f, virial, args_str)
   call asap_morse_stretch(this, at, e, local_e, f, virial, args_str)

   ! Unit conversion
   if (present(e)) e = e*HARTREE
   if (present(local_e)) local_e = local_e*HARTREE
   if (present(f)) f = f*(HARTREE/BOHR)
   if (present(virial)) virial = virial*HARTREE

   if (allocated(thefield)) deallocate(thefield)

end subroutine IPModel_ASAP2_Calc


subroutine IPModel_ASAP2_Print(this, file)
  type(IPModel_ASAP2), intent(in) :: this
  type(Inoutput), intent(inout),optional :: file

  integer :: ti, tj

  call Print("IPModel_ASAP2 : ASAP2 Potential", file=file)
  call Print("IPModel_ASAP2 : n_types = " // this%n_types, file=file)
  call Print("IPModel_ASAP2 : betapol = "//this%betapol//" maxipol = "//this%maxipol//" tolpol = "//this%tolpol//" pred_order = "//this%pred_order, file=file)
  call Print("IPModel_ASAP2 : yukalpha = "//this%yukalpha//" yuksmoothlength = "//this%yuksmoothlength, file=file)

  do ti=1, this%n_types
    call Print ("IPModel_ASAP2 : type " // ti // " atomic_num " // this%atomic_num(ti), file=file)
    call Print ("IPModel_ASAP2 : pol = "//this%pol(ti), file=file)
    call Print ("IPModel_ASAP2 : z   = "//this%z(ti), file=file)
   call verbosity_push_decrement()
    do tj =1,this%n_types
       call Print ("IPModel_ASAP2 : pair interaction ti tj " // ti // " " // tj // " Zi Zj " // this%atomic_num(ti) //&
            " " // this%atomic_num(tj), file=file)
       call Print ("IPModel_ASAP2 : pair " // this%D_ms(ti,tj) // " " // this%gamma_ms(ti,tj) // " " &
            // this%R_ms(ti,tj) // " " // this%B_pol(ti,tj) // " " // this%C_pol(ti, tj), file=file)
    end do
   call verbosity_pop()
  end do

end subroutine IPModel_ASAP2_Print

subroutine IPModel_ASAP2_read_params_xml(this, param_str)
  type(IPModel_ASAP2), intent(inout), target :: this
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
    call system_abort("IPModel_ASAP2_read_params_xml parsed file, but n_types = 0")
  endif

end subroutine IPModel_ASAP2_read_params_xml

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

  integer ti, tj, Zi, Zj

  if (name == 'ASAP_params') then ! new ASAP2 stanza

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
        call system_abort("Can't find n_types in ASAP_params")
      endif

      allocate(parse_ip%atomic_num(parse_ip%n_types))
      parse_ip%atomic_num = 0

      allocate(parse_ip%pol(parse_ip%n_types))
      parse_ip%pol = 0.0_dp
      allocate(parse_ip%z(parse_ip%n_types))

      allocate(parse_ip%D_ms(parse_ip%n_types,parse_ip%n_types))
      parse_ip%D_ms = 0.0_dp
      allocate(parse_ip%gamma_ms(parse_ip%n_types,parse_ip%n_types))
      parse_ip%gamma_ms = 0.0_dp
      allocate(parse_ip%R_ms(parse_ip%n_types,parse_ip%n_types))
      parse_ip%R_ms = 0.0_dp
      allocate(parse_ip%B_pol(parse_ip%n_types,parse_ip%n_types))
      parse_ip%B_pol = 0.0_dp
      allocate(parse_ip%C_pol(parse_ip%n_types,parse_ip%n_types))
      parse_ip%C_pol = 0.0_dp

      call QUIP_FoX_get_value(attributes, "cutoff_coulomb", value, status)
      if (status /= 0) call system_abort ("IPModel_ASAP2_read_params_xml cannot find cutoff_coulomb")
      read (value, *) parse_ip%cutoff_coulomb

      call QUIP_FoX_get_value(attributes, "cutoff_ms", value, status)
      if (status /= 0) call system_abort ("IPModel_ASAP2_read_params_xml cannot find cutoff_ms")
      read (value, *) parse_ip%cutoff_ms

      call QUIP_FoX_get_value(attributes, "betapol", value, status)
      if (status == 0) read (value, *) parse_ip%betapol

      call QUIP_FoX_get_value(attributes, "maxipol", value, status)
      if (status == 0) read (value, *) parse_ip%maxipol

      call QUIP_FoX_get_value(attributes, "tolpol", value, status)
      if (status == 0) read (value, *) parse_ip%tolpol

      call QUIP_FoX_get_value(attributes, "pred_order", value, status)
      if (status == 0) read (value, *) parse_ip%pred_order

      call QUIP_FoX_get_value(attributes, "yukalpha", value, status)
      if (status == 0) read (value, *) parse_ip%yukalpha

      call QUIP_FoX_get_value(attributes, "yuksmoothlength", value, status)
      if (status == 0) read (value, *) parse_ip%yuksmoothlength

    endif

  elseif (parse_in_ip .and. name == 'per_type_data') then

    call QUIP_FoX_get_value(attributes, "type", value, status)
    if (status /= 0) call system_abort ("IPModel_ASAP2_read_params_xml cannot find type")
    read (value, *) ti

    call QUIP_FoX_get_value(attributes, "atomic_num", value, status)
    if (status /= 0) call system_abort ("IPModel_ASAP2_read_params_xml cannot find atomic_num")
    read (value, *) parse_ip%atomic_num(ti)

    call QUIP_FoX_get_value(attributes, "pol", value, status)
    if (status /= 0) call system_abort ("IPModel_ASAP2_read_params_xml cannot find pol")
    read (value, *) parse_ip%pol(ti)

    call QUIP_FoX_get_value(attributes, "z", value, status)
    if (status /= 0) call system_abort ("IPModel_ASAP2_read_params_xml cannot find z")
    read (value, *) parse_ip%z(ti)

    if (allocated(parse_ip%type_of_atomic_num)) deallocate(parse_ip%type_of_atomic_num)
    allocate(parse_ip%type_of_atomic_num(maxval(parse_ip%atomic_num)))
    parse_ip%type_of_atomic_num = 0
    do ti=1, parse_ip%n_types
      if (parse_ip%atomic_num(ti) > 0) &
        parse_ip%type_of_atomic_num(parse_ip%atomic_num(ti)) = ti
    end do

  elseif (parse_in_ip .and. name == 'per_pair_data') then

    call QUIP_FoX_get_value(attributes, "atnum_i", value, status)
    if (status /= 0) call system_abort ("IPModel_SW_read_params_xml cannot find atnum_i")
    read (value, *) Zi
    call QUIP_FoX_get_value(attributes, "atnum_j", value, status)
    if (status /= 0) call system_abort ("IPModel_SW_read_params_xml cannot find atnum_j")
    read (value, *) Zj

    ti = get_type(parse_ip%type_of_atomic_num,Zi)
    tj = get_type(parse_ip%type_of_atomic_num,Zj)

    call QUIP_FoX_get_value(attributes, "D_ms", value, status)
    if (status /= 0) call system_abort ("IPModel_ASAP2_read_params_xml cannot find D_ms")
    read (value, *) parse_ip%D_ms(ti,tj)
    call QUIP_FoX_get_value(attributes, "gamma_ms", value, status)
    if (status /= 0) call system_abort ("IPModel_ASAP2_read_params_xml cannot find gamma_ms")
    read (value, *) parse_ip%gamma_ms(ti,tj)
    call QUIP_FoX_get_value(attributes, "R_ms", value, status)
    if (status /= 0) call system_abort ("IPModel_ASAP2_read_params_xml cannot find R_ms")
    read (value, *) parse_ip%R_ms(ti,tj)
    call QUIP_FoX_get_value(attributes, "B_pol", value, status)
    if (status /= 0) call system_abort ("IPModel_ASAP2_read_params_xml cannot find B_pol")
    read (value, *) parse_ip%B_pol(ti,tj)
    call QUIP_FoX_get_value(attributes, "C_pol", value, status)
    if (status /= 0) call system_abort ("IPModel_ASAP2_read_params_xml cannot find C_pol")
    read (value, *) parse_ip%C_pol(ti,tj)

    if (ti /= tj) then
      parse_ip%D_ms(tj,ti) = parse_ip%D_ms(ti,tj)
      parse_ip%gamma_ms(tj,ti) = parse_ip%gamma_ms(ti,tj)
      parse_ip%R_ms(tj,ti) = parse_ip%R_ms(ti,tj)
      parse_ip%B_pol(tj,ti) = parse_ip%B_pol(ti,tj)
      parse_ip%C_pol(tj,ti) = parse_ip%C_pol(ti,tj)
    endif

  endif

end subroutine IPModel_startElement_handler

subroutine IPModel_endElement_handler(URI, localname, name)
  character(len=*), intent(in)   :: URI
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name

  if (parse_in_ip) then
    if (name == 'ASAP2_params') then
      parse_in_ip = .false.
    end if
  endif

end subroutine IPModel_endElement_handler

end module IPModel_ASAP2_module
