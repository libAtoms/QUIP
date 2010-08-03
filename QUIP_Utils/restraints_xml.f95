module restraints_xml_module
use libatoms_module
use quip_common_module
implicit none
private

type(DynamicalSystem), pointer, private :: parse_ds
logical, private :: parse_in_restraints = .false.

public :: init_restraints

contains

   subroutine restraints_startElement_handler(URI, localname, name, attributes)
      character(len=*), intent(in)   :: URI  
      character(len=*), intent(in)   :: localname
      character(len=*), intent(in)   :: name 
      type(dictionary_t), intent(in) :: attributes

      character(len=1024) :: value
      integer :: status
      integer :: n
      integer :: atom_1, atom_2, atom_3
      real(dp) :: k
      real(dp) :: c, d

      if (name == 'restraints') then ! new restraints stanze
	 parse_in_restraints = .true.
	 call QUIP_FoX_get_value(attributes, "N", value, status)
	 if (status /= 0) call system_abort("restraint_startElement_handler failed to read N in restraints stanza")
	 read (value, *) n
	 if (allocated(parse_ds%restraint)) deallocate(parse_ds%restraint)
	 allocate(parse_ds%restraint(n))
      else if (parse_in_restraints) then
	 if (name == 'bond_length') then
	    call QUIP_FoX_get_value(attributes, "atom_1", value, status)
	    if (status /= 0) call system_abort("restraint_startElement_handler failed to read atom_1 in bond_length restraint")
	    read (value, *) atom_1
	    call QUIP_FoX_get_value(attributes, "atom_2", value, status)
	    if (status /= 0) call system_abort("restraint_startElement_handler failed to read atom_2 in bond_length restraint")
	    read (value, *) atom_2
	    call QUIP_FoX_get_value(attributes, "k", value, status)
	    if (status /= 0) call system_abort("restraint_startElement_handler failed to read k in bond_length restraint")
	    read (value, *) k
	    call QUIP_FoX_get_value(attributes, "d", value, status)
	    if (status == 0) then
	       read (value, *) d
	       call constrain_bondlength(parse_ds, atom_1, atom_2, d, restraint_k=k)
	    else
	       call constrain_bondlength(parse_ds, atom_1, atom_2, restraint_k=k)
	    endif
	 else if (name == 'bond_length_sq') then
	    call QUIP_FoX_get_value(attributes, "atom_1", value, status)
	    if (status /= 0) call system_abort("restraint_startElement_handler failed to read atom_1 in bond_length_sq restraint")
	    read (value, *) atom_1
	    call QUIP_FoX_get_value(attributes, "atom_2", value, status)
	    if (status /= 0) call system_abort("restraint_startElement_handler failed to read atom_2 in bond_length_sq restraint")
	    read (value, *) atom_2
	    call QUIP_FoX_get_value(attributes, "k", value, status)
	    if (status /= 0) call system_abort("restraint_startElement_handler failed to read k in bond_length_sq restraint")
	    read (value, *) k
	    d = -1.0_dp
	    call QUIP_FoX_get_value(attributes, "d", value, status)
	    if (status == 0) read (value, *) d

	    if (d >= 0.0_dp) then
	       call constrain_bondlength_sq(parse_ds, atom_1, atom_2, d, restraint_k=k)
	    else
	       call constrain_bondlength_sq(parse_ds, atom_1, atom_2, restraint_k=k)
	    endif
	 else if (name == 'bond_angle_cos') then
	    call QUIP_FoX_get_value(attributes, "atom_1", value, status)
	    if (status /= 0) call system_abort("restraint_startElement_handler failed to read atom_1 in bond_length_sq restraint")
	    read (value, *) atom_1
	    call QUIP_FoX_get_value(attributes, "atom_2", value, status)
	    if (status /= 0) call system_abort("restraint_startElement_handler failed to read atom_2 in bond_length_sq restraint")
	    read (value, *) atom_2
	    call QUIP_FoX_get_value(attributes, "atom_3", value, status)
	    if (status /= 0) call system_abort("restraint_startElement_handler failed to read atom_3 in bond_length_sq restraint")
	    read (value, *) atom_3
	    call QUIP_FoX_get_value(attributes, "k", value, status)
	    if (status /= 0) call system_abort("restraint_startElement_handler failed to read k in bond_length_sq restraint")
	    read (value, *) k
	    call QUIP_FoX_get_value(attributes, "c", value, status)
	    if (status == 0) then
	       read (value, *) c
	       call constrain_bondanglecos(parse_ds, atom_1, atom_2, atom_3, c, restraint_k=k)
	    else
	       call constrain_bondanglecos(parse_ds, atom_1, atom_2, atom_3, restraint_k=k)
	    endif
	 else
	    call system_abort("unknown restraint type '"//name//"' in restraints xml stanza")
	 endif
      endif

   end subroutine

   subroutine restraints_endElement_handler(URI, localname, name)
      character(len=*), intent(in)   :: URI  
      character(len=*), intent(in)   :: localname
      character(len=*), intent(in)   :: name 

      if (name == 'restraints') then ! end of restraints stanza
	 parse_in_restraints = .false.
      endif

   end subroutine

   subroutine init_restraints(this, param_str)
      type(DynamicalSystem), intent(inout), target :: this
      character(len=*), intent(in) :: param_str

      type(xml_t) :: fxml

      parse_ds => this

      call open_xml_string(fxml, param_str)

      call parse(fxml, startElement_handler = restraints_startElement_handler, &
	 endElement_handler = restraints_endElement_handler)

      call close_xml_t(fxml)

   end subroutine init_restraints

end module restraints_xml_module
