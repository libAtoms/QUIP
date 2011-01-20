module restraints_constraints_xml_module
use libatoms_module
use quip_common_module
implicit none
private

type(DynamicalSystem), pointer, private :: parse_ds
logical, private :: parse_in_restraints = .false., parse_in_constraints = .false.

public :: init_restraints_constraints

contains

   subroutine restraints_constraints_startElement_handler(URI, localname, name, attributes)
      character(len=*), intent(in)   :: URI  
      character(len=*), intent(in)   :: localname
      character(len=*), intent(in)   :: name 
      type(dictionary_t), intent(in) :: attributes

      character(len=1024) :: value
      integer :: status
      integer :: n
      integer :: atom_1, atom_2, atom_3
      real(dp) :: k
      character(len=1024) :: bound_str
      integer :: bound
      real(dp) :: c
      real(dp) :: d, plane_n(3), di, df, p
      real(dp) :: t0, tau
      real(dp) :: egap
      real(dp) :: tol
      character(len=20) :: type_str

      if (parse_in_restraints) type_str = "restraint"
      if (parse_in_constraints) type_str = "constraint"

      if (name == 'restraints' .or. name == 'constraints') then ! new restraints stanze
	 if (name == 'restraints') then
	    parse_in_restraints = .true.
	    parse_in_constraints = .false.
	    call QUIP_FoX_get_value(attributes, "N", value, status)
	    if (status /= 0) call system_abort("restraint_startElement_handler failed to read N in restraints stanza")
	    read (value, *) n
	    if (allocated(parse_ds%restraint)) deallocate(parse_ds%restraint)
	    allocate(parse_ds%restraint(n))
	 else ! constraints
	    parse_in_constraints = .true.
	    parse_in_restraints = .false.
	    call QUIP_FoX_get_value(attributes, "N", value, status)
	    if (status /= 0) call system_abort("restraint_startElement_handler failed to read N in constraints stanza")
	    read (value, *) n
	    if (allocated(parse_ds%constraint)) deallocate(parse_ds%constraint)
	    allocate(parse_ds%constraint(n))
	 endif

      else if (parse_in_restraints .or. parse_in_constraints) then

	 if (name == 'bond_length') then
	    call QUIP_FoX_get_value(attributes, "atom_1", value, status)
	    if (status /= 0) call system_abort("restraint_startElement_handler failed to read atom_1 in bond_length "//trim(type_str))
	    read (value, *) atom_1
	    call QUIP_FoX_get_value(attributes, "atom_2", value, status)
	    if (status /= 0) call system_abort("restraint_startElement_handler failed to read atom_2 in bond_length "//trim(type_str))
	    read (value, *) atom_2
	    if (parse_in_restraints) then
	       call QUIP_FoX_get_value(attributes, "k", value, status)
	       if (status /= 0) call system_abort("restraint_startElement_handler failed to read k in bond_length "//trim(type_str))
	       read (value, *) k
	       call QUIP_FoX_get_value(attributes, "bound", value, status)
	       if (status == 0) then
	          read (value, *) bound_str
                  call interpret_bound_string(trim(bound_str),bound)
               else
                  bound = BOTH_UPPER_AND_LOWER_BOUNDS
               endif
            else
               call QUIP_FoX_get_value(attributes, "tol", value, status) !constraint tolerance
	       if (status == 0) then
	          read (value, *) tol
               else
                  tol=-1._dp !default
               endif
	    endif
	    call QUIP_FoX_get_value(attributes, "d", value, status)
	    if (status == 0) then
	       read (value, *) d
	       if (parse_in_restraints) then
		  call constrain_bondlength(parse_ds, atom_1, atom_2, d, restraint_k=k, bound=bound)
	       else
		  call constrain_bondlength(parse_ds, atom_1, atom_2, d, tol=tol)
	       endif
	    else
	       if (parse_in_restraints) then
		  call constrain_bondlength(parse_ds, atom_1, atom_2, restraint_k=k, bound=bound)
	       else
		  call constrain_bondlength(parse_ds, atom_1, atom_2, tol=tol)
	       endif
	    endif

	 else if (name == 'relax_bond_length') then
	    call QUIP_FoX_get_value(attributes, "atom_1", value, status)
	    if (status /= 0) call system_abort("restraint_startElement_handler failed to read atom_1 in relax_bond_length "// &
	      trim(type_str))
	    read (value, *) atom_1
	    call QUIP_FoX_get_value(attributes, "atom_2", value, status)
	    if (status /= 0) call system_abort("restraint_startElement_handler failed to read atom_2 in relax_bond_length "// &
	      trim(type_str))
	    read (value, *) atom_2
	    if (parse_in_restraints) then
	       call QUIP_FoX_get_value(attributes, "k", value, status)
	       if (status /= 0) call system_abort("restraint_startElement_handler failed to read k in relax_bond_length "// &
		trim(type_str))
	       read (value, *) k
	       call QUIP_FoX_get_value(attributes, "bound", value, status)
	       if (status == 0) then
	          read (value, *) bound_str
                  call interpret_bound_string(trim(bound_str),bound)
               else
                  bound = BOTH_UPPER_AND_LOWER_BOUNDS
               endif
            else
               call QUIP_FoX_get_value(attributes, "tol", value, status) !constraint tolerance
	       if (status == 0) then
	          read (value, *) tol
               else
                  tol=-1._dp !default
               endif
	    endif
	    call QUIP_FoX_get_value(attributes, "t0", value, status)
	    if (status /= 0) then
	       call print("restraint_startElement_handler failed to read t0 in relax_bond_length "//trim(type_str)// &
			  ", defaulting to 0", PRINT_ALWAYS)
	       t0 = 0.0_dp
	    else
	       read (value, *) t0
	    endif
	    call QUIP_FoX_get_value(attributes, "tau", value, status)
	    if (status /= 0) call system_abort("restraint_startElement_handler failed to read tau in relax_bond_length "//trim(type_str))
	    read (value, *) tau
	    call QUIP_FoX_get_value(attributes, "df", value, status)
	    if (status /= 0) call system_abort("restraint_startElement_handler failed to read df in relax_bond_length "//trim(type_str))
	    read (value, *) df
	    call QUIP_FoX_get_value(attributes, "di", value, status)
	    if (status == 0) then
	       read (value, *) di
	       if (parse_in_restraints) then
		  call constrain_relax_bondlength(parse_ds, atom_1, atom_2, t0, tau, df, di, restraint_k=k, bound=bound)
	       else
		  call constrain_relax_bondlength(parse_ds, atom_1, atom_2, t0, tau, df, di, tol=tol)
	       endif
	    else
	       if (parse_in_restraints) then
		  call constrain_relax_bondlength(parse_ds, atom_1, atom_2, t0, tau, df, restraint_k=k, bound=bound)
	       else
		  call constrain_relax_bondlength(parse_ds, atom_1, atom_2, t0, tau, df, tol=tol)
	       endif
	    endif

	 else if (name == 'relax_bond_length_dev_pow') then
	    call QUIP_FoX_get_value(attributes, "atom_1", value, status)
	    if (status /= 0) call system_abort("restraint_startElement_handler failed to read atom_1 in relax_bond_length_dev_pow "// &
	      trim(type_str))
	    read (value, *) atom_1
	    call QUIP_FoX_get_value(attributes, "atom_2", value, status)
	    if (status /= 0) call system_abort("restraint_startElement_handler failed to read atom_2 in relax_bond_length_dev_pow "// &
	      trim(type_str))
	    read (value, *) atom_2
	    if (parse_in_restraints) then
	       call QUIP_FoX_get_value(attributes, "k", value, status)
	       if (status /= 0) call system_abort("restraint_startElement_handler failed to read k in relax_bond_length_dev_pow "// &
		trim(type_str))
	       read (value, *) k
	       call QUIP_FoX_get_value(attributes, "bound", value, status)
	       if (status == 0) then
	          read (value, *) bound_str
                  call interpret_bound_string(trim(bound_str),bound)
               else
                  bound = BOTH_UPPER_AND_LOWER_BOUNDS
               endif
            else
               call QUIP_FoX_get_value(attributes, "tol", value, status) !constraint tolerance
	       if (status == 0) then
	          read (value, *) tol
               else
                  tol=-1._dp !default
               endif
	    endif
	    call QUIP_FoX_get_value(attributes, "p", value, status)
	    if (status /= 0) call system_abort("restraint_startElement_handler failed to read p in relax_bond_length_dev_pow "//trim(type_str))
	    read (value, *) p
	    call QUIP_FoX_get_value(attributes, "t0", value, status)
	    if (status /= 0) then
	       call print("restraint_startElement_handler failed to read t0 in relax_bond_length_dev_pow "//trim(type_str)// &
			  ", defaulting to 0", PRINT_ALWAYS)
	       t0 = 0.0_dp
	    else
	       read (value, *) t0
	    endif
	    call QUIP_FoX_get_value(attributes, "tau", value, status)
	    if (status /= 0) call system_abort("restraint_startElement_handler failed to read tau in relax_bond_length_dev_pow "//trim(type_str))
	    read (value, *) tau
	    call QUIP_FoX_get_value(attributes, "df", value, status)
	    if (status /= 0) call system_abort("restraint_startElement_handler failed to read df in relax_bond_length_dev_pow "//trim(type_str))
	    read (value, *) df
	    call QUIP_FoX_get_value(attributes, "di", value, status)
	    if (status == 0) then
	       read (value, *) di
	       if (parse_in_restraints) then
		  call constrain_relax_bondlength_dev_pow(parse_ds, atom_1, atom_2, t0, tau, p, df, di, restraint_k=k, bound=bound)
	       else
		  call constrain_relax_bondlength_dev_pow(parse_ds, atom_1, atom_2, t0, tau, p, df, di, tol=tol)
	       endif
	    else
	       if (parse_in_restraints) then
		  call constrain_relax_bondlength_dev_pow(parse_ds, atom_1, atom_2, t0, tau, p, df, restraint_k=k, bound=bound)
	       else
		  call constrain_relax_bondlength_dev_pow(parse_ds, atom_1, atom_2, t0, tau, p, df, tol=tol)
	       endif
	    endif

	 else if (name == 'bond_length_sq') then

	    call QUIP_FoX_get_value(attributes, "atom_1", value, status)
	    if (status /= 0) call system_abort("restraint_startElement_handler failed to read atom_1 in bond_length_sq "//trim(type_str))
	    read (value, *) atom_1
	    call QUIP_FoX_get_value(attributes, "atom_2", value, status)
	    if (status /= 0) call system_abort("restraint_startElement_handler failed to read atom_2 in bond_length_sq "//trim(type_str))
	    read (value, *) atom_2
	    if (parse_in_restraints) then
	       call QUIP_FoX_get_value(attributes, "k", value, status)
	       if (status /= 0) call system_abort("restraint_startElement_handler failed to read k in bond_length_sq "//trim(type_str))
	       read (value, *) k
	       call QUIP_FoX_get_value(attributes, "bound", value, status)
	       if (status == 0) then
	          read (value, *) bound_str
                  call interpret_bound_string(trim(bound_str),bound)
               else
                  bound = BOTH_UPPER_AND_LOWER_BOUNDS
               endif
            else
               call QUIP_FoX_get_value(attributes, "tol", value, status) !constraint tolerance
	       if (status == 0) then
	          read (value, *) tol
               else
                  tol=-1._dp !default
               endif
	    endif
	    call QUIP_FoX_get_value(attributes, "d", value, status)
	    if (status == 0) then
	       read (value, *) d
	       if (parse_in_restraints) then
		  call constrain_bondlength_sq(parse_ds, atom_1, atom_2, d, restraint_k=k, bound=bound)
	       else
		  call constrain_bondlength_sq(parse_ds, atom_1, atom_2, d, tol=tol)
	       endif
	    else
	       if (parse_in_restraints) then
		  call constrain_bondlength_sq(parse_ds, atom_1, atom_2, restraint_k=k, bound=bound)
	       else
		  call constrain_bondlength_sq(parse_ds, atom_1, atom_2, tol=tol)
	       endif
	    endif

	 else if (name == 'bond_length_dev_pow') then

	    call QUIP_FoX_get_value(attributes, "atom_1", value, status)
	    if (status /= 0) call system_abort("restraint_startElement_handler failed to read atom_1 in bond_length_dev_pow "//trim(type_str))
	    read (value, *) atom_1
	    call QUIP_FoX_get_value(attributes, "atom_2", value, status)
	    if (status /= 0) call system_abort("restraint_startElement_handler failed to read atom_2 in bond_length_dev_pow "//trim(type_str))
	    read (value, *) atom_2
	    if (parse_in_restraints) then
	       call QUIP_FoX_get_value(attributes, "k", value, status)
	       if (status /= 0) call system_abort("restraint_startElement_handler failed to read k in bond_length_dev_pow "//trim(type_str))
	       read (value, *) k
	       call QUIP_FoX_get_value(attributes, "bound", value, status)
	       if (status == 0) then
	          read (value, *) bound_str
                  call interpret_bound_string(trim(bound_str),bound)
               else
                  bound = BOTH_UPPER_AND_LOWER_BOUNDS
               endif
            else
               call QUIP_FoX_get_value(attributes, "tol", value, status) !constraint tolerance
	       if (status == 0) then
	          read (value, *) tol
               else
                  tol=-1._dp !default
               endif
	    endif
	    call QUIP_FoX_get_value(attributes, "p", value, status)
	    if (status /= 0) call system_abort("restraint_startElement_handler failed to read p in bond_length_dev_pow "//trim(type_str))
	    read (value, *) p
	    call QUIP_FoX_get_value(attributes, "d", value, status)
	    if (status == 0) then
	       read (value, *) d
	       if (parse_in_restraints) then
		  call constrain_bondlength_dev_pow(parse_ds, atom_1, atom_2, p, d, restraint_k=k, bound=bound)
	       else
		  call constrain_bondlength_dev_pow(parse_ds, atom_1, atom_2, p, d, tol=tol)
	       endif
	    else
	       if (parse_in_restraints) then
		  call constrain_bondlength_dev_pow(parse_ds, atom_1, atom_2, p, restraint_k=k, bound=bound)
	       else
		  call constrain_bondlength_dev_pow(parse_ds, atom_1, atom_2, p, tol=tol)
	       endif
	    endif

	 else if (name == 'bond_angle_cos') then

	    call QUIP_FoX_get_value(attributes, "atom_1", value, status)
	    if (status /= 0) call system_abort("restraint_startElement_handler failed to read atom_1 in bond_angle_cos "//trim(type_str))
	    read (value, *) atom_1
	    call QUIP_FoX_get_value(attributes, "atom_2", value, status)
	    if (status /= 0) call system_abort("restraint_startElement_handler failed to read atom_2 in bond_angle_cos "//trim(type_str))
	    read (value, *) atom_2
	    call QUIP_FoX_get_value(attributes, "atom_3", value, status)
	    if (status /= 0) call system_abort("restraint_startElement_handler failed to read atom_3 in bond_angle_cos "//trim(type_str))
	    read (value, *) atom_3
	    if (parse_in_restraints) then
	       call QUIP_FoX_get_value(attributes, "k", value, status)
	       if (status /= 0) call system_abort("restraint_startElement_handler failed to read k in bond_angle_cos "//trim(type_str))
	       read (value, *) k
	       call QUIP_FoX_get_value(attributes, "bound", value, status)
	       if (status == 0) then
	          read (value, *) bound_str
                  call interpret_bound_string(trim(bound_str),bound)
               else
                  bound = BOTH_UPPER_AND_LOWER_BOUNDS
               endif
            else
               call QUIP_FoX_get_value(attributes, "tol", value, status) !constraint tolerance
	       if (status == 0) then
	          read (value, *) tol
               else
                  tol=-1._dp !default
               endif
	    endif
	    call QUIP_FoX_get_value(attributes, "c", value, status)
	    if (status == 0) then
	       read (value, *) c
	       if (parse_in_restraints) then
		  call constrain_bondanglecos(parse_ds, atom_1, atom_2, atom_3, c, restraint_k=k, bound=bound)
	       else
		  call constrain_bondanglecos(parse_ds, atom_1, atom_2, atom_3, c, tol=tol)
	       endif
	    else
	       if (parse_in_restraints) then
		  call constrain_bondanglecos(parse_ds, atom_1, atom_2, atom_3, restraint_k=k, bound=bound)
	       else
		  call constrain_bondanglecos(parse_ds, atom_1, atom_2, atom_3, tol=tol)
	       endif
	    endif

	 else if (name == 'bond_length_diff') then

	    call QUIP_FoX_get_value(attributes, "atom_1", value, status)
	    if (status /= 0) call system_abort("restraint_startElement_handler failed to read atom_1 in bond_length_diff "//trim(type_str))
	    read (value, *) atom_1
	    call QUIP_FoX_get_value(attributes, "atom_2", value, status)
	    if (status /= 0) call system_abort("restraint_startElement_handler failed to read atom_2 in bond_length_diff "//trim(type_str))
	    read (value, *) atom_2
	    call QUIP_FoX_get_value(attributes, "atom_3", value, status)
	    if (status /= 0) call system_abort("restraint_startElement_handler failed to read atom_3 in bond_length_diff "//trim(type_str))
	    read (value, *) atom_3
	    if (parse_in_restraints) then
	       call QUIP_FoX_get_value(attributes, "k", value, status)
	       if (status /= 0) call system_abort("restraint_startElement_handler failed to read k in bond_length_diff "//trim(type_str))
	       read (value, *) k
	       call QUIP_FoX_get_value(attributes, "bound", value, status)
	       if (status == 0) then
	          read (value, *) bound_str
                  call interpret_bound_string(trim(bound_str),bound)
               else
                  bound = BOTH_UPPER_AND_LOWER_BOUNDS
               endif
            else
               call QUIP_FoX_get_value(attributes, "tol", value, status) !constraint tolerance
	       if (status == 0) then
	          read (value, *) tol
               else
                  tol=-1._dp !default
               endif
	    endif
	    call QUIP_FoX_get_value(attributes, "d", value, status)
	    if (status == 0) then
	       read (value, *) d
	       if (parse_in_restraints) then
		  call constrain_bondlength_diff(parse_ds, atom_1, atom_2, atom_3, d, restraint_k=k, bound=bound)
	       else
		  call constrain_bondlength_diff(parse_ds, atom_1, atom_2, atom_3, d, tol=tol)
	       endif
	    else
	       if (parse_in_restraints) then
		  call constrain_bondlength_diff(parse_ds, atom_1, atom_2, atom_3, restraint_k=k, bound=bound)
	       else
		  call constrain_bondlength_diff(parse_ds, atom_1, atom_2, atom_3, tol=tol)
	       endif
	    endif

	 else if (name == 'gap_energy') then

	    if (parse_in_restraints) then
	       call QUIP_FoX_get_value(attributes, "k", value, status)
	       if (status /= 0) call system_abort("restraint_startElement_handler failed to read k in gap_energy "//trim(type_str))
	       read (value, *) k
	       call QUIP_FoX_get_value(attributes, "bound", value, status)
	       if (status == 0) then
	          read (value, *) bound_str
                  call interpret_bound_string(trim(bound_str),bound)
               else
                  bound = BOTH_UPPER_AND_LOWER_BOUNDS
               endif
            else
               call QUIP_FoX_get_value(attributes, "tol", value, status) !constraint tolerance
	       if (status == 0) then
	          read (value, *) tol
               else
                  tol=-1._dp !default
               endif
            endif
	    call QUIP_FoX_get_value(attributes, "egap", value, status)
	    if (status == 0) then
	       read (value, *) egap
	       if (parse_in_restraints) then
		  call constrain_gap_energy(parse_ds, egap, restraint_k=k, bound=bound)
	       else
		  call constrain_gap_energy(parse_ds, egap, tol=tol)
	       endif
	    else
	       call system_abort("restraint_startElement_handler failed to read egap in gap_energy "//trim(type_str))
	       !if (parse_in_restraints) then
	       !  call constrain_gap_energy(parse_ds, restraint_k=k)
	       !else
	       !  call constrain_gap_energy(parse_ds, tol=tol)
	       !endif
	    endif

	 else if (name == 'atom_plane') then

	    call QUIP_FoX_get_value(attributes, "atom", value, status)
	    if (status /= 0) call system_abort("restraint_startElement_handler failed to read atom in atom_plane "//trim(type_str))
	    read (value, *) atom_1
	    call QUIP_FoX_get_value(attributes, "normal", value, status)
	    if (status /= 0) call system_abort("restraint_startElement_handler failed to read normal in atom_plane "//trim(type_str))
	    read (value, *) plane_n
	    if (parse_in_restraints) then
	       call QUIP_FoX_get_value(attributes, "k", value, status)
	       if (status /= 0) call system_abort("restraint_startElement_handler failed to read k in atom_plane "//trim(type_str))
	       read (value, *) k
	       call QUIP_FoX_get_value(attributes, "bound", value, status)
	       if (status == 0) then
	          read (value, *) bound_str
                  call interpret_bound_string(trim(bound_str),bound)
               else
                  bound = BOTH_UPPER_AND_LOWER_BOUNDS
               endif
            else
               call QUIP_FoX_get_value(attributes, "tol", value, status) !constraint tolerance
	       if (status == 0) then
	          read (value, *) tol
               else
                  tol=-1._dp !default
               endif
	    endif
	    call QUIP_FoX_get_value(attributes, "d", value, status)
	    if (status == 0) then
	       read (value, *) d
	       if (parse_in_restraints) then
		  call constrain_atom_plane(parse_ds, atom_1, plane_n, d, restraint_k=k, bound=bound)
	       else
		  call constrain_atom_plane(parse_ds, atom_1, plane_n, d, tol=tol)
	       endif
	    else
	       if (parse_in_restraints) then
		  call constrain_atom_plane(parse_ds, atom_1, plane_n, restraint_k=k, bound=bound)
	       else
		  call constrain_atom_plane(parse_ds, atom_1, plane_n, tol=tol)
	       endif
	    endif
	 else
	    call system_abort("unknown " // trim(type_str) //" type '"//name//"' in " // trim(type_str) //"s xml stanza")
	 endif
      endif

   end subroutine

   subroutine restraints_constraints_endElement_handler(URI, localname, name)
      character(len=*), intent(in)   :: URI  
      character(len=*), intent(in)   :: localname
      character(len=*), intent(in)   :: name 

      if (name == 'restraints') then ! end of restraints stanza
	 parse_in_restraints = .false.
      else if (name == 'constraints') then ! end of constraints stanza
	 parse_in_constraints = .false.
      endif

   end subroutine

   !% Initialise restraints and constraints for a dynamical system from an XML containing string.
   !%
   !% XML format for constraints
   !%
   !% \begin{verbatim}
   !% <constraints N="n">
   !%   <bond_length atom_1="i" atom_2="j" d="r" />
   !%   <bond_length_sq atom_1="i" atom_2="j" d="r" />
   !%   <bond_length_diff atom_1="i" atom_2="j" atom_3="k" d="r" />
   !%   <bond_angle_cos atom_1="i" atom_2="j" atom_3="k" c="r" />
   !% </constraints N="n">
   !% \end{verbatim}
   !%
   !% In all types, constraint value (c or d) is optional.  If not 
   !% specified, defaults to value at initialization time, just like
   !% underlying constraint adding routines.
   !%
   !% Format for restraints is similar, except surrounding stanza name
   !% is `restraints', and the additional {\tt k} attribute
   !% specifying spring constant for each restraint
   !
   subroutine init_restraints_constraints(this, param_str)
      type(DynamicalSystem), intent(inout), target :: this
      character(len=*), intent(in) :: param_str

      type(xml_t) :: fxml

      parse_ds => this

      call open_xml_string(fxml, param_str)

      call parse(fxml, startElement_handler = restraints_constraints_startElement_handler, &
	 endElement_handler = restraints_constraints_endElement_handler)

      call close_xml_t(fxml)

   end subroutine init_restraints_constraints

   subroutine interpret_bound_string(bound_str,bound)
     character(len=*), intent(in) :: bound_str
     integer, intent(out) :: bound

     if ( trim(bound_str) == trim(BOUND_STRING(UPPER_BOUND)) ) then
        bound = UPPER_BOUND
     elseif ( trim(bound_str) == trim(BOUND_STRING(LOWER_BOUND)) ) then
        bound = LOWER_BOUND
     elseif ( trim(bound_str) == trim(BOUND_STRING(BOTH_UPPER_AND_LOWER_BOUNDS)) ) then
        bound = BOTH_UPPER_AND_LOWER_BOUNDS
     else
        call system_abort("bound must take a value of "//trim(BOUND_STRING(UPPER_BOUND))//", "// &
           trim(BOUND_STRING(LOWER_BOUND))//" and "//trim(BOUND_STRING(BOTH_UPPER_AND_LOWER_BOUNDS)))
     endif

   end subroutine interpret_bound_string

end module restraints_constraints_xml_module
