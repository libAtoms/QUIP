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
!X IPModel_BornMayer
!X
!% BornMayer potential of the form V_{ij} = A_{ij} exp(-b_{ij} r_{ij}) - c_{ij}/r_{ij}^6
!% Use together with an IPModel_Coulomb to include long range forces
!%
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#include "error.inc"

module IPModel_BornMayer_module

  use libatoms_module

  use mpi_context_module
  use QUIP_Common_module

  implicit none
  private

  include 'IPModel_interface.h'

  public :: IPModel_BornMayer
  type IPModel_BornMayer
     integer :: n_types = 0
     integer, allocatable :: atomic_num(:), type_of_atomic_num(:)
     real(dp), allocatable ::  A(:,:), b(:,:), c(:,:), cutoff_a(:,:), &
          energy_shift(:,:), linear_force_shift(:,:)

     real(dp) :: cutoff = 0.0_dp

     character(len=FIELD_LENGTH) :: label

  end type IPModel_BornMayer

  logical, private :: parse_in_ip, parse_matched_label
  type(IPModel_BornMayer), private, pointer :: parse_ip

  interface Initialise
     module procedure IPModel_BornMayer_Initialise_str
  end interface Initialise

  interface Finalise
     module procedure IPModel_BornMayer_Finalise
  end interface Finalise

  interface Print
     module procedure IPModel_BornMayer_Print
  end interface Print

  interface Calc
     module procedure IPModel_BornMayer_Calc
  end interface Calc

contains

  subroutine IPModel_BornMayer_Initialise_str(this, args_str, param_str)
    type(IPModel_BornMayer), intent(inout) :: this
    character(len=*), intent(in) :: args_str, param_str

    type(Dictionary) :: params

    call Finalise(this)

    call initialise(params)
    this%label=''
    call param_register(params, 'label', '', this%label, help_string="No help yet.  This source file was $LastChangedBy$")
    if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='IPModel_BornMayer_Initialise_str args_str')) then
       call system_abort("IPModel_BornMayer_Initialise_str failed to parse label from args_str="//trim(args_str))
    endif
    call finalise(params)

    call IPModel_BornMayer_read_params_xml(this, param_str)

    this%cutoff = maxval(this%cutoff_a)

  end subroutine IPModel_BornMayer_Initialise_str

  subroutine IPModel_BornMayer_Finalise(this)
    type(IPModel_BornMayer), intent(inout) :: this

    ! Add finalisation code here

    if (allocated(this%atomic_num)) deallocate(this%atomic_num)
    if (allocated(this%type_of_atomic_num)) deallocate(this%type_of_atomic_num)

    if (allocated(this%A)) deallocate(this%A)
    if (allocated(this%b)) deallocate(this%b)
    if (allocated(this%c)) deallocate(this%c)
    if (allocated(this%cutoff_a)) deallocate(this%cutoff_a)
    if (allocated(this%energy_shift)) deallocate(this%energy_shift)
    if (allocated(this%linear_force_shift)) deallocate(this%linear_force_shift)

    this%n_types = 0
    this%label = ''
  end subroutine IPModel_BornMayer_Finalise


  subroutine IPModel_BornMayer_Calc(this, at, e, local_e, f, virial, local_virial, args_str, mpi, error)
    type(IPModel_BornMayer), intent(inout):: this
    type(Atoms), intent(inout)      :: at
    real(dp), intent(out), optional :: e, local_e(:)
    real(dp), intent(out), optional :: f(:,:), local_virial(:,:)   !% Forces, dimensioned as \texttt{f(3,at%N)}, local virials, dimensioned as \texttt{local_virial(9,at%N)} 
    real(dp), intent(out), optional :: virial(3,3)
    character(len=*), optional      :: args_str
    type(MPI_Context), intent(in), optional :: mpi
    integer, intent(out), optional :: error

    type(Dictionary) :: params
    logical :: has_atom_mask_name
    character(FIELD_LENGTH) :: atom_mask_name
    logical, dimension(:), pointer :: atom_mask_pointer

    real(dp) :: r_scale, E_scale
    logical :: do_rescale_r, do_rescale_E
    integer i, ji, j, ti, tj, n_neigh_i
    real(dp) :: drij(3), drij_mag
    real(dp) :: drij_dri(3), drij_drj(3)
    real(dp) :: virial_i(3,3)
    real(dp) :: de, de_dr, de_drij


#ifdef _OPENMP
    real(dp) :: private_virial(3,3), private_e
    real(dp), allocatable :: private_f(:,:), private_local_e(:), private_local_virial(:,:)
#endif

    INIT_ERROR(error)

    if (present(e)) e = 0.0_dp
    if (present(local_e)) then
       call check_size('Local_E',local_e,(/at%N/),'IPModel_BornMayer_Calc', error)
       local_e = 0.0_dp
    endif
    if (present(f)) then
       call check_size('Force',f,(/3,at%Nbuffer/),'IPModel_BornMayer_Calc', error)
       f = 0.0_dp
    end if
    if (present(virial)) virial = 0.0_dp
    if (present(local_virial)) then
       call check_size('Local_virial',local_virial,(/9,at%Nbuffer/),'IPModel_BornMayer_Calc', error)
       local_virial = 0.0_dp
    endif

    atom_mask_pointer => null()
    if (present(args_str)) then
       call initialise(params)
       call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, help_string="Name of property containing logical mask for calculation")
       call param_register(params, 'r_scale', '1.0',r_scale, has_value_target=do_rescale_r, help_string="Recaling factor for distances. Default 1.0.")
       call param_register(params, 'E_scale', '1.0',E_scale, has_value_target=do_rescale_E, help_string="Recaling factor for energy. Default 1.0.")
       if(.not. param_read_line(params, args_str, ignore_unknown=.true.,task='IPModel_BornMayer_Calc args_str')) then
          RAISE_ERROR("IPModel_BornMayer_Calc failed to parse args_str='"//trim(args_str)//"'",error)
       endif
       call finalise(params)
       if(has_atom_mask_name) then
          if (.not. assign_pointer(at, trim(atom_mask_name) , atom_mask_pointer)) &
               call system_abort("IPModel_BornMayer_Calc did not find "//trim(atom_mask_name)//" property in the atoms object.")
       else
          atom_mask_pointer => null()
       endif
    else
       do_rescale_r = .false.
       do_rescale_E = .false.
    endif

    if (do_rescale_r) call print('IPModel_BornMayer_Calc: rescaling distances by factor '//r_scale, PRINT_VERBOSE)
    if (do_rescale_E) call print('IPModel_BornMayer_Calc: rescaling energy by factor '//E_scale, PRINT_VERBOSE)


#ifdef _OPENMP
    !$omp parallel default(none) private(i, ji, j, drij, drij_mag, ti, tj, drij_dri, drij_drj, de, de_dr, de_drij, private_virial, private_e, private_f, private_local_e, private_local_virial, n_neigh_i, virial_i) shared(this, at, e, virial, atom_mask_pointer, mpi, do_rescale_r, r_scale, f, local_e, local_virial)

    if (present(e)) private_e = 0.0_dp
    if (present(local_e)) then
       allocate(private_local_e(at%N))
       private_local_e = 0.0_dp
    endif
    if (present(f)) then
       allocate(private_f(3,at%N))
       private_f = 0.0_dp
    endif
    if (present(virial)) private_virial = 0.0_dp
    if (present(local_virial)) then
       allocate(private_local_virial(9,at%N))
       private_local_virial = 0.0_dp
    endif

    !$omp do
#endif
    do i=1, at%N
       if (present(mpi)) then
          if (mpi%active) then
             if (mod(i-1, mpi%n_procs) /= mpi%my_proc) cycle
          endif
       endif

       if(associated(atom_mask_pointer)) then
          if(.not. atom_mask_pointer(i)) cycle
       endif

       ti = get_type(this%type_of_atomic_num, at%Z(i))
       if (current_verbosity() >= PRINT_ANAL) call print ("IPModel_BornMayer_Calc i " // i // " " // atoms_n_neighbours(at,i), PRINT_ANAL)
       n_neigh_i = atoms_n_neighbours(at, i)
       do ji=1, n_neigh_i
          j = atoms_neighbour(at, i, ji, drij_mag, cosines=drij, max_dist=this%cutoff)
          if (j <= 0) cycle
          if (drij_mag .feq. 0.0_dp) cycle

          if (do_rescale_r) drij_mag = drij_mag*r_scale

          tj = get_type(this%type_of_atomic_num, at%Z(j))

          if (current_verbosity() >= PRINT_ANAL) call print ("IPModel_BornMayer_Calc i j " // i // " " // j, PRINT_ANAL)

          if (present(e) .or. present(local_e)) then

             de = IPModel_BornMayer_pairenergy(this, ti, tj, drij_mag)

             if (present(local_e)) then
#ifdef _OPENMP
                private_local_e(i) = private_local_e(i) + de
#else
                local_e(i) = local_e(i) + de
#endif
             endif
             if (present(e)) then
#ifdef _OPENMP
                private_e = private_e + de
#else
                e = e + de
#endif
             endif
          endif

          if (present(f) .or. present(virial) .or. present(local_virial)) then

             de_dr = IPModel_BornMayer_pairenergy_deriv(this, ti, tj, drij_mag)

             if (present(f)) then
#ifdef _OPENMP
                private_f(:,i) = private_f(:,i) + de_dr*drij
                private_f(:,j) = private_f(:,j) - de_dr*drij
#else
                f(:,i) = f(:,i) + de_dr*drij
                f(:,j) = f(:,j) - de_dr*drij
#endif
             endif

             if(present(virial) .or. present(local_virial)) virial_i = de_dr*(drij .outer. drij)*drij_mag

             if (present(virial)) then
#ifdef _OPENMP
                private_virial = private_virial - virial_i
#else
                virial = virial - virial_i
#endif
             endif
             if (present(local_virial)) then
#ifdef _OPENMP
                private_local_virial(:,i) = private_local_virial(:,i) - reshape(virial_i, (/9/))
#else
                local_virial(:,i) = local_virial(:,i) - reshape(virial_i, (/9/))
#endif
             endif

          endif
       end do
    end do

#ifdef _OPENMP
    !$omp critical
    if (present(e)) e = e + private_e
    if (present(f)) f = f + private_f
    if (present(local_e)) local_e = local_e + private_local_e
    if (present(virial)) virial = virial + private_virial
    if (present(local_virial)) local_virial = local_virial + private_local_virial
    !$omp end critical 

    if(allocated(private_f)) deallocate(private_f)
    if(allocated(private_local_e)) deallocate(private_local_e)
    if(allocated(private_local_virial)) deallocate(private_local_virial)

    !$omp end parallel
#endif

    if (present(mpi)) then
       if (present(e)) e = sum(mpi, e) 
       if (present(local_e)) call sum_in_place(mpi, local_e)
       if (present(f)) call sum_in_place(mpi, f)
       if (present(virial)) call sum_in_place(mpi, virial)
       if (present(local_virial)) call sum_in_place(mpi, local_virial)
    endif

    if (do_rescale_r) then
       if (present(f)) f = f*r_scale
    end if

    if (do_rescale_E) then
       if (present(e)) e = e*E_scale
       if (present(local_e)) local_e = local_e*E_scale
       if (present(f)) f = f*E_scale
       if (present(virial)) virial=virial*E_scale
       if (present(local_virial)) local_virial=local_virial*E_scale
    end if

    atom_mask_pointer => null()

  end subroutine IPModel_BornMayer_Calc

  function IPModel_BornMayer_pairenergy(this, ti, tj, r)
    type(IPModel_BornMayer), intent(in) :: this
    integer, intent(in) :: ti, tj   !% Atomic types.
    real(dp), intent(in) :: r       !% Distance.
    real(dp) :: IPModel_BornMayer_pairenergy

    if ((r .feq. 0.0_dp) .or. (r > this%cutoff_a(ti,tj))) then
       IPModel_BornMayer_pairenergy = 0.0_dp
       return
    endif

    IPModel_BornMayer_pairenergy = 0.5_dp*((this%A(ti,tj)*exp(-this%b(ti,tj)*r) - this%c(ti,tj)/(r**6.0_dp)) &
         - this%energy_shift(ti,tj) - this%linear_force_shift(ti,tj)*(r-this%cutoff_a(ti,tj)))

  end function IPModel_BornMayer_pairenergy

  function IPModel_BornMayer_pairenergy_deriv(this, ti, tj, r)
    type(IPModel_BornMayer), intent(in) :: this
    integer, intent(in) :: ti, tj   !% Atomic types.
    real(dp), intent(in) :: r       !% Distance.
    real(dp) :: IPModel_BornMayer_pairenergy_deriv

    if ((r .feq. 0.0_dp) .or. (r > this%cutoff_a(ti,tj))) then
       IPModel_BornMayer_pairenergy_deriv = 0.0_dp
       return
    endif

    IPModel_BornMayer_pairenergy_deriv = 0.5_dp*(-this%b(ti,tj)*this%A(ti,tj)*exp(-this%b(ti,tj)*r) &
         + 6.0_dp*this%c(ti,tj)/(r**7.0_dp) - this%linear_force_shift(ti,tj))

  end function IPModel_BornMayer_pairenergy_deriv


  subroutine IPModel_BornMayer_Print(this, file)
    type(IPModel_BornMayer), intent(in) :: this
    type(Inoutput), intent(inout),optional :: file

    integer :: ti, tj

    call Print("IPModel_BornMayer : BornMayer Potential", file=file)
    call Print("IPModel_BornMayer : n_types = " // this%n_types, file=file)
    do ti=1, this%n_types
       call Print ("IPModel_BornMayer : type " // ti // " atomic_num " // this%atomic_num(ti), file=file)
       do tj=1, this%n_types
          call Print ("IPModel_BornMayer : pair interaction ti tj " // ti // " " // tj // " Zi Zj " // this%atomic_num(ti) // &
               " " // this%atomic_num(tj), file=file)
          call Print ("IPModel_BornMayer : pair A=" // this%A(ti,tj) // " b=" // this%b(ti,tj) // " c=" //this%c(ti,tj), file=file)
          call Print ("                         cutoff="//this%cutoff_a(ti,tj)//" energy_shift="//this%energy_shift(ti,tj), file=file)
          call Print ("                         linear_force_shift="//this%linear_force_shift(ti,tj), file=file)
       end do
    end do

  end subroutine IPModel_BornMayer_Print

  subroutine IPModel_BornMayer_read_params_xml(this, param_str)
    type(IPModel_BornMayer), intent(inout), target :: this
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
       call system_abort("IPModel_BornMayer_read_params_xml parsed file, but n_types = 0")
    endif

  end subroutine IPModel_BornMayer_read_params_xml

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
    logical :: energy_shift, linear_force_shift

    integer :: status
    character(len=FIELD_LENGTH) :: value

    integer ti, tj, Zi, Zj

    if (name == 'BornMayer_params') then ! new BornMayer stanza

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
             call system_abort("Can't find n_types in BornMayer_params")
          endif

          allocate(parse_ip%atomic_num(parse_ip%n_types))
          parse_ip%atomic_num = 0

          allocate(parse_ip%A(parse_ip%n_types,parse_ip%n_types))
          parse_ip%A = 0.0_dp
          allocate(parse_ip%b(parse_ip%n_types,parse_ip%n_types))
          parse_ip%b = 0.0_dp
          allocate(parse_ip%c(parse_ip%n_types,parse_ip%n_types))
          parse_ip%c = 0.0_dp
          allocate(parse_ip%cutoff_a(parse_ip%n_types,parse_ip%n_types))
          parse_ip%cutoff_a = 0.0_dp
          allocate(parse_ip%energy_shift(parse_ip%n_types,parse_ip%n_types))
          allocate(parse_ip%linear_force_shift(parse_ip%n_types,parse_ip%n_types))
          parse_ip%energy_shift = 0.0_dp
          parse_ip%linear_force_shift = 0.0_dp
       endif


    elseif (parse_in_ip .and. name == 'per_type_data') then

       call QUIP_FoX_get_value(attributes, "type", value, status)
       if (status /= 0) call system_abort ("IPModel_BornMayer_read_params_xml cannot find type")
       read (value, *) ti

       call QUIP_FoX_get_value(attributes, "atomic_num", value, status)
       if (status /= 0) call system_abort ("IPModel_BornMayer_read_params_xml cannot find atomic_num")
       read (value, *) parse_ip%atomic_num(ti)

       if (allocated(parse_ip%type_of_atomic_num)) deallocate(parse_ip%type_of_atomic_num)
       allocate(parse_ip%type_of_atomic_num(maxval(parse_ip%atomic_num)))
       parse_ip%type_of_atomic_num = 0
       do ti=1, parse_ip%n_types
          if (parse_ip%atomic_num(ti) > 0) &
               parse_ip%type_of_atomic_num(parse_ip%atomic_num(ti)) = ti
       end do

       parse_ip%energy_shift = 0.0_dp
       parse_ip%linear_force_shift = 0.0_dp

    elseif (parse_in_ip .and. name == 'per_pair_data') then

       call QUIP_FoX_get_value(attributes, "atnum_i", value, status)
       if (status /= 0) call system_abort ("IPModel_BornMayer_read_params_xml cannot find atnum_i")
       read (value, *) Zi
       call QUIP_FoX_get_value(attributes, "atnum_j", value, status)
       if (status /= 0) call system_abort ("IPModel_BornMayer_read_params_xml cannot find atnum_j")
       read (value, *) Zj

       ti = get_type(parse_ip%type_of_atomic_num,Zi)
       tj = get_type(parse_ip%type_of_atomic_num,Zj)

       call QUIP_FoX_get_value(attributes, "A", value, status)
       if (status /= 0) call system_abort ("IPModel_BornMayer_read_params_xml cannot find A")
       read (value, *) parse_ip%A(ti,tj)
       call QUIP_FoX_get_value(attributes, "b", value, status)
       if (status /= 0) call system_abort ("IPModel_BornMayer_read_params_xml cannot find b")
       read (value, *) parse_ip%b(ti,tj)
       call QUIP_FoX_get_value(attributes, "c", value, status)
       if (status /= 0) call system_abort ("IPModel_BornMayer_read_params_xml cannot find c")
       read (value, *) parse_ip%c(ti,tj)
       call QUIP_FoX_get_value(attributes, "cutoff", value, status)
       if (status /= 0) call system_abort ("IPModel_LJ_read_params_xml cannot find cutoff")
       read (value, *) parse_ip%cutoff_a(ti,tj)

       call QUIP_FoX_get_value(attributes, "energy_shift", value, status)
       if (status /= 0) call system_abort ("IPModel_LJ_read_params_xml cannot find energy_shift")
       read (value, *) energy_shift
       if (energy_shift) parse_ip%energy_shift(ti,tj) = IPModel_BornMayer_pairenergy(parse_ip, ti, tj, parse_ip%cutoff_a(ti,tj))
       
       call QUIP_FoX_get_value(attributes, "linear_force_shift", value, status)
       if (status /= 0) call system_abort ("IPModel_LJ_read_params_xml cannot find linear_force_shift")
       read (value, *) linear_force_shift
       if (linear_force_shift) parse_ip%linear_force_shift(ti,tj) = IPModel_BornMayer_pairenergy_deriv(parse_ip, ti, tj, parse_ip%cutoff_a(ti,tj))

       if (ti /= tj) then
          parse_ip%A(tj,ti) = parse_ip%A(ti,tj)
          parse_ip%b(tj,ti) = parse_ip%b(ti,tj)
          parse_ip%c(tj,ti) = parse_ip%c(ti,tj)
          parse_ip%cutoff_a(tj,ti) = parse_ip%cutoff_a(ti,tj)
          parse_ip%energy_shift(tj,ti) = parse_ip%energy_shift(ti,tj)
          parse_ip%linear_force_shift(tj,ti) = parse_ip%linear_force_shift(ti,tj)
       endif
    end if

  end subroutine IPModel_startElement_handler

  subroutine IPModel_endElement_handler(URI, localname, name)
    character(len=*), intent(in)   :: URI
    character(len=*), intent(in)   :: localname
    character(len=*), intent(in)   :: name

    if (parse_in_ip) then
       if (name == 'BornMayer_params') then
          parse_in_ip = .false.
       end if
    endif

  end subroutine IPModel_endElement_handler

end module IPModel_BornMayer_module
