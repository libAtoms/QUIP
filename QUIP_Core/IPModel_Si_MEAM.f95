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
!X IPModel_Si_MEAM  module 
!X
!% Module for the Modified Embedded Atom Method potential for Si
!% The IPModel_Si_MEAM object contains all the parameters read 
!% from an 'Si_MEAN_params' XML stanza.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!
#include "error.inc"

module IPModel_Si_MEAM_module

use libatoms_module
use mpi_context_module
use QUIP_Common_module


implicit none
private

type IPModel_Si_MEAM
  integer :: n_types = 0
  integer, allocatable :: atomic_num(:), type_of_atomic_num(:)

  real(dp) :: cutoff = 0.0_dp

  type(spline), dimension(:,:), allocatable :: phi
  type(spline), dimension(:,:), allocatable :: rho
  type(spline), dimension(:,:), allocatable :: f
  type(spline), dimension(:), allocatable :: U
  type(spline), dimension(:,:,:), allocatable :: g

  real(dp), dimension(:,:), allocatable :: r_cut_phi
  real(dp), dimension(:,:), allocatable :: r_cut_rho
  real(dp), dimension(:,:), allocatable :: r_cut_f

  character(len=STRING_LENGTH) :: label
endtype IPModel_Si_MEAM

interface Initialise
  module procedure IPModel_Si_MEAM_Initialise_str
end interface Initialise

interface Finalise
  module procedure IPModel_Si_MEAM_Finalise
end interface Finalise

interface Print
  module procedure IPModel_Si_MEAM_Print
end interface Print

interface Calc
  module procedure IPModel_Si_MEAM_Calc
end interface Calc

logical, private :: parse_in_ip, parse_matched_label
type(IPModel_Si_MEAM), private, pointer :: parse_ip
integer :: parse_cur_type_i, parse_cur_type_j, parse_cur_type_k, parse_cur_point, n_spline

real(dp) :: yp1, ypn
real(dp), dimension(:), allocatable :: spline_x, spline_y
character(len=100) :: spline_function

public :: IPModel_Si_MEAM, Initialise, Finalise, Print, Calc

contains

subroutine IPModel_Si_MEAM_Initialise_str(this, args_str, param_str)
  type(IPModel_Si_MEAM), intent(inout) :: this
  character(len=*), intent(in) :: args_str, param_str

  type(Dictionary) :: params

  call Finalise(this)

  call initialise(params)
  this%label=''
  call param_register(params, 'label', '', this%label, help_string="No help yet.  This source file was $LastChangedBy$")
  if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='IPModel_Si_MEAM_Initialise_str args_str')) then
    call system_abort("IPModel_Si_MEAM_Initialise_str failed to parse label from args_str="//trim(args_str))
  endif
  call finalise(params)

  call IPModel_Si_MEAM_read_params_xml(this, param_str)

end subroutine IPModel_Si_MEAM_Initialise_str

subroutine IPModel_Si_MEAM_Finalise(this)
  type(IPModel_Si_MEAM), intent(inout) :: this

  integer :: ti, tj, tk

  do ti=1, this%n_types
    call finalise(this%U(ti))
    do tj=1, this%n_types
      call finalise(this%phi(ti,tj))
      call finalise(this%rho(ti,tj))
      call finalise(this%f(ti,tj))
      do tk=1, this%n_types
         call finalise(this%g(ti,tj,tk))
      enddo
    end do
  end do
    
  if(allocated(this%phi)) deallocate(this%phi)
  if(allocated(this%rho)) deallocate(this%rho)
  if(allocated(this%f)) deallocate(this%f)
  if(allocated(this%U)) deallocate(this%U)
  if(allocated(this%g)) deallocate(this%g)
  if(allocated(this%atomic_num)) deallocate(this%atomic_num)
  if(allocated(this%type_of_atomic_num)) deallocate(this%type_of_atomic_num)

  if(allocated(this%r_cut_phi)) deallocate(this%r_cut_phi)
  if(allocated(this%r_cut_rho)) deallocate(this%r_cut_rho)
  if(allocated(this%r_cut_f)) deallocate(this%r_cut_f)

  this%n_types = 0
  this%label = ''
end subroutine IPModel_Si_MEAM_Finalise

subroutine IPModel_Si_MEAM_Calc(this, at, e, local_e, f, virial, local_virial, args_str, mpi, error)
  type(IPModel_Si_MEAM), intent(inout) :: this
  type(Atoms), intent(in) :: at
  real(dp), intent(out), optional :: e !% \texttt{e} = System total energy
  real(dp), dimension(:), intent(out), optional :: local_e !% \texttt{local_e} = energy of each atom, vector dimensioned as \texttt{at%N}.
  real(dp), intent(out), optional :: f(:,:), local_virial(:,:)   !% Forces, dimensioned as \texttt{f(3,at%N)}, local virials, dimensioned as \texttt{local_virial(9,at%N)} 
  real(dp), dimension(3,3), intent(out), optional :: virial   !% Virial
   character(len=*), optional      :: args_str
  type(MPI_Context), intent(in), optional :: mpi
  integer, intent(out), optional :: error

  !local
  integer  :: i, j, k, n, nn, ti, tj, tk
  
  real(dp) :: r_ij, r_ik, n_i, f_ij, f_ik, theta_jik, g_jik, U_i, phi_ij, &
            & dphi_ij, df_ij, df_ik, dg_jik, dU_i, drho_ij
  real(dp), dimension(3)  :: u_ij, u_ik, dn_i, dn_j, dn_i_dr_ij
  real(dp), dimension(3,3) :: dn_i_drij_outer_rij

  type(Dictionary)                :: params
  logical :: has_atom_mask_name
  character(STRING_LENGTH) :: atom_mask_name
  real(dp) :: r_scale, E_scale
  logical :: do_rescale_r, do_rescale_E

  INIT_ERROR(error)

  if (present(e)) e = 0.0_dp
  if (present(local_e)) then
     call check_size('Local_E',local_e,(/at%N/),'IPModel_Si_MEAM_Calc', error)
     local_e = 0.0_dp
  endif
  if (present(f)) then
     call check_size('Force',f,(/3,at%Nbuffer/),'IPModel_Si_MEAM_Calc', error)
     f = 0.0_dp
  end if
  if (present(virial)) virial = 0.0_dp
  if (present(local_virial)) then
     call check_size('Local_virial',local_virial,(/9,at%Nbuffer/),'IPModel_Si_MEAM_Calc', error)
     local_virial = 0.0_dp
     RAISE_ERROR("IPModel_Si_MEAM_Calc: local_virial calculation requested but not supported yet.", error)
  endif

  if (present(args_str)) then
     call initialise(params)
     call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name, help_string="No help yet.  This source file was $LastChangedBy$")
     call param_register(params, 'r_scale', '1.0',r_scale, has_value_target=do_rescale_r, help_string="Recaling factor for distances. Default 1.0.")
     call param_register(params, 'E_scale', '1.0',E_scale, has_value_target=do_rescale_E, help_string="Recaling factor for energy. Default 1.0.")

     if(.not. param_read_line(params, args_str, ignore_unknown=.true.,task='IPModel_Si_MEAM_Calc args_str')) then
        RAISE_ERROR("IPModel_Si_MEAM_Calc failed to parse args_str='"//trim(args_str)//"'",error)
     endif
     call finalise(params)
     if(has_atom_mask_name) then
        RAISE_ERROR('IPModel_Si_MEAM_Calc: atom_mask_name found, but not supported', error)
     endif
     if (do_rescale_r .or. do_rescale_E) then
        RAISE_ERROR("IPModel_Si_MEAM_Calc: rescaling of potential with r_scale and E_scale not yet implemented!", error)
     end if
  endif

  !Loop over atoms
  do i = 1, at%N

     if (present(mpi)) then
	if(mpi%active) then
	   if(mod(i-1, mpi%n_procs) /= mpi%my_proc) cycle
	endif
     endif
     ti = get_type(this%type_of_atomic_num, at%Z(i))

     n_i = 0.0_dp
     dn_i = 0.0_dp
     dn_i_drij_outer_rij = 0.0_dp
     !Loop over neighbours
     do n = 1, atoms_n_neighbours(at,i)

        j = atoms_neighbour(at, i, n, distance=r_ij, cosines=u_ij) ! nth neighbour of atom i
        tj = get_type(this%type_of_atomic_num, at%Z(j))

        if( r_ij < this%r_cut_phi(ti,tj) ) then
            if(present(local_e) .or. present(e)) phi_ij = calc_y(this%phi(ti,tj),r_ij)
            if(present(local_e)) local_e(i) = local_e(i) + phi_ij*0.5_dp
            if(present(e)) e = e + phi_ij*0.5_dp

            if(present(f) .or. present(virial) ) dphi_ij = calc_dy(this%phi(ti,tj),r_ij)
            if(present(f)) then
               f(:,i) = f(:,i) + 0.5_dp*dphi_ij*u_ij
               f(:,j) = f(:,j) - 0.5_dp*dphi_ij*u_ij
            endif
            if(present(virial)) virial = virial - 0.5_dp * dphi_ij*(u_ij .outer. u_ij)*r_ij
        endif

        if( (r_ij >= this%r_cut_rho(ti,tj)).or.(r_ij >= this%r_cut_f(ti,tj)) ) cycle

        n_i = n_i + calc_y(this%rho(ti,tj),r_ij)
        if(present(f) .or. present(virial) ) then
           dn_i_dr_ij = calc_dy(this%rho(ti,tj),r_ij)*u_ij
           if(present(f)) dn_i = dn_i + dn_i_dr_ij
           if(present(virial)) dn_i_drij_outer_rij = dn_i_drij_outer_rij + (dn_i_dr_ij .outer. u_ij) * r_ij
        endif
 
        f_ij = calc_y(this%f(ti,tj),r_ij)
        if(present(f) .or. present(virial)) df_ij = calc_dy(this%f(ti,tj),r_ij)

        do nn = 1, atoms_n_neighbours(at,i)

           k = atoms_neighbour(at, i, nn, distance=r_ik, cosines=u_ik)
           tk = get_type(this%type_of_atomic_num, at%Z(k))

           if( j>=k ) cycle
           if( r_ik >= this%r_cut_f(ti,tk) ) cycle

           theta_jik = cosine( at, i, j, k )

           f_ik = calc_y(this%f(ti,tk),r_ik)
           if(present(f) .or. present(virial)) df_ik = calc_dy(this%f(ti,tk),r_ik)

           g_jik = calc_y(this%g(ti,tj,tk),theta_jik)
           if(present(f) .or. present(virial)) dg_jik = calc_dy(this%g(ti,tj,tk),theta_jik)

           n_i = n_i + f_ij * f_ik * g_jik

           if( present(f) ) &
               & dn_i = dn_i + g_jik * &
               & ( df_ij * f_ik * u_ij + &
               &   df_ik * f_ij * u_ik ) + &
               & f_ij * f_ik * dg_jik * &
               & ( u_ij / r_ik + u_ik / r_ij - &
               & theta_jik * ( u_ij / r_ij + u_ik / r_ik ) )

           if( present(virial) ) &
               & dn_i_drij_outer_rij = dn_i_drij_outer_rij + &
               & g_jik * ( df_ij * f_ik * (u_ij.outer.u_ij) * r_ij + &
                           df_ik * f_ij * (u_ik.outer.u_ik) * r_ik ) + &
               & f_ij * f_ik * dg_jik * &
               & ( (u_ij .outer. u_ik) - theta_jik * (u_ij .outer. u_ij) + &
               &   (u_ik .outer. u_ij) - theta_jik * (u_ik .outer. u_ik) )
        enddo

     enddo
     if(present(local_e) .or. present(e)) U_i = calc_y(this%U(ti),n_i)
     if(present(f) .or. present(virial)) dU_i = calc_dy(this%U(ti),n_i)

     if(present(local_e)) local_e(i) = local_e(i) + U_i
     if(present(e)) e = e + U_i
     if(present(f)) f(:,i) = f(:,i) + dU_i * dn_i
     if(present(virial)) virial = virial - dU_i * dn_i_drij_outer_rij

     if(present(f)) then
        do n = 1, atoms_n_neighbours(at,i)  !cross-terms

           j = atoms_neighbour(at, i, n, distance=r_ij, cosines=u_ij) ! nth neighbour of atom i

           tj = get_type(this%type_of_atomic_num, at%Z(j))
           if( (r_ij >= this%r_cut_rho(ti,tj)).or.(r_ij >= this%r_cut_f(ti,tj)) ) cycle

           drho_ij = calc_dy(this%rho(ti,tj),r_ij)

           f(:,j) = f(:,j) - dU_i * drho_ij * u_ij

           f_ij = calc_y(this%f(ti,tj),r_ij)
           df_ij = calc_dy(this%f(ti,tj),r_ij)

           dn_j = 0.0_dp

           do nn = 1, atoms_n_neighbours(at,i)

              k = atoms_neighbour(at, i, nn, distance=r_ik, cosines=u_ik)
              tk = get_type(this%type_of_atomic_num, at%Z(k))

              if( j==k ) cycle
              if( r_ik >= this%r_cut_f(ti,tk) ) cycle

              theta_jik = cosine( at, i, j, k )

              f_ik = calc_y(this%f(ti,tk),r_ik)
              df_ik = calc_dy(this%f(ti,tk),r_ik)

              g_jik = calc_y(this%g(ti,tj,tk),theta_jik)
              dg_jik = calc_dy(this%g(ti,tj,tk),theta_jik)

              dn_j = dn_j + df_ij * f_ik * g_jik * u_ij + &
                   & f_ij * f_ik * dg_jik * ( u_ik / r_ij - theta_jik * u_ij / r_ij )

           enddo

           f(:,j) = f(:,j) - dU_i * dn_j

        enddo
     endif
  enddo

  if (present(mpi)) then
     if (present(e)) e = sum(mpi, e)
     if (present(local_e)) call sum_in_place(mpi, local_e)
     if (present(f)) call sum_in_place(mpi, f)
     if (present(virial)) call sum_in_place(mpi, virial)
  endif

endsubroutine IPModel_Si_MEAM_Calc   

function calc_y(this,x)

   type(Spline), intent(in) :: this
   real(dp), intent(in)     :: x

   real(dp)                 :: calc_y

   if( x < min_knot(this) ) then
       calc_y = this%yp1 * (x - min_knot(this))
   elseif( x > max_knot(this) ) then
       calc_y = this%ypn * (x - max_knot(this))
   else
       calc_y = spline_value(this,x)
   endif

endfunction calc_y

function calc_dy(this,x)

   type(Spline), intent(in) :: this
   real(dp), intent(in)     :: x

   real(dp)                 :: calc_dy

   if( x < min_knot(this) ) then
       calc_dy = this%yp1
   elseif( x > max_knot(this) ) then
       calc_dy = this%ypn
   else
       calc_dy = spline_deriv(this,x)
   endif

endfunction calc_dy

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X 
!X XML param reader functions
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

subroutine IPModel_startElement_handler(URI, localname, name, attributes)
  character(len=*), intent(in)   :: URI  
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name 
  type(dictionary_t), intent(in) :: attributes

  integer status
  character(len=1024) :: value

  integer :: ti, Zi, Zj, Zk
  integer :: n_types
  real(dp) :: v

  if (name == 'Si_MEAM_params') then 

    if (parse_in_ip) &
      call system_abort("IPModel_startElement_handler entered Si_MEAM_params with parse_in true. Probably a bug in FoX (4.0.1, e.g.)")

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

      call QUIP_FoX_get_value(attributes, "n_types", value, status)
      if (status /= 0) call system_abort ("IPModel_Si_MEAM_read_params_xml cannot find n_types")
      read (value, *) parse_ip%n_types
      n_types = parse_ip%n_types

      allocate(parse_ip%atomic_num(parse_ip%n_types))
      parse_ip%atomic_num = 0

      allocate(parse_ip%r_cut_phi(n_types,n_types))
      parse_ip%r_cut_phi = 0.0_dp
      allocate(parse_ip%r_cut_rho(n_types,n_types))
      parse_ip%r_cut_rho = 0.0_dp
      allocate(parse_ip%r_cut_f(n_types,n_types))
      parse_ip%r_cut_f = 0.0_dp

      allocate( parse_ip%phi(n_types,n_types))
      allocate( parse_ip%rho(n_types,n_types))
      allocate( parse_ip%f(n_types,n_types))
      allocate( parse_ip%U(n_types))
      allocate( parse_ip%g(n_types,n_types,n_types))

    endif

  elseif (parse_in_ip .and. name == 'per_type_data') then

    call QUIP_FoX_get_value(attributes, "type", value, status)
    if (status /= 0) call system_abort ("IPModel_Si_MEAM_read_params_xml cannot find type")
    read (value, *) ti

    if (ti < 1 .or. ti > parse_ip%n_types) call system_abort("IPModel_Si_MEAM_read_params_xml got" // &
      " per_type_data type out of range " // ti // " n_types " // parse_ip%n_types)

    parse_cur_type_i = ti

    call QUIP_FoX_get_value(attributes, "atomic_num", value, status)
    if (status /= 0) call system_abort ("IPModel_Si_MEAM_read_params_xml cannot find atomic_num")
    read (value, *) parse_ip%atomic_num(ti)

    if (allocated(parse_ip%type_of_atomic_num)) deallocate(parse_ip%type_of_atomic_num)
    allocate(parse_ip%type_of_atomic_num(maxval(parse_ip%atomic_num)))
    parse_ip%type_of_atomic_num = 0
    do ti=1, parse_ip%n_types
      if (parse_ip%atomic_num(ti) > 0) &
	parse_ip%type_of_atomic_num(parse_ip%atomic_num(ti)) = ti
    end do

  elseif (parse_in_ip .and. name == 'per_pair_data') then

    call QUIP_FoX_get_value(attributes, "atomic_num_i", value, status)
    if (status /= 0) call system_abort ("IPModel_Si_MEAM_read_params_xml cannot find atomic_num_i")
    read (value, *) Zi
    call QUIP_FoX_get_value(attributes, "atomic_num_j", value, status)
    if (status /= 0) call system_abort ("IPModel_Si_MEAM_read_params_xml cannot find atomic_num_j")
    read (value, *) Zj

    parse_cur_type_i = get_type(parse_ip%type_of_atomic_num,Zi)
    parse_cur_type_j = get_type(parse_ip%type_of_atomic_num,Zj)

    call QUIP_FoX_get_value(attributes, "r_cut_phi", value, status)
    if (status /= 0) call system_abort ("IPModel_Si_MEAM_read_params_xml cannot find r_cut_phi")
    read (value, *) parse_ip%r_cut_phi(parse_cur_type_i, parse_cur_type_j)
    call QUIP_FoX_get_value(attributes, "r_cut_rho", value, status)
    if (status /= 0) call system_abort ("IPModel_Si_MEAM_read_params_xml cannot find r_cut_rho")
    read (value, *) parse_ip%r_cut_rho(parse_cur_type_i, parse_cur_type_j)
    call QUIP_FoX_get_value(attributes, "r_cut_f", value, status)
    if (status /= 0) call system_abort ("IPModel_Si_MEAM_read_params_xml cannot find r_cut_f")
    read (value, *) parse_ip%r_cut_f(parse_cur_type_i, parse_cur_type_j)

  elseif (parse_in_ip .and. name == 'per_triplet_data') then

    call QUIP_FoX_get_value(attributes, "atomic_num_i", value, status)
    if (status /= 0) call system_abort ("IPModel_Si_MEAM_read_params_xml cannot find atomic_num_i")
    read (value, *) Zi
    call QUIP_FoX_get_value(attributes, "atomic_num_j", value, status)
    if (status /= 0) call system_abort ("IPModel_Si_MEAM_read_params_xml cannot find atomic_num_j")
    read (value, *) Zj
    call QUIP_FoX_get_value(attributes, "atomic_num_k", value, status)
    if (status /= 0) call system_abort ("IPModel_Si_MEAM_read_params_xml cannot find atomic_num_k")
    read (value, *) Zk

    parse_cur_type_i = get_type(parse_ip%type_of_atomic_num,Zi)
    parse_cur_type_j = get_type(parse_ip%type_of_atomic_num,Zj)
    parse_cur_type_k = get_type(parse_ip%type_of_atomic_num,Zk)

  elseif (parse_in_ip .and. name == 'spline') then
    call QUIP_FoX_get_value(attributes, "spline_function", value, status)
    if (status /= 0) call system_abort ("IPModel_Si_MEAM_read_params_xml cannot find spline_function")
    read (value, *) spline_function
    call QUIP_FoX_get_value(attributes, "n_spline", value, status)
    if (status /= 0) call system_abort ("IPModel_Si_MEAM_read_params_xml cannot find n_spline")
    read (value, *) n_spline
    call QUIP_FoX_get_value(attributes, "yp1", value, status)
    if (status /= 0) call system_abort ("IPModel_Si_MEAM_read_params_xml cannot find yp1")
    read (value, *) yp1
    call QUIP_FoX_get_value(attributes, "ypn", value, status)
    if (status /= 0) call system_abort ("IPModel_Si_MEAM_read_params_xml cannot find ypn")
    read (value, *) ypn
    allocate( spline_x(n_spline), spline_y(n_spline) )
    parse_cur_point = 1
  elseif (parse_in_ip .and. name == 'point') then

    if (parse_cur_point > n_spline) call system_abort ("IPModel_Si_MEAM got too " // &
    &"many points " // parse_cur_point // " type " // parse_cur_type_i // " " // parse_cur_type_j // &
    & " in spline "// trim(spline_function))

    call QUIP_FoX_get_value(attributes, "x", value, status)
    if (status /= 0) call system_abort ("IPModel_EAM_Ercolessi_Adams_read_params_xml cannot find x")
    read (value, *) v
    spline_x(parse_cur_point) = v

    call QUIP_FoX_get_value(attributes, "y", value, status)
    if (status /= 0) call system_abort ("IPModel_EAM_Ercolessi_Adams_read_params_xml cannot find y")
    read (value, *) v
    spline_y(parse_cur_point) = v

    parse_cur_point = parse_cur_point + 1
  endif

end subroutine IPModel_startElement_handler

subroutine IPModel_endElement_handler(URI, localname, name)
  character(len=*), intent(in)   :: URI  
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name 

  if (parse_in_ip) then
    if (name == 'Si_MEAM_params') then
      parse_in_ip = .false.
    elseif (name == 'spline') then
      selectcase(trim(spline_function))
        case('phi')
          call spline_init(parse_ip%phi(parse_cur_type_i,parse_cur_type_j), spline_x, spline_y, yp1, ypn)
          if( parse_cur_type_i /= parse_cur_type_j ) &
          &   call spline_init(parse_ip%phi(parse_cur_type_j,parse_cur_type_i), spline_x, spline_y, yp1, ypn)
        case('rho')
          call spline_init(parse_ip%rho(parse_cur_type_i,parse_cur_type_j), spline_x, spline_y, yp1, ypn)
          if( parse_cur_type_i /= parse_cur_type_j ) &
          &   call spline_init(parse_ip%rho(parse_cur_type_j,parse_cur_type_i), spline_x, spline_y, yp1, ypn)
        case('f')
          call spline_init(parse_ip%f(parse_cur_type_i,parse_cur_type_j), spline_x, spline_y, yp1, ypn)
          if( parse_cur_type_i /= parse_cur_type_j ) &
          &   call spline_init(parse_ip%f(parse_cur_type_j,parse_cur_type_i), spline_x, spline_y, yp1, ypn)
        case('U')
          call spline_init(parse_ip%U(parse_cur_type_i), spline_x, spline_y, yp1, ypn)
        case('g')
          call spline_init(parse_ip%g(parse_cur_type_i,parse_cur_type_j,parse_cur_type_k), &
          & spline_x, spline_y, yp1, ypn)
          if( parse_cur_type_j /= parse_cur_type_k ) &
          &   call spline_init(parse_ip%g(parse_cur_type_i,parse_cur_type_k,parse_cur_type_j), &
          & spline_x, spline_y, yp1, ypn) 
        case default
          call system_abort('')
      endselect
      deallocate(spline_x, spline_y)
    endif
  endif

end subroutine IPModel_endElement_handler

subroutine IPModel_Si_MEAM_read_params_xml(this, param_str)
  type(IPModel_Si_MEAM), intent(inout), target :: this
  character(len=*), intent(in) :: param_str

  type (xml_t) :: fxml

  if (len(trim(param_str)) <= 0) return

  parse_ip => this
  parse_in_ip = .false.
  parse_matched_label = .false.

  call open_xml_string(fxml, param_str)

  call parse(fxml, &
    startElement_handler = IPModel_startElement_handler, &
    endElement_handler = IPModel_endElement_handler)

  call close_xml_t(fxml)

  if (this%n_types == 0) call system_abort("Si_MEAM Tried to parse, but n_types still 0")

  parse_ip%cutoff = max(maxval(parse_ip%r_cut_phi),maxval(parse_ip%r_cut_rho),maxval(parse_ip%r_cut_f))

end subroutine IPModel_Si_MEAM_read_params_xml


!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X 
!% Printing of potential parameters: number of different types, cutoff radius, atomic numbers, ect.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

subroutine IPModel_Si_MEAM_Print (this, file)
  type(IPModel_Si_MEAM), intent(in) :: this
  type(Inoutput), intent(inout),optional :: file

  integer :: ti, tj, tk

  call Print("IPModel_Si_MEAM : n_types = " // this%n_types // " cutoff = " // this%cutoff, file=file)

  do ti=1, this%n_types
    call Print("IPModel_Si_MEAM : type " // ti // " atomic_num " // this%atomic_num(ti), file=file)
    call verbosity_push_decrement()
    call Print("IPModel_Si_MEAM : U ", file=file)
    call Print(this%U(ti), file=file)
    do tj=1, this%n_types
      call Print("IPModel_Si_MEAM : pair "// ti // " " // tj // " r_cut_phi " // this%r_cut_phi(ti,tj), file=file)
      call Print("IPModel_Si_MEAM : pair phi ", file=file)
      call Print(this%phi(ti,tj), file=file)

      call Print("IPModel_Si_MEAM : pair "// ti // " " // tj // " r_cut_rho " // this%r_cut_rho(ti,tj), file=file)
      call Print("IPModel_Si_MEAM : pair rho ", file=file)
      call Print(this%rho(ti,tj), file=file)

      call Print("IPModel_Si_MEAM : pair "// ti // " " // tj // " r_cut_f " // this%r_cut_f(ti,tj), file=file)
      call Print("IPModel_Si_MEAM : pair f ", file=file)
      call Print(this%f(ti,tj), file=file)
      do tk=1, this%n_types
         call Print("IPModel_Si_MEAM : triplet g", file=file)
         call Print(this%g(ti,tj,tk), file=file)
      enddo
    end do
    call verbosity_pop()
  end do
    
end subroutine IPModel_Si_MEAM_Print

endmodule IPModel_Si_MEAM_module
!<Si_MEAM_params n_types="1" label="">
!<!-- Thomas J Lenosky, Babak Sadigh, Eduardo Alonso, Vasily V Bulatov, Tomas Diaz de la Rubia -->
!<!-- Jeongnim Kim, Arthur F Voter and Joel D Kress -->
!<!-- Highly optimized empirical potential model of silicon -->
!<!-- Modelling Simul. Mater. Sci. Eng. 8 (2000) 825-841 -->
!  <per_type_data atomic_num="14" type="1">
!    <spline spline_function="U" n_spline="8" yp1="0.73514" ypn="0.61652">
!      <point x="-1.7709300000" y="-1.0749300000" />
!      <point x="-0.3881514286" y="-0.2004500000" />
!      <point x=" 0.9946271429" y=" 0.4142200000" />
!      <point x=" 2.3774057143" y=" 0.8793900000" />
!      <point x=" 3.7601842857" y=" 1.2668900000" />
!      <point x=" 5.1429628571" y=" 1.6299800000" />
!      <point x=" 6.5257414286" y=" 1.9773800000" />
!      <point x=" 7.9085200000" y=" 2.3961800000" />
!    </spline>
!  </per_type_data>
!  <per_pair_data atomic_num_i="14" atomic_num_j="14" r_cut_phi="4.5" r_cut_rho="3.5" r_cut_f="3.5">
!    <spline spline_function="phi" n_spline="10" yp1="-42.66967" ypn="0.00000">
!      <point x="1.5000000000" y=" 6.9299400000" />
!      <point x="1.8333333333" y="-0.4399500000" />
!      <point x="2.1666666667" y="-1.7012300000" />
!      <point x="2.5000000000" y="-1.6247300000" />
!      <point x="2.8333333333" y="-0.9969600000" />
!      <point x="3.1666666667" y="-0.2739100000" />
!      <point x="3.5000000000" y="-0.0249900000" />
!      <point x="3.8333333333" y="-0.0178400000" />
!      <point x="4.1666666667" y="-0.0096100000" />
!      <point x="4.5000000000" y=" 0.0000000000" />
!    </spline>
!    <spline spline_function="rho" n_spline="11" yp1="-1.00000" ypn="0.00000">
!      <point x="1.5000000000" y=" 0.1374700000" />
!      <point x="1.7000000000" y="-0.1483100000" />
!      <point x="1.9000000000" y="-0.5597200000" />
!      <point x="2.1000000000" y="-0.7311000000" />
!      <point x="2.3000000000" y="-0.7628300000" />
!      <point x="2.5000000000" y="-0.7291800000" />
!      <point x="2.7000000000" y="-0.6662000000" />
!      <point x="2.9000000000" y="-0.5732800000" />
!      <point x="3.1000000000" y="-0.4069000000" />
!      <point x="3.3000000000" y="-0.1666200000" />
!      <point x="3.5000000000" y=" 0.0000000000" />
!    </spline>
!    <spline spline_function="f" n_spline="10" yp1="-3.61894" ypn="0.00000">
!      <point x="1.5000000000" y="1.2503100000" />
!      <point x="1.7222222222" y="0.8682100000" />
!      <point x="1.9444444444" y="0.6084600000" />
!      <point x="2.1666666667" y="0.4875600000" />
!      <point x="2.3888888889" y="0.4416300000" />
!      <point x="2.6111111111" y="0.3761000000" />
!      <point x="2.8333333333" y="0.2714500000" />
!      <point x="3.0555555556" y="0.1481400000" />
!      <point x="3.2777777778" y="0.0485500000" />
!      <point x="3.5000000000" y="0.0000000000" />
!    </spline>
!  </per_pair_data>
!  <per_triplet_data atomic_num_i="14" atomic_num_j="14" atomic_num_k="14">
!    <spline spline_function="g" n_spline="8" yp1="-13.95042" ypn="1.13462">
!      <point x="-1.0000000000" y="5.2541600000" />
!      <point x="-0.7428371429" y="2.3591500000" />
!      <point x="-0.4856742857" y="1.1959500000" />
!      <point x="-0.2285114286" y="1.2299500000" />
!      <point x=" 0.0286514286" y="2.0356500000" />
!      <point x=" 0.2858142857" y="3.4247400000" />
!      <point x=" 0.5429771429" y="4.9485900000" />
!      <point x=" 0.8001400000" y="5.6179900000" />
!    </spline>
!  </per_triplet_data>
!</Si_MEAM_params>
