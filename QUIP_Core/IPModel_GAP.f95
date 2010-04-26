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
!X IPModel_GAP module  
!X
!% Module for Gaussian Approximation Potential.
!%
!% The IPModel_GAP object contains all the parameters read from a
!% 'GAP_params' XML stanza.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
module IPModel_GAP_module

use libatoms_module

use mpi_context_module
use QUIP_Common_module
use IPEwald_module

#ifdef HAVE_GP
use bispectrum_module
use gp_sparse_module
#endif

implicit none

private 

include 'IPModel_interface.h'

! this stuff is here for now, but it should live somewhere else eventually
! lower down in the GP

public :: IPModel_GAP
type IPModel_GAP

  integer :: n_types = 0

  integer, allocatable :: atomic_num(:), type_of_atomic_num(:) !% Atomic number dimensioned as \texttt{n_types}.
  real(dp) :: cutoff = 0.0_dp                                  !% Cutoff for computing connection.

  ! bispectrum parameters
  integer :: j_max = 0
  real(dp) :: z0 = 0.0_dp
  integer :: n_species = 0                                       !% Number of atomic types.
  integer, dimension(:), allocatable :: Z
  logical :: do_ewald = .false.
  real(dp), dimension(116) :: z_eff = 0.0_dp
  real(dp), dimension(116) :: w_Z = 1.0_dp

  ! qw parameters
  integer :: qw_l_max = 0
  integer :: qw_f_n = 0
  logical :: qw_do_q = .false.
  logical :: qw_do_w = .false.
  real(dp), allocatable :: qw_cutoff(:)
  integer, allocatable :: qw_cutoff_f(:)
  real(dp), allocatable :: qw_cutoff_r1(:)

  character(len=256) :: datafile                               !% File name containing the GAP database
  character(len=value_len) :: datafile_coordinates             !% Coordinate system used in GAP database

  character(len=FIELD_LENGTH) :: label
  type(mpi_context) :: mpi

#ifdef HAVE_GP
  type(gp) :: my_gp
#endif

end type IPModel_GAP

logical :: parse_in_ip, parse_matched_label
type(IPModel_GAP), pointer :: parse_ip

interface Initialise
  module procedure IPModel_GAP_Initialise_str
end interface Initialise

interface Finalise
  module procedure IPModel_GAP_Finalise
end interface Finalise

interface Print
  module procedure IPModel_GAP_Print
end interface Print

interface Calc
  module procedure IPModel_GAP_Calc
end interface Calc

contains

subroutine IPModel_GAP_Initialise_str(this, args_str, param_str, mpi)
  type(IPModel_GAP), intent(inout) :: this
  character(len=*), intent(in) :: args_str, param_str
  type(mpi_context), intent(in), optional :: mpi

  type(Dictionary) :: my_dictionary
  integer :: i, n_species
  real(dp), dimension(:), allocatable :: w, z_eff

  call Finalise(this)

  if (present(mpi)) this%mpi = mpi

  ! now initialise the potential

  !call IPModel_GAP_read_params_xml(this, param_str)

#ifdef HAVE_GP
  !call gp_read_binary(this%my_gp, 'this%datafile')
  this%datafile = 'gp.dat'
  call gp_read_binary(this%my_gp, 'gp.dat')
  call read_string(my_dictionary, this%my_gp%comment)

  this%datafile_coordinates = ''
  if (.not. get_value(my_dictionary, 'coordinates', this%datafile_coordinates)) &
     this%datafile_coordinates = 'bispectrum'

  if (trim(this%datafile_coordinates) == 'bispectrum') then
     if( .not. ( get_value(my_dictionary,'cutoff',this%cutoff) .and. &
               & get_value(my_dictionary,'j_max',this%j_max) .and. &
               & get_value(my_dictionary,'z0',this%z0) ) ) &
     & call system_abort('Did not find bispectrum parameters in gp.dat file, &
     & might be old version or format not correct')
  elseif (trim(this%datafile_coordinates) == 'qw') then
     if (.not. get_value(my_dictionary, 'l_max', this%qw_l_max) .and. &
               get_value(my_dictionary, 'f_n', this%qw_f_n)) &
        call system_abort('Did not find qw parameters in gp.dat file')

     allocate(this%qw_cutoff(this%qw_f_n), this%qw_cutoff_f(this%qw_f_n), this%qw_cutoff_r1(this%qw_f_n))

     do i = 1, this%qw_f_n
        if (.not. (get_value(my_dictionary, 'cutoff_' // i, this%qw_cutoff(i)) .and. &
                   get_value(my_dictionary, 'cutoff_f_' // i, this%qw_cutoff_f(i)) .and. &
                   get_value(my_dictionary, 'cutoff_r1_' // i, this%qw_cutoff_r1(i)))) &
           call system_abort('Did not find qw parameters in gp.dat file')
     enddo

     if (.not. get_value(my_dictionary, 'do_q', this%qw_do_q)) this%qw_do_q = .true.
     if (.not. get_value(my_dictionary, 'do_w', this%qw_do_w)) this%qw_do_w = .true.

     this%cutoff = maxval(this%qw_cutoff)
  endif

  if( get_value(my_dictionary,'n_species',n_species) ) then
     allocate( this%Z(n_species), w(n_species), z_eff(n_species) )
     if( n_species == 1 ) then
        if( .not. ( get_value(my_dictionary,'Z',this%Z(1)) .and. get_value(my_dictionary,'w',w(1)) .and. &
        & get_value(my_dictionary,'z_eff',z_eff(1)) ) ) call system_abort('')
     else
        if( .not. ( get_value(my_dictionary,'Z',this%Z) .and. get_value(my_dictionary,'w',w) .and. &
        & get_value(my_dictionary,'z_eff',z_eff) ) ) call system_abort('')
     endif
     
     this%w_Z = 0.0_dp
     this%z_eff = 0.0_dp
     do i = 1, n_species
        this%w_Z(this%Z(i)) = w(i)
        this%z_eff(this%Z(i)) = z_eff(i)
     enddo
     deallocate(w, z_eff)
  endif
  if( .not. get_value(my_dictionary,'do_ewald',this%do_ewald) ) call system_abort('')
     
  call finalise(my_dictionary)
#endif  

end subroutine IPModel_GAP_Initialise_str

subroutine IPModel_GAP_Finalise(this)
  type(IPModel_GAP), intent(inout) :: this

  if (allocated(this%atomic_num)) deallocate(this%atomic_num)
  if (allocated(this%type_of_atomic_num)) deallocate(this%type_of_atomic_num)

  this%n_types = 0

  if (allocated(this%qw_cutoff)) deallocate(this%qw_cutoff)
  if (allocated(this%qw_cutoff_f)) deallocate(this%qw_cutoff_f)
  if (allocated(this%qw_cutoff_r1)) deallocate(this%qw_cutoff_r1)

  this%datafile = ''
  this%datafile_coordinates = ''

  this%label = ''

end subroutine IPModel_GAP_Finalise

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X 
!% The potential calculator: this routine computes energy, forces and the virial.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

subroutine IPModel_GAP_Calc(this, at, e, local_e, f, virial)
  type(IPModel_GAP), intent(in) :: this
  type(Atoms), intent(in) :: at
  real(dp), intent(out), optional :: e, local_e(:) !% \texttt{e} = System total energy, \texttt{local_e} = energy of each atom, vector dimensioned as \texttt{at%N}.  
  real(dp), intent(out), optional :: f(:,:)        !% Forces, dimensioned as \texttt{f(3,at%N)} 
  real(dp), intent(out), optional :: virial(3,3)   !% Virial

  real(dp), pointer :: w_e(:)
  real(dp) :: e_i, f_gp, f_gp_k
  real(dp), dimension(:), allocatable   :: local_e_in, w, charge
  real(dp), dimension(:,:,:), allocatable   :: virial_in
  real(dp), dimension(:,:), allocatable   :: vec, f_ewald
  real(dp), dimension(:,:,:), allocatable   :: jack
  integer :: d, i, j, k, n, nei_max, jn

  real(dp) :: e_ewald
  real(dp), dimension(3,3) :: virial_ewald

  integer, dimension(3) :: shift

#ifdef HAVE_GP
  type(fourier_so4), save :: f_hat
  type(grad_fourier_so4), save :: df_hat
  type(bispectrum_so4), save :: bis
  type(grad_bispectrum_so4), save :: dbis

  type(fourier_so3), save :: f3_hat
  type(grad_fourier_so3), save :: df3_hat
  type(qw_so3), save :: qw
  type(grad_qw_so3), save :: dqw

  !$omp threadprivate(f_hat,df_hat,bis,dbis)  
  !$omp threadprivate(f3_hat,df3_hat,qw,dqw)  
#endif  

  if (present(e)) then
     e = 0.0_dp
     e_ewald = 0.0_dp
  endif

  if (present(local_e)) then
     call check_size('Local_E',local_e,(/at%N/),'IPModel_GAP_Calc')
     local_e = 0.0_dp
  endif

  if (present(f)) then 
     call check_size('Force',f,(/3,at%N/),'IPModel_GAP_Calc')
     f = 0.0_dp
  end if
  if (present(virial)) then
     virial = 0.0_dp
     virial_ewald = 0.0_dp
  endif

  if(present(e) .or. present(local_e) ) then
     allocate(local_e_in(at%N))
     local_e_in = 0.0_dp
  endif

  if (present(virial)) then
     allocate(virial_in(3,3,at%N))
     virial_in = 0.0_dp
  endif

  if (.not. assign_pointer(at, "weight", w_e)) nullify(w_e)

  allocate(w(at%N))
  do i = 1, at%N
     w(i) = this%w_Z(at%Z(i))
  enddo

#ifdef HAVE_GP
  if (trim(this%datafile_coordinates) == 'bispectrum') then
     d = j_max2d(this%j_max)
     call cg_initialise(this%j_max,2)
  elseif (trim(this%datafile_coordinates) == 'qw') then
     d = ((this%qw_l_max / 2) + 1) * this%qw_f_n
     if (this%qw_do_q .and. this%qw_do_w) d = d * 2
  endif
#endif

  nei_max = 0
  do i = 1, at%N
     if( nei_max < (atoms_n_neighbours(at,i)+1) ) nei_max = atoms_n_neighbours(at,i)+1
  enddo

  allocate(vec(d,at%N),jack(d,3*nei_max,at%N))

!$omp parallel 
#ifdef HAVE_GP
  if (trim(this%datafile_coordinates) == 'bispectrum') then
     call initialise(f_hat,this%j_max,this%z0,this%cutoff)
     if(present(f).or.present(virial)) call initialise(df_hat,this%j_max,this%z0,this%cutoff)
  elseif (trim(this%datafile_coordinates) == 'qw') then
     call initialise(f3_hat, this%qw_l_max, this%qw_cutoff, this%qw_cutoff_f, this%qw_cutoff_r1)
     call initialise(qw, this%qw_l_max, this%qw_f_n, do_q = this%qw_do_q, do_w = this%qw_do_w)
     if (present(f) .or. present(virial)) then
        call initialise(df3_hat, this%qw_l_max, this%qw_cutoff, this%qw_cutoff_f, this%qw_cutoff_r1)
        call initialise(dqw, this%qw_l_max, this%qw_f_n, do_q = this%qw_do_q, do_w = this%qw_do_w)
     endif
  endif
#endif
!$omp do private(n)
  do i = 1, at%N

#ifdef HAVE_GP
     if (trim(this%datafile_coordinates) == 'bispectrum') then
        call fourier_transform(f_hat,at,i,w)
        call calc_bispectrum(bis,f_hat)
        call bispectrum2vec(bis,vec(:,i))
        if(present(f).or.present(virial)) then
           do n = 0, atoms_n_neighbours(at,i)
              call fourier_transform(df_hat,at,i,n,w)
              call calc_bispectrum(dbis,f_hat,df_hat)
              call bispectrum2vec(dbis,jack(:,3*n+1:3*(n+1),i))
           enddo
        endif
     elseif (trim(this%datafile_coordinates) == 'qw') then
        call fourier_transform(f3_hat, at, i)
        call calc_qw(qw, f3_hat)
        call qw2vec(qw, vec(:,i))
        if (present(f) .or. present(virial)) then
           do n = 0, atoms_n_neighbours(at, i)
              call fourier_transform(df3_hat, at, i, n)
              call calc_qw(dqw, f3_hat, df3_hat)
              call qw2vec(dqw, jack(:,3*n+1:3*(n+1),i))
           enddo
        endif
     endif
#endif
  enddo
!$omp end do 
#ifdef HAVE_GP
  if (trim(this%datafile_coordinates) == 'bispectrum') then
     call finalise(f_hat)
     call finalise(df_hat)
     call finalise(bis)
     call finalise(dbis)
  elseif (trim(this%datafile_coordinates) == 'qw') then
     call finalise(f3_hat)
     call finalise(df3_hat)
     call finalise(qw)
     call finalise(dqw)
  endif
#endif
!$omp end parallel
    
!$omp parallel do private(k,f_gp,f_gp_k,n,j,jn,shift)
  do i = 1, at%N
     if(present(e) .or. present(local_e)) then
#ifdef HAVE_GP
        call gp_predict(gp_data=this%my_gp, mean=local_e_in(i),x_star=vec(:,i),Z=at%Z(i))
#endif
     endif

     if(present(f).or.present(virial)) then
        do k = 1, 3
           f_gp = 0.0_dp
           !kk = (i-1)*3 + k
       
           !call gp_predict(gp_data=this%my_gp, mean=f_gp_k,x_star=vec(:,i),x_prime_star=jack(:,kk,i))
#ifdef HAVE_GP
           call gp_predict(gp_data=this%my_gp, mean=f_gp_k,x_star=vec(:,i),x_prime_star=jack(:,k,i),Z=at%Z(i))
#endif
           f_gp = f_gp - f_gp_k
       
           if( present(virial) ) virial_in(:,k,i) = virial_in(:,k,i) - f_gp_k*at%pos(:,i)

           do n = 1, atoms_n_neighbours(at,i)
              j = atoms_neighbour(at,i,n,jn=jn,shift=shift)
       
!              if(jn==i)cycle              
#ifdef HAVE_GP
              call gp_predict(gp_data=this%my_gp,mean=f_gp_k,x_star=vec(:,j),x_prime_star=jack(:,jn*3+k,j),Z=at%Z(j))
#endif
              !call gp_predict(gp_data=this%my_gp,mean=f_gp_k,x_star=vec(:,j),x_prime_star=jack(:,kk,j))
              f_gp = f_gp - f_gp_k
              !if( present(virial) ) virial_in(:,k,i) = virial_in(:,k,i) - f_gp_k*( at%pos(:,i) - matmul(shift,at%lattice) )
              if( present(virial) ) virial_in(:,k,i) = virial_in(:,k,i) - f_gp_k*( at%pos(:,i) - matmul(at%lattice,shift) )
           enddo
       
           if(present(f)) f(k,i) = f_gp
        enddo
     endif
  enddo

  deallocate(vec,jack)

  if( this%do_ewald ) then
     allocate(charge(at%N))
     if(present(f)) allocate(f_ewald(3,at%N))
     do i = 1, at%N
        charge(i) = this%z_eff(at%Z(i))
     enddo

     if( present(e) .and. .not.present(f) .and. .not.present(virial)) call Ewald_calc(at,e=e_ewald,charge=charge)
     if( present(f) .and. .not.present(e) .and. .not.present(virial)) call Ewald_calc(at,f=f_ewald,charge=charge)
     if( present(virial) .and. .not.present(e) .and. .not.present(f)) call Ewald_calc(at,virial=virial_ewald,charge=charge)

     if( present(e) .and. present(f) .and. .not.present(virial)) call Ewald_calc(at,e=e_ewald,f=f_ewald,charge=charge)
     if( present(e) .and. present(virial) .and. .not.present(f)) call Ewald_calc(at,e=e_ewald,virial=virial_ewald,charge=charge)
     if( present(f) .and. present(virial) .and. .not.present(e)) call Ewald_calc(at,f=f_ewald,virial=virial_ewald,charge=charge)

     if( present(e) .and. present(f) .and. present(virial)) call Ewald_calc(at,e=e_ewald,f=f_ewald,virial=virial_ewald,charge=charge)
     
     if(present(f)) f = f + f_ewald
     deallocate(charge)
     if(allocated(f_ewald)) deallocate(f_ewald)
  endif

  if(present(e)) e = sum(local_e_in) + e_ewald
  if(present(local_e)) local_e = local_e_in
  if(present(virial)) virial = sum(virial_in,dim=3) + virial_ewald

  if(allocated(local_e_in)) deallocate(local_e_in)
  if(allocated(virial_in)) deallocate(virial_in)
  if(allocated(w)) deallocate(w)

end subroutine IPModel_GAP_Calc

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X 
!% XML param reader functions.
!% An example for XML stanza is given below, please notice that
!% they are simply dummy parameters for testing purposes, with no physical meaning.
!%
!%> <GAP_params datafile="file" n_types="1" label="default">
!%> <per_type_data type="1" atomic_num="14" />
!%> </GAP_params>
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

subroutine IPModel_startElement_handler(URI, localname, name, attributes)
  character(len=*), intent(in)   :: URI  
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name 
  type(dictionary_t), intent(in) :: attributes

  integer :: status
  character(len=1024) :: value

  logical shifted
  integer ti, tj

  if (name == 'GAP_params') then ! new GAP stanza
     
    if (parse_in_ip) &
          call system_abort("IPModel_startElement_handler entered GAP_params with parse_in true. Probably a bug in FoX (4.0.1, e.g.)")
     
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
	call system_abort("Can't find n_types in GAP_params")
      endif

      allocate(parse_ip%atomic_num(parse_ip%n_types))
      parse_ip%atomic_num = 0

      parse_ip%datafile = ''
      call QUIP_FoX_get_value(attributes, 'datafile', value, status)
      if (status == 0) then
         parse_ip%datafile = trim(value)
      else
	call system_abort("Can't find datafile in GAP_params")
      endif

    endif

  elseif (parse_in_ip .and. name == 'per_type_data') then

    call QUIP_FoX_get_value(attributes, "type", value, status)
    if (status /= 0) call system_abort ("IPModel_GAP_read_params_xml cannot find type")
    read (value, *) ti

    call QUIP_FoX_get_value(attributes, "atomic_num", value, status)
    if (status /= 0) call system_abort ("IPModel_GAP_read_params_xml cannot find atomic_num")
    read (value, *) parse_ip%atomic_num(ti)

    ! build reverse lookup table to atomic types
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
    if (name == 'GAP_params') then
      parse_in_ip = .false.
    end if
  endif

end subroutine IPModel_endElement_handler

subroutine IPModel_GAP_read_params_xml(this, param_str)
  type(IPModel_GAP), intent(inout), target :: this
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
    call system_abort("IPModel_GAP_read_params_xml parsed file, but n_types = 0")
  endif

end subroutine IPModel_GAP_read_params_xml

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X 
!% Printing of GAP parameters: number of different types, cutoff radius, atomic numbers, etc.
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

subroutine IPModel_GAP_Print (this, file)
  type(IPModel_GAP), intent(in) :: this
  type(Inoutput), intent(inout),optional :: file

  integer :: ti, tj

  call Print("IPModel_GAP : Gaussian Approximation Potential", file=file)
  call Print("IPModel_GAP : datafile = "//trim(this%datafile)//" n_types = " // this%n_types // " cutoff = " // this%cutoff, file=file)

  do ti=1, this%n_types
    call Print ("IPModel_GAP : type " // ti // " atomic_num " // this%atomic_num(ti), file=file)
  end do

end subroutine IPModel_GAP_Print

end module IPModel_GAP_module
