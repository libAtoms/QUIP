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

#ifdef HAVE_GP
use bispectrum_module
use qw_continous_module
use gp_sparse_module
#endif

implicit none

private 

include 'IPModel_interface.h'

! this stuff is here for now, but it should live somewhere else eventually
! lower down in the GP

public :: IPModel_GAP
type IPModel_GAP
  integer :: n_types = 0                                       !% Number of atomic types.
  integer, allocatable :: atomic_num(:), type_of_atomic_num(:) !% Atomic number dimensioned as \texttt{n_types}.
  real(dp) :: cutoff = 0.0_dp                                  !% Cutoff for computing connection.

  ! bispectrum parameters
  integer :: j_max = 0
  real(dp) :: z0 = 0.0_dp
  real(dp), dimension(:), allocatable :: w_Z

  ! qw parameters
  integer :: qw_dim = 0
  integer, allocatable :: qw_type(:)
  integer, allocatable :: qw_l(:)
  real(dp), allocatable :: qw_cutoff(:)
  real(dp), allocatable :: qw_wf_char_length(:)
  integer, allocatable :: qw_wf_type(:)
  logical, allocatable :: qw_do_weight(:)

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
  integer, dimension(:), allocatable :: Z
  real(dp), dimension(:), allocatable :: w
  integer :: i

  call Finalise(this)

  if (present(mpi)) this%mpi = mpi

  ! now initialise the potential

  !call IPModel_GAP_read_params_xml(this, param_str)

#ifdef HAVE_GP
  !call gp_read_binary(this%my_gp, 'this%datafile')
  this%datafile = 'gp.dat'
  call gp_read_binary(this%my_gp, 'gp.dat')
  call read_string(my_dictionary, this%my_gp%comment)
#endif

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
     if (.not. get_value(my_dictionary, 'qw_dim', this%qw_dim)) &
        call system_abort('Did not find qw dimensionality')

     allocate(this%qw_type(this%qw_dim), this%qw_l(this%qw_dim), this%qw_cutoff(this%qw_dim), &
              this%qw_wf_char_length(this%qw_dim), this%qw_wf_type(this%qw_dim), this%qw_do_weight(this%qw_dim))

     do i = 1, this%qw_dim
        if (.not. (get_value(my_dictionary, 'qw_type_' // i, this%qw_type(i)) .and. &
                   get_value(my_dictionary, 'qw_l_' // i, this%qw_l(i)) .and. &
                   get_value(my_dictionary, 'qw_cutoff_' // i, this%qw_cutoff(i)) .and. &
                   get_value(my_dictionary, 'qw_wf_char_length_' // i, this%qw_wf_char_length(i)) .and. &
                   get_value(my_dictionary, 'qw_wf_type_' // i, this%qw_wf_type(i)) .and. &
                   get_value(my_dictionary, 'qw_do_weight_' // i, this%qw_do_weight(i)))) &
           call system_abort('Did not find qw parameters in')
     enddo

     this%cutoff = maxval(this%qw_cutoff)
  endif


  if( get_value(my_dictionary,'n_species',n_species) ) then
     allocate( Z(n_species), w(n_species) )
     if( .not. ( get_value(my_dictionary,'Z',Z) .and. get_value(my_dictionary,'w',w) ) ) &
     & call system_abort('')
     allocate(this%w_Z(maxval(Z)))
     do i = 1, n_species
        this%w_Z(Z(i)) = w(i)
     enddo
  endif
     
  call finalise(my_dictionary)
#endif  

  if( allocated(Z) ) deallocate(Z)
  if( allocated(w) ) deallocate(w)

end subroutine IPModel_GAP_Initialise_str

subroutine IPModel_GAP_Finalise(this)
  type(IPModel_GAP), intent(inout) :: this

  if (allocated(this%atomic_num)) deallocate(this%atomic_num)
  if (allocated(this%type_of_atomic_num)) deallocate(this%type_of_atomic_num)

  this%n_types = 0

  if (allocated(this%qw_type)) deallocate(this%qw_type)
  if (allocated(this%qw_l)) deallocate(this%qw_l)
  if (allocated(this%qw_cutoff)) deallocate(this%qw_cutoff)
  if (allocated(this%qw_wf_char_length)) deallocate(this%qw_wf_char_length)
  if (allocated(this%qw_wf_type)) deallocate(this%qw_wf_type)
  if (allocated(this%qw_do_weight)) deallocate(this%qw_do_weight)

  this%qw_dim = 0

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
  real(dp), dimension(:), allocatable   :: local_e_in, w
  real(dp), dimension(:,:,:), allocatable   :: virial_in
  real(dp), dimension(:,:), allocatable   :: vec
  real(dp), dimension(:,:,:), allocatable   :: jack
  integer :: d, i, j, k, l, n, nei_max, jn

  real(dp) :: qw_in(this%qw_dim), qw_prime_in(this%qw_dim,3)

  integer, dimension(3) :: shift

#ifdef HAVE_GP
  type(fourier_so4), save :: f_hat
  type(grad_fourier_so4), save :: df_hat
  type(bispectrum_so4), save :: bis
  type(grad_bispectrum_so4), save :: dbis
#endif  

!$omp threadprivate(f_hat,df_hat,bis,dbis)  

  if (present(e)) e = 0.0_dp
  if (present(local_e)) then
     call check_size('Local_E',local_e,(/at%N/),'IPModel_GAP_Calc')
     local_e = 0.0_dp
  endif
  if (present(f)) then 
     call check_size('Force',f,(/3,at%N/),'IPModel_GAP_Calc')
     f = 0.0_dp
  end if
  if (present(virial)) virial = 0.0_dp

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
  if( allocated(this%w_Z) ) then
     do i = 1, at%N
        w(i) = this%w_Z(at%Z(i))
     enddo
  else
     w = 1.0_dp
  endif

  if (trim(this%datafile_coordinates) == 'bispectrum') then

#ifdef HAVE_GP
  d = j_max2d(this%j_max)
  call cg_initialise(this%j_max,2)
#endif

  nei_max = 0
  do i = 1, at%N
     if( nei_max < (atoms_n_neighbours(at,i)+1) ) nei_max = atoms_n_neighbours(at,i)+1
  enddo

  allocate(vec(d,at%N),jack(d,3*nei_max,at%N))

!$omp parallel 
#ifdef HAVE_GP
  call initialise(f_hat,this%j_max,this%z0,this%cutoff)
  if(present(f).or.present(virial)) call initialise(df_hat,this%j_max,this%z0,this%cutoff)
#endif
!$omp do private(n)
  do i = 1, at%N

#ifdef HAVE_GP
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
#endif
  enddo
!$omp end do 
#ifdef HAVE_GP
  call finalise(f_hat)
  call finalise(df_hat)
  call finalise(bis)
  call finalise(dbis)
#endif
!$omp end parallel
    
!$omp parallel do private(k,f_gp,f_gp_k,n,j,jn,shift)
  do i = 1, at%N
     if(present(e) .or. present(local_e)) then
#ifdef HAVE_GP
        call gp_predict(gp_data=this%my_gp, mean=local_e_in(i),x_star=vec(:,i))
#endif
     endif

     if(present(f).or.present(virial)) then
        do k = 1, 3
           f_gp = 0.0_dp
           !kk = (i-1)*3 + k
       
           !call gp_predict(gp_data=this%my_gp, mean=f_gp_k,x_star=vec(:,i),x_prime_star=jack(:,kk,i))
#ifdef HAVE_GP
           call gp_predict(gp_data=this%my_gp, mean=f_gp_k,x_star=vec(:,i),x_prime_star=jack(:,k,i))
#endif
           f_gp = f_gp - f_gp_k
       
           if( present(virial) ) virial_in(:,k,i) = virial_in(:,k,i) - f_gp_k*at%pos(:,i)

           do n = 1, atoms_n_neighbours(at,i)
              j = atoms_neighbour(at,i,n,jn=jn,shift=shift)
       
!              if(jn==i)cycle              
#ifdef HAVE_GP
              call gp_predict(gp_data=this%my_gp,mean=f_gp_k,x_star=vec(:,j),x_prime_star=jack(:,jn*3+k,j))
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

  elseif (trim(this%datafile_coordinates) == 'qw') then

#ifdef HAVE_GP

  do i = 1, at%N
     if (present(f) .or. present(virial)) then
        do j = 1, this%qw_dim
           if (this%qw_type(j) == 1) then
              call calc_qw(at, i, this%qw_l(j), q = qw_in(j), grad_q = qw_prime_in(j,:), d = i, cutoff = this%qw_cutoff(j), &
                           wf_char_length = this%qw_wf_char_length(j), wf_type = this%qw_wf_type(j), do_weight = this%qw_do_weight(j))
           elseif (this%qw_type(j) == 2) then
              call calc_qw(at, i, this%qw_l(j), w = qw_in(j), grad_w = qw_prime_in(j,:), d = i, cutoff = this%qw_cutoff(j), &
                           wf_char_length = this%qw_wf_char_length(j), wf_type = this%qw_wf_type(j), do_weight = this%qw_do_weight(j))
           endif
        enddo
     elseif ((.not. (present(f) .or. present(virial))) .and. (present(e) .or. present(local_e))) then 
        do j = 1, this%qw_dim
           if (this%qw_type(j) == 1) then
              call calc_qw(at, i, this%qw_l(j), q = qw_in(j), cutoff = this%qw_cutoff(j), &
                           wf_char_length = this%qw_wf_char_length(j), wf_type = this%qw_wf_type(j), do_weight = this%qw_do_weight(j))
           elseif (this%qw_type(j) == 2) then
              call calc_qw(at, i, this%qw_l(j), w = qw_in(j), cutoff = this%qw_cutoff(j), &
                           wf_char_length = this%qw_wf_char_length(j), wf_type = this%qw_wf_type(j), do_weight = this%qw_do_weight(j))
           endif
        enddo
     endif

     if (present(e) .or. present(local_e)) local_e_in(i) = gp_mean(this%my_gp, qw_in)

     if (present(f) .or. present(virial)) then
        do j = 1, 3
           f_gp_k = - gp_mean(this%my_gp, qw_in, qw_prime_in(:,j))

           if (present(f)) f(j,i) = f_gp_k
           if (present(virial)) virial_in(:,j,i) = - f_gp_k * at%pos(:,i)
        enddo

        do k = 1, atoms_n_neighbours(at, i)
           l = atoms_neighbour(at, i, k)

           do j = 1, this%qw_dim
              if (this%qw_type(j) == 1) then
                 call calc_qw(at, l, this%qw_l(j), q = qw_in(j), grad_q = qw_prime_in(j,:), d = i, cutoff = this%qw_cutoff(j), &
                              wf_char_length = this%qw_wf_char_length(j), wf_type = this%qw_wf_type(j), do_weight = this%qw_do_weight(j))
              elseif (this%qw_type(j) == 2) then
                 call calc_qw(at, l, this%qw_l(j), w = qw_in(j), grad_w = qw_prime_in(j,:), d = i, cutoff = this%qw_cutoff(j), &
                              wf_char_length = this%qw_wf_char_length(j), wf_type = this%qw_wf_type(j), do_weight = this%qw_do_weight(j))
              endif
           enddo

           do j = 1, 3
              f_gp_k = - gp_mean(this%my_gp, qw_in, qw_prime_in(:,j))

              if (present(f)) f(j,i) = f(j,i) + f_gp_k
              if (present(virial)) virial_in(:,j,i) = virial_in(:,j,i) - f_gp_k * at%pos(:,i)
           enddo
        enddo
     endif
  enddo

#endif

  endif

  if(present(e)) e = sum(local_e_in)
  if(present(local_e)) local_e = local_e_in
  if(present(virial)) virial = sum(virial_in,dim=3)

  if(allocated(local_e_in)) deallocate(local_e_in)
  if(allocated(virial_in)) deallocate(virial_in)
  deallocate(w)

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
