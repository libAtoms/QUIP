! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   libAtoms+QUIP: atomistic simulation library
! HND X
! HND X   Portions of this code were written by
! HND X     Albert Bartok-Partay, Silvia Cereda, Gabor Csanyi, James Kermode,
! HND X     Ivan Solt, Wojciech Szlachta, Csilla Varnai, Steven Winfield.
! HND X
! HND X   Copyright 2006-2010.
! HND X
! HND X   Not for distribution
! HND X
! HND X   Portions of this code were written by Noam Bernstein as part of
! HND X   his employment for the U.S. Government, and are not subject
! HND X   to copyright in the USA.
! HND X
! HND X   When using this software, please cite the following reference:
! HND X
! HND X   http://www.libatoms.org
! HND X
! HND X  Additional contributions by
! HND X    Alessio Comisso, Chiara Gattinoni, and Gianpietro Moras
! HND X
! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

module teach_sparse_mod

  use libatoms_module
  use descriptors_module
  use gp_predict_module
  use gp_teach_module
  use fox_wxml
  use potential_module

  implicit none

  integer, parameter :: SPARSE_LENGTH = 10000
  integer, parameter :: THETA_LENGTH = 10000

  type teach_sparse
     type(Atoms), dimension(:), allocatable :: at
     character(len=STRING_LENGTH) :: at_file='', ip_args = '', &
     energy_parameter_name, force_parameter_name, virial_parameter_name, config_type_parameter_name, config_type_sigma_string
     character(len=10240) :: command_line = ''
     real(dp) :: e0
     real(dp) :: max_cutoff
     real(dp), dimension(3) :: default_sigma
     real(dp) :: sparse_jitter
     logical :: do_core = .false., do_sparse, has_config_type_sigma, do_e0_avg = .true.

     integer :: n_frame, n_coordinate, n_ener, n_force, n_virial, min_save, n_species
     type(extendable_str) :: quip_string
     type(Potential) :: core_pot

     type(gpFull) :: my_gp
     type(gpSparse) :: gp_sp

     type(descriptor), dimension(:), allocatable :: my_descriptor
     character(len=STRING_LENGTH), dimension(99) :: descriptor_str

     real(dp), dimension(:), allocatable :: delta, f0, theta, theta_uniform, zeta
     real(dp), dimension(:,:), allocatable :: sigma
     integer, dimension(:), allocatable :: n_sparseX, sparse_method, target_type, n_cross, n_descriptors, species_Z, covariance_type
     integer, dimension(:,:), allocatable :: config_type_n_sparseX
     character(len=STRING_LENGTH), dimension(:), allocatable :: theta_file, sparse_file, theta_fac_string, config_type, config_type_n_sparseX_string
     logical, dimension(:), allocatable :: mark_sparse_atoms, add_species

     logical :: sparseX_separate_file

  endtype teach_sparse
     
  private

  public :: teach_n_from_xyz
  public :: teach_data_from_xyz
  public :: e0_avg_from_xyz
  public :: w_Z_from_xyz
  public :: teach_sparse
  public :: teach_sparse_print_xml
  public :: file_print_xml
!  public :: print_sparse
  public :: parse_config_type_sigma
  public :: parse_config_type_n_sparseX
  public :: read_teach_xyz
  public :: read_descriptors
  public :: get_species_xyz
  public :: add_multispecies_descriptors

contains

  subroutine read_teach_xyz(this)

    type(teach_sparse), intent(inout) :: this

    type(cinoutput) :: xyzfile
    integer :: n_con

    if( allocated(this%at) ) then
       do n_con = 1, this%n_frame
          call finalise(this%at(n_con))
       enddo
       deallocate(this%at)
       this%n_frame = 0
    endif

    call initialise(xyzfile,this%at_file)
    this%n_frame = xyzfile%n_frame

    allocate(this%at(this%n_frame))

    do n_con = 1, this%n_frame
       call read(xyzfile,this%at(n_con),frame=n_con-1)
       call set_cutoff(this%at(n_con), this%max_cutoff)
       call calc_connect(this%at(n_con))
    enddo

    call finalise(xyzfile)

  endsubroutine read_teach_xyz

  subroutine read_descriptors(this)

    type(teach_sparse), intent(inout) :: this

    integer :: i

    this%max_cutoff = 0.0_dp

    if(allocated(this%my_descriptor)) then
       do i = 1, size(this%my_descriptor)
          call finalise(this%my_descriptor(i))
       enddo
       deallocate(this%my_descriptor)
    endif

    allocate(this%my_descriptor(this%n_coordinate))
    do i = 1, this%n_coordinate
       call initialise(this%my_descriptor(i),this%descriptor_str(i))
       if( this%max_cutoff < cutoff(this%my_descriptor(i)) ) this%max_cutoff = cutoff(this%my_descriptor(i))
    enddo

  endsubroutine read_descriptors

  subroutine teach_n_from_xyz(this)

    type(teach_sparse), intent(inout) :: this

    integer :: n_con
    logical :: has_ener, has_force, has_virial
    real(dp) :: ener, virial(3,3)
    real(dp), pointer :: f(:,:)
    integer :: i
    integer :: n_descriptors, n_cross

    allocate(this%n_cross(this%n_coordinate))
    allocate(this%n_descriptors(this%n_coordinate))

    this%n_cross = 0
    this%n_descriptors = 0
    this%n_ener = 0
    this%n_force = 0
    this%n_virial = 0

    do n_con = 1, this%n_frame

       has_ener = get_value(this%at(n_con)%params,this%energy_parameter_name,ener)
       has_force = assign_pointer(this%at(n_con),this%force_parameter_name, f)
       has_virial = get_value(this%at(n_con)%params,this%virial_parameter_name,virial)

       if( has_ener ) then
          this%n_ener = this%n_ener + 1
       endif

       if( has_force ) then
          this%n_force = this%n_force + this%at(n_con)%N*3
       endif

       if( has_virial ) then
          this%n_virial = this%n_virial + 6
       endif

       do i = 1, this%n_coordinate
          if( has_ener .or. has_force .or. has_virial ) then
             call descriptor_sizes(this%my_descriptor(i),this%at(n_con),n_descriptors,n_cross)
          endif

          if( has_force ) then
             this%n_cross(i) = this%n_cross(i) + n_cross*3
          endif

          if( has_virial ) then
             this%n_cross(i) = this%n_cross(i) + n_cross*6
          endif

          if( has_ener .or. has_force .or. has_virial ) this%n_descriptors(i) = this%n_descriptors(i) + n_descriptors
       enddo

    enddo

  end subroutine teach_n_from_xyz

  subroutine teach_data_from_xyz(this)

    type(teach_sparse), intent(inout) :: this

    type(inoutput) :: theta_inout
    type(descriptor_data) :: my_descriptor_data

    integer :: d
    integer :: n_con
    logical :: has_ener, has_force, has_virial, has_config_type, has_energy_sigma, has_force_sigma, has_virial_sigma
    real(dp) :: ener, ener_core, my_cutoff, energy_sigma, force_sigma, virial_sigma
    real(dp), dimension(3) :: pos
    real(dp), dimension(3,3) :: virial, virial_core
    real(dp), dimension(:), allocatable :: theta, theta_fac
    real(dp), dimension(:,:), pointer :: f
    real(dp), dimension(:,:), allocatable :: f_core
    integer, dimension(:,:), allocatable :: force_loc, permutations
    integer :: ie, i, j, n, k, l, i_coordinate
    integer, dimension(:), allocatable :: xloc
    integer, dimension(3,3) :: virial_loc

    integer :: i_config_type, n_config_type, n_theta_fac
    character(len=STRING_LENGTH) :: config_type
    character(len=THETA_LENGTH) :: theta_string
    character(len=STRING_LENGTH), dimension(:), allocatable :: theta_string_array

    my_cutoff = 0.0_dp
    call gp_setParameters(this%my_gp,this%n_coordinate,this%n_ener,this%n_force+this%n_virial,this%sparse_jitter)

    do i_coordinate = 1, this%n_coordinate
       d = descriptor_dimensions(this%my_descriptor(i_coordinate))

       call gp_setParameters(this%my_gp,i_coordinate, d, this%n_descriptors(i_coordinate), this%n_cross(i_coordinate), this%delta(i_coordinate), this%f0(i_coordinate), &
                             covariance_type=this%covariance_type(i_coordinate) )
       call gp_addDescriptor(this%my_gp,i_coordinate,trim(this%descriptor_str(i_coordinate)))

       allocate(permutations(d,descriptor_n_permutations(this%my_descriptor(i_coordinate))))
       call descriptor_permutations(this%my_descriptor(i_coordinate),permutations)
       call gp_setPermutations(this%my_gp,i_coordinate,permutations)
       deallocate(permutations)

       my_cutoff = max(my_cutoff,cutoff(this%my_descriptor(i_coordinate)))
    enddo

    call print("Number of target energies (property name: "//trim(this%energy_parameter_name)//") found: "//this%n_ener)
    call print("Number of target forces (property name: "//trim(this%force_parameter_name)//") found: "//this%n_force)
    call print("Number of target virials (property name: "//trim(this%virial_parameter_name)//") found: "//this%n_virial)

    if( this%do_core ) call Initialise(this%core_pot, args_str=this%ip_args, param_str=string(this%quip_string))

    do n_con = 1, this%n_frame

       has_ener = get_value(this%at(n_con)%params,this%energy_parameter_name,ener)
       has_force = assign_pointer(this%at(n_con),this%force_parameter_name, f)
       has_virial = get_value(this%at(n_con)%params,this%virial_parameter_name,virial)
       has_config_type = get_value(this%at(n_con)%params,this%config_type_parameter_name,config_type)

       has_energy_sigma = get_value(this%at(n_con)%params,'energy_sigma',energy_sigma)
       has_force_sigma = get_value(this%at(n_con)%params,'force_sigma',force_sigma)
       has_virial_sigma = get_value(this%at(n_con)%params,'virial_sigma',virial_sigma)

       if( has_config_type ) then
          config_type = trim(config_type)
       else
          config_type = "default"
       endif

       if( .not. allocated(this%config_type) ) call system_abort('config_type not allocated')
       n_config_type = 0
       do i_config_type = 1, size(this%config_type)
          if( trim(this%config_type(i_config_type)) == trim(config_type) ) n_config_type = i_config_type
       enddo

       if( n_config_type == 0 ) then ! get the number of the "default" type as default
          do i_config_type = 1, size(this%config_type)
             if( trim(this%config_type(i_config_type)) == "default" ) n_config_type = i_config_type
          enddo
       endif

       if( this%do_core ) then
          allocate(f_core(3,this%at(n_con)%N))
          ener_core = 0.0_dp
          f_core = 0.0_dp
          virial_core = 0.0_dp

          if( this%at(n_con)%cutoff < max(cutoff(this%core_pot),my_cutoff) ) then
             call set_cutoff(this%at(n_con), max(cutoff(this%core_pot),my_cutoff))
             call calc_connect(this%at(n_con))
          endif

          if(has_virial .and. has_force) then
             call calc(this%core_pot,this%at(n_con),energy=ener_core,force=f_core,virial=virial_core)
          elseif(has_force) then
             call calc(this%core_pot,this%at(n_con),energy=ener_core,force=f_core)
          else
             call calc(this%core_pot,this%at(n_con),energy=ener_core)
          end if

          if(has_ener) ener = ener - ener_core
          if(has_force) f = f - f_core
          if(has_virial) virial = virial - virial_core
          deallocate(f_core)
       endif

       if(has_ener) then
          ener = ener - this%at(n_con)%N*this%e0     
       endif

       if( this%at(n_con)%cutoff < my_cutoff ) then
          call set_cutoff(this%at(n_con),my_cutoff)
          call calc_connect(this%at(n_con))
       endif

       if( .not. has_energy_sigma ) energy_sigma = this%sigma(1,n_config_type)*this%at(n_con)%N
       if( .not. has_force_sigma ) force_sigma = this%sigma(2,n_config_type)
       if( .not. has_virial_sigma ) virial_sigma = this%sigma(3,n_config_type)*this%at(n_con)%N

       if( has_ener ) ie = gp_addFunctionValue(this%my_gp,ener, energy_sigma)
       if(has_force) then
          allocate(force_loc(3,this%at(n_con)%N))
          do i = 1, this%at(n_con)%N
             do k = 1, 3
                force_loc(k,i) = gp_addFunctionDerivative(this%my_gp,-f(k,i),force_sigma)
             enddo
          enddo
       endif
       if(has_virial) then
          ! check if virial is symmetric
          if( sum((virial - transpose(virial))**2) .fne. 0.0_dp ) &
          call print_warning('virial not symmetric, now symmetrised')


          ! Now symmetrise matrix
          virial = ( virial + transpose(virial) ) / 2.0_dp

          do k = 1, 3
             do l = k, 3
                virial_loc(l,k) = gp_addFunctionDerivative(this%my_gp,-virial(l,k),virial_sigma)
             enddo
          enddo
       endif

       if(has_ener .or. has_force .or. has_virial ) then
          do i_coordinate = 1, this%n_coordinate

             call calc(this%my_descriptor(i_coordinate),this%at(n_con),my_descriptor_data, &
             do_descriptor=.true.,do_grad_descriptor=has_force .or. has_virial)

             allocate(xloc(size(my_descriptor_data%x)))

             if( has_ener ) then
                do i = 1, size(my_descriptor_data%x)
                   if( .not. my_descriptor_data%x(i)%has_data) cycle
                   xloc(i) = gp_addCoordinates(this%my_gp,my_descriptor_data%x(i)%data(:),i_coordinate, &
                   cutoff_in=my_descriptor_data%x(i)%covariance_cutoff, current_y=ie,config_type=n_config_type)
                enddo
             else
                do i = 1, size(my_descriptor_data%x)
                   if( .not. my_descriptor_data%x(i)%has_data) cycle
                   xloc(i) = gp_addCoordinates(this%my_gp,my_descriptor_data%x(i)%data(:),i_coordinate, &
                   cutoff_in=my_descriptor_data%x(i)%covariance_cutoff, config_type=n_config_type)
                enddo
             endif


             if(has_force) then
                do i = 1, size(my_descriptor_data%x)
                   do n = lbound(my_descriptor_data%x(i)%ii,1), ubound(my_descriptor_data%x(i)%ii,1)
                      if( .not. my_descriptor_data%x(i)%has_grad_data(n)) cycle
                      j = my_descriptor_data%x(i)%ii(n)

                      do k = 1, 3
                         call gp_addCoordinateDerivatives(this%my_gp,my_descriptor_data%x(i)%grad_data(:,k,n),i_coordinate, &
                         force_loc(k,j), xloc(i), dcutoff_in=my_descriptor_data%x(i)%grad_covariance_cutoff(k,n) )
                      enddo
                   enddo
                enddo

             endif

             if(has_virial) then
                do k = 1, 3
                   do l = k, 3

                      do i = 1, size(my_descriptor_data%x)
                         do n = lbound(my_descriptor_data%x(i)%ii,1), ubound(my_descriptor_data%x(i)%ii,1)
                            if( .not. my_descriptor_data%x(i)%has_grad_data(n)) cycle
                            j = my_descriptor_data%x(i)%ii(n)
                            pos = my_descriptor_data%x(i)%pos(:,n)
                            call gp_addCoordinateDerivatives(this%my_gp,my_descriptor_data%x(i)%grad_data(:,k,n)*pos(l), i_coordinate, &
                            virial_loc(l,k), xloc(i), dcutoff_in=my_descriptor_data%x(i)%grad_covariance_cutoff(k,n)*pos(l))
                         enddo
                      enddo

                   enddo
                enddo
             endif

             if(allocated(xloc)) deallocate(xloc)
             call finalise(my_descriptor_data)
          enddo
       endif

       if(allocated(force_loc)) deallocate(force_loc)
    enddo

    do i_coordinate = 1, this%n_coordinate
       if(len_trim(this%theta_file(i_coordinate)) > 0) then
          allocate(theta_string_array(this%my_gp%coordinate(i_coordinate)%d))
          allocate(theta(this%my_gp%coordinate(i_coordinate)%d))

          call initialise(theta_inout,trim(this%theta_file(i_coordinate)))
          read(theta_inout%unit,'(a)') theta_string
          call split_string(theta_string,' :','{}',theta_string_array,d,matching=.true.)
          if(this%my_gp%coordinate(i_coordinate)%d /= d) call system_abort('File '//trim(this%theta_file(i_coordinate))//' does not contain the right number of hyperparameters')
          do i = 1, d
             theta(i) = string_to_real(trim(theta_string_array(i)))
          enddo
          call gp_setTheta(this%my_gp,i_coordinate,theta)
          deallocate(theta_string_array)
          deallocate(theta)
          call finalise(theta_inout)
       elseif(this%theta_uniform(i_coordinate) .fne. 0.0_dp) then
          allocate(theta(this%my_gp%coordinate(i_coordinate)%d))
          theta = this%theta_uniform(i_coordinate)
          call gp_setTheta(this%my_gp,i_coordinate,theta)
          deallocate(theta)
       else
          if(this%covariance_type(i_coordinate) == COVARIANCE_BOND_REAL_SPACE) then
             call gp_setTheta(this%my_gp,i_coordinate,(/ this%theta(i_coordinate) /))
          elseif(this%covariance_type(i_coordinate) == COVARIANCE_DOT_PRODUCT) then
             call gp_setTheta(this%my_gp,i_coordinate,(/ this%zeta(i_coordinate) /))
          elseif( this%covariance_type(i_coordinate) == COVARIANCE_ARD_SE ) then
             allocate(theta_string_array(this%my_gp%coordinate(i_coordinate)%d))
             allocate(theta_fac(this%my_gp%coordinate(i_coordinate)%d))
             call split_string(trim(this%theta_fac_string(i_coordinate))," ",'{}',theta_string_array,n_theta_fac,matching=.true.)

             if(n_theta_fac == 1) then
                theta_fac = string_to_real(theta_string_array(1))
             elseif(n_theta_fac == this%my_gp%coordinate(i_coordinate)%d) then
                do i = 1, this%my_gp%coordinate(i_coordinate)%d
                   theta_fac(i) = string_to_real(theta_string_array(i))
                enddo
             else
                call system_abort("theta_fac can only contain one value or as many as dimensions the descriptor is")
             endif
             call gp_setThetaFactor(this%my_gp,i_coordinate,theta_fac,useSparseX=.false.)
          
             deallocate(theta_fac)
             deallocate(theta_string_array)
          endif
       endif


    enddo

    if( this%do_core ) call Finalise(this%core_pot)

    if( this%do_sparse ) call gp_sparsify(this%my_gp,n_sparseX=this%config_type_n_sparseX,default_all=(this%n_sparseX/=0),sparseMethod=this%sparse_method)

  end subroutine teach_data_from_xyz

  subroutine e0_avg_from_xyz(this)

    type(teach_sparse), intent(inout) :: this

    integer :: n_con
    integer :: n_ener
    logical :: has_ener
    real(dp) :: ener, ener_core

    if( this%do_core ) call Initialise(this%core_pot, this%ip_args, param_str=string(this%quip_string))

    n_ener = 0

    if(this%do_e0_avg) then
       this%e0 = 0.0_dp
    else
       this%e0 = huge(1.0_dp)
    endif

    do n_con = 1, this%n_frame

       has_ener = get_value(this%at(n_con)%params,trim(this%energy_parameter_name),ener)

       if( has_ener ) then

          ener_core = 0.0_dp
          if( this%do_core ) then
             if( this%at(n_con)%cutoff < cutoff(this%core_pot) ) then
                call set_cutoff(this%at(n_con), cutoff(this%core_pot))
                call calc_connect(this%at(n_con))
             endif
             call calc(this%core_pot,this%at(n_con),energy=ener_core)
          endif

          if(this%do_e0_avg) then
             this%e0 = this%e0 + (ener-ener_core) / this%at(n_con)%N
          else
             if( this%e0 < (ener-ener_core) / this%at(n_con)%N ) this%e0 = (ener-ener_core) / this%at(n_con)%N
          endif

          n_ener = n_ener + 1
       endif
    enddo

    if( n_ener > 0 ) then
       if(this%do_e0_avg) this%e0 = this%e0 / n_ener
    else
       this%e0 = 0.0_dp
    endif

    if( this%do_core ) call Finalise(this%core_pot)

  end subroutine e0_avg_from_xyz

  subroutine w_Z_from_xyz(this)

    type(teach_sparse), intent(inout) :: this

    type(cinoutput) :: xyzfile
    type(atoms) :: at

    call initialise(xyzfile,this%at_file)

    call read(xyzfile,at,frame=0)
    !call get_weights(at,this%w_Z)
    call finalise(at)

    call finalise(xyzfile)

  end subroutine w_Z_from_xyz

  subroutine teach_sparse_print_xml(this,filename,sparseX_separate_file)
     type(teach_sparse), intent(in) :: this
     character(len=*), intent(in) :: filename
     logical, intent(in), optional :: sparseX_separate_file

     type(xmlf_t) :: xf
     type(extendable_str) :: gap_string
     character(len=STRING_LENGTH) :: gp_tmp_file, gp_label
     type(inoutput) :: gp_inout
     integer, dimension(8) :: values
     logical :: my_sparseX_separate_file

     call date_and_time(values=values)
     ! Get totally unique label for GAP. This will be used at various places.
     write(gp_label,'("GAP_"7(i0,"_")i0)') values

     ! Unique temporary file
     gp_tmp_file = 'tmp_'//trim(gp_label)//'.xml'

     ! Print GAP part of the potential into the temporary file.
     call xml_OpenFile(gp_tmp_file,xf,addDecl=.false.)

     call xml_NewElement(xf,"GAP_params")
     call xml_AddAttribute(xf,"label",trim(gp_label))
     call xml_AddAttribute(xf,"svn_version",""//current_version())

     call xml_NewElement(xf,"GAP_data")
     call xml_AddAttribute(xf,"do_core",""//this%do_core)
     call xml_AddAttribute(xf,"e0",""//this%e0)

     call xml_EndElement(xf,"GAP_data")

     my_sparseX_separate_file = optional_default(.false., sparseX_separate_file)

     ! Print GP bit of the potential
     if (my_sparseX_separate_file) then
	call gp_printXML(this%gp_sp,xf,label=gp_label,sparseX_base_filename=trim(filename)//".sparseX")
     else
	call gp_printXML(this%gp_sp,xf,label=gp_label)
     endif

     ! Print the command line used for the teaching
     if(len(trim(this%command_line))> 0 ) then
        call xml_NewElement(xf,"command_line")
        call xml_AddCharacters(xf,trim(this%command_line),parsed=.false.)
        call xml_EndElement(xf,"command_line")
     endif

     ! Print the teaching configurations used for this potential.
     if(len(trim(this%at_file)) > 0 ) call file_print_xml(this%at_file,xf)

     call xml_EndElement(xf,"GAP_params")
     call xml_Close(xf)

     ! Now read back into an extendable string what we have just printed out.
     call read(gap_string, trim(gp_tmp_file), keep_lf=.true.)

     ! Initialise the final file
     call initialise(gp_inout,trim(filename),action=OUTPUT)

     ! Open a unique root element for the xml
     call print('<'//trim(gp_label)//'>',file=gp_inout)
     !call system_command('echo "<'//trim(gp_label)//'>" >>'//trim(filename))

     if(this%do_core) then
        ! Create the sum potential xml entry (by hand)
        call print('<Potential label="'//trim(gp_label)//'" init_args="Sum init_args_pot1={'//trim(this%ip_args)//'} init_args_pot2={IP GAP label='//trim(gp_label)//'}"/>',file=gp_inout)
        !call system_command('echo "<Potential label=\"'//trim(gp_label)//'\" init_args=\"Sum init_args_pot1={'//trim(this%ip_args)//'} init_args_pot2={IP GAP label='//trim(gp_label)//'}\"/>" >>'//trim(filename))

        ! Now add the core potential that was used.
        call print(string(this%quip_string),file=gp_inout)
        !call system_command('echo "'//string(this%quip_string)//' >>'//trim(filename))

     endif

     ! Add the GAP potential
     call print(string(gap_string),file=gp_inout)
     !call system_command('cat '//trim(gp_tmp_file)//' >>'//trim(filename))

     ! Close the root element
     call print('</'//trim(gp_label)//'>',file=gp_inout)
     !call system_command('echo "</'//trim(gp_label)//'>" >>'//trim(filename))

     call finalise(gp_inout)
     call finalise(gap_string)

     ! Delete the temporary file
     call system_command('rm -f '//trim(gp_tmp_file))

  endsubroutine teach_sparse_print_xml

  subroutine file_print_xml(this,xf)
     character(len=*), intent(in) :: this
     type(xmlf_t), intent(inout) :: xf

     type(inoutput) :: atfile
     character(len=10240) :: line
     integer :: iostat

     call initialise(atfile,trim(this))
     call xml_NewElement(xf,"XYZ_data")
     call xml_AddNewLine(xf)

     do
        read(atfile%unit,'(a)',iostat=iostat) line
        if(iostat < 0) then
           exit
        elseif(iostat > 0) then
           call system_abort('file_print_xml: unkown error ('//iostat//') while reading '//trim(this))
        endif
        call xml_AddCharacters(xf,trim(line),parsed=.false.)
        call xml_AddNewLine(xf)
     enddo
     call xml_EndElement(xf,"XYZ_data")
     call finalise(atfile)

  endsubroutine file_print_xml

!  subroutine print_sparse(this)
!    type(teach_sparse), intent(in) :: this
!    type(cinoutput) :: xyzfile, xyzfile_out
!    type(atoms) :: at, at_out
!
!    integer :: li, ui, n_con
!    logical, dimension(:), allocatable :: x
!    logical, dimension(:), pointer :: sparse
!
!    if(this%do_mark_sparse_atoms) then
!
!       allocate(x(this%n_descriptors))
!       x = .false.
!       x(this%r) = .true.
!
!       call initialise(xyzfile,this%at_file)
!       call initialise(xyzfile_out,this%mark_sparse_atoms,action=OUTPUT)
!
!       li = 0
!       ui = 0
!       do n_con = 1, xyzfile%n_frame
!          call read(xyzfile,at,frame=n_con-1)
!          at_out = at
!
!          call add_property(at_out,'sparse',.false.,ptr=sparse)
!
!          li = ui + 1
!          ui = ui + at%N
!          if(any( x(li:ui) )) sparse(find_indices(x(li:ui))) = .true.
!
!          call write(at_out,xyzfile_out,properties="species:pos:sparse")
!       enddo
!       call finalise(xyzfile)
!       call finalise(xyzfile_out)
!       deallocate(x)
!
!    endif
!
!  endsubroutine print_sparse

  subroutine parse_config_type_sigma(this)
    type(teach_sparse), intent(inout) :: this
    character(len=STRING_LENGTH), dimension(99) :: config_type_sigma_fields
    integer :: config_type_sigma_num_fields, i_default, i, n_config_type

    if( this%has_config_type_sigma ) then
       call split_string(this%config_type_sigma_string,':','{}',config_type_sigma_fields,config_type_sigma_num_fields,matching=.true.)

       n_config_type = config_type_sigma_num_fields / 4

       ! find "default" if present
       i_default = 0
       do i = 1, config_type_sigma_num_fields, 4
          if( trim(config_type_sigma_fields(i)) == "default" ) i_default = i
       enddo

       if( i_default == 0 ) then
          ! no default present in the string, we add it, and it'll be the last one
          n_config_type = n_config_type + 1
          i_default = n_config_type
          config_type_sigma_fields(config_type_sigma_num_fields+1) = "default"
          config_type_sigma_fields(config_type_sigma_num_fields+2) = ""//this%default_sigma(1)
          config_type_sigma_fields(config_type_sigma_num_fields+3) = ""//this%default_sigma(2)
          config_type_sigma_fields(config_type_sigma_num_fields+4) = ""//this%default_sigma(3)
          config_type_sigma_num_fields = config_type_sigma_num_fields + 4
       endif

       allocate(this%config_type(n_config_type))
       allocate(this%sigma(3,n_config_type))

       do i = 1, n_config_type 
          this%config_type(i) = trim(config_type_sigma_fields(4*(i-1)+1))
          this%sigma(1,i) = string_to_real(config_type_sigma_fields(4*(i-1)+2))
          this%sigma(2,i) = string_to_real(config_type_sigma_fields(4*(i-1)+3))
          this%sigma(3,i) = string_to_real(config_type_sigma_fields(4*(i-1)+4))
       enddo

       call print('Sparse points and target errors per pre-defined types of configurations')
       do i = 1, n_config_type
          call print(""//trim(this%config_type(i))//"  "//this%sigma(:,i))
       enddo
    else
       allocate(this%config_type(1))
       allocate(this%sigma(3,1))
       this%config_type(1)= "default"
       this%sigma(:,1) = this%default_sigma
    endif

  endsubroutine parse_config_type_sigma

  subroutine parse_config_type_n_sparseX(this)
    type(teach_sparse), intent(inout) :: this

    integer :: i, j, i_default, i_coordinate, i_config_type, config_type_n_sparseX_num_fields, n_config_type, new_config_types
    character(len=STRING_LENGTH), dimension(99) :: config_type_n_sparseX_fields
    logical :: config_type_present

    if( .not. allocated(this%config_type) ) call system_abort('config_type not allocated, call parse_config_type_sigma first')

    do i = 1, size(this%config_type)
       if( trim(this%config_type(i)) == "default" ) i_default = i
    enddo

    ! Check first if we have more new config types than we had from config_type_sigma
    do i_coordinate = 1, this%n_coordinate
       if( this%n_sparseX(i_coordinate) == 0 .and. len_trim(this%config_type_n_sparseX_string(i_coordinate)) > 0) then
          call split_string(this%config_type_n_sparseX_string(i_coordinate),':','{}',config_type_n_sparseX_fields,config_type_n_sparseX_num_fields,matching=.true.)

          n_config_type = size(this%config_type)
          new_config_types = 0 ! Assume there are no new config_types
          do j = 1, config_type_n_sparseX_num_fields, 2 ! loop over config_types in the descriptor string
             config_type_present = .false.
             do i = 1, n_config_type ! loop over config_types previously set
                if( trim(this%config_type(i)) == trim(config_type_n_sparseX_fields(j)) ) config_type_present = .true. ! Found config_type among old ones
             enddo
             if(.not.config_type_present) new_config_types = new_config_types + 1 ! Increment as it's a genuine new config_type
          enddo
          if( new_config_types > 0 ) then
             call reallocate(this%config_type, n_config_type + new_config_types, copy=.true.)
             call reallocate(this%sigma,3,n_config_type + new_config_types, copy=.true.)

             i_config_type = n_config_type
             do j = 1, config_type_n_sparseX_num_fields, 2 ! loop over config_types in the descriptor string
                config_type_present = .false.
                do i = 1, n_config_type ! loop over config_types previously set
                   if( trim(this%config_type(i)) == trim(config_type_n_sparseX_fields(j)) ) config_type_present = .true. ! Found config_type among old ones
                enddo
                if(.not.config_type_present) then ! it's a genuine new config_type
                   i_config_type = i_config_type + 1
                   this%config_type(i_config_type) = trim(config_type_n_sparseX_fields(j))
                   this%sigma(:,i_config_type) = this%sigma(:,i_default)
                endif
             enddo
          endif

       elseif(this%n_sparseX(i_coordinate) > 0 .and. len_trim(this%config_type_n_sparseX_string(i_coordinate)) > 0) then
          call system_abort('Confused: cannot specify both n_sparseX and config_type_n_sparseX')


       elseif(this%n_sparseX(i_coordinate) == 0 .and. len_trim(this%config_type_n_sparseX_string(i_coordinate)) == 0) then
          call system_abort('Confused: either n_sparseX or config_type_n_sparse has to be specified')
       endif

    enddo

    n_config_type = size(this%config_type)
    allocate(this%config_type_n_sparseX(n_config_type,this%n_coordinate))
    this%config_type_n_sparseX = 0

    do i_coordinate = 1, this%n_coordinate
       if( this%n_sparseX(i_coordinate) == 0 .and. len_trim(this%config_type_n_sparseX_string(i_coordinate)) > 0) then
          call split_string(this%config_type_n_sparseX_string(i_coordinate),':','{}',config_type_n_sparseX_fields,config_type_n_sparseX_num_fields,matching=.true.)

          do j = 1, config_type_n_sparseX_num_fields, 2 ! loop over config_types in the descriptor string
             do i = 1, n_config_type ! loop over config_types previously set
                if( trim(this%config_type(i)) == trim(config_type_n_sparseX_fields(j)) ) &
                   this%config_type_n_sparseX(i,i_coordinate) = string_to_int( config_type_n_sparseX_fields(j+1) )
             enddo
          enddo
          !this%n_sparseX(i_coordinate) = sum( this%config_type_n_sparseX(:,i_coordinate) )

       elseif( this%n_sparseX(i_coordinate) > 0 .and. len_trim(this%config_type_n_sparseX_string(i_coordinate)) == 0) then
          this%config_type_n_sparseX(i_default,i_coordinate) = this%n_sparseX(i_coordinate)
       endif
    enddo

  endsubroutine parse_config_type_n_sparseX

  subroutine get_species_xyz(this)
    type(teach_sparse), intent(inout) :: this

    integer :: n_con, i
    integer, dimension(116) :: species_present

    this%n_species = 0
    species_present = 0

    do n_con = 1, this%n_frame
       do i = 1, this%at(n_con)%N
          if( all(this%at(n_con)%Z(i) /= species_present) ) then
             this%n_species = this%n_species + 1
             species_present(this%n_species) = this%at(n_con)%Z(i)
          endif
       enddo
    enddo

    allocate(this%species_Z(this%n_species))
    this%species_Z = species_present(1:this%n_species)
    
  endsubroutine get_species_xyz

  subroutine add_multispecies_descriptors(this)
    type(teach_sparse), intent(inout) :: this

    integer :: i_coordinate, i, n_descriptor_str, i_add_species
    character(STRING_LENGTH), dimension(:), allocatable :: descriptor_str_i, new_descriptor_str

    ! temporary arrays
    real(dp), dimension(:), allocatable :: delta, f0, theta_uniform, theta, zeta
    integer, dimension(:), allocatable :: n_sparseX, sparse_method, covariance_type
    character(len=STRING_LENGTH), dimension(:), allocatable :: theta_file, sparse_file, theta_fac_string, config_type_n_sparseX_string
    logical, dimension(:), allocatable :: mark_sparse_atoms

    n_descriptor_str = 0
    do i_coordinate = 1, this%n_coordinate
       if( this%add_species(i_coordinate) ) then

          call print('Old descriptor: {'//trim(this%descriptor_str(i_coordinate))//'}')
          call descriptor_str_add_species(this%descriptor_str(i_coordinate),this%species_Z,descriptor_str_i)
          call reallocate(new_descriptor_str, n_descriptor_str+size(descriptor_str_i),copy=.true.)

          call reallocate(delta, n_descriptor_str+size(descriptor_str_i),copy=.true.)
          call reallocate(f0, n_descriptor_str+size(descriptor_str_i),copy=.true.)
          call reallocate(n_sparseX, n_descriptor_str+size(descriptor_str_i),copy=.true.)
          call reallocate(config_type_n_sparseX_string, n_descriptor_str+size(descriptor_str_i),copy=.true.)
          call reallocate(sparse_method, n_descriptor_str+size(descriptor_str_i),copy=.true.)
          call reallocate(theta_fac_string, n_descriptor_str+size(descriptor_str_i),copy=.true.)
          call reallocate(theta_uniform, n_descriptor_str+size(descriptor_str_i),copy=.true.)
          call reallocate(theta_file, n_descriptor_str+size(descriptor_str_i),copy=.true.)
          call reallocate(sparse_file, n_descriptor_str+size(descriptor_str_i),copy=.true.)
          call reallocate(mark_sparse_atoms, n_descriptor_str+size(descriptor_str_i),copy=.true.)
          call reallocate(covariance_type, n_descriptor_str+size(descriptor_str_i),copy=.true.)
          call reallocate(theta, n_descriptor_str+size(descriptor_str_i),copy=.true.)
          call reallocate(zeta, n_descriptor_str+size(descriptor_str_i),copy=.true.)

          do i = 1, size(descriptor_str_i)
             i_add_species = index(descriptor_str_i(i),'add_species')
             if(i_add_species /= 0) descriptor_str_i(i)(i_add_species:i_add_species+len('add_species')-1) = '           '

             new_descriptor_str(i+n_descriptor_str) = trim(descriptor_str_i(i))
             call print('New descriptor: {'//trim(descriptor_str_i(i))//'}')

             delta(i+n_descriptor_str) = this%delta(i_coordinate)
             f0(i+n_descriptor_str) = this%f0(i_coordinate)
             n_sparseX(i+n_descriptor_str) = this%n_sparseX(i_coordinate)
             config_type_n_sparseX_string(i+n_descriptor_str) = this%config_type_n_sparseX_string(i_coordinate)
             sparse_method(i+n_descriptor_str) = this%sparse_method(i_coordinate)
             theta_fac_string(i+n_descriptor_str) = this%theta_fac_string(i_coordinate)
             theta_uniform(i+n_descriptor_str) = this%theta_uniform(i_coordinate)
             theta_file(i+n_descriptor_str) = this%theta_file(i_coordinate)
             sparse_file(i+n_descriptor_str) = this%sparse_file(i_coordinate)
             mark_sparse_atoms(i+n_descriptor_str) = this%mark_sparse_atoms(i_coordinate)
             covariance_type(i+n_descriptor_str) = this%covariance_type(i_coordinate)
             theta(i+n_descriptor_str) = this%theta(i_coordinate)
             zeta(i+n_descriptor_str) = this%zeta(i_coordinate)

          enddo
          n_descriptor_str = n_descriptor_str + size(descriptor_str_i)
          deallocate(descriptor_str_i)

       else
          n_descriptor_str = n_descriptor_str + 1

          call reallocate(new_descriptor_str, n_descriptor_str,copy=.true.)
          call reallocate(delta, n_descriptor_str,copy=.true.)
          call reallocate(f0, n_descriptor_str,copy=.true.)
          call reallocate(n_sparseX, n_descriptor_str,copy=.true.)
          call reallocate(config_type_n_sparseX_string, n_descriptor_str,copy=.true.)
          call reallocate(sparse_method, n_descriptor_str,copy=.true.)
          call reallocate(theta_fac_string, n_descriptor_str,copy=.true.)
          call reallocate(theta_uniform, n_descriptor_str,copy=.true.)
          call reallocate(theta_file, n_descriptor_str,copy=.true.)
          call reallocate(sparse_file, n_descriptor_str,copy=.true.)
          call reallocate(mark_sparse_atoms, n_descriptor_str,copy=.true.)
          call reallocate(covariance_type, n_descriptor_str,copy=.true.)
          call reallocate(theta, n_descriptor_str,copy=.true.)
          call reallocate(zeta, n_descriptor_str,copy=.true.)

          new_descriptor_str(n_descriptor_str) = trim(this%descriptor_str(i_coordinate))
          delta(n_descriptor_str) = this%delta(i_coordinate)
          f0(n_descriptor_str) = this%f0(i_coordinate)
          n_sparseX(n_descriptor_str) = this%n_sparseX(i_coordinate)
          config_type_n_sparseX_string(n_descriptor_str) = this%config_type_n_sparseX_string(i_coordinate)
          sparse_method(n_descriptor_str) = this%sparse_method(i_coordinate)
          theta_fac_string(n_descriptor_str) = this%theta_fac_string(i_coordinate)
          theta_uniform(n_descriptor_str) = this%theta_uniform(i_coordinate)
          theta_file(n_descriptor_str) = this%theta_file(i_coordinate)
          sparse_file(n_descriptor_str) = this%sparse_file(i_coordinate)
          mark_sparse_atoms(n_descriptor_str) = this%mark_sparse_atoms(i_coordinate)
          covariance_type(n_descriptor_str) = this%covariance_type(i_coordinate)
          theta(n_descriptor_str) = this%theta(i_coordinate)
          zeta(n_descriptor_str) = this%zeta(i_coordinate)

          call print('Unchanged descriptor: {'//trim(this%descriptor_str(i_coordinate))//'}')
       endif

    enddo
    call reallocate(this%delta, n_descriptor_str)
    call reallocate(this%f0, n_descriptor_str)
    call reallocate(this%n_sparseX, n_descriptor_str)
    call reallocate(this%config_type_n_sparseX_string, n_descriptor_str)
    call reallocate(this%sparse_method, n_descriptor_str)
    call reallocate(this%theta_fac_string, n_descriptor_str)
    call reallocate(this%theta_uniform, n_descriptor_str)
    call reallocate(this%theta_file, n_descriptor_str)
    call reallocate(this%sparse_file, n_descriptor_str)
    call reallocate(this%mark_sparse_atoms, n_descriptor_str)
    call reallocate(this%covariance_type, n_descriptor_str)
    call reallocate(this%theta, n_descriptor_str)
    call reallocate(this%zeta, n_descriptor_str)

    this%descriptor_str(1:n_descriptor_str) = new_descriptor_str
    this%delta = delta
    this%f0 = f0
    this%n_sparseX = n_sparseX
    this%config_type_n_sparseX_string = config_type_n_sparseX_string
    this%sparse_method = sparse_method
    this%theta_fac_string = theta_fac_string
    this%theta_uniform = theta_uniform
    this%theta_file = theta_file
    this%sparse_file = sparse_file
    this%mark_sparse_atoms = mark_sparse_atoms
    this%covariance_type = covariance_type
    this%theta = theta
    this%zeta = zeta

    this%n_coordinate = n_descriptor_str

    if(allocated(delta)) deallocate(delta)
    if(allocated(f0)) deallocate(f0)
    if(allocated(n_sparseX)) deallocate(n_sparseX)
    if(allocated(config_type_n_sparseX_string)) deallocate(config_type_n_sparseX_string)
    if(allocated(sparse_method)) deallocate(sparse_method)
    if(allocated(theta_fac_string)) deallocate(theta_fac_string)
    if(allocated(theta_uniform)) deallocate(theta_uniform)
    if(allocated(theta_file)) deallocate(theta_file)
    if(allocated(sparse_file)) deallocate(sparse_file)
    if(allocated(mark_sparse_atoms)) deallocate(mark_sparse_atoms)
    if(allocated(covariance_type)) deallocate(covariance_type)
    if(allocated(theta)) deallocate(theta)
    if(allocated(zeta)) deallocate(zeta)

  endsubroutine add_multispecies_descriptors

end module teach_sparse_mod
