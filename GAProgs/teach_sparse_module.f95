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

  type sparse_types
     character(len=FIELD_LENGTH) :: type
     integer :: m
  endtype sparse_types

  type teach_sparse
     character(len=FIELD_LENGTH) :: at_file='', ip_args = '', &
     energy_property_name, force_property_name, virial_property_name, coordinates
     character(len=10240) :: command_line = ''
     real(dp) :: r_cut, e0, z0, f0, dlt, theta_fac
     real(dp), dimension(3) :: sgm
     logical :: do_core = .false., &
     qw_no_q, qw_no_w, do_sigma, do_delta, do_theta, do_sparx, do_f0, &
     do_theta_fac, do_test_gp_gradient, do_cluster, do_pivot, do_sparse, &
     has_m_sparse_in_type, do_pca

     integer :: d, m, j_max, qw_l_max, n, nn, ne, n_ener, n_force, n_virial, min_steps, min_save, n_species, &
     qw_f_n
     type(extendable_str) :: quip_string
     type(gp) :: my_gp

     real(dp), dimension(99) :: qw_cutoff, qw_cutoff_r1
     real(dp), dimension(:), allocatable :: w_Z, yf, ydf, dlta, pca_mean
     real(dp), dimension(:,:), allocatable :: x, xd, theta, pca_matrix
     integer, dimension(:), allocatable :: lf, ldf, xf, xdf, xz, target_type, r, species_Z
     integer, dimension(99) :: qw_cutoff_f
     type(sparse_types), dimension(:), allocatable :: m_sparse_in_type

  endtype teach_sparse
     
  private

  public :: teach_n_from_xyz
  public :: teach_data_from_xyz
  public :: e0_avg_from_xyz
  public :: w_Z_from_xyz
  public :: teach_sparse
  public :: teach_sparse_print_xml
  public :: file_print_xml

contains

  subroutine teach_n_from_xyz(this)

    type(teach_sparse), intent(inout) :: this

    type(cinoutput) :: xyzfile
    type(atoms) :: at
    integer :: n_max, n_con
    logical :: has_ener, has_force, has_virial
    real(dp) :: ener, virial(3,3)
    real(dp), pointer :: f(:,:)
    integer :: i
    integer, dimension(116) :: species_present

    call initialise(xyzfile,this%at_file)

    n_max = xyzfile%n_frame

    this%n = 0
    this%nn = 0
    this%ne = 0
    this%n_ener = 0
    this%n_force = 0
    this%n_virial = 0
    species_present = 0
    this%n_species = 0

    do n_con = 1, n_max
       call read(xyzfile,at,frame=n_con-1)

       has_ener = get_value(at%params,this%energy_property_name,ener)
       has_force = assign_pointer(at,this%force_property_name, f)
       has_virial = get_value(at%params,this%virial_property_name,virial)

       if( has_ener .or. has_force .or. has_virial ) then
          call set_cutoff(at,this%r_cut)
          call calc_connect(at)
       endif

       if( has_ener ) then
          this%n_ener = this%n_ener + 1
          this%ne = this%ne + at%N
       endif

       if( has_force ) then
          this%n_force = this%n_force + at%N*3
          do i = 1, at%N
             this%n = this%n + 3*(atoms_n_neighbours(at,i)+1)
          enddo
       endif

       if( has_virial ) then
          this%n_virial = this%n_virial + 6
          do i = 1, at%N
             this%n = this%n + 6*(atoms_n_neighbours(at,i)+1)
          enddo
       endif

       if( has_ener .or. has_force .or. has_virial ) this%nn = this%nn + at%N

       do i = 1, at%N
          if( all(at%Z(i) /= species_present) ) then
             this%n_species = this%n_species + 1
             species_present(this%n_species) = at%Z(i)
          endif
       enddo

       call finalise(at)
    enddo

    call finalise(xyzfile)

    select case(trim(this%coordinates))
    case('water_monomer')
       this%n_species = 1
       allocate(this%species_Z(this%n_species))
       this%species_Z = 8
       if(3*this%n_ener /= this%nn) call system_abort('coordinates type water_monomer, but 3*n_ener /= ne')
       this%nn = this%n_ener
       this%ne = this%n_ener
    case('water_dimer')
       this%n_species = 1
       allocate(this%species_Z(this%n_species))
       this%species_Z = 8
       if(6*this%n_ener /= this%nn) call system_abort('coordinates type water_dimer, but 6*n_ener /= ne')
       this%nn = this%n_ener
       this%ne = this%n_ener
    case default
       allocate(this%species_Z(this%n_species))
       this%species_Z = species_present(1:this%n_species)
    endselect

  end subroutine teach_n_from_xyz

  subroutine teach_data_from_xyz(this)

    type(teach_sparse), intent(inout) :: this

    type(cinoutput) :: xyzfile
    type(atoms) :: at
    type(fourier_so4) :: f_hat
    type(grad_fourier_so4) :: df_hat
    type(bispectrum_so4) :: bis
    type(grad_bispectrum_so4) :: dbis
    type(fourier_so3) :: f3_hat
    type(grad_fourier_so3) :: df3_hat
    type(qw_so3) :: qw
    type(grad_qw_so3) :: dqw

    type(Potential) :: core_pot
    
    integer :: d
    integer :: n_max, n_con
    logical :: has_ener, has_force, has_virial
    real(dp) :: ener, ener_core
    real(dp), dimension(3,3) :: virial, virial_core
    real(dp), pointer :: f(:,:)
    real(dp), dimension(:,:), allocatable :: f_core
    real(dp), allocatable :: vec(:,:), jack(:,:,:), w(:)
    integer :: shift(3)
    integer, dimension(3,1) :: water_monomer_index
    integer, dimension(3,2) :: water_dimer_index
    real(dp), dimension(3) :: pos, water_monomer_v
    real(dp), dimension(WATER_DIMER_D) :: water_dimer_v
    integer :: li, ui, nn, ix, ie, i_con, i_ener, w_con
    integer :: nei_max
    integer :: i, j, n, k, l, jx, jn

    select case(trim(this%coordinates))
    case('qw')

       call initialise(f3_hat, this%qw_l_max, &
       this%qw_cutoff(1:this%qw_f_n), this%qw_cutoff_f(1:this%qw_f_n), &
       this%qw_cutoff_r1(1:this%qw_f_n))
       call initialise(df3_hat, this%qw_l_max, &
       this%qw_cutoff(1:this%qw_f_n), this%qw_cutoff_f(1:this%qw_f_n), &
       this%qw_cutoff_r1(1:this%qw_f_n))
       call initialise(qw, this%qw_l_max, &
       this%qw_f_n, do_q = (.not. this%qw_no_q), do_w = (.not. this%qw_no_w))
       call initialise(dqw, this%qw_l_max, &
       this%qw_f_n, do_q = (.not. this%qw_no_q), do_w = (.not. this%qw_no_w))

       this%d = qw2d(qw)
    case('bispectrum')
       call initialise(f_hat,this%j_max, this%z0,this%r_cut)
       call initialise(df_hat,this%j_max, this%z0,this%r_cut)

       this%d = j_max2d(f_hat%j_max)

    case('water_monomer')
       this%d = 3
    case('water_dimer')
       this%d = WATER_DIMER_D
    case default
       call system_abort('Unknown coordinates '//trim(this%coordinates))
    endselect
    d = this%d

    allocate(this%x(d,this%nn),this%xd(d,this%n), &
    this%yf(this%n_ener), this%ydf(this%n_force+this%n_virial), &
    this%lf(this%n_ener),this%ldf(this%n_force+this%n_virial), &
    this%xf(this%ne),this%xdf(this%n))
    allocate(this%xz(this%nn))
    allocate(this%target_type(this%n_ener+this%n_force+this%n_virial))

    if( this%do_core ) call Initialise(core_pot, args_str=this%ip_args, param_str=string(this%quip_string))

    call initialise(xyzfile,this%at_file)

    n_max = xyzfile%n_frame

    li = 0
    ui = 0
    nn = 0
    ix = 0
    ie = 0
    i_con = 0
    i_ener = 0
    w_con = 0

    do n_con = 1, n_max
       call read(xyzfile,at,frame=n_con-1)

       has_ener = get_value(at%params,this%energy_property_name,ener)
       has_force = assign_pointer(at,this%force_property_name, f)
       has_virial = get_value(at%params,this%virial_property_name,virial)

       if( this%do_core ) then
          allocate(f_core(3,at%N))
          ener_core = 0.0_dp
          f_core = 0.0_dp
          virial_core = 0.0_dp

          call set_cutoff(at, max(cutoff(core_pot),this%r_cut))
          call calc_connect(at)

          call calc(core_pot,at,energy=ener_core,force=f_core,virial=virial_core)

          if(has_ener) ener = ener - ener_core
          if(has_force) f = f - f_core
          if(has_virial) virial = virial - virial_core
          deallocate(f_core)
       endif

       if(has_ener) then
          select case(trim(this%coordinates))
          case('water_monomer','water_dimer')
             ener = ener - this%e0     
          case default
             ener = ener - at%N*this%e0     
          endselect
       endif

       call set_cutoff(at,this%r_cut)
       call calc_connect(at)

       nei_max = 0
       do i = 1, at%N
          if(nei_max < atoms_n_neighbours(at,i)+1) nei_max = atoms_n_neighbours(at,i)+1
       enddo

       if(has_ener .or. has_force .or. has_virial ) allocate(vec(d,at%N))
       if(has_force .or. has_virial) then
          allocate(jack(d,3*nei_max,at%N))
          jack = 0.0_dp
       endif
       allocate(w(at%N))

       do i = 1, at%N
          w(i) = this%w_Z(at%Z(i))
       enddo

       if(has_ener .or. has_force .or. has_virial ) then
          select case(trim(this%coordinates))
          case('water_monomer')
             if( at%N /= 3 ) call system_abort('Number of atoms is '//at%N//', not a single water molecule')
             call find_water_monomer(at,water_monomer_index)
             water_monomer_v = water_monomer(at,water_monomer_index(:,1))

             w_con = w_con + 1
             this%x(:,w_con) = water_monomer_v
             this%xz(w_con) = 8
             this%xf(w_con) = w_con
             this%yf(w_con) = ener
             this%lf(w_con) = w_con
             this%target_type(w_con) = 1

          case('water_dimer')
             if( at%N /= 6 ) call system_abort('Number of atoms is '//at%N//', not two water molecules')
             call find_water_monomer(at,water_dimer_index)
             water_dimer_v = water_dimer(at,water_dimer_index(:,1),water_dimer_index(:,2),this%r_cut)

             w_con = w_con + 1
             this%x(:,w_con) = water_dimer_v
             this%xz(w_con) = 8
             this%xf(w_con) = w_con
             this%yf(w_con) = ener
             this%lf(w_con) = w_con
             this%target_type(w_con) = 1

          case('qw')
             do i = 1, at%N
                call fourier_transform(f3_hat,at,i)
                call calc_qw(qw,f3_hat)
                call qw2vec(qw,vec(:,i))
                if(has_force .or. has_virial) then
                   do n = 0, atoms_n_neighbours(at,i)
                      call fourier_transform(df3_hat,at,i,n)
                      call calc_qw(dqw,f3_hat,df3_hat)
                      call qw2vec(dqw,jack(:,3*n+1:3*(n+1),i))
                   enddo
                endif
             enddo
          case('bispectrum')
             do i = 1, at%N
                call fourier_transform(f_hat,at,i,w)
                call calc_bispectrum(bis,f_hat)
                call bispectrum2vec(bis,vec(:,i))
                if(has_force .or. has_virial) then
                   do n = 0, atoms_n_neighbours(at,i)
                      call fourier_transform(df_hat,at,i,n,w)
                      call calc_bispectrum(dbis,f_hat,df_hat)
                      call bispectrum2vec(dbis,jack(:,3*n+1:3*(n+1),i))
                   enddo
                endif
             enddo
          case default
             call system_abort('Unknown coordinates '//trim(this%coordinates))
          endselect

          select case(trim(this%coordinates))
          case('water_monomer','water_dimer')

          case default
             do i = 1, at%N
                ix = i_con + i
                this%x(:,ix) = vec(:,i)
                this%xz(ix) = at%Z(i)

                if( has_ener ) then
                   ie = ie + 1
                   this%xf(ie) = ix
                endif
           
                if(has_force) then
                   do k = 1, 3
                      li = ui + 1
                      ui = ui + 1
                      nn = nn+1
                      this%xdf(ui) = ix
                      this%xd(:,ui) = jack(:,k,i)
                      do n = 1, atoms_n_neighbours(at,i)
                         j = atoms_neighbour(at,i,n,jn=jn)
                         ui = ui + 1
                         jx = i_con + j
                         this%xdf(ui) = jx
                         this%xd(:,ui) = jack(:,jn*3+k,j)
                      enddo

                      this%ydf(nn) = -f(k,i)
                      this%ldf(nn) = ui
                      !sigma(nn+n_ener) = sgm(2)
                      this%target_type(nn+this%n_ener) = 2
                   enddo
                endif
             enddo

             if(has_virial) then
                ! check if virial is symmetric
                if( sum((virial - transpose(virial))**2) .fne. 0.0_dp ) &
                call print_warning('virial not symmetric, now symmetrised')


                ! Now symmetrise matrix
                virial = ( virial + transpose(virial) ) / 2.0_dp

                do k = 1, 3
                   do l = k, 3
                      nn = nn+1

                      do i = 1, at%N
                         ix = i_con + i

                         ui = ui + 1
                         this%xdf(ui) = ix
                         this%xd(:,ui) = jack(:,k,i)*at%pos(l,i)
                         do n = 1, atoms_n_neighbours(at,i)
                            j = atoms_neighbour(at,i,n,jn=jn,shift=shift)
                            ui = ui + 1
                            jx = i_con + j
                            this%xdf(ui) = jx
                            pos = at%pos(:,i) - matmul(at%lattice,shift)
                            this%xd(:,ui) = jack(:,jn*3+k,j)*pos(l)
                         enddo
                      enddo

                      this%ydf(nn) = -virial(l,k) 
                      this%ldf(nn) = ui
                      !sigma(nn+n_ener) = sgm(3)
                      this%target_type(nn+this%n_ener) = 3
                   enddo
                enddo
             endif

             i_con = i_con + at%N

             if(has_ener) then
                i_ener = i_ener + 1
                this%yf(i_ener) = ener
                this%lf(i_ener) = ie
                !sigma(i_ener) = sgm(1)
                this%target_type(i_ener) = 1
             endif
          endselect
       endif

       if(allocated(vec)) deallocate(vec)
       if(allocated(jack)) deallocate(jack)
       if(allocated(w)) deallocate(w)

       call finalise(at)
    enddo

    call finalise(xyzfile)

    if( this%do_core ) then
       call Finalise(core_pot)
    endif

  end subroutine teach_data_from_xyz

  subroutine e0_avg_from_xyz(this)

    type(teach_sparse), intent(inout) :: this

    type(Potential) :: core_pot
    type(cinoutput) :: xyzfile
    type(atoms) :: at
    integer :: n_max, n_con
    integer :: n_ener
    logical :: has_ener
    real(dp) :: ener, ener_core
    integer :: i

    if( this%do_core ) call Initialise(core_pot, this%ip_args, param_str=string(this%quip_string))

    call initialise(xyzfile,this%at_file)

    n_max = xyzfile%n_frame

    n_ener = 0
    this%e0 = 0.0_dp

    do n_con = 1, n_max
       call read(xyzfile,at,frame=n_con-1)

       has_ener = get_value(at%params,'Energy',ener)

       if( has_ener ) then

          ener_core = 0.0_dp
          if( this%do_core ) then
             call set_cutoff(at, cutoff(core_pot))
             call calc_connect(at)
             call calc(core_pot,at,energy=ener_core)
          endif

          select case(trim(this%coordinates))
          case('water_monomer','water_dimer')
             this%e0 = this%e0 + (ener-ener_core)
          case default
             this%e0 = this%e0 + (ener-ener_core) / at%N
          endselect

          n_ener = n_ener + 1
       endif

       call finalise(at)
    enddo

    if( n_ener > 0 ) this%e0 = this%e0 / n_ener

    call finalise(xyzfile)

    if( this%do_core ) call Finalise(core_pot)

  end subroutine e0_avg_from_xyz

  subroutine w_Z_from_xyz(this)

    type(teach_sparse), intent(inout) :: this

    type(cinoutput) :: xyzfile
    type(atoms) :: at

    allocate(this%w_Z(maxval(this%species_Z)))
    select case(trim(this%coordinates))
    case('water_monomer','water_dimer')
       this%w_Z(8) = 1.0_dp
    case default
       call initialise(xyzfile,this%at_file)

       call read(xyzfile,at,frame=0)
       call get_weights(at,this%w_Z)
       call finalise(at)

       call finalise(xyzfile)
    endselect

  end subroutine w_Z_from_xyz

  subroutine teach_sparse_print_xml(this,filename)
     type(teach_sparse), intent(in) :: this
     character(len=*), intent(in) :: filename
     type(xmlf_t) :: xf
     integer :: i
     type(extendable_str) :: gap_string
     character(len=FIELD_LENGTH) :: gp_tmp_file, gp_label
     type(inoutput) :: gp_inout
     integer, dimension(8) :: values

     call date_and_time(values=values)
     ! Get totally unique label for GAP. This will be used at various places.
     write(gp_label,'("GAP_"7(i0,"_")i0)') values

     ! Unique temporary file
     gp_tmp_file = 'tmp_'//trim(gp_label)//'.xml'

     ! Print GAP part of the potential into the temporary file.
     call xml_OpenFile(gp_tmp_file,xf,addDecl=.false.)

     call xml_NewElement(xf,"GAP_params")
     call xml_AddAttribute(xf,"label",trim(gp_label))

     call xml_NewElement(xf,"GAP_data")
     call xml_AddAttribute(xf,"n_species",""//this%n_species)
     call xml_AddAttribute(xf,"do_core",""//this%do_core)
     call xml_AddAttribute(xf,"e0",""//this%e0)
     call xml_AddAttribute(xf,"f0",""//this%f0)
     call xml_AddAttribute(xf,"do_pca",""//this%do_pca)

     select case(trim(this%coordinates))
     case('water_monomer')
        call xml_AddAttribute(xf,"coordinates","water_monomer")
        call xml_NewElement(xf,"water_monomer_params")
        call xml_AddAttribute(xf,"cutoff",""//this%r_cut)
        call xml_EndElement(xf,"water_monomer_params")
     case('water_dimer')
        call xml_AddAttribute(xf,"coordinates","water_dimer")
        call xml_NewElement(xf,"water_dimer_params")
        call xml_AddAttribute(xf,"cutoff",""//this%r_cut)
        call xml_EndElement(xf,"water_dimer_params")
     case('qw')
        call xml_AddAttribute(xf,"coordinates","qw")

        call xml_NewElement(xf,"qw_so3_params")
        call xml_AddAttribute(xf,"l_max",""//this%qw_l_max)
        call xml_AddAttribute(xf,"n_radial",""//this%qw_f_n)
        call xml_AddAttribute(xf,"do_q",""//(.not.this%qw_no_q))
        call xml_AddAttribute(xf,"do_w",""//(.not.this%qw_no_w))

        do i = 1, this%qw_f_n
           call xml_NewElement(xf,"radial_function")
           call xml_AddAttribute(xf,"i",""//i)
           call xml_AddAttribute(xf,"cutoff",""//this%qw_cutoff(i))
           call xml_AddAttribute(xf,"cutoff_type",""//this%qw_cutoff_f(i))
           call xml_AddAttribute(xf,"cutoff_r1",""//this%qw_cutoff_r1(i))
           call xml_EndElement(xf,"radial_function")
        enddo

        call xml_EndElement(xf,"qw_so3_params")
     case('bispectrum')
        call xml_AddAttribute(xf,"coordinates","bispectrum")
        
        call xml_NewElement(xf,"bispectrum_so4_params")
        call xml_AddAttribute(xf,"cutoff",""//this%r_cut)
        call xml_AddAttribute(xf,"j_max",""//this%j_max)
        call xml_AddAttribute(xf,"z0",""//this%z0)
        call xml_EndElement(xf,"bispectrum_so4_params")
     case default
        call system_abort('Unknown coordinates '//trim(this%coordinates))
     endselect

     if(this%do_pca) then
        call xml_NewElement(xf,"PCA_matrix")
        call xml_AddAttribute(xf,"n",""//this%d)
        call xml_NewElement(xf,"PCA_mean")
        call xml_AddCharacters(xf, ""//this%pca_mean)
        call xml_EndElement(xf,"PCA_mean")
        do i = 1, this%d
           call xml_NewElement(xf,"row")
           call xml_AddAttribute(xf,"i",""//i)
           call xml_AddCharacters(xf, ""//this%pca_matrix(:,i))
           call xml_EndElement(xf,"row")
        enddo
        call xml_EndElement(xf,"PCA_matrix")
     endif

     do i = 1, this%n_species
        call xml_NewElement(xf,"per_type_data")
        call xml_AddAttribute(xf,"i",""//i)
        call xml_AddAttribute(xf,"atomic_num",""//this%species_Z(i))
        call xml_AddAttribute(xf,"weight",""//this%w_Z(this%species_Z(i)))
        call xml_EndElement(xf,"per_type_data")
     enddo
     call xml_EndElement(xf,"GAP_data")

     ! Print GP bit of the potential
     call gp_print_xml(this%my_gp,xf,label=gp_label)

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

     if(this%do_core) then
        ! Create the sum potential xml entry (by hand)
        call print('<Potential label="'//trim(gp_label)//'" init_args="Sum init_args_pot1={'//trim(this%ip_args)//'} init_args_pot2={IP GAP label='//trim(gp_label)//'}"/>',file=gp_inout)

        ! Now add the core potential that was used.
        call print(string(this%quip_string),file=gp_inout)
     endif

     ! Add the GAP potential
     call print(string(gap_string),file=gp_inout)

     ! Close the root element
     call print('</'//trim(gp_label)//'>',file=gp_inout)

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

end module teach_sparse_mod
