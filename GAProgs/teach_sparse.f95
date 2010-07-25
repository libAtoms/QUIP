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

module teach_sparse_module

  use libatoms_module
  use bispectrum_module
  use fox_wxml
  use IPEwald_module
  use QUIP_module

  implicit none

  type teach_sparse
     character(len=FIELD_LENGTH) :: at_file, ip_args = '', &
     energy_property_name, force_property_name, virial_property_name
     real(dp) :: r_cut, e0, z0, f0, dlt, theta_fac
     real(dp), dimension(3) :: sgm
     logical :: do_core = .false., do_ewald = .false. , do_ewald_corr = .false., &
     do_qw_so3, qw_no_q, qw_no_w, do_sigma, do_delta, do_theta, do_sparx, do_f0, &
     do_theta_fac, do_test_gp_gradient, do_cluster, do_pivot

     integer :: d, m, j_max, qw_l_max, n, nn, ne, n_ener, n_force, n_virial, min_steps, min_save, n_species, &
     qw_f_n
     type(extendable_str) :: quip_string

     real(dp), dimension(116) :: z_eff
     real(dp), dimension(99) :: qw_cutoff, qw_cutoff_r1
     real(dp), dimension(:), allocatable :: w_Z, yf, ydf, dlta
     real(dp), dimension(:,:), allocatable :: x, xd, theta
     integer, dimension(:), allocatable :: lf, ldf, xf, xdf, xz, target_type, r, species_Z
     integer, dimension(99) :: qw_cutoff_f

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
    call query(xyzfile)

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
       call add_property(at,'charge',0.0_dp,n_cols=1)

       has_ener = get_value(at%params,this%energy_property_name,ener)
       has_force = assign_pointer(at,this%force_property_name, f)
       has_virial = get_value(at%params,this%virial_property_name,virial)

       if( has_ener .or. has_force .or. has_virial ) then
          call atoms_set_cutoff(at,this%r_cut)
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
          this%n_virial = this%n_virial + 9
          do i = 1, at%N
             this%n = this%n + 9*(atoms_n_neighbours(at,i)+1)
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

    allocate(this%species_Z(this%n_species))
    this%species_Z = species_present(1:this%n_species)

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
    
    logical :: do_qw_so3
    integer :: d
    integer :: n_max, n_con
    logical :: has_ener, has_force, has_virial
    real(dp) :: ener, ener_ewald, ener_ewald_corr, ener_core
    real(dp), dimension(3,3) :: virial, virial_ewald, virial_ewald_corr, virial_core
    real(dp), pointer :: charge(:), f(:,:)
    real(dp), dimension(:,:), allocatable :: f_ewald, f_ewald_corr, f_core
    real(dp), allocatable :: vec(:,:), jack(:,:,:), w(:)
    integer :: shift(3)
    real(dp) :: pos(3)
    integer :: li, ui, nn, ix, ie, i_con, i_ener
    integer :: nei_max
    integer :: i, j, n, k, l, jx, jn

    if(this%do_qw_so3) then

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
    else 
       call initialise(f_hat,this%j_max, &
       this%z0,this%r_cut)
       call initialise(df_hat,this%j_max, &
       this%z0,this%r_cut)

       this%d = j_max2d(f_hat%j_max)
    endif
    d = this%d

    allocate(this%x(d,this%nn),this%xd(d,this%n), &
    this%yf(this%n_ener), this%ydf(this%n_force+this%n_virial), &
    this%lf(this%n_ener),this%ldf(this%n_force+this%n_virial), &
    this%xf(this%ne),this%xdf(this%n))
    allocate(this%xz(this%nn))
    allocate(this%target_type(this%n_ener+this%n_force+this%n_virial))

    if( this%do_core ) call Initialise(core_pot, args_str=this%ip_args, param_str=string(this%quip_string))

    call initialise(xyzfile,this%at_file)
    call query(xyzfile)

    n_max = xyzfile%n_frame

    li = 0
    ui = 0
    nn = 0
    ix = 0
    ie = 0
    i_con = 0
    i_ener = 0

    do n_con = 1, n_max
       call read(xyzfile,at,frame=n_con-1)
       call add_property(at, 'charge',0.0_dp,n_cols=1)

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

          call calc(core_pot,at,e=ener_core,f=f_core,virial=virial_core)

          if(has_ener) ener = ener - ener_core
          if(has_force) f = f - f_core
          if(has_virial) virial = virial - virial_core
          deallocate(f_core)
       endif

       if( this%do_ewald ) then
          allocate(f_ewald(3,at%N))
          if( .not. assign_pointer(at, 'charge', charge) ) call system_abort('Could not assign pointer')
          do i = 1, at%N
             charge(i) = this%z_eff(at%Z(i))
          enddo
          call Ewald_calc(at,e=ener_ewald,f=f_ewald,virial=virial_ewald)
          if(has_ener) ener = ener - ener_ewald
          if(has_force) f = f - f_ewald
          if(has_virial) virial = virial - virial_ewald
          deallocate(f_ewald)
          if( this%do_ewald_corr ) then
             allocate(f_ewald_corr(3,at%N))
             call Ewald_corr_calc(at,e=ener_ewald_corr,f=f_ewald_corr,virial=virial_ewald_corr,cutoff=this%r_cut)
             if(has_ener) ener = ener + ener_ewald_corr
             if(has_force) f = f + f_ewald_corr
             if(has_virial) virial = virial + virial_ewald_corr
             deallocate(f_ewald_corr)
          endif
       endif

       if(has_ener) ener = ener - at%N*this%e0     

       if( at%cutoff /= this%r_cut ) then
          call atoms_set_cutoff(at,this%r_cut)
          call calc_connect(at)
       endif

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
          do i = 1, at%N
             if (do_qw_so3) then
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
             else
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
             endif
          enddo

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
             do k = 1, 3
                do l = 1, 3
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
    real(dp) :: ener, ener_ewald, ener_ewald_corr, ener_core
    real(dp), pointer :: charge(:)
    integer :: i

    if( this%do_core ) call Initialise(core_pot, this%ip_args, param_str=string(this%quip_string))

    call initialise(xyzfile,this%at_file)
    call query(xyzfile)

    n_max = xyzfile%n_frame

    n_ener = 0
    this%e0 = 0.0_dp

    do n_con = 1, n_max
       call read(xyzfile,at,frame=n_con-1)
       call add_property(at,'charge',0.0_dp,n_cols=1)

       has_ener = get_value(at%params,'Energy',ener)

       if( has_ener ) then

          ener_core = 0.0_dp
          if( this%do_core ) then
             call set_cutoff(at, cutoff(core_pot))
             call calc_connect(at)
             call calc(core_pot,at,e=ener_core)
          endif

          ener_ewald = 0.0_dp
          ener_ewald_corr = 0.0_dp
          if( this%do_ewald ) then
             if( .not. assign_pointer(at, 'charge', charge) ) call system_abort('Could not assign pointer')
             do i = 1, at%N
                charge(i) = this%z_eff(at%Z(i))
             enddo
             call Ewald_calc(at,e=ener_ewald)
             if( this%do_ewald_corr ) then
                call set_cutoff(at, this%r_cut)
                call calc_connect(at)
                call Ewald_corr_calc(at,e=ener_ewald_corr,cutoff=this%r_cut)
             endif
          endif

          this%e0 = this%e0 + (ener-ener_ewald+ener_ewald_corr-ener_core) / at%N

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
    call initialise(xyzfile,this%at_file)

    call read(xyzfile,at,frame=0)
    call get_weights(at,this%w_Z)
    call finalise(at)

    call finalise(xyzfile)

  end subroutine w_Z_from_xyz

  subroutine teach_sparse_print_xml(this,xf)
     type(teach_sparse), intent(in) :: this
     type(xmlf_t), intent(inout) :: xf
     integer :: i

     call xml_NewElement(xf,"GAP_data")
     call xml_AddAttribute(xf,"n_species",""//this%n_species)
     call xml_AddAttribute(xf,"do_ewald",""//this%do_ewald)
     call xml_AddAttribute(xf,"do_ewald_corr",""//this%do_ewald_corr)
     call xml_AddAttribute(xf,"do_core",""//this%do_core)
     call xml_AddAttribute(xf,"e0",""//this%e0)
     call xml_AddAttribute(xf,"f0",""//this%f0)

     if(this%do_qw_so3) then
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
     else
        call xml_AddAttribute(xf,"coordinates","bispectrum")
        
        call xml_NewElement(xf,"bispectrum_so4_params")
        call xml_AddAttribute(xf,"cutoff",""//this%r_cut)
        call xml_AddAttribute(xf,"j_max",""//this%j_max)
        call xml_AddAttribute(xf,"z0",""//this%z0)
        call xml_EndElement(xf,"bispectrum_so4_params")
     endif

     do i = 1, this%n_species
        call xml_NewElement(xf,"per_type_data")
        call xml_AddAttribute(xf,"i",""//i)
        call xml_AddAttribute(xf,"atomic_num",""//this%species_Z(i))
        call xml_AddAttribute(xf,"weight",""//this%w_Z(this%species_Z(i)))
        call xml_AddAttribute(xf,"charge",""//this%z_eff(this%species_Z(i)))
        call xml_EndElement(xf,"per_type_data")
     enddo

     if(this%do_core) then
        call xml_NewElement(xf,"core_params")
        call xml_AddAttribute(xf,"ip_args",trim(this%ip_args))
        call xml_AddCharacters(xf,string(this%quip_string),parsed=.true.)
        call xml_EndElement(xf,"core_params")
     endif

     call xml_EndElement(xf,"GAP_data")
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

end module teach_sparse_module

program teach_sparse_program

  use teach_sparse_module
  use libatoms_module
  use bispectrum_module
  use gp_sparse_module
  use fox_wxml
  use clustering_module

  implicit none

  integer, parameter :: SPARSE_LENGTH = 10000
  integer, parameter :: SPARSE_N_FIELDS = 2000
  integer, parameter :: THETA_LENGTH = 10000
  real(dp), parameter :: THETA_MIN = 0.000000001

  type(teach_sparse) :: main_teach_sparse
  type(inoutput) :: bispectrum_inout, theta_inout, sparse_inout
  type(Dictionary) :: params
  type(gp) :: my_gp
  type(gp_sparse) :: gp_sp

  character(len=FIELD_LENGTH) :: qw_cutoff_string, qw_cutoff_f_string, qw_cutoff_r1_string, &
  theta_file, sparse_file, z_eff_string, bispectrum_file
  real(dp) :: mem_required, mem_total, mem_free
  logical :: has_e0, has_f0, has_theta_file, has_sparse_file, has_bispectrum_file, test_gp_gradient_result

  character(len=FIELD_LENGTH), dimension(232) :: z_eff_fields
  integer :: num_z_eff_fields
  character(len=FIELD_LENGTH), dimension(99) :: qw_cutoff_fields, qw_cutoff_f_fields, qw_cutoff_r1_fields
  character(len=SPARSE_LENGTH) :: sparse_string
  character(len=FIELD_LENGTH), dimension(:), allocatable :: sparse_string_array
  character(len=THETA_LENGTH) :: theta_string
  character(len=FIELD_LENGTH), dimension(:), allocatable :: theta_string_array
  integer :: i, j, k, l, o, dd, dt
  character(len=FIELD_LENGTH) :: gp_file
  character(len=10240) :: command_line

  type(xmlf_t) :: xml_file

  call system_initialise(verbosity=PRINT_NORMAL)
  call initialise(params)
  call param_register(params, 'at_file', PARAM_MANDATORY, main_teach_sparse%at_file)
  call param_register(params, 'm', '50', main_teach_sparse%m)
  call param_register(params, 'r_cut', '2.75', main_teach_sparse%r_cut)
  call param_register(params, 'j_max', '4', main_teach_sparse%j_max)
  call param_register(params, 'z0', '0.0', main_teach_sparse%z0)
  call param_register(params, 'qw_so3', 'F', main_teach_sparse%do_qw_so3)
  call param_register(params, 'l_max', '6', main_teach_sparse%qw_l_max)
  call param_register(params, 'cutoff', '', qw_cutoff_string)
  call param_register(params, 'cutoff_f', '', qw_cutoff_f_string)
  call param_register(params, 'cutoff_r1', '', qw_cutoff_r1_string)
  call param_register(params, 'no_q', 'F', main_teach_sparse%qw_no_q)
  call param_register(params, 'no_w', 'F', main_teach_sparse%qw_no_w)
  call param_register(params, 'e0', '0.0', main_teach_sparse%e0, has_e0)
  call param_register(params, 'f0', '0.0', main_teach_sparse%f0, has_f0)
  call param_register(params, 'sgm', '0.1 0.1 0.1', main_teach_sparse%sgm)
  call param_register(params, 'dlt', '1.0', main_teach_sparse%dlt)
  call param_register(params, 'theta_file', '', theta_file, has_theta_file)
  call param_register(params, 'sparse_file', '', sparse_file, has_sparse_file)
  call param_register(params, 'theta_fac', '1.5', main_teach_sparse%theta_fac)
  call param_register(params, 'do_sigma', 'F', main_teach_sparse%do_sigma)
  call param_register(params, 'do_delta', 'F', main_teach_sparse%do_delta)
  call param_register(params, 'do_theta', 'F', main_teach_sparse%do_theta)
  call param_register(params, 'do_sparx', 'F', main_teach_sparse%do_sparx)
  call param_register(params, 'do_f0', 'F', main_teach_sparse%do_f0)
  call param_register(params, 'do_theta_fac', 'F', main_teach_sparse%do_theta_fac)
  call param_register(params, 'do_cluster', 'F', main_teach_sparse%do_cluster)
  call param_register(params, 'do_pivot', 'F', main_teach_sparse%do_pivot)
  call param_register(params, 'min_steps', '10', main_teach_sparse%min_steps)
  call param_register(params, 'min_save', '0', main_teach_sparse%min_save)
  call param_register(params, 'z_eff', '', z_eff_string,main_teach_sparse%do_ewald)
  call param_register(params, 'do_test_gp_gradient', 'F', main_teach_sparse%do_test_gp_gradient)
  call param_register(params, 'bispectrum_file', '', bispectrum_file, has_bispectrum_file)
  call param_register(params, 'ip_args', '', main_teach_sparse%ip_args, main_teach_sparse%do_core)
  call param_register(params, 'do_ewald_corr', 'T', main_teach_sparse%do_ewald_corr)
  call param_register(params, 'energy_property_name', 'energy', main_teach_sparse%energy_property_name)
  call param_register(params, 'force_property_name', 'force', main_teach_sparse%force_property_name)
  call param_register(params, 'virial_property_name', 'virial', main_teach_sparse%virial_property_name)

  if (.not. param_read_args(params, do_check = .true., command_line=command_line)) then
     call print("Usage: teach_sparse [at_file=file] [m=50] &
     & [r_cut=2.75] [j_max=4] [z0=0.0] [qw_so3] [l_max=6] [cutoff={:}] [cutoff_f={:}] [cutoff_r1={:}] [no_q] [no_w] &
     & [e0=0.0] [f0=avg] [sgm={0.1 0.1 0.1}] [dlt=1.0] [theta_file=file] [sparse_file=file] [theta_fac=3.0] &
     & [do_sigma=F] [do_delta=F] [do_theta=F] [do_sparx=F] [do_f0=F] [do_theta_fac=F] &
     & [do_cluster=F] [do_pivot=F] [min_steps=10] [min_save=0] [z_eff={Ga:1.0:N:-1.0}] &
     & [do_test_gp_gradient=F] [bispectrum_file=file] [ip_args={}] [do_ewald_corr=F] &
     & [energy_property_name=energy] [force_property_name=force] [virial_property_name=virial]")
     call system_abort('Exit: Mandatory argument(s) missing...')
  endif
  call finalise(params)

  if( count( (/has_e0,has_f0/) ) > 1 ) &
  & call print('Warning - you have specified both e0 and f0 - careful!')

  if( count( (/has_sparse_file,main_teach_sparse%do_cluster,main_teach_sparse%do_pivot/) ) > 1 ) &
  & call system_abort('There has been more than one method specified for sparsification.')
     
  main_teach_sparse%z_eff = 0.0_dp
  if(main_teach_sparse%do_ewald) then
     call parse_string(z_eff_string,':',z_eff_fields,num_z_eff_fields)
     do i = 1, num_z_eff_fields, 2
        j = atomic_number_from_symbol(z_eff_fields(i))
        if(j < 1 .or. j > 116) call system_abort("Invalid atomic number "//j//" parsed from "//z_eff_fields(i))
        main_teach_sparse%z_eff(j) = string_to_real(z_eff_fields(i+1))
     enddo
  endif
  main_teach_sparse%do_ewald_corr = main_teach_sparse%do_ewald .and. main_teach_sparse%do_ewald_corr

  if (main_teach_sparse%do_qw_so3) then
     main_teach_sparse%qw_cutoff = 0.0_dp
     main_teach_sparse%qw_cutoff_f = 0
     main_teach_sparse%qw_cutoff_r1 = 0.0_dp
     call parse_string(qw_cutoff_string, ':', qw_cutoff_fields, main_teach_sparse%qw_f_n)
     call parse_string(qw_cutoff_f_string, ':', qw_cutoff_f_fields, main_teach_sparse%qw_f_n)
     call parse_string(qw_cutoff_r1_string, ':', qw_cutoff_r1_fields, main_teach_sparse%qw_f_n)
     do i = 1, main_teach_sparse%qw_f_n
        main_teach_sparse%qw_cutoff(i) = string_to_real(qw_cutoff_fields(i))
        main_teach_sparse%qw_cutoff_f(i) = string_to_int(qw_cutoff_f_fields(i))
        main_teach_sparse%qw_cutoff_r1(i) = string_to_real(qw_cutoff_r1_fields(i))
     enddo
     main_teach_sparse%r_cut = maxval(main_teach_sparse%qw_cutoff(1:main_teach_sparse%qw_f_n))
  else
     main_teach_sparse%z0 = max(1.0_dp,main_teach_sparse%z0) * main_teach_sparse%r_cut/(PI-0.02_dp)
     main_teach_sparse%d = j_max2d(main_teach_sparse%j_max)
  endif

  if(main_teach_sparse%do_core) &
     call read(main_teach_sparse%quip_string, "core_params.xml")

  call teach_n_from_xyz(main_teach_sparse)

  if (.not. has_f0) then
     call e0_avg_from_xyz(main_teach_sparse)
  end if

  call w_Z_from_xyz(main_teach_sparse)

  call teach_data_from_xyz(main_teach_sparse)

  if( has_bispectrum_file ) then
     call initialise(bispectrum_inout,bispectrum_file,action=OUTPUT)
     do i = 1, main_teach_sparse%nn/3
        write(bispectrum_inout%unit,"("//main_teach_sparse%d//"f16.8)") main_teach_sparse%x(:,i)
     enddo
     call finalise(bispectrum_inout)
  endif

  allocate(main_teach_sparse%dlta(main_teach_sparse%n_species))

  main_teach_sparse%dlta = main_teach_sparse%dlt 

  if( has_sparse_file ) then
     allocate(sparse_string_array(SPARSE_N_FIELDS))
     call initialise(sparse_inout,sparse_file)
     read(sparse_inout%unit,'(a)') sparse_string
     call parse_string(sparse_string,' ',sparse_string_array,main_teach_sparse%m)
     allocate(main_teach_sparse%r(main_teach_sparse%m))
     do i = 1, main_teach_sparse%m
        main_teach_sparse%r(i) = string_to_int(sparse_string_array(i))
     enddo
     deallocate(sparse_string_array)
     call finalise(sparse_inout)
  elseif(main_teach_sparse%do_cluster) then
     allocate(main_teach_sparse%r(main_teach_sparse%m))
     call bisect_kmedoids(main_teach_sparse%x,main_teach_sparse%m,med=main_teach_sparse%r,theta_fac=main_teach_sparse%theta_fac)
  elseif(main_teach_sparse%do_pivot) then
     allocate(main_teach_sparse%r(main_teach_sparse%m))
     call pivot(main_teach_sparse%x, main_teach_sparse%r,theta_fac=main_teach_sparse%theta_fac)
  else
     allocate(main_teach_sparse%r(main_teach_sparse%m))
     call fill_random_integer(main_teach_sparse%r,size(main_teach_sparse%x,2))
  endif
  call sort_array(main_teach_sparse%r)

  call print('')
  call print('Atomic environments used in sparsification')
  call print(main_teach_sparse%r)
  call print('')

  allocate(main_teach_sparse%theta(main_teach_sparse%d,main_teach_sparse%n_species))

  if( has_theta_file ) then
     allocate(theta_string_array(main_teach_sparse%d))
     call initialise(theta_inout,theta_file)
     read(theta_inout%unit,'(a)') theta_string
     call parse_string(theta_string,' ',theta_string_array,dt)
     if(main_teach_sparse%d /= dt) call system_abort('File '//trim(theta_file)//' does not contain the right number of hyperparameters')
     do k = 1, main_teach_sparse%n_species  
        do dd = 1, main_teach_sparse%d
           main_teach_sparse%theta(dd,k) = string_to_real(theta_string_array(dd+(k-1)*main_teach_sparse%d))
        enddo
     enddo
     deallocate(theta_string_array)
     call finalise(theta_inout)
  else
     do k = 1, main_teach_sparse%n_species
        do dd = 1, main_teach_sparse%d
!           theta(dd,k) = ( maxval(x(dd,:),mask=(xz(:)==species_Z(k))) - minval(x(dd,:),mask=(xz(:)==species_Z(k))) )
           main_teach_sparse%theta(dd,k) = ( maxval(main_teach_sparse%x(dd,main_teach_sparse%r),&
           mask=(main_teach_sparse%xz(main_teach_sparse%r)==main_teach_sparse%species_Z(k))) &
           - minval(main_teach_sparse%x(dd,main_teach_sparse%r),&
           mask=(main_teach_sparse%xz(main_teach_sparse%r)==main_teach_sparse%species_Z(k))) )
!           theta(dd) = sqrt( & !take square root
!                          & sum( x(dd,:)**2 ) / size(x(dd,:)) - &
!                          & (sum( x(dd,:) ) / size(x(dd,:)))**2 )

           if( main_teach_sparse%theta(dd,k) >= THETA_MIN ) then
              main_teach_sparse%theta(dd,k) = main_teach_sparse%theta_fac*main_teach_sparse%theta(dd,k)
           else
              main_teach_sparse%theta(dd,k) = 1.0_dp
           endif
        enddo
     enddo
  endif

  ! Stop execution if required memory is greater than the available memory. 
  ! The biggest arrays allocated are 2*sr*(nx+nxd), where sr is the
  ! number of sparse points, nx and nxd are the number of bispectra and partial
  ! derivatives.
  mem_required = 2.0_dp * real(size(main_teach_sparse%r),dp) * (real(size(main_teach_sparse%xf),dp) &
  + real(size(main_teach_sparse%xdf),dp)) * real(dp,dp) / (1024.0_dp**3)
  call mem_info(mem_total,mem_free)
  mem_total = mem_total / (1024.0_dp**3)
  mem_free = mem_free / (1024.0_dp**3)

  call print('Memory required (approx.): '//mem_required//' GB')
  if( mem_required > mem_free ) call system_abort('Required memory ('//mem_required//' GB) exceeds available memory ('//mem_free//' GB).')

  call gp_sparsify(gp_sp,main_teach_sparse%r,&
  main_teach_sparse%sgm,main_teach_sparse%dlta,main_teach_sparse%theta,&
  main_teach_sparse%yf,main_teach_sparse%ydf,main_teach_sparse%x,main_teach_sparse%xd,&
  main_teach_sparse%xf,main_teach_sparse%xdf,main_teach_sparse%lf,main_teach_sparse%ldf,&
  main_teach_sparse%xz,main_teach_sparse%species_Z,(/(main_teach_sparse%f0,i=1,main_teach_sparse%n_species)/),&
  main_teach_sparse%target_type)

  !deallocate(x,xd,xf,xdf,yf,ydf,lf,ldf)

  call print('')
  call print('theta')
  do l = 1, size(gp_sp%theta, 2)
     do o = 1, size(gp_sp%theta, 1)
        call print(real(gp_sp%theta(o,l),kind=dp))
     enddo
  enddo
  call print('')

  call enable_timing()
  
  if( main_teach_sparse%do_test_gp_gradient ) then
     call verbosity_push(PRINT_NERD)
     test_gp_gradient_result = test_gp_gradient(gp_sp,&
     sigma=main_teach_sparse%do_sigma,delta=main_teach_sparse%do_delta,&
     theta=main_teach_sparse%do_theta,sparx=main_teach_sparse%do_sparx,&
     f0=main_teach_sparse%do_f0,theta_fac=main_teach_sparse%do_theta_fac)
     call verbosity_pop()
  endif

  ! Conjugate gradient minimiser's counter starts at 1, and stops when it reaches min_steps,
  ! so if min_steps is equal to 1, no iterations are made!
 
  if (main_teach_sparse%min_save == 0) main_teach_sparse%min_save = main_teach_sparse%min_steps

  k = 0
  do i = 1, ((main_teach_sparse%min_steps / main_teach_sparse%min_save) + 1)
     if (k == main_teach_sparse%min_steps) exit

     if ((main_teach_sparse%min_steps - k) >= main_teach_sparse%min_save) then
        j = minimise_gp_gradient(gp_sp,max_steps=(main_teach_sparse%min_save + 1),&
        sigma=main_teach_sparse%do_sigma,delta=main_teach_sparse%do_delta,&
        theta=main_teach_sparse%do_theta,sparx=main_teach_sparse%do_sparx,&
        f0=main_teach_sparse%do_f0,theta_fac=main_teach_sparse%do_theta_fac)
        k = k + main_teach_sparse%min_save
     elseif ((main_teach_sparse%min_steps - k) < main_teach_sparse%min_save) then
        j = minimise_gp_gradient(gp_sp,max_steps=(main_teach_sparse%min_steps - k + 1),&
        sigma=main_teach_sparse%do_sigma,delta=main_teach_sparse%do_delta,&
        theta=main_teach_sparse%do_theta,sparx=main_teach_sparse%do_sparx,&
        f0=main_teach_sparse%do_f0,theta_fac=main_teach_sparse%do_theta_fac)

        k = main_teach_sparse%min_steps
     endif

     call print('')
     call print(k // ' iterations completed:')
     call print('delta')
     call print(real(gp_sp%delta,kind=dp))
     call print('theta')
     do l = 1, size(gp_sp%theta, 2)
        do o = 1, size(gp_sp%theta, 1)
           call print(real(gp_sp%theta(o,l),kind=dp))
        enddo
     enddo
     call print('f0')
     call print(real(gp_sp%f0,kind=dp))
     call print('')

     call initialise(my_gp,gp_sp)

     gp_file = 'gp_'//main_teach_sparse%m//'_'//k//'.xml'

     call xml_OpenFile(gp_file,xml_file)
     call xml_NewElement(xml_file,"GAP_params")
     call teach_sparse_print_xml(main_teach_sparse,xml_file)
     call gp_print_xml(my_gp,xml_file)
     call xml_NewElement(xml_file,"commmmand_line")
     call xml_AddCharacters(xml_file,trim(command_line),parsed=.false.)
     call xml_EndElement(xml_file,"commmmand_line")
     call file_print_xml(main_teach_sparse%at_file,xml_file)
     call xml_EndElement(xml_file,"GAP_params")
     call xml_Close(xml_file)


     call gp_print_binary(my_gp,trim(gp_file))
     call system_command('ln -fs '//trim(gp_file)//' gp.xml')

     call finalise(my_gp)
  enddo

  call print("model parameters:")
  call print("r_cut     = "//main_teach_sparse%r_cut)
  if (main_teach_sparse%do_qw_so3) then
     call print("l_max     = "//main_teach_sparse%qw_l_max)
     call print("cutoff    = "//qw_cutoff_string)
     call print("cutoff_f  = "//qw_cutoff_f_string)
     call print("cutoff_r1 = "//qw_cutoff_r1_string)
     call print("q         = "//(.not. main_teach_sparse%qw_no_q))
     call print("w         = "//(.not. main_teach_sparse%qw_no_w))
  else
     call print("j_max     = "//main_teach_sparse%j_max)
     call print("z0        = "//main_teach_sparse%z0)
  endif
  call print("n_species = "//main_teach_sparse%n_species)
  call print("species_Z = "//main_teach_sparse%species_Z)
  call print("w         = "//main_teach_sparse%w_Z(main_teach_sparse%species_Z))
  call print("z_eff     = "//main_teach_sparse%z_eff(main_teach_sparse%species_Z))
  call print("do_ewald  = "//main_teach_sparse%do_ewald)
  call print("do_ewald_corr  = "//main_teach_sparse%do_ewald_corr)
  call print("e0        = "//main_teach_sparse%e0)

  call finalise(gp_sp)

  call system_finalise()

end program teach_sparse_program
