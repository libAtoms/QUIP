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
     character(len=STRING_LENGTH) :: type
     integer :: m
     real(dp), dimension(3) :: sgm
  endtype sparse_types

  type teach_sparse
     character(len=STRING_LENGTH) :: at_file='', ip_args = '', &
     energy_parameter_name, force_parameter_name, virial_parameter_name, coordinates, config_type_parameter_name, mark_sparse_atoms = '', mask_name = ''
     character(len=10240) :: command_line = ''
     real(dp) :: r_cut, e0, z0, f0, dlt, theta_fac
     real(dp), dimension(3) :: sgm
     logical :: do_core = .false., &
     qw_no_q, qw_no_w, do_sigma, do_delta, do_theta, do_sparx, do_f0, &
     do_theta_fac, do_test_gp_gradient, do_cluster, do_pivot, do_sparse, &
     has_config_type_hypers, do_pca, do_mark_sparse_atoms, use_rdf

     integer :: d, m, j_max, qw_l_max, n, nn, ne, n_ener, n_force, n_virial, min_steps, min_save, n_species, &
     qw_f_n, cosnx_l_max, cosnx_n_max
     type(extendable_str) :: quip_string
     type(gp) :: my_gp

     real(dp), dimension(99) :: qw_cutoff, qw_cutoff_r1
     real(dp), dimension(:), allocatable :: w_Z, yf, ydf, dlta, pca_mean, sigma, NormFunction
     real(dp), dimension(:,:), allocatable :: x, xd, theta, pca_matrix, RadialTransform, rdf
     integer, dimension(:), allocatable :: lf, ldf, xf, xdf, xz, target_type, r, species_Z, config_type
     integer, dimension(99) :: qw_cutoff_f
     type(sparse_types), dimension(:), allocatable :: config_type_hypers

  endtype teach_sparse
     
  private

  public :: teach_n_from_xyz
  public :: teach_data_from_xyz
  public :: e0_avg_from_xyz
  public :: w_Z_from_xyz
  public :: teach_sparse
  public :: teach_sparse_print_xml
  public :: file_print_xml
  public :: print_sparse

contains

  subroutine teach_n_from_xyz(this)

    type(teach_sparse), intent(inout) :: this

    type(cinoutput) :: xyzfile
    type(atoms) :: at
    integer :: n_max, n_con, n_at
    logical :: has_ener, has_force, has_virial, has_mask
    real(dp) :: ener, virial(3,3), d, bin_width
    real(dp), pointer :: f(:,:)
    integer :: i, j, n, num_bins, n_neighb
    logical, dimension(:), pointer :: mask
    integer, dimension(116) :: species_present

    type(Table) :: distances

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

    if ( this%use_rdf ) then
       call allocate(distances,0,1,0,0,1000000)
       call set_increment(distances,1000000)
    endif

    do n_con = 1, n_max
       call read(xyzfile,at,frame=n_con-1)

       has_mask = .false.
       if (len_trim(this%mask_name) > 0) has_mask = assign_pointer(at, this%mask_name, mask)
       if (.not. has_mask) then
          allocate(mask(at%n))
          mask(:) = .true.
       end if

       has_ener = get_value(at%params,this%energy_parameter_name,ener)
       has_force = assign_pointer(at,this%force_parameter_name, f)
       has_virial = get_value(at%params,this%virial_parameter_name,virial)

       if (has_mask .and. (has_ener .or. has_virial)) &
            call system_abort('Has mask and energy or virial')

       if (has_mask .and. trim(this%coordinates) /= 'bispectrum') &
            call system_abort('Masked teaching only implemented for bispectrum descriptor')

       n_at = count(mask)

       if( has_ener .or. has_force .or. has_virial .or. this%use_rdf ) then
          call set_cutoff(at,this%r_cut)
          call calc_connect(at)
       endif

       if ( this%use_rdf ) then
          do i = 1, at%N
             do n = 1, atoms_n_neighbours(at,i)
                j = atoms_neighbour(at,i,n,distance=d)
                call append(distances, d)
             enddo
          enddo
       endif

       if( has_ener ) then
          this%n_ener = this%n_ener + 1
          this%ne = this%ne + at%N
       endif

       if( has_force ) then
          this%n_force = this%n_force + 3*count(mask)
          ! count number of neighbour pairs (i,j) where both i and j are in mask
          n_neighb = 0
          do i = 1, at%N
             if (.not. mask(i)) cycle
             do n = 1, atoms_n_neighbours(at, i)
                j = atoms_neighbour(at, i, n)
                if (.not. mask(j)) cycle
                n_neighb = n_neighb + 1
             end do
          enddo
          this%n = this%n + 3*n_at + 3*n_neighb
       endif

       if( has_virial ) then
          this%n_virial = this%n_virial + 6
          do i = 1, at%N
             this%n = this%n + 6*(atoms_n_neighbours(at,i)+1)
          enddo
       endif

       do i = 1, at%N
          if (.not. mask(i)) cycle
          if( all(at%Z(i) /= species_present) ) then
             this%n_species = this%n_species + 1
             species_present(this%n_species) = at%Z(i)
          endif
       enddo

       if( has_ener .or. has_force .or. has_virial ) this%nn = this%nn + n_at

       call finalise(at)
       if (.not. has_mask) deallocate(mask)
    enddo

    call finalise(xyzfile)

    if ( this%use_rdf ) then
       num_bins = int(sqrt(real(distances%N,dp)))
       allocate(this%rdf(num_bins,2))
       
       bin_width = this%r_cut / real(num_bins,dp)
       do i = 1, num_bins
          this%rdf(i,1) = (real(i,dp) - 0.5_dp) * bin_width
       enddo
    
       this%rdf(:,2) = histogram(distances%real(1,1:distances%N), 0.0_dp, this%r_cut, num_bins)
       call finalise(distances)
    endif

    select case(trim(this%coordinates))
    case('hf_dimer')
       this%n_species = 1
       allocate(this%species_Z(this%n_species))
       this%species_Z = 9
       this%nn = this%nn / 4
       this%ne = this%n_ener
       this%n = this%n_force
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
       this%nn = this%nn / 6
       this%ne = this%n_ener
       this%n = this%n_force
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
    type(cosnx) :: my_cosnx
    type(fourier_so3) :: f3_hat
    type(grad_fourier_so3) :: df3_hat
    type(qw_so3) :: qw
    type(grad_qw_so3) :: dqw

    type(Potential) :: core_pot
    
    integer :: d
    integer :: n_max, n_con, n_at
    logical :: has_ener, has_force, has_virial, has_config_type, has_mask
    real(dp) :: ener, ener_core
    real(dp), dimension(3,3) :: virial, virial_core
    real(dp), pointer :: f(:,:)
    real(dp), dimension(:,:), allocatable :: f_core
    real(dp), allocatable :: vec(:,:), jack(:,:,:), w(:), dvec(:,:,:,:)
    integer :: shift(3)
    integer, dimension(3,1) :: water_monomer_index
    integer, dimension(3,2) :: water_dimer_index
    real(dp), dimension(3) :: pos, water_monomer_v
    real(dp), dimension(WATER_DIMER_D) :: water_dimer_v
    logical, dimension(:), pointer :: mask
    integer :: li, ui, nn, ix, ie, i_con, i_ener, w_con
    integer :: nei_max
    integer :: i, j, n, k, l, jx, jn, ii, jj

    integer :: it, n_config_type
    character(len=STRING_LENGTH) :: config_type

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

    case('cosnx')
    
       if( this%use_rdf ) then
          call initialise(my_cosnx,this%cosnx_l_max, this%cosnx_n_max, this%r_cut, w = this%rdf )
       else
          call initialise(my_cosnx,this%cosnx_l_max, this%cosnx_n_max, this%r_cut )
       endif
       allocate(this%NormFunction(this%cosnx_n_max),this%RadialTransform(this%cosnx_n_max,this%cosnx_n_max))
       this%NormFunction = my_cosnx%NormFunction
       this%RadialTransform = my_cosnx%RadialTransform
    
       this%d = (this%cosnx_l_max+1)*this%cosnx_n_max

    case('hf_dimer')
       this%d = 6
    case('water_monomer')
       this%d = 3
    case('water_dimer')
       this%d = WATER_DIMER_D
    case default
       call system_abort('Unknown coordinates '//trim(this%coordinates))
    endselect
    d = this%d

    call print("Number of target energies (property name: "//trim(this%energy_parameter_name)//") found: "//this%n_ener)
    call print("Number of target forces (property name: "//trim(this%force_parameter_name)//") found: "//this%n_force)
    call print("Number of target virials (property name: "//trim(this%virial_parameter_name)//") found: "//this%n_virial)

    allocate(this%x(d,this%nn),this%xd(d,this%n), &
    this%yf(this%n_ener), this%ydf(this%n_force+this%n_virial), &
    this%lf(this%n_ener),this%ldf(this%n_force+this%n_virial), &
    this%xf(this%ne),this%xdf(this%n))
    allocate(this%xz(this%nn))
    allocate(this%target_type(this%n_ener+this%n_force+this%n_virial),this%config_type(this%nn))

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

       has_mask = .false.
       if (len_trim(this%mask_name) > 0) has_mask = assign_pointer(at, this%mask_name, mask)
       if (.not. has_mask) then
          allocate(mask(at%n))
          mask(:) = .true.
       end if
       n_at = count(mask)

       call print('Config #'//n_con//' has mask with '//count(mask)//' selected atoms')
       
       has_ener = get_value(at%params,this%energy_parameter_name,ener)
       has_force = assign_pointer(at,this%force_parameter_name, f)
       has_virial = get_value(at%params,this%virial_parameter_name,virial)
       has_config_type = get_value(at%params,this%config_type_parameter_name,config_type)

       if (has_mask .and. (has_ener .or. has_virial)) &
            call system_abort('Has mask and energy or virial')

       if (has_mask .and. trim(this%coordinates) /= 'bispectrum') &
            call system_abort('Masked teaching only implemented for bispectrum descriptor')

       if( has_config_type ) then
          config_type = lower_case(trim(config_type))
       else
          config_type = "default"
       endif

       if( .not. allocated(this%config_type_hypers) ) call system_abort('config_type_hypers not allocated')
       n_config_type = 0
       do it = lbound(this%config_type_hypers,dim=1), ubound(this%config_type_hypers,dim=1)
          if( trim(this%config_type_hypers(it)%type) == trim(config_type) ) n_config_type = it
       enddo

       if( this%do_core ) then
          allocate(f_core(3,at%N))
          ener_core = 0.0_dp
          f_core = 0.0_dp
          virial_core = 0.0_dp

          call set_cutoff(at, max(cutoff(core_pot),this%r_cut))
          call calc_connect(at)

          if(has_virial .and. has_force) then
             call calc(core_pot,at,energy=ener_core,force=f_core,virial=virial_core)
          elseif(has_force) then
             call calc(core_pot,at,energy=ener_core,force=f_core)
          else
             call calc(core_pot,at,energy=ener_core)
          endif

          if(has_ener) ener = ener - ener_core
          if(has_force) f = f - f_core
          if(has_virial) virial = virial - virial_core
          deallocate(f_core)
       endif

       if(has_ener) then
          select case(trim(this%coordinates))
          case('water_monomer','water_dimer','hf_dimer')
             ener = ener - this%e0     
          case default
             ener = ener - at%n*this%e0     
          endselect
       endif

       call set_cutoff(at,this%r_cut)
       call calc_connect(at)

       nei_max = 0
       do i = 1, at%N
          if(nei_max < atoms_n_neighbours(at,i)+1) nei_max = atoms_n_neighbours(at,i)+1
       enddo

       if(has_ener .or. has_force .or. has_virial ) allocate(vec(d,n_at))
       if(has_force .or. has_virial) then
          select case(trim(this%coordinates))
          case('water_dimer')
             allocate(dvec(3,6,d,1))
             dvec = 0.0_dp
          case('hf_dimer')
             allocate(dvec(3,4,d,1))
             dvec = 0.0_dp
          case default
             allocate(jack(d,3*nei_max,n_at))
             jack = 0.0_dp
          endselect
       endif
       allocate(w(at%N))

       do i = 1, at%N
          w(i) = this%w_Z(at%Z(i))
       enddo

       if(has_ener .or. has_force .or. has_virial ) then
          select case(trim(this%coordinates))
          case('hf_dimer')
             if( at%N /= 4 ) call system_abort('Number of atoms is '//at%N//', not two HF molecules')

             w_con = w_con + 1

             if(has_ener) then
                ie = ie + 1
                this%xf(ie) = w_con
                this%yf(ie) = ener
                this%lf(ie) = ie
                this%target_type(ie) = 1
             endif

             if(has_force) then
                li = ui + 1
                ui = ui + 12
                dvec(:,:,:,1) = hf_dimer_grad(at)
                this%xd(:,li:ui) = transpose(reshape(dvec(:,:,:,1), (/12,d/)))
                this%ydf(li:ui) = -reshape(f(:,:),(/12/))
                this%xdf(li:ui) = w_con
                this%ldf(li:ui) = (/(i, i=li,ui)/)
                this%target_type(this%n_ener+li:this%n_ener+ui) = 2
             endif

             this%x(:,w_con) = hf_dimer(at)
             this%xz(w_con) = 9
             this%config_type(w_con) = n_config_type

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
             if( at%n /= 6 ) call system_abort('Number of atoms is '//at%n//', not two water molecules')

             call find_water_monomer(at,water_dimer_index)
             call water_dimer(at,water_dimer_index(:,1),water_dimer_index(:,2),this%r_cut, vec = water_dimer_v,dvec=dvec)

             w_con = w_con + 1

             if(has_ener) then
                ie = ie + 1
                this%xf(ie) = w_con
                this%yf(ie) = ener
                this%lf(ie) = ie
                this%target_type(ie) = 1
             endif

             if(has_force) then
                li = ui + 1
                ui = ui + 18
                this%xd(:,li:ui) = transpose(reshape(dvec(:,:,:,1), (/18,d/)))
                !this%ydf(li:ui) = -reshape(f(:,(/water_dimer_index(:,1),water_dimer_index(:,2)/)),(/18/))
                this%ydf(li:ui) = -reshape(f(:,:),(/18/))                
                this%xdf(li:ui) = w_con
                this%ldf(li:ui) = (/(i, i=li,ui)/)
                this%target_type(this%n_ener+li:this%n_ener+ui) = 2
             endif

             this%x(:,w_con) = water_dimer_v
             this%xz(w_con) = 8
             this%config_type(w_con) = n_config_type

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
             ii = 0
             do i = 1, at%N
                if (.not. mask(i)) cycle
                ii = ii + 1
                call fourier_transform(f_hat,at,i,w)
                call calc_bispectrum(bis,f_hat)
                call bispectrum2vec(bis,vec(:,ii))
                if(has_force .or. has_virial) then
                   do n = 0, atoms_n_neighbours(at,i)
                      call fourier_transform(df_hat,at,i,n,w)
                      call calc_bispectrum(dbis,f_hat,df_hat)
                      call bispectrum2vec(dbis,jack(:,3*n+1:3*(n+1),ii))
                   enddo
                endif
             enddo
          case('cosnx')
             do i = 1, at%N
                call calc_cosnx(my_cosnx,at,vec(:,i),i,w)
                if(has_force .or. has_virial) then
                   do n = 0, atoms_n_neighbours(at,i)
                      call calc_grad_cosnx(my_cosnx,at,jack(:,3*n+1:3*(n+1),i),i,n,w)
                   enddo
                endif
             enddo
          case default
             call system_abort('Unknown coordinates '//trim(this%coordinates))
          endselect

          select case(trim(this%coordinates))
          case('water_monomer','water_dimer','hf_dimer')

          case default
             ii = 0
             do i = 1, at%N
                if (.not. mask(i)) cycle
                ii = ii + 1

                ix = i_con + ii
                this%x(:,ix) = vec(:,ii)
                this%config_type(ix) = n_config_type
                this%xz(ix) = at%Z(i)

                if( has_ener ) then
                   ie = ie + 1
                   this%xf(ie) = ix
                endif
           
                if(has_force) then
                   do k = 1, 3
                      nn = nn+1
                      ui = ui + 1
                      this%xdf(ui) = ix
                      this%xd(:,ui) = jack(:,k,ii)

                      do n = 1, atoms_n_neighbours(at,i)
                         j = atoms_neighbour(at,i,n,jn=jn)

                         if (.not. mask(j)) cycle
                         jj = count(mask(1:j)) ! find reduced index of neighbour
                         ui = ui + 1
                         jx = i_con + jj
                         this%xdf(ui) = jx
                         this%xd(:,ui) = jack(:,jn*3+k,jj)
                      enddo

                      this%ydf(nn) = -f(k,i)
                      this%ldf(nn) = ui
                      !sigma(nn+n_ener) = sgm(2)
                      this%target_type(nn+this%n_ener) = 3*n_config_type + 2
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
                      this%target_type(nn+this%n_ener) = 3*n_config_type + 3
                   enddo
                enddo
             endif

             i_con = i_con + n_at

             if(has_ener) then
                i_ener = i_ener + 1
                this%yf(i_ener) = ener
                this%lf(i_ener) = ie
                !sigma(i_ener) = sgm(1)
                this%target_type(i_ener) = 3*n_config_type + 1
             endif
          endselect
       endif

       if(allocated(vec)) deallocate(vec)
       if(allocated(dvec)) deallocate(dvec)
       if(allocated(jack)) deallocate(jack)
       if(allocated(w)) deallocate(w)
       if(.not. has_mask) deallocate(mask)

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

    if( this%do_core ) call Initialise(core_pot, this%ip_args, param_str=string(this%quip_string))

    call initialise(xyzfile,this%at_file)

    n_max = xyzfile%n_frame

    n_ener = 0
    this%e0 = 0.0_dp

    do n_con = 1, n_max
       call read(xyzfile,at,frame=n_con-1)

       has_ener = get_value(at%params,trim(this%energy_parameter_name),ener)

       if( has_ener ) then

          ener_core = 0.0_dp
          if( this%do_core ) then
             call set_cutoff(at, cutoff(core_pot))
             call calc_connect(at)
             call calc(core_pot,at,energy=ener_core)
          endif

          select case(trim(this%coordinates))
          case('water_monomer','water_dimer','hf_dimer')
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
    this%w_Z = 0.0_dp
    select case(trim(this%coordinates))
    case('water_monomer','water_dimer')
       this%w_Z(8) = 1.0_dp
    case('hf_dimer')
       this%w_Z(9) = 1.0_dp
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
     character(len=STRING_LENGTH) :: gp_tmp_file, gp_label
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
     call xml_AddAttribute(xf,"svn_version",""//current_version())

     call xml_NewElement(xf,"GAP_data")
     call xml_AddAttribute(xf,"n_species",""//this%n_species)
     call xml_AddAttribute(xf,"do_core",""//this%do_core)
     call xml_AddAttribute(xf,"e0",""//this%e0)
     call xml_AddAttribute(xf,"f0",""//this%f0)
     call xml_AddAttribute(xf,"do_pca",""//this%do_pca)

     select case(trim(this%coordinates))
     case('hf_dimer')
        call xml_AddAttribute(xf,"coordinates","hf_dimer")
        call xml_NewElement(xf,"hf_dimer_params")
        call xml_AddAttribute(xf,"cutoff",""//this%r_cut)
        call xml_EndElement(xf,"hf_dimer_params")
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
     case('cosnx')
        call xml_AddAttribute(xf,"coordinates","cosnx")
     
        call xml_NewElement(xf,"cosnx_params")
        call xml_AddAttribute(xf,"cutoff",""//this%r_cut)
        call xml_AddAttribute(xf,"l_max",""//this%cosnx_l_max)
        call xml_AddAttribute(xf,"n_max",""//this%cosnx_n_max)
     
        call xml_NewElement(xf,"NormFunction")
        call xml_AddCharacters(xf, ""//this%NormFunction)
        call xml_EndElement(xf,"NormFunction")
     
        call xml_NewElement(xf,"RadialTransform")
        do i = 1, this%cosnx_n_max
           call xml_NewElement(xf,"RadialTransform_row")
           call xml_AddAttribute(xf,"i",""//i)
           call xml_AddCharacters(xf, ""//this%RadialTransform(:,i))
           call xml_EndElement(xf,"RadialTransform_row")
        enddo
        call xml_EndElement(xf,"RadialTransform")
     
        call xml_EndElement(xf,"cosnx_params")

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

  subroutine print_sparse(this)
    type(teach_sparse), intent(in) :: this
    type(cinoutput) :: xyzfile, xyzfile_out
    type(atoms) :: at, at_out

    integer :: li, ui, n_con
    logical, dimension(:), allocatable :: x
    logical, dimension(:), pointer :: sparse

    if(this%do_mark_sparse_atoms) then

       allocate(x(this%nn))
       x = .false.
       x(this%r) = .true.

       call initialise(xyzfile,this%at_file)
       call initialise(xyzfile_out,this%mark_sparse_atoms,action=OUTPUT)

       li = 0
       ui = 0
       do n_con = 1, xyzfile%n_frame
          call read(xyzfile,at,frame=n_con-1)
          at_out = at

          call add_property(at_out,'sparse',.false.,ptr=sparse)

          li = ui + 1
          ui = ui + at%N
          if(any( x(li:ui) )) sparse(find_indices(x(li:ui))) = .true.

          call write(at_out,xyzfile_out,properties="species:pos:sparse")
       enddo
       call finalise(xyzfile)
       call finalise(xyzfile_out)
       deallocate(x)

    endif

  endsubroutine print_sparse

end module teach_sparse_mod
