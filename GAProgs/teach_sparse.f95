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
  use IPEwald_module
  use QUIP_module

  implicit none

  type ip_core
     logical :: do_core = .false.
     character(len=10000) :: quip_string = ''
     character(len=FIELD_LENGTH) :: ip_args = ''
  endtype ip_core
     
     
  interface teach_data_from_xyz
     module procedure xyzfile_teach_data_from_xyz_so4, xyzfile_teach_data_from_xyz_qw_so3
  endinterface teach_data_from_xyz

  private

  public :: teach_n_from_xyz
  public :: teach_data_from_xyz
  public :: e0_avg_from_xyz
  public :: w_Z_from_xyz
  public :: ip_core

  character(len=FIELD_LENGTH)::energy_property_name, force_property_name, virial_property_name
  public :: energy_property_name, force_property_name, virial_property_name

contains

  subroutine teach_n_from_xyz(at_file, r_cut, n, nn, ne, n_ener, n_force, n_virial, species_present, n_species)

    character(len=FIELD_LENGTH), intent(in) :: at_file
    real(dp), intent(in) :: r_cut
    integer, intent(out) :: n, nn, ne, n_ener, n_force, n_virial, species_present(116), n_species

    type(cinoutput) :: xyzfile
    type(atoms) :: at
    integer :: n_max, n_con
    logical :: has_ener, has_force, has_virial
    real(dp) :: ener, virial(3,3)
    real(dp), pointer :: f(:,:)
    integer :: i

    call initialise(xyzfile,at_file)
    call query(xyzfile)

    n_max = xyzfile%n_frame

    n = 0
    nn = 0
    ne = 0
    n_ener = 0
    n_force = 0
    n_virial = 0
    species_present = 0
    n_species = 0

    do n_con = 1, n_max
       call read(xyzfile,at,frame=n_con-1)
       call add_property(at,'charge',0.0_dp,n_cols=1)

       has_ener = get_value(at%params,energy_property_name,ener)
       has_force = assign_pointer(at,force_property_name, f)
       has_virial = get_value(at%params,virial_property_name,virial)

       if( has_ener .or. has_force .or. has_virial ) then
          call atoms_set_cutoff(at,r_cut)
          call calc_connect(at)
       endif

       if( has_ener ) then
          n_ener = n_ener + 1
          ne = ne + at%N
       endif

       if( has_force ) then
          n_force = n_force + at%N*3
          do i = 1, at%N
             n = n + 3*(atoms_n_neighbours(at,i)+1)
          enddo
       endif

       if( has_virial ) then
          n_virial = n_virial + 9
          do i = 1, at%N
             n = n + 9*(atoms_n_neighbours(at,i)+1)
          enddo
       endif

       if( has_ener .or. has_force .or. has_virial ) nn = nn + at%N

       do i = 1, at%N
          if( all(at%Z(i) /= species_present) ) then
             n_species = n_species + 1
             species_present(n_species) = at%Z(i)
          endif
       enddo

       call finalise(at)
    enddo

    call finalise(xyzfile)

  end subroutine teach_n_from_xyz

  subroutine xyzfile_teach_data_from_xyz_so4(at_file, f_hat, df_hat, r_cut, do_ewald, do_ewald_corr, z_eff, core, e0, w_Z, &
                                     x, xd, yf, ydf, lf, ldf, xf, xdf, xz, target_type )

    character(len=FIELD_LENGTH), intent(in) :: at_file
    type(fourier_so4), intent(inout) :: f_hat
    type(grad_fourier_so4), intent(inout) :: df_hat
    real(dp), intent(in) :: r_cut, z_eff(116), e0, w_Z(:)
    logical, intent(in) :: do_ewald, do_ewald_corr
    type(ip_core), intent(in) :: core
    real(dp), intent(out) :: x(:,:), xd(:,:), yf(:), ydf(:)
    integer, intent(out) :: lf(:), ldf(:), xf(:), xdf(:), xz(:), target_type(:)

    call xyzfile_teach_data_from_xyz(at_file, r_cut, do_ewald, do_ewald_corr, z_eff, core, e0, w_Z, &
                                     x, xd, yf, ydf, lf, ldf, xf, xdf, xz, target_type, &
                                     f_hat = f_hat, df_hat = df_hat)

  end subroutine xyzfile_teach_data_from_xyz_so4

  subroutine xyzfile_teach_data_from_xyz_qw_so3(at_file, f3_hat, df3_hat, qw, dqw, r_cut, do_ewald, do_ewald_corr, z_eff, core, e0, w_Z, &
                                                x, xd, yf, ydf, lf, ldf, xf, xdf, xz , target_type)

    character(len=FIELD_LENGTH), intent(in) :: at_file
    type(fourier_so3), intent(inout) :: f3_hat
    type(grad_fourier_so3), intent(inout) :: df3_hat
    type(qw_so3), intent(inout) :: qw
    type(grad_qw_so3), intent(inout) :: dqw
    real(dp), intent(in) :: r_cut, z_eff(116), e0, w_Z(:)
    logical, intent(in) :: do_ewald, do_ewald_corr
    type(ip_core), intent(in) :: core
    real(dp), intent(out) :: x(:,:), xd(:,:), yf(:), ydf(:)
    integer, intent(out) :: lf(:), ldf(:), xf(:), xdf(:), xz(:), target_type(:)

    call xyzfile_teach_data_from_xyz(at_file, r_cut, do_ewald, do_ewald_corr, z_eff, core, e0, w_Z, &
                                     x, xd, yf, ydf, lf, ldf, xf, xdf, xz, target_type, &
                                     f3_hat = f3_hat, df3_hat = df3_hat, qw = qw, dqw = dqw)

  end subroutine xyzfile_teach_data_from_xyz_qw_so3

  subroutine xyzfile_teach_data_from_xyz(at_file, r_cut, do_ewald, do_ewald_corr, z_eff, core, e0, w_Z, &
                                           x, xd, yf, ydf, lf, ldf, xf, xdf, xz, target_type, &
                                           f_hat, df_hat, f3_hat, df3_hat, qw, dqw)

    character(len=FIELD_LENGTH), intent(in) :: at_file
    real(dp), intent(in) :: r_cut, z_eff(116), e0, w_Z(:)
    logical, intent(in) :: do_ewald, do_ewald_corr
    type(ip_core), intent(in) :: core
    real(dp), intent(out) :: x(:,:), xd(:,:), yf(:), ydf(:)
    integer, intent(out) :: lf(:), ldf(:), xf(:), xdf(:), xz(:), target_type(:)
    type(fourier_so4), intent(inout), optional :: f_hat
    type(grad_fourier_so4), intent(inout), optional :: df_hat
    type(fourier_so3), intent(inout), optional :: f3_hat
    type(grad_fourier_so3), intent(inout), optional :: df3_hat
    type(qw_so3), intent(inout), optional :: qw
    type(grad_qw_so3), intent(inout), optional :: dqw

    type(cinoutput) :: xyzfile
    type(atoms) :: at
    type(bispectrum_so4) :: bis
    type(grad_bispectrum_so4) :: dbis

    type(Potential) :: core_pot
    
    logical :: do_qw_so3
    integer :: d
    integer :: n_max, n_con
    integer :: n_ener
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

    if (present(f_hat)) then
       do_qw_so3 = .false.
       d = j_max2d(f_hat%j_max)
    elseif (present(f3_hat) .and. present(qw)) then
       do_qw_so3 = .true.
       d = qw2d(qw)
    endif

    if( core%do_core ) call Initialise(core_pot, core%ip_args, core%quip_string)

    call initialise(xyzfile,at_file)
    call query(xyzfile)

    n_max = xyzfile%n_frame

    n_ener = size(yf)

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

       has_ener = get_value(at%params,energy_property_name,ener)
       has_force = assign_pointer(at,force_property_name, f)
       has_virial = get_value(at%params,virial_property_name,virial)

       if( core%do_core ) then
          allocate(f_core(3,at%N))
          ener_core = 0.0_dp
          f_core = 0.0_dp
          virial_core = 0.0_dp

          call set_cutoff(at, max(cutoff(core_pot),r_cut))
          call calc_connect(at)

          call calc(core_pot,at,e=ener_core,f=f_core,virial=virial_core)

          if(has_ener) ener = ener - ener_core
          if(has_force) f = f - f_core
          if(has_virial) virial = virial - virial_core
          deallocate(f_core)
       endif

       if( do_ewald ) then
          allocate(f_ewald(3,at%N))
          if( .not. assign_pointer(at, 'charge', charge) ) call system_abort('Could not assign pointer')
          do i = 1, at%N
             charge(i) = z_eff(at%Z(i))
          enddo
          call Ewald_calc(at,e=ener_ewald,f=f_ewald,virial=virial_ewald)
          if(has_ener) ener = ener - ener_ewald
          if(has_force) f = f - f_ewald
          if(has_virial) virial = virial - virial_ewald
          deallocate(f_ewald)
          if( do_ewald_corr ) then
             allocate(f_ewald_corr(3,at%N))
             call Ewald_corr_calc(at,e=ener_ewald_corr,f=f_ewald_corr,virial=virial_ewald_corr)
             if(has_ener) ener = ener + ener_ewald_corr
             if(has_force) f = f + f_ewald_corr
             if(has_virial) virial = virial + virial_ewald_corr
             deallocate(f_ewald_corr)
          endif
       endif

       if(has_ener) ener = ener - at%N*e0     

       if( at%cutoff /= r_cut ) then
          call atoms_set_cutoff(at,r_cut)
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
          w(i) = w_Z(at%Z(i))
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
             x(:,ix) = vec(:,i)
             xz(ix) = at%Z(i)

             if( has_ener ) then
                ie = ie + 1
                xf(ie) = ix
             endif
        
             if(has_force) then
                do k = 1, 3
                   li = ui + 1
                   ui = ui + 1
                   nn = nn+1
                   xdf(ui) = ix
                   xd(:,ui) = jack(:,k,i)
                   do n = 1, atoms_n_neighbours(at,i)
                      j = atoms_neighbour(at,i,n,jn=jn)
                      ui = ui + 1
                      jx = i_con + j
                      xdf(ui) = jx
                      xd(:,ui) = jack(:,jn*3+k,j)
                   enddo

                   ydf(nn) = -f(k,i)
                   ldf(nn) = ui
                   !sigma(nn+n_ener) = sgm(2)
                   target_type(nn+n_ener) = 2
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
                      xdf(ui) = ix
                      xd(:,ui) = jack(:,k,i)*at%pos(l,i)
                      do n = 1, atoms_n_neighbours(at,i)
                         j = atoms_neighbour(at,i,n,jn=jn,shift=shift)
                         ui = ui + 1
                         jx = i_con + j
                         xdf(ui) = jx
                         pos = at%pos(:,i) - matmul(at%lattice,shift)
                         xd(:,ui) = jack(:,jn*3+k,j)*pos(l)
                      enddo
                   enddo

                   ydf(nn) = -virial(l,k) 
                   ldf(nn) = ui
                   !sigma(nn+n_ener) = sgm(3)
                   target_type(nn+n_ener) = 3
                enddo
             enddo
          endif

          i_con = i_con + at%N

          if(has_ener) then
             i_ener = i_ener + 1
             yf(i_ener) = ener
             lf(i_ener) = ie
             !sigma(i_ener) = sgm(1)
             target_type(i_ener) = 1
          endif
       endif

       if(allocated(vec)) deallocate(vec)
       if(allocated(jack)) deallocate(jack)
       if(allocated(w)) deallocate(w)

       call finalise(at)
    enddo

    call finalise(xyzfile)

    if( core%do_core ) then
       call Finalise(core_pot)
    endif

  end subroutine xyzfile_teach_data_from_xyz

  subroutine e0_avg_from_xyz(at_file, do_ewald, do_ewald_corr, z_eff, core, e0)

    character(len=FIELD_LENGTH), intent(in) :: at_file
    real(dp), intent(in) :: z_eff(116)
    logical, intent(in) :: do_ewald, do_ewald_corr
    type(ip_core), intent(in) :: core
    real(dp), intent(out) :: e0

    type(Potential) :: core_pot
    type(cinoutput) :: xyzfile
    type(atoms) :: at
    integer :: n_max, n_con
    integer :: n_ener
    logical :: has_ener
    real(dp) :: ener, ener_ewald, ener_ewald_corr, ener_core
    real(dp), pointer :: charge(:)
    integer :: i

    if( core%do_core ) call Initialise(core_pot, core%ip_args, core%quip_string)

    call initialise(xyzfile,at_file)
    call query(xyzfile)

    n_max = xyzfile%n_frame

    n_ener = 0
    e0 = 0.0_dp

    do n_con = 1, n_max
       call read(xyzfile,at,frame=n_con-1)
       call add_property(at,'charge',0.0_dp,n_cols=1)

       has_ener = get_value(at%params,'Energy',ener)

       if( has_ener ) then

          ener_core = 0.0_dp
          if( core%do_core ) then
             call set_cutoff(at, cutoff(core_pot))
             call calc_connect(at)
             call calc(core_pot,at,e=ener_core)
          endif

          ener_ewald = 0.0_dp
          ener_ewald_corr = 0.0_dp
          if( do_ewald ) then
             if( .not. assign_pointer(at, 'charge', charge) ) call system_abort('Could not assign pointer')
             do i = 1, at%N
                charge(i) = z_eff(at%Z(i))
             enddo
             call Ewald_calc(at,e=ener_ewald)
             if( do_ewald_corr ) call Ewald_corr_calc(at,e=ener_ewald_corr)
          endif

          e0 = e0 + (ener-ener_ewald+ener_ewald_corr-ener_core) / at%N

          n_ener = n_ener + 1
       endif

       call finalise(at)
    enddo

    if( n_ener > 0 ) e0 = e0 / n_ener

    call finalise(xyzfile)

    if( core%do_core ) call Finalise(core_pot)

  end subroutine e0_avg_from_xyz

  subroutine w_Z_from_xyz(at_file, w_Z)

    character(len=FIELD_LENGTH), intent(in) :: at_file
    real(dp), intent(out) :: w_Z(:)

    type(cinoutput) :: xyzfile
    type(atoms) :: at

    call initialise(xyzfile,at_file)

    call read(xyzfile,at,frame=0)
    call get_weights(at,w_Z)
    call finalise(at)

    call finalise(xyzfile)

  end subroutine w_Z_from_xyz

end module teach_sparse_module

program teach_sparse

  use teach_sparse_module
  use libatoms_module
  use bispectrum_module
  use gp_sparse_module
  use clustering_module

  implicit none

  integer, parameter :: SPARSE_LENGTH = 10000
  integer, parameter :: SPARSE_N_FIELDS = 2000
  integer, parameter :: THETA_LENGTH = 10000
  real(dp), parameter :: THETA_MIN = 0.000000001

  type(inoutput) :: bispectrum_inout, theta_inout, sparse_inout
  type(Dictionary) :: params
  type(fourier_so4) :: f_hat
  type(grad_fourier_so4) :: df_hat
  type(fourier_so3) :: f3_hat
  type(grad_fourier_so3) :: df3_hat
  type(qw_so3) :: qw
  type(grad_qw_so3) :: dqw
  type(gp) :: my_gp
  type(gp_sparse) :: gp_sp

  character(len=FIELD_LENGTH) :: at_file, qw_cutoff_string, qw_cutoff_f_string, qw_cutoff_r1_string, theta_file, sparse_file, z_eff_string, bispectrum_file, ip_args
  integer :: m, j_max, qw_l_max, min_steps, min_save
  real(dp) :: r_cut, z0, e0, f0, sgm(3), dlt, theta_fac, mem_required, mem_available
  logical :: do_qw_so3, qw_no_q, qw_no_w, has_e0, has_f0, has_theta_file, has_sparse_file, &
  & do_sigma, do_delta, do_theta, do_sparx, do_f0, do_theta_fac, do_test_gp_gradient, has_bispectrum_file, do_cluster, do_pivot, &
  & do_core, do_ewald, do_ewald_corr, test_gp_gradient_result

  real(dp), dimension(116) :: z_eff
  character(len=FIELD_LENGTH), dimension(232) :: z_eff_fields
  integer :: num_z_eff_fields
  integer :: d
  character(len=FIELD_LENGTH), dimension(99) :: qw_cutoff_fields, qw_cutoff_f_fields, qw_cutoff_r1_fields
  integer :: qw_f_n
  real(dp), dimension(99) :: qw_cutoff, qw_cutoff_r1
  integer, dimension(99) :: qw_cutoff_f
  integer :: n, nn, ne, n_ener, n_force, n_virial
  integer, dimension(116) :: species_present
  integer :: n_species
  integer, dimension(:), allocatable :: species_Z
  real(dp), dimension(:), allocatable :: w_Z
  real(dp), dimension(:,:), allocatable :: x, xd
  real(dp), dimension(:), allocatable :: yf, ydf
  integer, dimension(:), allocatable :: lf, ldf, xf, xdf, xz, target_type 
  real(dp), dimension(:), allocatable :: dlta
  real(dp), dimension(:,:), allocatable :: theta
  integer, dimension(:), allocatable :: r
  character(len=SPARSE_LENGTH) :: sparse_string
  character(len=FIELD_LENGTH), dimension(:), allocatable :: sparse_string_array
  character(len=THETA_LENGTH) :: theta_string
  character(len=FIELD_LENGTH), dimension(:), allocatable :: theta_string_array
  integer :: i, j, k, l, o, dd, dt
  character(len=FIELD_LENGTH) :: gp_file

  type(extendable_str)  :: quip_params_str
  type(ip_core) :: core

  call system_initialise(verbosity=NORMAL)
  call initialise(params)
  call param_register(params, 'at_file', PARAM_MANDATORY, at_file)
  call param_register(params, 'm', '50', m)
  call param_register(params, 'r_cut', '2.75', r_cut)
  call param_register(params, 'j_max', '4', j_max)
  call param_register(params, 'z0', '0.0', z0)
  call param_register(params, 'qw_so3', 'F', do_qw_so3)
  call param_register(params, 'l_max', '6', qw_l_max)
  call param_register(params, 'cutoff', '', qw_cutoff_string)
  call param_register(params, 'cutoff_f', '', qw_cutoff_f_string)
  call param_register(params, 'cutoff_r1', '', qw_cutoff_r1_string)
  call param_register(params, 'no_q', 'F', qw_no_q)
  call param_register(params, 'no_w', 'F', qw_no_w)
  call param_register(params, 'e0', '0.0', e0, has_e0)
  call param_register(params, 'f0', '0.0', f0, has_f0)
  call param_register(params, 'sgm', '0.1 0.1 0.1', sgm)
  call param_register(params, 'dlt', '1.0', dlt)
  call param_register(params, 'theta_file', '', theta_file, has_theta_file)
  call param_register(params, 'sparse_file', '', sparse_file, has_sparse_file)
  call param_register(params, 'theta_fac', '3.0', theta_fac)
  call param_register(params, 'do_sigma', 'F', do_sigma)
  call param_register(params, 'do_delta', 'F', do_delta)
  call param_register(params, 'do_theta', 'F', do_theta)
  call param_register(params, 'do_sparx', 'F', do_sparx)
  call param_register(params, 'do_f0', 'F', do_f0)
  call param_register(params, 'do_theta_fac', 'F', do_theta_fac)
  call param_register(params, 'do_cluster', 'F', do_cluster)
  call param_register(params, 'do_pivot', 'F', do_pivot)
  call param_register(params, 'min_steps', '10', min_steps)
  call param_register(params, 'min_save', '0', min_save)
  call param_register(params, 'z_eff', '', z_eff_string,do_ewald)
  call param_register(params, 'do_test_gp_gradient', 'F', do_test_gp_gradient)
  call param_register(params, 'bispectrum_file', '', bispectrum_file, has_bispectrum_file)
  call param_register(params, 'ip_args', '', ip_args, do_core)
  call param_register(params, 'do_ewald_corr', 'F', do_ewald_corr)
  call param_register(params, 'energy_property_name', 'energy', energy_property_name)
  call param_register(params, 'force_property_name', 'force', force_property_name)
  call param_register(params, 'virial_property_name', 'virial', virial_property_name)

  if (.not. param_read_args(params, do_check = .true.)) then
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

  if( count( (/has_sparse_file,do_cluster,do_pivot/) ) > 1 ) &
  & call system_abort('There has been more than one method specified for sparsification.')
     
  z_eff = 0.0_dp
  if(do_ewald) then
     call parse_string(z_eff_string,':',z_eff_fields,num_z_eff_fields)
     do i = 1, num_z_eff_fields, 2
        j = atomic_number_from_symbol(z_eff_fields(i))
        if(j < 1 .or. j > 116) call system_abort("Invalid atomic number "//j//" parsed from "//z_eff_fields(i))
        z_eff(j) = string_to_real(z_eff_fields(i+1))
     enddo
  endif
  do_ewald_corr = do_ewald .and. do_ewald_corr

  if (do_qw_so3) then
     qw_cutoff = 0.0_dp
     qw_cutoff_f = 0
     qw_cutoff_r1 = 0.0_dp
     call parse_string(qw_cutoff_string, ':', qw_cutoff_fields, qw_f_n)
     call parse_string(qw_cutoff_f_string, ':', qw_cutoff_f_fields, qw_f_n)
     call parse_string(qw_cutoff_r1_string, ':', qw_cutoff_r1_fields, qw_f_n)
     do i = 1, qw_f_n
        qw_cutoff(i) = string_to_real(qw_cutoff_fields(i))
        qw_cutoff_f(i) = string_to_int(qw_cutoff_f_fields(i))
        qw_cutoff_r1(i) = string_to_real(qw_cutoff_r1_fields(i))
     enddo

     call initialise(f3_hat, qw_l_max, qw_cutoff(1:qw_f_n), qw_cutoff_f(1:qw_f_n), qw_cutoff_r1(1:qw_f_n))
     call initialise(df3_hat, qw_l_max, qw_cutoff(1:qw_f_n), qw_cutoff_f(1:qw_f_n), qw_cutoff_r1(1:qw_f_n))
     call initialise(qw, qw_l_max, qw_f_n, do_q = (.not. qw_no_q), do_w = (.not. qw_no_w))
     call initialise(dqw, qw_l_max, qw_f_n, do_q = (.not. qw_no_q), do_w = (.not. qw_no_w))

     d = qw2d(qw)

     r_cut = maxval(qw_cutoff(1:qw_f_n))
  else
     z0 = max(1.0_dp,z0) * r_cut/(PI-0.02_dp)

     call initialise(f_hat,j_max,z0,r_cut)
     call initialise(df_hat,j_max,z0,r_cut)

     d = j_max2d(j_max)
  endif

  if(do_core) then
     call read(quip_params_str, "quip_params.xml")
     core%ip_args = ip_args
     core%do_core = do_core
     core%quip_string = string(quip_params_str)
  endif

  call teach_n_from_xyz(at_file, r_cut, n, nn, ne, n_ener, n_force, n_virial, species_present, n_species)

  allocate(species_Z(n_species))
  species_Z = species_present(1:n_species)

  if (.not. has_f0) then
     call e0_avg_from_xyz(at_file, do_ewald, do_ewald_corr, z_eff, core, f0)
  end if

  allocate(w_Z(maxval(species_Z)))

  call w_Z_from_xyz(at_file, w_Z)

  allocate(x(d,nn),xd(d,n),yf(n_ener),ydf(n_force+n_virial),lf(n_ener),ldf(n_force+n_virial),xf(ne),xdf(n))
  allocate(xz(nn))
  allocate(target_type(n_ener+n_force+n_virial))

  if (do_qw_so3) then
     call teach_data_from_xyz(at_file, f3_hat, df3_hat, qw, dqw, r_cut, do_ewald, do_ewald_corr, z_eff, core, e0, w_Z, &
                              x, xd, yf, ydf, lf, ldf, xf, xdf, xz, target_type)
  else
     call teach_data_from_xyz(at_file, f_hat, df_hat, r_cut, do_ewald, do_ewald_corr, z_eff, core, e0, w_Z, &
                              x, xd, yf, ydf, lf, ldf, xf, xdf, xz, target_type)
  endif

  if( has_bispectrum_file ) then
     call initialise(bispectrum_inout,bispectrum_file,action=OUTPUT)
     do i = 1, nn/3
        write(bispectrum_inout%unit,"("//d//"f16.8)") x(:,i)
     enddo
     call finalise(bispectrum_inout)
  endif

  allocate(dlta(n_species))

  dlta = dlt 

  if( has_sparse_file ) then
     allocate(sparse_string_array(SPARSE_N_FIELDS))
     call initialise(sparse_inout,sparse_file)
     read(sparse_inout%unit,'(a)') sparse_string
     call parse_string(sparse_string,' ',sparse_string_array,m)
     allocate(r(m))
     do i = 1, m
        r(i) = string_to_int(sparse_string_array(i))
     enddo
     deallocate(sparse_string_array)
     call finalise(sparse_inout)
  elseif(do_cluster) then
     allocate(r(m))
     call bisect_kmedoids(x,m,med=r,theta_fac=theta_fac)
  elseif(do_pivot) then
     allocate(r(m))
     call pivot(x, r,theta_fac=theta_fac)
  else
     allocate(r(m))
     call fill_random_integer(r,size(x,2))
  endif
  call sort_array(r)

  call print('')
  call print('Atomic environments used in sparsification')
  call print(r)
  call print('')

  allocate(theta(d,n_species))

  if( has_theta_file ) then
     allocate(theta_string_array(d))
     call initialise(theta_inout,theta_file)
     read(theta_inout%unit,'(a)') theta_string
     call parse_string(theta_string,' ',theta_string_array,dt)
     if(d /= dt) call system_abort('File '//theta_file//'does not contain the right number of hyperparameters')
     do k = 1, n_species  
        do dd = 1, d
           theta(dd,k) = string_to_real(theta_string_array(dd+(k-1)*d))
        enddo
     enddo
     deallocate(theta_string_array)
     call finalise(theta_inout)
  else
     do k = 1, n_species
        do dd = 1, d
!           theta(dd,k) = ( maxval(x(dd,:),mask=(xz(:)==species_Z(k))) - minval(x(dd,:),mask=(xz(:)==species_Z(k))) )
           theta(dd,k) = ( maxval(x(dd,r),mask=(xz(r)==species_Z(k))) - minval(x(dd,r),mask=(xz(r)==species_Z(k))) )
!           theta(dd) = sqrt( & !take square root
!                          & sum( x(dd,:)**2 ) / size(x(dd,:)) - &
!                          & (sum( x(dd,:) ) / size(x(dd,:)))**2 )

           if( theta(dd,k) >= THETA_MIN ) then
              theta(dd,k) = theta_fac*theta(dd,k)
           else
              theta(dd,k) = 1.0_dp
           endif
        enddo
     enddo
  endif

  ! Stop execution if required memory is greater than the available memory. 
  ! The biggest arrays allocated are 2*sr*(nx+nxd), where sr is the
  ! number of sparse points, nx and nxd are the number of bispectra and partial
  ! derivatives.
  mem_required = 2.0_dp * real(size(r),dp) * (real(size(xf),dp) + real(size(xdf),dp)) * real(dp,dp) / (1024.0_dp**3)
  call mem_info(mem_available)
  mem_available = mem_available / (1024.0_dp**3)

  if( mem_required > mem_available ) call system_abort('Required memory ('//mem_required//' GB) exceeds available memory ('//mem_available//' GB).')

  call gp_sparsify(gp_sp,r,sgm,dlta,theta,yf,ydf,x,xd,xf,xdf,lf,ldf,xz,species_Z,(/(f0,i=1,n_species)/),target_type)

  deallocate(x,xd,xf,xdf,yf,ydf,lf,ldf)

  call print('')
  call print('theta')
  do l = 1, size(gp_sp%theta, 2)
     do o = 1, size(gp_sp%theta, 1)
        call print(real(gp_sp%theta(o,l),kind=dp))
     enddo
  enddo
  call print('')

  call enable_timing()
  
  if( do_test_gp_gradient ) then
     call verbosity_push(NERD)
     test_gp_gradient_result = test_gp_gradient(gp_sp,sigma=do_sigma,delta=do_delta,theta=do_theta,sparx=do_sparx,f0=do_f0,theta_fac=do_theta_fac)
     call verbosity_pop()
  endif

  ! Conjugate gradient minimiser's counter starts at 1, and stops when it reaches min_steps,
  ! so if min_steps is equal to 1, no iterations are made!
 
  if (min_save == 0) min_save = min_steps

  k = 0
  do i = 1, ((min_steps / min_save) + 1)
     if (k == min_steps) exit

     if ((min_steps - k) >= min_save) then
        j = minimise_gp_gradient(gp_sp,max_steps=(min_save + 1),sigma=do_sigma,delta=do_delta,theta=do_theta,sparx=do_sparx,f0=do_f0,theta_fac=do_theta_fac)
        k = k + min_save
     elseif ((min_steps - k) < min_save) then
        j = minimise_gp_gradient(gp_sp,max_steps=(min_steps - k + 1),sigma=do_sigma,delta=do_delta,theta=do_theta,sparx=do_sparx,f0=do_f0,theta_fac=do_theta_fac)
        k = min_steps
     endif

     call print('')
     call print(k // ' iterations completed:')
     !call print('sigma')
     !call print(real(gp_sp%sigma,kind=dp))
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

     if (do_qw_so3) then
        my_gp%comment = "coordinates=qw l_max="//qw_l_max//" f_n="//qw_f_n
        do l = 1, qw_f_n
           my_gp%comment = trim(my_gp%comment)//" cutoff_"//l//"="//qw_cutoff(l)
           my_gp%comment = trim(my_gp%comment)//" cutoff_f_"//l//"="//qw_cutoff_f(l)
           my_gp%comment = trim(my_gp%comment)//" cutoff_r1_"//l//"="//qw_cutoff_r1(l)
        enddo
        if (.not. qw_no_q) then
           my_gp%comment = trim(my_gp%comment)//" do_q=T"
        else
           my_gp%comment = trim(my_gp%comment)//" do_q=F"
        endif
        if (.not. qw_no_w) then
           my_gp%comment = trim(my_gp%comment)//" do_w=T"
        else
           my_gp%comment = trim(my_gp%comment)//" do_w=F"
        endif
     else
        my_gp%comment = "coordinates=bispectrum cutoff="//r_cut//" j_max="//j_max//" z0="//z0//" n_species="//n_species//" Z={"//species_Z//&
        & "} w={"//w_Z(species_Z)//"} do_ewald="//do_ewald//" do_ewald_corr="//do_ewald_corr//" z_eff={"//z_eff(species_Z)//"}"
     endif
  
     if( core%do_core ) then
        my_gp%comment = trim(my_gp%comment)//" do_core=T ip_args={"//trim(core%ip_args)//"} quip_string={"//trim(core%quip_string)//"}"
     endif
     my_gp%comment = trim(my_gp%comment)//" sgm={"//sgm//"}"
     my_gp%comment = trim(my_gp%comment)//" e0="//e0//" f0="//f0

     gp_file = 'gp_'//m//'_'//k//'.dat'

     call gp_print_binary(my_gp,trim(gp_file))
     call system_command('ln -fs '//trim(gp_file)//' gp.dat')

     call finalise(my_gp)
  enddo

  call print("model parameters:")
  call print("r_cut     = "//r_cut)
  if (do_qw_so3) then
     call print("l_max     = "//qw_l_max)
     call print("cutoff    = "//qw_cutoff_string)
     call print("cutoff_f  = "//qw_cutoff_f_string)
     call print("cutoff_r1 = "//qw_cutoff_r1_string)
     call print("q         = "//(.not. qw_no_q))
     call print("w         = "//(.not. qw_no_w))
  else
     call print("j_max     = "//j_max)
     call print("z0        = "//z0)
  endif
  call print("n_species = "//n_species)
  call print("species_Z = "//species_Z)
  call print("w         = "//w_Z(species_Z))
  call print("z_eff     = "//z_eff(species_Z))
  call print("do_ewald  = "//do_ewald)
  call print("do_ewald_corr  = "//do_ewald_corr)
  call print("e0        = "//e0)

  call finalise(gp_sp)

  call system_finalise()

end program teach_sparse
