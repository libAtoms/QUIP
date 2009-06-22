program md

  use libatoms_module
  use bispectrum_module
  use gp_module

  implicit none

  type(atoms) :: at
  type(inoutput) :: xyzfile

  type(fourier_so4) :: f_hat

  type(bispectrum_so4) :: bis

  real(dp) :: ener, e_corr, r_cut, z0, e
  real(dp), dimension(:), allocatable :: theta, y
  real(dp), dimension(:,:), allocatable :: vec, x
  integer  :: i, n_con, stat, d, n_max, n_flds, j_max, dd

  type(gp) :: my_gp

  character(len=1024) :: ener_string
  character(len=1024), dimension(2) :: ener_string_parse
  
  call system_initialise(verbosity=NORMAL, enable_timing=.true.)

  if ( NUM_COMMAND_ARGS == 3 ) then
     call initialise(xyzfile,trim(COMMAND_ARG(1)))
     r_cut = string_to_real(trim(COMMAND_ARG(2)))
  else
     call system_abort('Usage: teach <xyzfile1> <r_cut>')
  endif
  j_max = 2
  z0 = r_cut / (PI-0.02_dp)

  call initialise(f_hat,j_max,z0)

  d = j_max2d(j_max)
  allocate( theta(2+d) )

  n_max = 0
  do 
     call read_xyz(at,xyzfile,status=stat)
     if( stat /= 0) exit
     n_max = n_max + 1
     
     call finalise(at)
  enddo
  call finalise(xyzfile)
  
  call initialise(xyzfile,trim(COMMAND_ARG(1)))

  allocate(x(d,n_max),y(n_max))

  do n_con = 1, n_max

     call read_xyz(at,xyzfile,status=stat,comment=ener_string)
     if( stat /= 0) exit

     call parse_string(ener_string,'=',ener_string_parse,n_flds)
     ener = string_to_real(ener_string_parse(2))

     call atoms_set_cutoff(at,r_cut)
     call calc_connect(at)

     allocate(vec(d,at%N))

     call calc_c(at,e=e_corr)
     ener = ener - e_corr
     
     do i = 1, at%N
        call fourier_transform(f_hat,at,i)
        call calc_bispectrum(bis,f_hat)
        call bispectrum2vec(bis,vec(:,i))
     enddo
     i = 1
     x(:,n_con) = vec(:,1)
     y(n_con) = ener / at%N

     deallocate(vec)
     call finalise(at)
  enddo

  call finalise(xyzfile)

  theta(1) = 0.1_dp
  theta(2) = 5.0_dp
  do dd = 1, d
!     theta(dd) = ( maxval(x(dd,:)) - minval(x(dd,:)) ) * 2
     theta(dd+2) = sqrt( & !take square root
                       & sum( x(dd,:)**2 ) / size(x(dd,:)) - &
                       & (sum( x(dd,:) ) / size(x(dd,:)))**2 )
     if( theta(dd+2) .fne. 0.0_dp ) then
        theta(dd+2) = theta(dd+2) * 10.0_dp
     else
        theta(dd+2) = 1.0_dp
     endif
  enddo
     
  call initialise(my_gp,'ard',theta,y,x)
  deallocate(x,y)

  call print('theta_init')
  call print(my_gp%theta)

  print*, minimise_gp_gradient(my_gp)

  call print('theta')
  call print(my_gp%theta)

  call initialise(xyzfile,trim(COMMAND_ARG(3)))
  do 

     call read_xyz(at,xyzfile,status=stat,comment=ener_string)
     if( stat /= 0) exit

     call parse_string(ener_string,'=',ener_string_parse,n_flds)
     ener = string_to_real(ener_string_parse(2))

     call atoms_set_cutoff(at,r_cut)
     call calc_connect(at)

     allocate(vec(d,at%N))

     call calc_c(at,e=e_corr)
    
     e = 0.0_dp 
     do i = 1, at%N
        call fourier_transform(f_hat,at,i)
        call calc_bispectrum(bis,f_hat)
        call bispectrum2vec(bis,vec(:,i))
        e = gp_mean(my_gp,vec(:,i)) + e
     enddo
     e = e + e_corr
     call print(ener//"  "//e//"  "//(e-ener))

     deallocate(vec)
     call finalise(at)
  enddo
  call system_finalise()

  contains

  function ff(xx,data)
     real(dp), dimension(:) :: xx
     real(dp) :: ff
     character(*),optional :: data(:)

     at%pos = reshape(xx,(/3,at%N/))
     call calc_c(at,e=ff)
  endfunction ff

  function dff(xx,data)
     real(dp), dimension(:) :: xx
     real(dp), dimension(size(xx)) :: dff
     character(*),optional :: data(:)

     real(dp), dimension(3,at%N) :: f

     at%pos = reshape(xx,(/3,at%N/))
     call calc_c(at,f=f)
     dff = - reshape(f,(/at%N*3/))
  endfunction dff

  subroutine calc_c(at,e,f)
     type(atoms), intent(in) :: at
     real(dp), intent(out), optional :: e
     real(dp), dimension(3,at%N), intent(inout), optional :: f
     integer :: i, j, n

     real(dp), dimension(19), parameter :: x_V = 0.2_dp*(/(i, i=0, 18)/) + 2.2_dp
     real(dp), dimension(19), parameter :: y_V = (/ 0.0_dp, 0.0470336_dp, 0.0766591_dp, 0.055763_dp, 0.022853_dp, &
     & 0.00556713_dp, -0.00119112_dp, -0.00407387_dp, -0.00502886_dp, -0.00483256_dp, -0.00376998_dp, &
     & -0.00247485_dp, -0.00135794_dp, -0.000915401_dp, -0.000781631_dp, -0.00053046_dp, &
     & -0.000272229_dp, -0.0000355333_dp, 0.0_dp /)
     real(dp), dimension(6), parameter :: x_rho =  0.2_dp*(/(i, i=0, 5)/) + 0.8_dp
     real(dp), dimension(6), parameter :: y_rho = (/-0.615676_dp, -0.314382_dp, &
     & -0.159363_dp, -0.0443928_dp, -0.00726638_dp, 0.0_dp /)
     real(dp) :: r_ij
     real(dp), dimension(3) :: u_ij, drho_i
     real(dp), dimension(:), allocatable :: rho

     type(spline) :: my_spline_V, my_spline_rho
     type(atoms) :: my_at

     call initialise(my_spline_V,x_V,y_V,0.0_dp,0.0_dp)
     call initialise(my_spline_rho,x_rho,y_rho,2.10618_dp,0.0_dp)

     my_at = at
     call set_cutoff(my_at,5.8_dp)
     call calc_connect(my_at)

     if(present(e)) e = 0.0_dp
     if(present(f)) f = 0.0_dp
     allocate(rho(my_at%N))
     rho = 0.0_dp

     do i = 1, my_at%N
        drho_i = 0.0_dp
        do n = 1, atoms_n_neighbours(my_at,i)
           j = atoms_neighbour(my_at,i,n,distance=r_ij,cosines=u_ij)

           if( r_ij < 1.8_dp ) then
              rho(i) = rho(i) + spline_value(my_spline_rho,r_ij)
              !rho(i) = rho(i) + 7.534_dp*exp(-3.43944_dp*r_ij)
              !if(i==1) print*,r_ij
              if(present(f)) drho_i = drho_i + spline_deriv(my_spline_rho,r_ij) * u_ij
           endif
           
           if( r_ij < 2.2_dp ) cycle
           if(present(e)) e = e + spline_value(my_spline_V,r_ij)*0.5_dp
           if(present(f)) f(:,i) = f(:,i) + spline_deriv(my_spline_V,r_ij) * u_ij
        enddo
        if(present(e)) e = e + u(rho(i))
        if(present(f)) f(:,i) = f(:,i) + du(rho(i)) * drho_i

!        if(i==1) print*,rho(i)
     enddo
           

     if(present(f)) then
        do i = 1, my_at%N
           do n = 1, atoms_n_neighbours(my_at,i)
              j = atoms_neighbour(my_at,i,n,distance=r_ij,cosines=u_ij)
              if( r_ij < 1.8_dp ) f(:,i) = f(:,i) + du(rho(j)) * spline_deriv(my_spline_rho,r_ij) * u_ij
           enddo
        enddo
     endif

     call finalise(my_spline_V)
     call finalise(my_spline_rho)
     call finalise(my_at)
     deallocate(rho)

  endsubroutine calc_c

  function u(x)
     real(dp), intent(in) :: x
     real(dp) :: u

     u = 10.3054*x + 39.7293*x**2 + 15.5667*x**3 + 1.96467*x**4
  endfunction u

  function du(x)
     real(dp), intent(in) :: x
     real(dp) :: du

     du = 10.3054 + 2*39.7293*x + 3*15.5667*x**2 + 4*1.96467*x**3
  endfunction du

end program md
