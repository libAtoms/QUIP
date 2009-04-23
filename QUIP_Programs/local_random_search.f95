!%%%%%%%%%%%%%%%%%%  custom fast LJ
module mylj
use libAtoms_module

implicit none

contains

pure function lj(x, data)
  real(dp), intent(in) :: x(:)
  character,optional, intent(in)::data(:)
  real(dp) :: lj
  real(dp) r2, tmp
  integer i,j

  lj = 0.0_dp

  do i=1,size(x)/3
     do j=i+1,size(x)/3
        r2 = sum((x((i-1)*3+1:(i-1)*3+3)-x((j-1)*3+1:(j-1)*3+3))**2)
        tmp = (1/r2)**3
        lj = lj+tmp*tmp-tmp
     end do
  end do
end function lj

pure function dlj(x, data)
  real(dp), intent(in) :: x(:)
  character,optional, intent(in)::data(:)
  real(dp) :: dlj(size(x))
  real(dp) r2, tmp, dr(3), r
  integer i,j


  dlj = 0.0_dp

  do i=1,size(x)/3
     do j=i+1,size(x)/3
        r2 = sum((x((i-1)*3+1:(i-1)*3+3)-x((j-1)*3+1:(j-1)*3+3))**2)
        r = sqrt(r2)
        dr = (x((i-1)*3+1:(i-1)*3+3)-x((j-1)*3+1:(j-1)*3+3))/r
        tmp = (1/r2)**3
        dlj((i-1)*3+1:(i-1)*3+3) = dlj((i-1)*3+1:(i-1)*3+3)+(-12.0_dp*tmp*tmp+6.0_dp*tmp)/r*dr
        dlj((j-1)*3+1:(j-1)*3+3) = dlj((j-1)*3+1:(j-1)*3+3)-(-12.0_dp*tmp*tmp+6.0_dp*tmp)/r*dr
     end do
  end do
  
end function dlj

pure function elj(x)
  real(dp), intent(in) :: x(:)
  real(dp) :: elj(size(x)/3)
  real(dp) r2, tmp
  integer i,j

  elj = 0.0_dp

  do i=1,size(x)/3
     do j=1,size(x)/3
        if(i==j) cycle
        r2 = sum((x((i-1)*3+1:(i-1)*3+3)-x((j-1)*3+1:(j-1)*3+3))**2)
        tmp = (1/r2)**3
        elj(i) = elj(i)+0.5*(tmp*tmp-tmp)
     end do
  end do
end function elj



end module mylj
!%%%%%%%%%%%%%%%%%% end of custom fast LJ

program local_random_search

use libAtoms_module
use Potential_module
use MetaPotential_module
use mylj

implicit none

  type(Potential) pot
  type(MetaPotential) mpot
  type(inoutput) params, movie
  type(Atoms), target:: at


  integer, parameter :: N = 38

  real(dp) :: r(3), r0(3), etmp(N), e(N), v, eold, rr
  real(dp), pointer :: atlocale(:)
  real(dp), allocatable :: x(:)
  logical :: allok, status
  integer :: i, jj, j, k, niter, count, nrand, idx(N)



  call system_initialise()
  
!  call Initialise(params, "quip_params.xml")
!  call Initialise(pot, 'IP LJ', params)
!  call Initialise(pot, 'IP LJ', '<LJ_params n_types="1" label="default"><per_type_data type="1" atomic_num="1" /><per_pair_data type1="1" type2="1" sigma="1.0" eps6="1.0" eps12="1.0" cutoff="10.0" shifted="T" /></LJ_params>')

!  call Initialise(mpot, 'Simple', pot)
  
  call initialise(at, N, reshape((/100.0_dp,0.0_dp,0.0_dp,0.0_dp,100.0_dp,0.0_dp,0.0_dp,0.0_dp, 100.0_dp/), (/3,3/)))
  call set_cutoff(at, cutoff(pot)+0.2)
  call add_property(at, "local_e", 0.0_dp);
  call initialise(movie, 'movie.xyz', OUTPUT)
  call set_atoms(at, 1)
  status = assign_pointer(at, 'local_e', atlocale)

  allocate(x(N*3))

  do k=1,1000
     at%pos = 0.0_dp
     call randomise(at%pos, 6.0_dp)
     at%pos = at%pos+50.0_dp
     
     allok = .false.
     do while(allok .eqv. .false.)
        x = reshape(at%pos, (/at%N*3/))
        atlocale = elj(x)
!        call calc(mpot, at, local_e=e)
        call print (atlocale)
        allok = .true.
        do i=1,N
           if(atlocale(i) > 0.0_dp) then 
              call randomise(at%pos(:,i), 2.0_dp)
              allok = .false.
           end if
        end do
     end do


     x = reshape(at%pos, (/at%N*3/))
     atlocale = elj(x)
     
     do i=1,20
!        niter =  minim(mpot, at, 'cg', 1.0e-10_dp, 5000, &
!             do_pos=.true., do_print=.false., print_inoutput=movie)
!        call calc_connect(at)
!        call calc(mpot, at, local_e=atlocale)

        niter = minim(x, lj, dlj, 'cg', 1.0e-10_dp, 5000)
        call print('Relaxation finished after '//niter//' iterations')

        at%pos = reshape(x, (/3,at%N/))
        atlocale = elj(x)
        call print_xyz(at, movie, ('Energy='//sum(atlocale)), all_properties=.true.)
        
        !r = 0
        !call randomise(r, 0.7_dp)
        !at%pos(:,maxloc(e)) = at%pos(:,maxloc(e)) + r
        
        
        !do j=1,N
        !   if(e(j) > minv+(maxv-minv)*0.7) then 
        !      r = 0
        !      call randomise(r, 1.0_dp)
        !      at%pos(:,maxloc(e)) = at%pos(:,maxloc(e)) + r
        !   end if
        !end do     
        
        ! find out nrand-th highest energy atom
        etmp = elj(x)
        idx = (/ (i, i=1,size(etmp)) /)
        call sort_array(etmp, idx)

        rr = 10.0_dp ! randomize by this amount

        do j=1,1
           jj = idx(N-j+1)
           call print('Randomsing '//jj)
           r0 = x((jj-1)*3+1:(jj-1)*3+3)
           atlocale(jj) = 1.0_dp
           do while (atlocale(jj) > 0.0_dp)
              r = 0
              call randomise(r, rr)
              x((jj-1)*3+1:(jj-1)*3+3) = r0 + r
              atlocale = elj(x)
           end do
        enddo
        
     end do
  end do

  call system_finalise()

contains

end program local_random_search


