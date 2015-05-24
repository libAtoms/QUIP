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
use mylj

implicit none

  type(Potential) mpot
  type(inoutput) params
  type(cinoutput) movie
  type(Atoms), target:: at


  integer, parameter :: N = 13
  real, parameter :: SIGMA = 2.35_dp

  real(dp) :: r(3), r0(3), etmp(N), e(N), v, eold, rr
  real(dp), pointer :: atlocale(:)
  real(dp), allocatable :: x(:), x0(:), rrr(:)
  logical :: allok, status
  integer :: i, jj, j, k, niter, count, nrand, idx(N)



  call system_initialise()
  
!  call Initialise(params, "quip_params.xml")
!  call Initialise(mpot, 'IP LJ', params)
!  call Initialise(mpot, 'IP LJ', '<LJ_params n_types="1" label="default"><per_type_data type="1" atomic_num="1" /><per_pair_data type1="1" type2="1" sigma="1.0" eps6="1.0" eps12="1.0" cutoff="10.0" shifted="T" /></LJ_params>')

  call Initialise(mpot, 'FilePot command=./castep_driver.sh')
  
  call initialise(at, N, reshape((/10.0_dp,0.0_dp,0.0_dp,0.0_dp,10.0_dp,0.0_dp,0.0_dp,0.0_dp, 10.0_dp/), (/3,3/)))
  call set_cutoff(at, cutoff(mpot)+0.2)
  call initialise(movie, 'movie.xyz', OUTPUT)
  call set_atoms(at, 14)

  call add_property(at, "local_e", 0.0_dp);
  status = assign_pointer(at, 'local_e', atlocale)


  allocate(x(N*3), x0(N*3), rrr(N*3))

  do k=1,1000

     call random_initial_pos_chris(at, 5.0_dp)

     ! TIMING TEST
     
     !x = reshape(at%pos, (/at%N*3/))
     !eold = 0.0_dp
     !do i=1,1000000
     !   !call random_number(rrr)
     !   !call randomise(x, 1e-6_dp)
     !   x = x + 1.0e-10_dp
     !   eold = eold+lj(x)
     !end do
     !call print(eold)
     !stop


     !x0 = reshape(at%pos, (/at%N*3/))
     !atlocale = elj(x0)
     !call print_xyz(at, movie, ('ljenergy='//sum(atlocale)), all_properties=.true.)
     
     do i=1,10
!        niter =  minim(mpot, at, 'cg', 1.0e-10_dp, 5000, &
!             do_pos=.true., do_print=.false., print_inoutput=movie)
!        call calc_connect(at)
!        call calc(mpot, at, local_energy=atlocale)

        !x = x0
        !niter = minim(x, lj, dlj, 'cg', 1.0e-10_dp, 5000)
        !call print('Relaxation finished after '//niter//' iterations')
        !atlocale = elj(x)
        !at%pos = reshape(x, (/3,at%N/))
        !call print_xyz(at, movie, ('comment="quipcg" Energy='//sum(atlocale)), all_properties=.true.)
 

        !at%pos = reshape(x0, (/3,at%N/))

        call Calc(mpot, at, energy=eold)
        call read(at, 'filepot.0.out')
        call zero_sum(at%pos)
        call add_property(at, "local_e", 0.0_dp);
        status = assign_pointer(at, 'local_e', atlocale)
        x = reshape(at%pos, (/at%N*3/))
        atlocale = elj(x)
	call set_value(at%properties, "comment", "castepbfgs")
	call set_value(at%properties, "ljenergy", sum(atlocale))
        call write(at, movie)
        
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
        
        !find highest energy atoms
        etmp = elj(x)
        idx = (/ (i, i=1,size(etmp)) /)
        call sort_array(etmp, idx)
        rr = 3.0_dp ! randomize by max this amount
        
        do j=1,1
           jj = idx(N-j+1)
           call print('Randomising '//jj)
           r0 = x((jj-1)*3+1:(jj-1)*3+3)
           atlocale(jj) = 1.0_dp
           do while (atlocale(jj) > 0.0_dp)
              r = 0
              call randomise(r, rr)
              x((jj-1)*3+1:(jj-1)*3+3) = r0 + r
              atlocale = elj(x)
           end do
        enddo
        
        at%pos = reshape(x, (/3,at%N/))
        

     end do
  end do

  call system_finalise()

contains

subroutine random_initial_pos(at, rad)
  type(Atoms), intent(inout) :: at
  real(dp), intent(in) :: rad
  logical :: bad
  real(dp), allocatable :: x(:)
  real(dp), pointer :: atlocale(:)
  integer :: i

  allocate(x(at%N*3))

  at%pos = 0.0_dp
  call randomise(at%pos, rad)
  !     at%pos = at%pos+50.0_dp
  
  bad = .true.
  do while(bad)
     x = reshape(at%pos, (/at%N*3/))
     atlocale = elj(x)
     !        call calc(mpot, at, local_energy=e)
     call print (atlocale)
     bad = .false.
     do i=1,at%N
        if(atlocale(i) > 10.0_dp) then 
           call randomise(at%pos(:,i), 2.0_dp)
           bad = .true.
        end if
     end do
  end do
  
end subroutine random_initial_pos

subroutine random_initial_pos_chris(at, rad)
  type(Atoms), intent(inout) :: at
  real(dp), intent(in) :: rad
  logical :: bad
  real(dp) :: cut

  cut = rad/(real(at%N)**0.3333333_dp)*0.6_dp

  at%pos = 0.0_dp
  do i=1,at%N
     bad = .true.
     do while(bad)
        r = 0.0_dp
        call randomise(r, rad)
        if(norm(r) > rad) then
           bad = .true.
           cycle
        end if

        bad = .false.
        do j=1,i-1
           if(norm(at%pos(:,j)-r) < cut) bad = .true.
        end do
     end do
     at%pos(:,i) = r
  end do

end subroutine random_initial_pos_chris

end program local_random_search


