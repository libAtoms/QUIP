program test_QC_QUIP_Wrapper
use QC_QUIP_Wrapper_module
  integer N
  double precision, allocatable :: pos(:,:), weight(:), local_e(:), f(:,:)
  integer, allocatable :: Z(:)
  double precision Lz

  integer err
  integer i

  N = 3
  allocate(pos(3,N), weight(N), local_E(N), f(3,N), Z(N))
  do i=1, N
    Z(i) = 47
    pos(1,i) = i*2.0D0
    pos(2:3,i) = 0.0D0
    if (i <= 2) then
      weight(i) = 1.0D0
    else
      weight(i) = 0.0D0
    end if
  end do

  call QC_QUIP_Initialise("EP LJ", err)
  if (err /= 0) then
    print *, "Error ", err, "calling Initialise"
    stop
  endif

  Lz = 2.0D0
  call QC_QUIP_calc(Lz, pos, Z, weight, local_e, f, err)
  if (err /= 0) then
    print *, "Error ", err, "calling calc"
    stop
  endif

  print *, "local_e " , local_e
  print *, "total_e ", sum(weight*local_e)
  do i=1, N
    write (*,'("atoms ", 3F10.5," ",I0," ",F5.3,3F10.5)') pos(:,i),Z(i),weight(i),f(:,i)
  end do

end program
