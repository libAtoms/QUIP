program parse_fc
use libatoms_module
implicit none
  character(len=100) :: f_s, f_fc, cmd
  integer :: N_at, N_fc, s(3), i, j, fc_i_a(1), fc_i, Zi, Zj
  integer, allocatable :: Z(:), bond_list_i(:), bond_list(:)
  double precision :: m, p(3), dr(3)
  double precision, allocatable :: phi2(:), phi3(:), phi4(:), r0(:)
  double precision, allocatable :: phi2_table(:,:,:), phi3_table(:,:,:), phi4_table(:,:,:), r0_table(:,:,:)
  integer :: n_types, max_n_fcs

  call system_initialise()

  if (cmd_arg_count() /= 2) then
    call get_cmd_arg(0, cmd)
    call system_abort("Usage: "// trim(cmd)//" struct_file fc_file")
  endif

  call get_cmd_arg(1, f_s)
  call get_cmd_arg(2, f_fc)

  open (unit=10, file=f_s, status="OLD")
  read (10,*) N_at
  ! print *, "reading ",N_at," atoms"
  read (10,*) p
  read (10,*) p
  read (10,*) p
  allocate(Z(N_at))
  do i=1, N_at
    read (10,*) p, Z(i), m
    if (abs(m-55.85) < 0.1) Z(i) = 26
    if (abs(m-121.75) < 0.1) Z(i) = 51
    if (abs(m-138.92) < 0.1) Z(i) = 57
  end do
  close(unit=10)

  open (unit=10, file=f_fc, status="OLD")
  read (10,*) N_fc
  ! print *, "reading ",N_fc," force constants"
  allocate(phi2(N_fc), phi3(N_fc), phi4(N_fc), r0(N_fc))
  do i=1, N_fc
    read (10,*) phi2(i), phi3(i), phi4(i), s, dr, r0(i)
  end do
  allocate(bond_list_i(N_at+1))
  allocate(bond_list(N_fc))
  read (10,*) bond_list_i
  read (10,*) N_fc
  read (10,*) bond_list

  allocate(phi2_table(130,130,30))
  allocate(phi3_table(130,130,30))
  allocate(phi4_table(130,130,30))
  allocate(r0_table(130,130,30))
  r0_table = -1.0
  phi2_table = 0.0
  phi3_table = 0.0
  phi4_table = 0.0

  do i=1, N_at
    ! print '("atom ", I6)', i
    do j=bond_list_i(i), bond_list_i(i+1)-1
      ! print '("bond ",I3," Zs ",2I3," phis ",3F25.15," r0 ",F25.15)', j-bond_list_i(i)+1, &
	! Z(i), Z(bond_list(j)), phi2(j), phi3(j), phi4(j), r0(j)
      Zi = Z(i)
      Zj = Z(bond_list(j))
      if (minval(abs(r0_table(Zi,Zj,:)-r0(j))) > 1.0e-4_dp) then
	fc_i_a = minloc(r0_table(Zi,Zj,:))
	fc_i = fc_i_a(1)
	phi2_table(Zi,Zj,fc_i) = phi2(j)
	phi3_table(Zi,Zj,fc_i) = phi3(j)
	phi4_table(Zi,Zj,fc_i) = phi4(j)
	r0_table(Zi,Zj,fc_i) = r0(j)
      endif
    end do
  end do

  n_types = count(count(r0_table(:,:,1) > 0.0, 2) > 0)
  max_n_fcs = maxval(count(r0_table > 0.0, 3))
  call print('<FC_params n_types="'// n_types //'" ideal_struct_file="" max_n_fcs="'//max_n_fcs//'" cutoff="'//maxval(r0_table)+1.0//'" label="" >')
  do Zi=1, size(phi2_table,1)
    do Zj=1, size(phi2_table,2)
      if (maxval(r0_table(Zi,Zj,:)) > 0.0) then
	do fc_i=1, count(r0_table(Zi,Zj,:) > 0.0)
	  call print('  <FC atnum_i="'//Zi//'" atnum_j="'//Zj//'" fc_i="'//fc_i//'" r0="'//r0_table(Zi,Zj,fc_i)//'" phi2="'//phi2_table(Zi,Zj,fc_i)//&
	    '" phi3="'//phi3_table(Zi,Zj,fc_i)//'" phi4="'//phi4_table(Zi,Zj,fc_i)//'" />')
	end do
      endif
    end do
  end do
  call print ("</FC_params>")
end program
