module test_SE_mod

use System_module
use Atoms_module
use QUIP_Module
implicit none

contains

subroutine print_G_dev(at, tbsys, G, Gref, prefix)
  type(Atoms), intent(in) :: at
  type(TBsystem), intent(in) :: tbsys
  complex(dp), intent(in) :: G(:,:), Gref(:,:)
  character(len=*) prefix

  integer i, orb_min, orb_max
  real(dp) :: lp(3)

  do i=1, at%N
    orb_min = tbsys%first_orb_of_atom(i)
    orb_max = orb_min + n_orbs_of_Z(tbsys%tbmodel, at%Z(i)) - 1
    lp = at%g .mult. at%pos(:,i)
    call print (prefix // " G dev " // (lp(3)*norm(at%lattice(:,3))) // " " // maxval(abs(G(orb_min:orb_max,:)-Gref(orb_min:orb_max,:))))
  end do
end subroutine

subroutine print_local_e_dev (ref, at, ref_e, local_e, prefix)
  type(Atoms), intent(in) :: ref, at
  real(dp), intent(in) :: ref_e(:), local_e(:)
  character(len=*) :: prefix

  integer i, ii, min_i
  real(dp) :: min_dist

  real(dp) :: lp (3)

  do i=1, at%N
    min_dist = 1e38_dp
    min_i = -1
    do ii=1, ref%N
      if (norm(ref%pos(:,ii)-at%pos(:,i)) < min_dist) then
	min_i = ii
	min_dist = norm(ref%pos(:,ii)-at%pos(:,i))
      endif
    end do
    if (min_i < 1 .or. min_i > ref%N) call system_abort("garbage index matching to bulk " // min_i)
    if (min_dist > 1e-1_dp) then
      call print ("bulk configuration")
      call print(ref)
      call system_abort("Couldn't find match to cluster atom " // i // " at " // at%pos(:,i))
    endif
    lp = at%g .mult.at%pos(:,i)
    call print (prefix // " local_e dev " // (lp(3)*norm(at%lattice(:,3))) // " " // (local_e(i) - ref_e(min_i)))
  end do

end subroutine

subroutine print_struct(at, tbc, m, prefix)
implicit none
  type(Atoms), intent(in) :: at
  type(TBCalculate), intent(in) :: tbc
  type(TBMatrix), intent(in) :: m
  character(len=*) :: prefix

  integer i, j, i_orb_min, i_orb_max, j_orb_min, j_orb_max

  do i=1, at%N
    i_orb_min = tbc%tbsys%first_orb_of_atom(i)
    i_orb_max = i_orb_min + n_orbs_of_Z(tbc%tbsys%tbmodel, at%Z(i)) - 1
    if (m%is_complex) then
      if (maxval(abs(m%data_z(1)%data(i_orb_min:i_orb_max,:))) == 0.0_dp) then
	  call Print (prefix // " SE struct " // i // " (" // atoms_n_neighbours(at,i) // ") all 0's")
      else
	do j=1, at%N
	j_orb_min = tbc%tbsys%first_orb_of_atom(j)
	j_orb_max = j_orb_min + n_orbs_of_Z(tbc%tbsys%tbmodel, at%Z(j)) - 1
	  if (maxval(abs(m%data_z(1)%data(i_orb_min:i_orb_max,j_orb_min:j_orb_max))) /= 0.0_dp) then
	    call print (prefix // " SE struct " // i // " (" // atoms_n_neighbours(at,i) // ") " // j // " " // &
			atoms_n_neighbours(at,j) // " " // &
			maxval(abs(m%data_z(1)%data(i_orb_min:i_orb_max,j_orb_min:j_orb_max))))
	  endif
	end do
      endif
    else
      if (maxval(abs(m%data_d(1)%data(i_orb_min:i_orb_max,:))) == 0.0_dp) then
	  call Print (prefix // " SE struct " // i // " (" // atoms_n_neighbours(at,i) // ") all 0's")
      else
	do j=1, at%N
	j_orb_min = tbc%tbsys%first_orb_of_atom(j)
	j_orb_max = j_orb_min + n_orbs_of_Z(tbc%tbsys%tbmodel, at%Z(j)) - 1
	  if (maxval(abs(m%data_d(1)%data(i_orb_min:i_orb_max,j_orb_min:j_orb_max))) /= 0.0_dp) then
	    call print (prefix // " SE struct " // i // " (" // atoms_n_neighbours(at,i) // ") " // j // " " // &
			atoms_n_neighbours(at,j) // " " // &
			maxval(abs(m%data_d(1)%data(i_orb_min:i_orb_max,j_orb_min:j_orb_max))))
	  endif
	end do
      endif
    endif
  end do
end subroutine print_struct

subroutine make_dia_111_si(at, n1, n2, n3)
  type(Atoms), intent(inout) :: at

  integer n1, n2, n3
  integer i1, i2, i3, ii
  real(dp) :: a1(3), a2(3), a3(3), o(3)
  real(dp) :: lc = 5.43_dp
  real(dp) :: basis(3,6)
  real(dp) :: lat(3,3)

  a1 = (/0.5_dp, -0.5_dp, 0.0_dp/)
  a2 = (/0.0_dp, 0.5_dp, -0.5_dp/)
  a3 = (/1.0_dp, 1.0_dp, 1.0_dp/)
  basis(:,1) = (/0.0_dp, 0.0_dp, 0.0_dp/)
  basis(:,2) = (/0.5_dp, 0.5_dp, 0.0_dp/)
  basis(:,3) = (/1.0_dp, 0.5_dp, 0.5_dp/)
  basis(:,4) = (/0.25_dp, 0.25_dp, 0.25_dp/)
  basis(:,5) = (/0.75_dp, 0.75_dp, 0.25_dp/)
  basis(:,6) = (/1.25_dp, 0.75_dp, 0.75_dp/)

  lat(:,1) = n1*a1*lc
  lat(:,2) = n2*a2*lc
  lat(:,3) = n3*a3*lc
  call atoms_initialise(at,0,lat)

  a1 = a1 * lc
  a2 = a2 * lc
  a3 = a3 * lc
  basis = basis*lc

  do i1 = 1, n1
  do i2 = 1, n2
  do i3 = 1, n3
    o = a1*i1 + a2*i2 + a3*i3
    do ii=1, 6
      call add_atoms(at, o+basis(:,ii), 14)
    end do
  end do
  end do
  end do

end subroutine make_dia_111_si

subroutine make_fcc_111_ag(at, n1, n2, n3)
  type(Atoms), intent(inout) :: at

  integer n1, n2, n3
  integer i1, i2, i3, ii
  real(dp) :: a1(3), a2(3), a3(3), o(3)
  real(dp) :: lc = 4.09_dp
  real(dp) :: basis(3,3)
  real(dp) :: lat(3,3)

  a1 = (/0.5_dp, -0.5_dp, 0.0_dp/)
  a2 = (/0.0_dp, 0.5_dp, -0.5_dp/)
  a3 = (/1.0_dp, 1.0_dp, 1.0_dp/)
  basis(:,1) = (/0.0_dp, 0.0_dp, 0.0_dp/)
  basis(:,2) = (/0.5_dp, 0.5_dp, 0.0_dp/)
  basis(:,3) = (/1.0_dp, 0.5_dp, 0.5_dp/)

  lat(:,1) = n1*a1*lc
  lat(:,2) = n2*a2*lc
  lat(:,3) = n3*a3*lc
  call atoms_initialise(at,0,lat)

  a1 = a1 * lc
  a2 = a2 * lc
  a3 = a3 * lc
  basis = basis*lc

  do i1 = 1, n1
  do i2 = 1, n2
  do i3 = 1, n3
    o = a1*i1 + a2*i2 + a3*i3
    do ii=1, 3
      call add_atoms(at, o+basis(:,ii), 47)
    end do
  end do
  end do
  end do

end subroutine make_fcc_111_ag

subroutine make_sys_si(bulk, at_i, at_o, rv)
  type(Atoms), intent(inout) :: bulk, at_i, at_o
  real(dp), intent(in), optional :: rv

  real(dp) :: my_rv
  integer i

  logical, allocatable :: mask(:)
  real(dp) :: lc(3)

  my_rv = 0.0_dp
  if (present(rv)) my_rv = rv

  call make_dia_111_si(bulk, 3, 3, 8)
  do i=1, bulk%N
    bulk%pos(:,i) = bulk%pos(:,i) - 0.5_dp*bulk%lattice(:,3)/norm(bulk%lattice(:,3))
  end do
  call randomise(bulk%pos, my_rv)
  call map_into_cell(bulk)

  allocate(mask(bulk%N))

  mask = .false.
  do i=1, bulk%N
    lc = bulk%g .mult. bulk%pos(:,i)
    if (lc(3) < my_rv*2.0_dp/norm(bulk%lattice(:,3))) then
      mask(i) = .true.
    endif
  end do
  call select(at_i, bulk, mask)

  mask = .false.
  do i=1, bulk%N
    lc = bulk%g .mult. bulk%pos(:,i)
    if (lc(3) >= my_rv*2.0_dp/norm(bulk%lattice(:,3))) then
      mask(i) = .true.
    endif
  end do
  call select(at_o, bulk, mask)

  call Finalise(bulk)
  bulk = at_i
  call add_atoms(bulk, data=at_o%data)

end subroutine

subroutine make_sys_ag(bulk, at_i, at_o, rv)
  type(Atoms), intent(inout) :: bulk, at_i, at_o
  real(dp), intent(in), optional :: rv

  real(dp) :: my_rv
  integer i

  logical, allocatable :: mask(:)
  real(dp) :: lc(3)

  my_rv = 0.0_dp
  if (present(rv)) my_rv = rv

  call make_fcc_111_ag(bulk, 3, 3, 8)
  do i=1, bulk%N
    bulk%pos(:,i) = bulk%pos(:,i) - 0.5_dp*bulk%lattice(:,3)/norm(bulk%lattice(:,3))
  end do
  call randomise(bulk%pos, my_rv)
  call map_into_cell(bulk)

  allocate(mask(bulk%N))

  mask = .false.
  do i=1, bulk%N
    lc = bulk%g .mult. bulk%pos(:,i)
    if (lc(3) < my_rv*2.0_dp/norm(bulk%lattice(:,3))) then
      mask(i) = .true.
    endif
  end do
  call select(at_i, bulk, mask)

  mask = .false.
  do i=1, bulk%N
    lc = bulk%g .mult. bulk%pos(:,i)
    if (lc(3) >= my_rv*2.0_dp/norm(bulk%lattice(:,3))) then
      mask(i) = .true.
    endif
  end do
  call select(at_o, bulk, mask)

  call Finalise(bulk)
  bulk = at_i
  call add_atoms(bulk, data=at_o%data)

end subroutine

end module test_SE_mod

program test_SE
use System_module
use Extendable_Str_module
use Atoms_module
use TBCalculate_module
use ApproxFermi_module
use test_SE_mod
implicit none

  type(inoutput) :: params, out
  type(extendable_str) :: params_str
  type(TBCalculate) :: tbc_bulk, tbc_bulk_id, tbc_i, tbc_SE_i, tbc_o, tbc_i_passivated
  type(Atoms) :: bulk, bulk_id, at_i, at_o, at_i_passivated
  type(ApproxFermi) :: af
  real(dp) :: E_bulk, E_i, E_o
  real(dp) :: GF_E_bulk, GF_E_i, GF_E_o
  type(table) :: cluster_list
  real(dp) :: band_width, exact_band_width

  real(dp), allocatable :: local_e(:), local_e_i(:), local_E_i_passivated(:), local_e_i_SE(:)

  complex(dp), allocatable :: Aio(:,:), Aoi(:,:), Aoo(:,:), Goo(:,:), tt(:,:), Gii_exact(:,:), Aoo_diag_inv_vec(:), Aoo_offdiag(:,:)
  type(TBMatrix), allocatable :: SEii(:)
  integer i_min, i_max, o_min, o_max
  integer i_orb_min, i_orb_max, j_orb_min, j_orb_max

  integer ip, i, j

  call system_initialise()
  use_intrinsic_blas = .false.
  call Initialise(out, 'stdout', OUTPUT)

  call Initialise(params, "quip_params.xml", INPUT)
  call Initialise(params_str)
  call read(params_str, params%unit)
#ifdef SI
  call Initialise(tbc_bulk, "DFTB", string(params_str))
  call Initialise(tbc_bulk_id, "DFTB", string(params_str))
  call Initialise(tbc_i, "DFTB", string(params_str))
  call Initialise(tbc_i_passivated, "DFTB", string(params_str))
  call Initialise(tbc_SE_i, "DFTB", string(params_str))
  call Initialise(tbc_o, "DFTB", string(params_str))
#else 
#ifdef AG
  call Initialise(tbc_bulk, "NRL-TB Silver", string(params_str))
  call Initialise(tbc_bulk_id, "NRL-TB Silver", string(params_str))
  call Initialise(tbc_i, "NRL-TB Silver", string(params_str))
  call Initialise(tbc_i_passivated, "NRL-TB Silver", string(params_str))
  call Initialise(tbc_SE_i, "NRL-TB Silver", string(params_str))
  call Initialise(tbc_o, "NRL-TB Silver", string(params_str))
#else
DIE
#endif
#endif

#ifdef SI
  call make_sys_si(bulk_id, at_i, at_o)
#else
#ifdef AG
  call make_sys_ag(bulk_id, at_i, at_o)
#else
DIE
#endif
#endif
  call calc_connect(bulk_id)

#ifdef SI
  call make_sys_si(bulk, at_i, at_o, 0.01_dp)
#else
#ifdef AG
  call make_sys_ag(bulk, at_i, at_o, 0.01_dp)
#else
DIE
#endif
#endif
  call calc_connect(bulk)
  call calc_connect(at_i)
  call calc_connect(at_o)

  out%prefix = "+AT_I+"
  call print_xyz(at_i, out)
  out%prefix = ""
  out%prefix = "+AT_O+"
  call print_xyz(at_o, out)
  out%prefix = ""
  out%prefix = "+BULK+"
  call print_xyz(bulk, out)
  out%prefix = ""

  call table_allocate(cluster_list, 4, 0, 0, 0)
  do i=1, at_i%N
    call append(cluster_list, (/i, 0, 0, 0/))
  end do
  at_i_passivated = create_cluster(bulk, cluster_list, .true., same_lattice = .true., in_out_in = .false.)
  out%prefix = "+AT_I_PASSIVATED+"
  call print_xyz(at_i_passivated, out)
  out%prefix = ""

#ifdef SI
  band_width = 18.0_dp
#else
#ifdef AG
  band_width = 10.0_dp
#else
DIE
#endif
#endif

  call Setup_atoms(tbc_bulk_id, bulk_id)
  E_bulk = calc_diag(tbc_bulk_id, bulk_id, Fermi_T = 0.1_dp)
  exact_band_width = tbc_bulk_id%fermi_E - minval(tbc_bulk_id%evals%data_d)
  call print ("True band width " // exact_band_width)
  if (exact_band_width > band_width) call system_abort ("band width too small")

  call Initialise(af, 0.0_dp, 0.1_dp, band_width)
  E_bulk = calc_diag(tbc_bulk_id, bulk_id, Fermi_T = 0.1_dp, AF=af)

  call Setup_atoms(tbc_bulk, bulk)
  E_bulk = calc_diag(tbc_bulk, bulk, Fermi_T = 0.1_dp, AF=af)

  call print("got AF from bulk")
  call print(af)

  call Setup_atoms(tbc_i, at_i)
  E_i = calc_diag(tbc_i, at_i, use_Fermi_E = .true., AF=af)
  call Setup_atoms(tbc_o, at_o)
  E_o = calc_diag(tbc_o, at_o, use_Fermi_E = .true., AF=af)

  call Print ("E " // E_bulk // " " // E_i // " " // E_o)

  allocate(local_e(bulk%N))
  allocate(local_e_i(at_i%N))
  allocate(local_e_i_passivated(at_i_passivated%N))
  allocate(local_e_i_SE(at_i%N))

  GF_E_bulk= calc_GF(tbc_bulk, bulk, .true., af%Fermi_E, 0.1_dp, band_width, local_e = local_e)
  GF_E_i = calc_GF(tbc_i, at_i, .true., af%Fermi_E, 0.1_dp, band_width, local_e = local_e_i)
  GF_E_o = calc_GF(tbc_o, at_o, .true., af%Fermi_E, 0.1_dp, band_width)

  call print ("bulk local_e")
  call print (local_e)

  call print("tbc_bulk%fermi_E " // tbc_bulk%fermi_E // " " // af%fermi_E)

  call Print ("GF_E " // GF_E_bulk // " " // GF_E_i // " " // GF_E_o)
  allocate(Aio(tbc_i%tbsys%N,tbc_o%tbsys%N))
  allocate(tt(tbc_i%tbsys%N,tbc_o%tbsys%N))
  allocate(Aoi(tbc_o%tbsys%N,tbc_i%tbsys%N))
  allocate(Aoo(tbc_o%tbsys%N,tbc_o%tbsys%N))
  allocate(Goo(tbc_o%tbsys%N,tbc_o%tbsys%N))
  allocate(Gii_exact(tbc_i%tbsys%N,tbc_i%tbsys%N))

  i_min = 1
  i_max = tbc_i%tbsys%N
  o_min = tbc_i%tbsys%N+1
  o_max = tbc_bulk%tbsys%N

  call print ("i min max " // i_min // " " // i_max // " o min max " // o_min // " " // o_max)

  allocate(SEii(af%n_poles))

  do ip=1, af%n_poles
    call Initialise(SEii(ip), tbc_i%tbsys%N, tbc_i%tbsys%n_matrices, .true.)
  end do

  Gii_exact = tbc_bulk%GF%G(1)%data_z(1)%data(i_min:i_max,i_min:i_max)

 call print ("after exact calc")
 call system("ps auxw | egrep 'MEM|test_SE' | grep -v grep")
! not terminated
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call print_G_dev(at_i, tbc_i%tbsys, tbc_i%GF%G(1)%data_z(1)%data, Gii_exact, "not terminated")
  call print_local_e_dev(bulk, at_i, local_e, local_e_i, "not terminated")

  call Setup_atoms(tbc_SE_i, at_i)

call Print("")
! good approx Aoo^-1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  allocate(Aoo_diag_inv_vec(tbc_o%tbsys%N))
  allocate(Aoo_offdiag(tbc_o%tbsys%N, tbc_o%tbsys%N))

  do ip=1, af%n_poles
    Aio = af%z(ip)*tbc_bulk%tbsys%S%data_d(1)%data(i_min:i_max,o_min:o_max) - tbc_bulk%tbsys%H%data_d(1)%data(i_min:i_max,o_min:o_max)
    Aoi = af%z(ip)*tbc_bulk%tbsys%S%data_d(1)%data(o_min:o_max,i_min:i_max) - tbc_bulk%tbsys%H%data_d(1)%data(o_min:o_max,i_min:i_max)
    Aoo = af%z(ip)*tbc_bulk%tbsys%S%data_d(1)%data(o_min:o_max,o_min:o_max) - tbc_bulk%tbsys%H%data_d(1)%data(o_min:o_max,o_min:o_max)

    !! call inverse(Aoo, Goo, positive_in = .false.)
    ! Goo = Aoo_diag^-1 - Aoo_diag^-1 Aoo_offdiag Aoo_diag^-1

    Aoo_diag_inv_vec = diag(Aoo)
    Aoo_diag_inv_vec = 1.0_dp/Aoo_diag_inv_vec

    Aoo_offdiag = Aoo
    do i=1, size(Aoo,1)
      Aoo_offdiag(i,i) = 0.0_dp
    end do
    call matrix_product_vect_asdiagonal_sub(tt,Aoo_offdiag,Aoo_diag_inv_vec)
    call matrix_product_vect_asdiagonal_sub(Goo,Aoo_diag_inv_vec,tt)
    Goo = -Goo
    do i=1, size(Goo,1)
      Goo(i,i) = Goo(i,i) + Aoo_diag_inv_vec(i)
    end do

 call print ("starting matmul approx2 pole " // ip)
    call matrix_product_sub(tt, Aio, Goo)
    call matrix_product_sub(SEii(ip)%data_z(1)%data, tt, Aoi)
 call print ("done matmul approx2 pole " // ip)
    do i=1, size(SEii(ip)%data_z(1)%data,1)
    do j=1, size(SEii(ip)%data_z(1)%data,2)
      SEii(ip)%data_z(1)%data(i,j) = -SEii(ip)%data_z(1)%data(i,j)
    end do
    end do
    ! call Print ("got SE " // ip)
 call print ("done approx2 pole " // ip)
  end do
  deallocate(Aoo_diag_inv_vec)
  deallocate(Aoo_offdiag)

  GF_E_i = calc_GF(tbc_SE_i, at_i, .true., af%Fermi_E, 0.1_dp, band_width, SelfEnergy=SEii, local_e = local_e_i_SE)
  call print_G_dev(at_i, tbc_SE_i%tbsys, tbc_SE_i%GF%G(1)%data_z(1)%data, Gii_exact, "approx2 Goo")
  call print_local_e_dev(bulk, at_i, local_e, local_e_i_SE, "approx2 Goo")


call Print("")
! exact SE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do ip=1, af%n_poles
    Aio = af%z(ip)*tbc_bulk%tbsys%S%data_d(1)%data(i_min:i_max,o_min:o_max) - tbc_bulk%tbsys%H%data_d(1)%data(i_min:i_max,o_min:o_max)
    Aoi = af%z(ip)*tbc_bulk%tbsys%S%data_d(1)%data(o_min:o_max,i_min:i_max) - tbc_bulk%tbsys%H%data_d(1)%data(o_min:o_max,i_min:i_max)
    Aoo = af%z(ip)*tbc_bulk%tbsys%S%data_d(1)%data(o_min:o_max,o_min:o_max) - tbc_bulk%tbsys%H%data_d(1)%data(o_min:o_max,o_min:o_max)

    call inverse(Aoo, Goo, positive_in = .false.)

    call matrix_product_sub(tt, Aio, Goo)
    call matrix_product_sub(SEii(ip)%data_z(1)%data, tt, Aoi)
    SEii(ip)%data_z(1)%data = - SEii(ip)%data_z(1)%data
    ! call Print ("got SE " // ip)
  end do
 call print("onsite SE(1)")
 call print(SEii(1)%data_z(1)%data(1:4,1:4))
  GF_E_i = calc_GF(tbc_SE_i, at_i, .true., af%Fermi_E, 0.1_dp, band_width, SelfEnergy=SEii, local_e = local_e_i_SE)
  call print_G_dev(at_i, tbc_SE_i%tbsys, tbc_SE_i%GF%G(1)%data_z(1)%data, Gii_exact, "exact SE")
  call print_local_e_dev(bulk, at_i, local_e, local_e_i_SE, "exact SE")
!!  call print_struct(at_i, tbc_i, tbc_i%tbsys%H, "H")
!!  call print_struct(at_i, tbc_SE_i, SEii(1), "SE")

call Print("")
! truncate SE, spherically symmetric part
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do ip=1, af%n_poles
    Aio = af%z(ip)*tbc_bulk_id%tbsys%S%data_d(1)%data(i_min:i_max,o_min:o_max) - tbc_bulk_id%tbsys%H%data_d(1)%data(i_min:i_max,o_min:o_max)
    Aoi = af%z(ip)*tbc_bulk_id%tbsys%S%data_d(1)%data(o_min:o_max,i_min:i_max) - tbc_bulk_id%tbsys%H%data_d(1)%data(o_min:o_max,i_min:i_max)
    Aoo = af%z(ip)*tbc_bulk_id%tbsys%S%data_d(1)%data(o_min:o_max,o_min:o_max) - tbc_bulk_id%tbsys%H%data_d(1)%data(o_min:o_max,o_min:o_max)

    call inverse(Aoo, Goo, positive_in = .false.)

    call matrix_product_sub(tt, Aio,Goo)
    call matrix_product_sub(SEii(ip)%data_z(1)%data, tt, Aoi)
    SEii(ip)%data_z(1)%data = -SEii(ip)%data_z(1)%data

    do i=1, at_i%N
      if (atoms_n_neighbours(at_i,i) == 4) then
	i_orb_min = tbc_i%tbsys%first_orb_of_atom(i)
	i_orb_max = i_orb_min + n_orbs_of_Z(tbc_i%tbsys%tbmodel, at_i%Z(i)) - 1
!! if (ip == 1) call print ("maxval of zeroed SE " // maxval(abs(SEii(ip)%data_z(1)%data(:,i_orb_min:i_orb_max))))
	SEii(ip)%data_z(1)%data(:,i_orb_min:i_orb_max) = 0.0_dp
	SEii(ip)%data_z(1)%data(i_orb_min:i_orb_max,:) = 0.0_dp
      else
	i_orb_min = tbc_i%tbsys%first_orb_of_atom(i)
	i_orb_max = i_orb_min + n_orbs_of_Z(tbc_i%tbsys%tbmodel, at_i%Z(i)) - 1
	SEii(ip)%data_z(1)%data(:,i_orb_min+1:i_orb_max) = 0.0_dp
	SEii(ip)%data_z(1)%data(i_orb_min+1:i_orb_max,:) = 0.0_dp
      endif
    end do
    ! call Print ("got SE " // ip)
    SEii(ip)%data_z(1)%data = SEii(ip)%data_z(1)%data
  end do
  GF_E_i = calc_GF(tbc_SE_i, at_i, .true., af%Fermi_E, 0.1_dp, band_width, SelfEnergy=SEii, local_e = local_e_i_SE)
  call print_G_dev(at_i, tbc_SE_i%tbsys, tbc_SE_i%GF%G(1)%data_z(1)%data, Gii_exact, "trunc s_SE")
  call print_local_e_dev(bulk, at_i, local_e, local_e_i_SE, "trunc s_SE")

 call Print("")
#ifdef SI
 ! H terminated
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 
   call Setup_atoms(tbc_i_passivated, at_i_passivated)
   GF_E_bulk= calc_GF(tbc_i_passivated, at_i_passivated, .true., af%Fermi_E, 0.1_dp, band_width, local_e = local_e_i_passivated)
   call print_G_dev(at_i, tbc_i%tbsys, tbc_i_passivated%GF%G(1)%data_z(1)%data(i_min:i_max,i_min:i_max), Gii_exact, "H terminated")
   call print_local_e_dev(bulk, at_i, local_e, local_e_i_passivated, "H terminated")
 
 call Print("")
#endif
! truncate Goo, nearest neighbor only
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do ip=1, af%n_poles
    Aio = af%z(ip)*tbc_bulk_id%tbsys%S%data_d(1)%data(i_min:i_max,o_min:o_max) - tbc_bulk_id%tbsys%H%data_d(1)%data(i_min:i_max,o_min:o_max)
    Aoi = af%z(ip)*tbc_bulk_id%tbsys%S%data_d(1)%data(o_min:o_max,i_min:i_max) - tbc_bulk_id%tbsys%H%data_d(1)%data(o_min:o_max,i_min:i_max)
    Aoo = af%z(ip)*tbc_bulk_id%tbsys%S%data_d(1)%data(o_min:o_max,o_min:o_max) - tbc_bulk_id%tbsys%H%data_d(1)%data(o_min:o_max,o_min:o_max)

    call inverse(Aoo, Goo, positive_in = .false.)
    do i=1, at_o%N
      if (atoms_n_neighbours(at_o,i) == 4) then
	i_orb_min = tbc_o%tbsys%first_orb_of_atom(i)
	i_orb_max = i_orb_min + n_orbs_of_Z(tbc_o%tbsys%tbmodel, at_o%Z(i)) - 1
!! if (ip == 1) call print ("maxval of zeroed Goo " // maxval(abs(Goo(:,i_orb_min:i_orb_max))))
	Goo(:,i_orb_min:i_orb_max) = 0.0_dp
	Goo(i_orb_min:i_orb_max,:) = 0.0_dp
      endif
    end do
    call matrix_product_sub(tt, Aio,Goo)
    call matrix_product_sub(SEii(ip)%data_z(1)%data, tt, Aoi)
    SEii(ip)%data_z(1)%data = -SEii(ip)%data_z(1)%data
    ! call Print ("got SE " // ip)
  end do
  GF_E_i = calc_GF(tbc_SE_i, at_i, .true., af%Fermi_E, 0.1_dp, band_width, SelfEnergy=SEii, local_e = local_e_i_SE)
  call print_G_dev(at_i, tbc_SE_i%tbsys, tbc_SE_i%GF%G(1)%data_z(1)%data, Gii_exact, "nn Goo")
  call print_local_e_dev(bulk, at_i, local_e, local_e_i_SE, "nn Goo")

call Print("")
! truncate Goo, nearest neighbor onsite only
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do ip=1, af%n_poles
    Aio = af%z(ip)*tbc_bulk_id%tbsys%S%data_d(1)%data(i_min:i_max,o_min:o_max) - tbc_bulk_id%tbsys%H%data_d(1)%data(i_min:i_max,o_min:o_max)
    Aoi = af%z(ip)*tbc_bulk_id%tbsys%S%data_d(1)%data(o_min:o_max,i_min:i_max) - tbc_bulk_id%tbsys%H%data_d(1)%data(o_min:o_max,i_min:i_max)
    Aoo = af%z(ip)*tbc_bulk_id%tbsys%S%data_d(1)%data(o_min:o_max,o_min:o_max) - tbc_bulk_id%tbsys%H%data_d(1)%data(o_min:o_max,o_min:o_max)

    call inverse(Aoo, Goo, positive_in = .false.)
    do i=1, at_o%N
      if (atoms_n_neighbours(at_o,i) == 4) then
	i_orb_min = tbc_o%tbsys%first_orb_of_atom(i)
	i_orb_max = i_orb_min + n_orbs_of_Z(tbc_o%tbsys%tbmodel, at_o%Z(i)) - 1
	Goo(:,i_orb_min:i_orb_max) = 0.0_dp
	Goo(i_orb_min:i_orb_max,:) = 0.0_dp
      else
	i_orb_min = tbc_o%tbsys%first_orb_of_atom(i)
	i_orb_max = i_orb_min + n_orbs_of_Z(tbc_o%tbsys%tbmodel, at_o%Z(i)) - 1
	if (i_orb_min > 1) then
	  Goo(:,1:i_orb_min-1) = 0.0_dp
	  Goo(1:i_orb_min-1,:) = 0.0_dp
	endif
	if (i_orb_max < tbc_o%tbsys%N) then
	  Goo(:,i_orb_max+1:tbc_o%tbsys%N) = 0.0_dp
	  Goo(i_orb_max+1:tbc_o%tbsys%N,:) = 0.0_dp
	endif
      endif
    end do
    call matrix_product_sub(tt, Aio,Goo)
    call matrix_product_sub(SEii(ip)%data_z(1)%data, tt, Aoi)
    SEii(ip)%data_z(1)%data = -SEii(ip)%data_z(1)%data
    ! call Print ("got SE " // ip)
  end do
  GF_E_i = calc_GF(tbc_SE_i, at_i, .true., af%Fermi_E, 0.1_dp, band_width, SelfEnergy=SEii, local_e = local_e_i_SE)
  call print_G_dev(at_i, tbc_SE_i%tbsys, tbc_SE_i%GF%G(1)%data_z(1)%data, Gii_exact, "onsite Goo")
  call print_local_e_dev(bulk, at_i, local_e, local_e_i_SE, "onsite Goo")

call Print("")
! bad approx Aoo^-1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do ip=1, af%n_poles
    Aio = af%z(ip)*tbc_bulk%tbsys%S%data_d(1)%data(i_min:i_max,o_min:o_max) - tbc_bulk%tbsys%H%data_d(1)%data(i_min:i_max,o_min:o_max)
    Aoi = af%z(ip)*tbc_bulk%tbsys%S%data_d(1)%data(o_min:o_max,i_min:i_max) - tbc_bulk%tbsys%H%data_d(1)%data(o_min:o_max,i_min:i_max)
    Aoo = af%z(ip)*tbc_bulk%tbsys%S%data_d(1)%data(o_min:o_max,o_min:o_max) - tbc_bulk%tbsys%H%data_d(1)%data(o_min:o_max,o_min:o_max)

    !! call inverse(Aoo, Goo, positive_in = .false.)
    ! Goo = 2-Aoo
    Goo = -Aoo
    call add_identity(Goo)
    call add_identity(Goo)

    call matrix_product_sub(tt, Aio,Goo)
    call matrix_product_sub(SEii(ip)%data_z(1)%data, tt, Aoi)
    SEii(ip)%data_z(1)%data = -SEii(ip)%data_z(1)%data
    ! call Print ("got SE " // ip)
  end do
  GF_E_i = calc_GF(tbc_SE_i, at_i, .true., af%Fermi_E, 0.1_dp, band_width, SelfEnergy=SEii, local_e = local_e_i_SE)
  call print_G_dev(at_i, tbc_SE_i%tbsys, tbc_SE_i%GF%G(1)%data_z(1)%data, Gii_exact, "approx1 Goo")
  call print_local_e_dev(bulk, at_i, local_e, local_e_i_SE, "approx1 Goo")


call Print("")
! trunc Goo, nearest neighbor only, s part only
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do ip=1, af%n_poles
    Aio = af%z(ip)*tbc_bulk_id%tbsys%S%data_d(1)%data(i_min:i_max,o_min:o_max) - tbc_bulk_id%tbsys%H%data_d(1)%data(i_min:i_max,o_min:o_max)
    Aoi = af%z(ip)*tbc_bulk_id%tbsys%S%data_d(1)%data(o_min:o_max,i_min:i_max) - tbc_bulk_id%tbsys%H%data_d(1)%data(o_min:o_max,i_min:i_max)
    Aoo = af%z(ip)*tbc_bulk_id%tbsys%S%data_d(1)%data(o_min:o_max,o_min:o_max) - tbc_bulk_id%tbsys%H%data_d(1)%data(o_min:o_max,o_min:o_max)

    call inverse(Aoo, Goo, positive_in = .false.)
    do i=1, at_o%N
      if (atoms_n_neighbours(at_o,i) == 4) then
	i_orb_min = tbc_o%tbsys%first_orb_of_atom(i)
	i_orb_max = i_orb_min + n_orbs_of_Z(tbc_o%tbsys%tbmodel, at_o%Z(i)) - 1
!! if (ip == 1) call print ("maxval of zeroed Goo " // maxval(abs(Goo(:,i_orb_min:i_orb_max))))
	Goo(:,i_orb_min:i_orb_max) = 0.0_dp
	Goo(i_orb_min:i_orb_max,:) = 0.0_dp
      else
	i_orb_min = tbc_o%tbsys%first_orb_of_atom(i)
	i_orb_max = i_orb_min + n_orbs_of_Z(tbc_o%tbsys%tbmodel, at_o%Z(i)) - 1
	Goo(:,i_orb_min+1:i_orb_max) = 0.0_dp
	Goo(i_orb_min+1:i_orb_max,:) = 0.0_dp
      endif
    end do
    call matrix_product_sub(tt, Aio,Goo)
    call matrix_product_sub(SEii(ip)%data_z(1)%data, tt, Aoi)
    SEii(ip)%data_z(1)%data = -SEii(ip)%data_z(1)%data
    ! call Print ("got SE " // ip)
  end do
  GF_E_i = calc_GF(tbc_SE_i, at_i, .true., af%Fermi_E, 0.1_dp, band_width, SelfEnergy=SEii, local_e = local_e_i_SE)
  call print_G_dev(at_i, tbc_SE_i%tbsys, tbc_SE_i%GF%G(1)%data_z(1)%data, Gii_exact, "nn s_Goo")
  call print_local_e_dev(bulk, at_i, local_e, local_e_i_SE, "trunc s_Goo")

call Print("")
! trunc SE, nn only
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do ip=1, af%n_poles
    Aio = af%z(ip)*tbc_bulk_id%tbsys%S%data_d(1)%data(i_min:i_max,o_min:o_max) - tbc_bulk_id%tbsys%H%data_d(1)%data(i_min:i_max,o_min:o_max)
    Aoi = af%z(ip)*tbc_bulk_id%tbsys%S%data_d(1)%data(o_min:o_max,i_min:i_max) - tbc_bulk_id%tbsys%H%data_d(1)%data(o_min:o_max,i_min:i_max)
    Aoo = af%z(ip)*tbc_bulk_id%tbsys%S%data_d(1)%data(o_min:o_max,o_min:o_max) - tbc_bulk_id%tbsys%H%data_d(1)%data(o_min:o_max,o_min:o_max)

    call inverse(Aoo, Goo, positive_in = .false.)

    call matrix_product_sub(tt, Aio,Goo)
    call matrix_product_sub(SEii(ip)%data_z(1)%data, tt, Aoi)
    SEii(ip)%data_z(1)%data = -SEii(ip)%data_z(1)%data

    do i=1, at_i%N
      if (atoms_n_neighbours(at_i,i) == 4) then
	i_orb_min = tbc_i%tbsys%first_orb_of_atom(i)
	i_orb_max = i_orb_min + n_orbs_of_Z(tbc_i%tbsys%tbmodel, at_i%Z(i)) - 1
!! if (ip == 1) call print ("maxval of zeroed SE " // maxval(abs(SEii(ip)%data_z(1)%data(:,i_orb_min:i_orb_max))))
	SEii(ip)%data_z(1)%data(:,i_orb_min:i_orb_max) = 0.0_dp
	SEii(ip)%data_z(1)%data(i_orb_min:i_orb_max,:) = 0.0_dp
      endif
    end do
    ! call Print ("got SE " // ip)
  end do
  GF_E_i = calc_GF(tbc_SE_i, at_i, .true., af%Fermi_E, 0.1_dp, band_width, SelfEnergy=SEii, local_e = local_e_i_SE)
  call print_G_dev(at_i, tbc_SE_i%tbsys, tbc_SE_i%GF%G(1)%data_z(1)%data, Gii_exact, "nn SE")
  call print_local_e_dev(bulk, at_i, local_e, local_e_i_SE, "nn SE")

call Print("")
! trunc SE, onsite nn only
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do ip=1, af%n_poles
    Aio = af%z(ip)*tbc_bulk_id%tbsys%S%data_d(1)%data(i_min:i_max,o_min:o_max) - tbc_bulk_id%tbsys%H%data_d(1)%data(i_min:i_max,o_min:o_max)
    Aoi = af%z(ip)*tbc_bulk_id%tbsys%S%data_d(1)%data(o_min:o_max,i_min:i_max) - tbc_bulk_id%tbsys%H%data_d(1)%data(o_min:o_max,i_min:i_max)
    Aoo = af%z(ip)*tbc_bulk_id%tbsys%S%data_d(1)%data(o_min:o_max,o_min:o_max) - tbc_bulk_id%tbsys%H%data_d(1)%data(o_min:o_max,o_min:o_max)

    call inverse(Aoo, Goo, positive_in = .false.)

    call matrix_product_sub(tt, Aio,Goo)
    call matrix_product_sub(SEii(ip)%data_z(1)%data, tt, Aoi)
    SEii(ip)%data_z(1)%data = -SEii(ip)%data_z(1)%data

    do i=1, at_i%N
      if (atoms_n_neighbours(at_i,i) == 4) then
	i_orb_min = tbc_i%tbsys%first_orb_of_atom(i)
	i_orb_max = i_orb_min + n_orbs_of_Z(tbc_i%tbsys%tbmodel, at_i%Z(i)) - 1
	SEii(ip)%data_z(1)%data(:,i_orb_min:i_orb_max) = 0.0_dp
	SEii(ip)%data_z(1)%data(i_orb_min:i_orb_max,:) = 0.0_dp
      else
	i_orb_min = tbc_i%tbsys%first_orb_of_atom(i)
	i_orb_max = i_orb_min + n_orbs_of_Z(tbc_i%tbsys%tbmodel, at_i%Z(i)) - 1
      	if (i_orb_min > 1) then
	  SEii(ip)%data_z(1)%data(:,1:i_orb_min-1) = 0.0_dp
	  SEii(ip)%data_z(1)%data(1:i_orb_min-1,:) = 0.0_dp
	endif
	if (i_orb_max < tbc_i%tbsys%N) then
	  SEii(ip)%data_z(1)%data(:,i_orb_max+1:tbc_i%tbsys%N) = 0.0_dp
	  SEii(ip)%data_z(1)%data(i_orb_max+1:tbc_i%tbsys%N,:) = 0.0_dp
	endif
      endif
    end do
  end do
  GF_E_i = calc_GF(tbc_SE_i, at_i, .true., af%Fermi_E, 0.1_dp, band_width, SelfEnergy=SEii, local_e = local_e_i_SE)
  call print_G_dev(at_i, tbc_SE_i%tbsys, tbc_SE_i%GF%G(1)%data_z(1)%data, Gii_exact, "nn SE")
  call print_local_e_dev(bulk, at_i, local_e, local_e_i_SE, "nn SE")

call Print("")
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call system_finalise()

end program
