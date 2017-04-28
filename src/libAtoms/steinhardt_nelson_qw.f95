#include "error.inc"

module steinhardt_nelson_qw_module

use error_module
use system_module
use units_module
use dictionary_module
use atoms_types_module
use atoms_module
use linearalgebra_module
use angular_functions_module

implicit none
private

public :: calc_qw, calc_qw_grad

contains

subroutine calc_qw(this,l,do_q,do_w,cutoff,cutoff_transition_width,mask,error)
   type(atoms), intent(inout) :: this
   integer, intent(in) :: l
   logical, intent(in), optional :: do_q, do_w
   real(dp), intent(in), optional :: cutoff, cutoff_transition_width
   logical, dimension(:), intent(in), optional :: mask
   integer, intent(out), optional :: error

   real(dp), dimension(:), pointer :: ql, wl
   real(dp), dimension(3) :: diff_ij
   real(dp) :: ql_crystal, wl_crystal, my_cutoff, my_cutoff_transition_width, r_ij, f_cut
   logical :: my_do_q, my_do_w, do_smooth_cutoff

   complex(dp), dimension(:), allocatable :: spherical_c
   complex(dp), dimension(:), allocatable :: spherical_c_crystal
   integer :: i, j, m, n, m1, m2, m3
   real(dp) :: n_bonds, n_bonds_crystal

   INIT_ERROR(error)

   my_do_q = optional_default(.false.,do_q)
   my_do_w = optional_default(.false.,do_w)
   my_cutoff = optional_default(this%cutoff,cutoff)
   my_cutoff_transition_width = optional_default(0.0_dp,cutoff_transition_width)

   do_smooth_cutoff = present(cutoff_transition_width)

   if(my_cutoff > this%cutoff) then
      call set_cutoff(this,my_cutoff)
      call calc_connect(this)
   endif
   
   if(my_do_q) call add_property(this, 'q'//l, 0.0_dp, ptr=ql)
   if(my_do_w) call add_property(this, 'w'//l, 0.0_dp, ptr=wl)

   if(present(mask)) then
      call check_size("mask",mask,this%N,"calc_qw",error)
   endif

   allocate(spherical_c(-l:l), spherical_c_crystal(-l:l))
   spherical_c = CPLX_ZERO
   spherical_c_crystal = CPLX_ZERO


   f_cut = 1.0_dp
   n_bonds_crystal = 0.0_dp

   do i = 1, this%N
      if(present(mask)) then
         if( .not. mask(i) ) cycle
      endif
      n_bonds = 0.0_dp
      do n = 1, n_neighbours(this,i)
	 j = neighbour(this,i,n,diff=diff_ij,distance=r_ij)

	 if(r_ij > my_cutoff) cycle

	 if(do_smooth_cutoff) then
	    f_cut = coordination_function(r_ij,my_cutoff,my_cutoff_transition_width)
	 endif
	 n_bonds = n_bonds + f_cut

	 do m = -l, l
	    spherical_c(m) = spherical_c(m) + SphericalYCartesian(l,m,diff_ij) * f_cut
	 enddo

      enddo

      ! spherical_c_crystal is defined as in Steinhardt and Nelson, from ( \sum_i spherical_c(i) ) / ( \sum_i n_bonds(i))
      ! An altenative would be \sum_i ( spherical_c(i) / n_bonds(i) )
      spherical_c_crystal = spherical_c_crystal + spherical_c
      n_bonds_crystal = n_bonds_crystal + n_bonds

      if(n_bonds > 0.0_dp) spherical_c = spherical_c / n_bonds

      if(my_do_q) ql(i) = sqrt( dot_product(spherical_c,spherical_c) * 4.0_dp * PI / (2.0_dp * l + 1.0_dp) )

      if(my_do_w) then
	 wl(i) = 0.0_dp
	 do m1 = -l, l
	    do m2 = -l, l
	       m3 = -m1-m2
	       if( m3 >= -l .and. m3 <= l ) wl(i) = wl(i) + spherical_c(m1) * spherical_c(m2) * spherical_c(m3) * wigner3j(l,m1,l,m2,l,m3)
	    enddo
	 enddo
	 
	 wl(i) = wl(i) / sqrt( dot_product(spherical_c,spherical_c)**3 )
      endif

   enddo

   if (n_bonds_crystal > 0.0_dp ) spherical_c_crystal = spherical_c_crystal / n_bonds_crystal
   if (my_do_q) then
      ql_crystal = sqrt( dot_product(spherical_c_crystal,spherical_c_crystal) * 4.0_dp * PI / (2.0_dp * l + 1.0_dp) )
      call set_value(this%params, "q"//l//"_global",ql_crystal)
   endif

   if(my_do_w) then
      wl_crystal = 0.0_dp
      do m1 = -l, l
	 do m2 = -l, l
	    m3 = -m1-m2
	    if( m3 >= -l .and. m3 <= l ) wl_crystal = wl_crystal + spherical_c_crystal(m1) * spherical_c_crystal(m2) * spherical_c_crystal(m3) * wigner3j(l,m1,l,m2,l,m3)
	 enddo
      enddo
      
      wl_crystal = wl_crystal / sqrt( dot_product(spherical_c_crystal(:),spherical_c_crystal(:))**3 )
      call set_value(this%params, "w"//l//"_global",wl_crystal)
   endif

   if(allocated(spherical_c)) deallocate(spherical_c)
   if(allocated(spherical_c_crystal)) deallocate(spherical_c_crystal)

endsubroutine calc_qw

subroutine calc_qw_grad(this,grad_ind,l,do_q,do_w,cutoff,cutoff_transition_width)
   type(atoms), intent(inout) :: this
   integer, intent(in) :: grad_ind
   integer, intent(in) :: l
   logical, intent(in), optional :: do_q, do_w
   real(dp), optional :: cutoff, cutoff_transition_width

   integer :: i, j, n, m
   logical :: my_do_q, my_do_w, do_smooth_cutoff
   real(dp) :: my_cutoff, my_cutoff_transition_width, ql
   real(dp) :: f_cut, df_cut(3), r_ij, diff_ij(3), n_bonds
   complex(dp), allocatable :: spherical_c(:), dspherical_c(:,:,:)
   real(dp), allocatable :: dn_bonds(:,:)
   real(dp), pointer :: ql_grad(:,:), wl_grad(:,:)

   my_do_q = optional_default(.false.,do_q)
   my_do_w = optional_default(.false.,do_w)
   my_cutoff = optional_default(this%cutoff,cutoff)
   my_cutoff_transition_width = optional_default(0.0_dp,cutoff_transition_width)

   do_smooth_cutoff = present(cutoff_transition_width)

   if(my_cutoff > this%cutoff) then
      call set_cutoff(this,my_cutoff)
      call calc_connect(this)
   endif

   if(my_do_q) call add_property(this, 'q'//l//'_grad', 0.0_dp, n_cols=3, ptr2=ql_grad)
   if(my_do_w) call add_property(this, 'w'//l//'_grad', 0.0_dp, n_cols=3, ptr2=wl_grad)

   if (my_do_q) ql_grad = 0.0_dp
   if (my_do_w) wl_grad = 0.0_dp

   allocate(spherical_c(-l:l), dspherical_c(3,-l:l,this%N), dn_bonds(3,this%N))
   spherical_c = CPLX_ZERO
   dspherical_c = CPLX_ZERO

   f_cut = 1.0_dp
   df_cut = 0.0_dp

   n_bonds = 0.0_dp
   dn_bonds = 0.0_dp

   i = grad_ind
   do n = 1, n_neighbours(this,i)
      j = neighbour(this,i,n,diff=diff_ij,distance=r_ij)

      if(r_ij > my_cutoff) cycle

      if(do_smooth_cutoff) then
	 f_cut = coordination_function(r_ij,my_cutoff,my_cutoff_transition_width)
	 df_cut(:) = -dcoordination_function(r_ij,my_cutoff,my_cutoff_transition_width)*diff_ij/r_ij
      endif
      n_bonds = n_bonds + f_cut
      dn_bonds(:,i) = dn_bonds(:,i) + df_cut(:)
      dn_bonds(:,j) = dn_bonds(:,j) - df_cut(:)

      do m = -l, l
	 spherical_c(m) = spherical_c(m) + SphericalYCartesian(l,m,diff_ij) * f_cut
	 dspherical_c(:,m,i) = dspherical_c(:,m,i) - GradSphericalYCartesian(l,m,diff_ij) * f_cut + &
						     SphericalYCartesian(l,m,diff_ij)*df_cut
	 dspherical_c(:,m,j) = dspherical_c(:,m,j) + GradSphericalYCartesian(l,m,diff_ij) * f_cut - &
						     SphericalYCartesian(l,m,diff_ij)*df_cut
      enddo
   enddo

   if(n_bonds > 0.0_dp) then
      do m=-l,l
	 dspherical_c(:,m,i) = ((n_bonds*dspherical_c(:,m,i))-(spherical_c(m)*dn_bonds(:,i))) / (n_bonds**2)
	 do n = 1, n_neighbours(this,i)
	    j = neighbour(this,i,n)
	    dspherical_c(:,m,j) = ((n_bonds*dspherical_c(:,m,j))-(spherical_c(m)*dn_bonds(:,j))) / (n_bonds**2)
	 end do
      end do
      spherical_c(:) = spherical_c(:) / n_bonds
   end if

   if(my_do_q) then
      ql = sqrt( dot_product(spherical_c(:),spherical_c(:)) * 4.0_dp * PI / (2.0_dp * l + 1.0_dp) )
      do m=-l,l
	 ql_grad(:,i) = ql_grad(:,i) + (1.0_dp/ql)*spherical_c(m)*conjg(dspherical_c(:,m,i)) * (4.0_dp * PI / (2.0_dp * l + 1.0_dp) )
	 do n = 1, n_neighbours(this,i)
	    j = neighbour(this,i,n)
	    ql_grad(:,j) = ql_grad(:,j) + (1.0_dp/ql)*spherical_c(m)*conjg(dspherical_c(:,m,j)) * (4.0_dp * PI / (2.0_dp * l + 1.0_dp) )
	 end do
      end do
   endif

   if(my_do_w) then
      call system_abort("calc_qw_grad does not support grad of w yet")
   endif

end subroutine calc_qw_grad

end module steinhardt_nelson_qw_module


