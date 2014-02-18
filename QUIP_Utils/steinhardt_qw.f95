module qw_module

use libatoms_module
#ifdef HAVE_GAP
use descriptors_module, only: calc_qw
#endif

implicit none

private

#ifdef HAVE_GAP
public :: calc_steinhardt_bond_order
#endif

contains

#ifdef HAVE_GAP
  subroutine calc_steinhardt_bond_order(this,l,do_q,do_w,cutoff,cutoff_transition_width)
     type(atoms), intent(inout) :: this
     integer, intent(in) :: l
     logical, intent(in), optional :: do_q, do_w
     real(dp), optional :: cutoff, cutoff_transition_width

     call calc_qw(this,l,do_q,do_w,cutoff,cutoff_transition_width)
  endsubroutine calc_steinhardt_bond_order
#endif

endmodule qw_module
