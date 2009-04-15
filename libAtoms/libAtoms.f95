!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X     libAtoms: atomistic simulation library
!X     
!X     Copyright 2006-2007.
!X
!X     Authors: Gabor Csanyi, Steven Winfield, James Kermode
!X     Contributors: Noam Bernstein, Alessio Comisso
!X
!X     The source code is released under the GNU General Public License,
!X     version 2, http://www.gnu.org/copyleft/gpl.html
!X
!X     If you would like to license the source code under different terms,
!X     please contact Gabor Csanyi, gabor@csanyi.net
!X
!X     When using this software, please cite the following reference:
!X
!X     http://www.libatoms.org
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X  libAtoms module
!X
!% This is a container module, "use" it if you want to use libatoms stuff:
!%> 	  use libAtoms_module
!% 
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! $Id: libAtoms.f95,v 1.2 2007-09-21 14:31:00 nb326 Exp $
! $Log: not supported by cvs2svn $
! Revision 1.1  2007/07/09 14:49:54  gc121
! container module for the whole package
!

module libAtoms_module
  use system_module
  use units_module
  use linearalgebra_module
  use minimization_module
  use table_module
  use spline_module
  use sparse_module
  use periodictable_module
  use paramreader_module
  use rigidbody_module
  use quaternions_module
  use extendable_str_module
  use dictionary_module
  use constraints_module
  use atoms_module
  use thermostat_module
  use dynamicalsystem_module
  use structures_module
  use frametools_module
  use nye_tensor_module
  use clusters_module
  use Topology_module
  use cinoutput_module

#ifdef HAVE_NETCDF_F
  use netcdf_module
#endif HAVE_NETCDF_F

end module libAtoms_module
