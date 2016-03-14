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

!X
!X  libAtoms module
!X
!% This is a container module, "use" it if you want to use libatoms stuff:
!%> 	  use libAtoms_module
!% 
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

module libAtoms_module
  use error_module
  use system_module
  use linkedlist_module
  use MPI_context_module
  use units_module
  use linearalgebra_module
  use minimization_module
  use extendable_str_module
  use dictionary_module
  use table_module
  use spline_module
  use sparse_module
  use periodictable_module
  use cinoutput_module
  use paramreader_module
  use atoms_types_module
  use connection_module
  use atoms_module
  use quaternions_module
  use rigidbody_module
  use steinhardt_nelson_qw_module
  use constraints_module
  use thermostat_module
  use barostat_module
  use dynamicalsystem_module
  use clusters_module
  use structures_module
  use frametools_module
  use nye_tensor_module
  use Topology_module
  use atoms_ll_module
  use ringstat_module
  use histogram1d_module
  use domaindecomposition_module
  use k_means_clustering_module
  use SocketTools_module
  use partition_module
end module libAtoms_module
