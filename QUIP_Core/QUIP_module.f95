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
!X QUIP Module 
!X
!% This is a container module. ``use" it if you want to use QUIP stuff
!%>       use QUIP_module 
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

module QUIP_internals_module
  use functions_module
  use quip_common_module
  use ewald_module
  use IPModel_lj_module
  use IPModel_sw_module
  use IPModel_tersoff_module
#ifdef HAVE_TB
  use tb_mixing_module
  use tb_common_module
  use tbmodel_bowler_module
  use tbmodel_dftb_module
  use tbmodel_nrl_tb_module
#endif
end module QUIP_internals_module

!% This is a container module. ``use" it if you want to use QUIP stuff
!%>       use QUIP_module 
module QUIP_module
  use mpi_context_module
  use scalapack_module
  use matrix_module
  use rs_sparsematrix_module
#ifdef HAVE_TB
  use approxfermi_module
  use tb_kpoints_module
  use tbmodel_module
  use tbmatrix_module
  use tbsystem_module
  use tb_greensfunctions_module
  use tb_module
#endif
  use ip_module
  use potential_module
end module QUIP_module
