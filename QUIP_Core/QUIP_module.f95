!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X     QUIP: quantum mechanical and interatomic potential simulation package
!X     
!X     Portions written by Noam Bernstein, while working at the
!X     Naval Research Laboratory, Washington DC. 
!X
!X     Portions written by Gabor Csanyi, Copyright 2006-2007.   
!X
!X     When using this software,  please cite the following reference:
!X
!X     reference
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X QUIP Module 
!X
!% This is a container module. ``use" it if you want to use QUIP stuff
!%>       use QUIP_module 
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

module QUIP_internals_module
  use IPModel_lj_module
  use IPModel_sw_module
  use IPModel_tersoff_module
  use functions_module
  use ewald_module
  use tb_mixing_module
  use quip_common_module
  use tb_common_module
  use tbmodel_bowler_module
  use tbmodel_dftb_module
  use tbmodel_nrl_tb_module
end module QUIP_internals_module

!% This is a container module. ``use" it if you want to use QUIP stuff
!%>       use QUIP_module 
module QUIP_module
  use tbsystem_module
  use tb_module
  use tbmodel_module
  use tb_kpoints_module
  use matrix_module
  use approxfermi_module
  use ip_module
  use tb_greensfunctions_module
  use mpi_context_module
  use potential_module
  use rs_sparsematrix_module
  use scalapack_module
  use tbmatrix_module
  use metapotential_module
end module QUIP_module
