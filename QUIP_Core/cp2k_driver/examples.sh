#!/bin/sh

#../test_CP2K Run_Type=MM coord_file=diala_zw.xyz Topology_Print=-1 Connect_Cutoff=1.8 > water_test_CP2K_MM_aala.out
#../test_CP2K Topology_Print=-2 coord_file=diala_zw.xyz Connect_Cutoff=1.8 > water_test_CP2K_MM_aala.out
#../test_CP2K Topology_Print=-2 coord_file=aala.xyz Connect_Cutoff=1.8 > water_test_CP2K_MM_water.out
#../test_CP2K coord_file=aala.xyz Topology_Print=-1 > water_test_CP2K_MM_aala.out
#../test_CP2K coord_file=diala_zw.xyz Topology_Print=-1 > water_test_CP2K_MM_aala.out
#../test_CP2K coord_file=coord_MCL2.xyz > water_test_CP2K_MM_MCL2.out
#../test_CP2K Run_Type=QS QM_flag=-1 coord_file=H2O32_300K.xyz serial=F > water_test_CP2K_QS.out
#../test_CP2K Run_Type=QMMM QM_flag=1 QM_flag_MM=1 coord_file=H2O32_300K.xyz Topology_Print=-1 qm_list_filename=qmlist.dat serial=T > water_test_CP2K_QMMM_H2O.out
#../test_CP2K Run_Type=QMMM QM_flag=1 QM_flag_MM=0 coord_file=H2O32_300K.xyz qm_list_filename=qmlist.dat serial=F > water_test_CP2K_QMMM_H2O.out
#../test_CP2K Run_Type=QMMM QM_flag=2 QM_flag_MM=2 Inner_QM_Radius=2.5 Outer_QM_Radius=3.0 Equilib_Time=0.0 Run_Time=0. coord_file=H2O-512.xyz new_coord_file=water_test_CP2K_QMMM_H2O512.xyz qm_list_filename=qmlist.dat.ext serial=F > water_test_CP2K_QMMM_H2O512.out
#../test_CP2K Run_Type=QMMM QM_flag=2 QM_flag_MM=2 Inner_QM_Radius=3.0 Outer_QM_Radius=3.0 Equilib_Time=0.0 Run_Time=0. coord_file=H2O-512.xyz new_coord_file=water_test_CP2K_QMMM_H2O512_ext.xyz qm_list_filename=qmlist.dat.ext serial=F > water_test_CP2K_QMMM_H2O512_ext.out
#../test_CP2K Run_Type=QMMM coord_file=coord_MCL2.xyz QM_flag=1 serial=F > water_test_CP2K_QMMM_MCL2.out
