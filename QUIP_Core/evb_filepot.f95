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

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X EVB_filepot_template program
!X
!% filepot program for EVB energy and force calculation through CP2K
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

program EVB_filepot_template

  use libatoms_module
  use potential_module

  implicit none

    type(Atoms)                   :: at1,at2
    real(dp), pointer             :: at1_pos(:,:), at2_pos(:,:)
    type(Potential)               :: cp2k_fast_pot, pot
    character(STRING_LENGTH)       :: args_str
    integer                       :: stat, error
    real(dp)                      :: energy
    real(dp), dimension(:,:), allocatable :: f1

    !ARGS_STR INPUT PARAMETERS
    type(Dictionary)        :: cli_params
    type(Dictionary)        :: params
    character(STRING_LENGTH) :: evb_params_line
    character(FIELD_LENGTH) :: posfilename
    character(FIELD_LENGTH) :: xyz_template_filename
    character(FIELD_LENGTH) :: outfilename
    logical                 :: only_GAP
    character(FIELD_LENGTH) :: psf_file1
    character(FIELD_LENGTH) :: psf_file2
    character(FIELD_LENGTH) :: filepot_program
    character(FIELD_LENGTH) :: cp2k_calc_args
    character(FIELD_LENGTH) :: evb_params_filename
    character(FIELD_LENGTH) :: mm_args_str
    character(FIELD_LENGTH) :: topology_suffix1
    character(FIELD_LENGTH) :: topology_suffix2
    integer                 :: form_bond(2), break_bond(2)
    real(dp)                :: diagonal_dE2, &
                               offdiagonal_A12, &
                               offdiagonal_mu12, &
                               offdiagonal_mu12_square, &
                               offdiagonal_r0
    logical                 :: save_forces, &
                               save_energy
    type(Inoutput)          :: evb_params_file
    logical                 :: exists
    character(len=FIELD_LENGTH) :: filename
    integer                        :: tmp_run_dir_i


    call system_initialise(verbosity=PRINT_SILENT,enable_timing=.true.)
    call verbosity_push(PRINT_NORMAL)
    mainlog%prefix="EVB_FILEPOT"
    call system_timer('evb_filepot_template')

    !INPUT PARAMETERS
    call initialise(cli_params)
    call param_register(cli_params,"pos","", posfilename, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(cli_params,"xyz_template","", xyz_template_filename, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(cli_params,"out","", outfilename, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(cli_params,"only_GAP","F", only_gap, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(cli_params,"psf_file1","", psf_file1, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(cli_params,"psf_file2","", psf_file2, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(cli_params,"filepot_program","", filepot_program, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(cli_params,"cp2k_calc_args","", cp2k_calc_args, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(cli_params,"param_file","evb_params.dat", evb_params_filename, help_string="No help yet.  This source file was $LastChangedBy$")
    if (.not. param_read_args(cli_params)) then
      !call print("Usage: pos xyz_template out only_GAP param_file",PRINT_ALWAYS)
      call system_abort('EVB_filepot: could not parse argument line')
    endif
    call finalise(cli_params)


    !EVB params file with EVB parameters and with params that can overwrite the defaults
    inquire(file=trim(evb_params_filename), exist=exists)
    if (.not.exists) call system_abort("WARNING: no "//trim(evb_params_filename)//" file found.")

    call initialise(evb_params_file,trim(evb_params_filename),ACTION=INPUT)
    evb_params_line=read_line(evb_params_file,stat)
    if (len_trim(evb_params_line)>0) then
       call initialise(params)
       call param_register(params,"pos",""//trim(posfilename), posfilename, help_string="No help yet.  This source file was $LastChangedBy$")
       call param_register(params,"xyz_template",""//trim(xyz_template_filename), xyz_template_filename, help_string="No help yet.  This source file was $LastChangedBy$")
       call param_register(params,"out",""//trim(outfilename), outfilename, help_string="No help yet.  This source file was $LastChangedBy$")
       call param_register(params,"only_GAP",""//only_gap, only_gap, help_string="No help yet.  This source file was $LastChangedBy$")
       call param_register(params,"psf_file1",""//trim(psf_file1), psf_file1, help_string="No help yet.  This source file was $LastChangedBy$")
       call param_register(params,"psf_file2",""//trim(psf_file2), psf_file2, help_string="No help yet.  This source file was $LastChangedBy$")
       call param_register(params,"filepot_program",""//trim(filepot_program), filepot_program, help_string="No help yet.  This source file was $LastChangedBy$")
       call param_register(params,"cp2k_calc_args",""//trim(cp2k_calc_args), cp2k_calc_args, help_string="No help yet.  This source file was $LastChangedBy$")
       call param_register(params, 'mm_args_str', "", mm_args_str, help_string="No help yet.  This source file was $LastChangedBy$")
       call param_register(params, 'topology_suffix1', "_EVB1", topology_suffix1, help_string="No help yet.  This source file was $LastChangedBy$")
       call param_register(params, 'topology_suffix2', "_EVB2", topology_suffix2, help_string="No help yet.  This source file was $LastChangedBy$")
       call param_register(params, 'form_bond', '0 0', form_bond, help_string="No help yet.  This source file was $LastChangedBy$")
       call param_register(params, 'break_bond', '0 0', break_bond, help_string="No help yet.  This source file was $LastChangedBy$")
       call param_register(params, 'diagonal_dE2', '0.0', diagonal_dE2, help_string="No help yet.  This source file was $LastChangedBy$")
       call param_register(params, 'offdiagonal_A12', '0.0', offdiagonal_A12, help_string="No help yet.  This source file was $LastChangedBy$")
       call param_register(params, 'offdiagonal_mu12', '0.0', offdiagonal_mu12, help_string="No help yet.  This source file was $LastChangedBy$")
       call param_register(params, 'offdiagonal_mu12_square', '0.0', offdiagonal_mu12_square, help_string="No help yet.  This source file was $LastChangedBy$")
       call param_register(params, 'offdiagonal_r0', '0.0', offdiagonal_r0, help_string="No help yet.  This source file was $LastChangedBy$")
       call param_register(params, 'save_forces', "T", save_forces, help_string="No help yet.  This source file was $LastChangedBy$")
       call param_register(params, 'save_energy', 'T', save_energy, help_string="No help yet.  This source file was $LastChangedBy$")
       call param_register(params, 'tmp_run_dir_i', '-1', tmp_run_dir_i, help_string="No help yet.  This source file was $LastChangedBy$")
       if (.not. param_read_line(params, evb_params_line, ignore_unknown=.true.,task='EVB_filepot_template args_str')) then
          call system_abort('EVB_filepot_template failed to parse args_str="'//trim(evb_params_line)//'"')
       endif
       call finalise(params)
    else
       call system_abort("WARNING: file with EVB params "//trim(evb_params_filename)//" is empty.  It must contain 1 only line with all the params.")
    endif
    evb_params_line=read_line(evb_params_file,stat)
    if (stat==0 .and. len_trim(evb_params_line)>0) then
       call system_abort("File with EVB params "//trim(evb_params_filename)//" has more than 1 lines and the second line is not empty.  It must contain 1 only line with all the params.")
    endif
    call finalise(evb_params_file)


    !remove old output file
    call system_command("rm -f "//trim(outfilename),stat)


    !read positions and template XYZ with lattice
    if (len_trim(posfilename)==0) call system_abort("No posfilename given.")
    inquire(file=trim(posfilename), exist=exists)
    if (.not.exists) call system_abort("No "//trim(posfilename)//" file found.")
    if (len_trim(xyz_template_filename)==0) call system_abort("No xyz_template_filename given.")
    inquire(file=trim(xyz_template_filename), exist=exists)
    if (.not.exists) call system_abort("No "//trim(xyz_template_filename)//" file found.")
    !call print("Reading at1")
    call read(at1,posfilename)
    !call print("Reading at2")
    call read(at2,xyz_template_filename)

    if (at1%N/=at2%N) call system_abort("EVB_filepot_template: at1 and at2 has different number of atoms.")
    if (.not.assign_pointer(at1,"pos",at1_pos)) call system_abort("EVB_filepot_template: Could not find pos property in at1.")
    if (.not.assign_pointer(at2,"pos",at2_pos)) call system_abort("EVB_filepot_template: Could not find pos property in at2.")

    at2_pos(1:3,1:at2%N) = at1_pos(1:3,1:at2%N)

    !possible i/o on /tmp
    if (tmp_run_dir_i>0) then
       filename="/tmp/cp2k_run_"//tmp_run_dir_i//"/filepot"
    else
       filename="filepot"
    endif

    !call print("Writing filepot.0.xyz")
    !!call write(at2,"filepot.0.xyz", properties="species:pos:mol_id:atom_res_number",real_format='f17.10',error=error)
    !call write(at2,trim(filename)//".0.xyz", properties="species:pos",real_format='%17.10f',error=error)
         !HANDLE_ERROR(error)
    call finalise(at1)


    !call print("Initialise potentials")
    !setup EVB potential (pot and args_str)
    if (len_trim(filepot_program)==0) call system_abort("EVB_filepot_template: Filepot not given.")
    !!call initialise(cp2k_fast_pot,'FilePot command='//trim(filepot_program)//' property_list=species:pos:avgpos:mol_id:atom_res_number min_cutoff=0.0')
    call initialise(cp2k_fast_pot,'FilePot filename='//trim(filename)//' command='//trim(filepot_program)//' property_list=species:pos min_cutoff=0.0')

    args_str=""
    !if (len_trim(cp2k_calc_args)==0) call print("WARNING: EVB_filepot_template: cp2k_calc_args not given.",PRINT_ALWAYS)
    if (len_trim(cp2k_calc_args)>0) mm_args_str=trim(mm_args_str)//" "//trim(cp2k_calc_args)
    if (len_trim(mm_args_str)>0) &
       args_str=trim(args_str)//" mm_args_str="//'"'//trim(mm_args_str)//'"'
    if (len_trim(topology_suffix1)>0) &
       args_str=trim(args_str)//" topology_suffix1="//trim(topology_suffix1)
    if (len_trim(topology_suffix2)>0) &
       args_str=trim(args_str)//" topology_suffix2="//trim(topology_suffix2)
    if (any(form_bond/=0)) &
       args_str=trim(args_str)//' form_bond={'//form_bond//'}'
    if (any(break_bond/=0)) &
       args_str=trim(args_str)//' break_bond={'//break_bond//'}'
    if (.not.(diagonal_dE2.feq.0.0_dp)) &
       args_str=trim(args_str)//" diagonal_dE2="//diagonal_dE2
    if (.not.(offdiagonal_A12.feq.0.0_dp)) &
       args_str=trim(args_str)//" offdiagonal_A12="//offdiagonal_A12
    if (.not.(offdiagonal_mu12.feq.0.0_dp)) &
       args_str=trim(args_str)//" offdiagonal_mu12="//offdiagonal_mu12
    if (.not.(offdiagonal_mu12_square.feq.0.0_dp)) &
       args_str=trim(args_str)//" offdiagonal_mu12_square="//offdiagonal_mu12_square
    if (.not.(offdiagonal_r0.feq.0.0_dp)) &
       args_str=trim(args_str)//" offdiagonal_r0="//offdiagonal_r0
    args_str=trim(args_str)//" save_forces="//save_forces
    args_str=trim(args_str)//" save_energy="//save_energy

    call initialise(pot,args_str='EVB=T '//trim(args_str), &
       pot1=cp2k_fast_pot)
    !call print(pot)


    !copy psf files if needed
    if (len_trim(psf_file1)>0) then
       !call print("cp "//trim(psf_file1)//" quip_cp2k"//trim(topology_suffix1)//".psf")
       call system_command("cp "//trim(psf_file1)//" quip_cp2k"//trim(topology_suffix1)//".psf")
    endif
    if (len_trim(psf_file2)>0) then
       !call print("cp "//trim(psf_file2)//" quip_cp2k"//trim(topology_suffix2)//".psf")
       call system_command("cp "//trim(psf_file2)//" quip_cp2k"//trim(topology_suffix2)//".psf")
    endif


    !call print("Calculate energy and forces")
    !calc E and F with EVB potential
    allocate(f1(3,at2%N))
    call calc(pot,at2,energy=energy,force=f1,args_str=" EVB_gap=gap calc_force=forces calc_energy=energy")
    deallocate(f1)


    !call print("Write output")
    if (only_GAP) then
       call write(at2,filename=trim(outfilename),properties="gap_force")
    else
       call write(at2,filename=trim(outfilename))
    endif
       call finalise(at2)


    call system_timer('evb_filepot_template')
    call system_finalise

end program EVB_filepot_template
