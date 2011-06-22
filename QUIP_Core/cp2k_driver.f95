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
!X cp2k_driver_template_module
!X
!% Driver for CP2K code using a template input file
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#include "error.inc"
module cp2k_driver_template_module
use libatoms_module
use topology_module
implicit none

private
integer, parameter, private :: CP2K_LINE_LENGTH = 1024 !Max line length to be printed into the CP2K input files

public :: do_cp2k_calc


contains

  subroutine do_cp2k_calc_fake(at, f, e, args_str)
    type(Atoms), intent(inout) :: at
    real(dp), intent(out) :: f(:,:), e
    character(len=*), intent(in) :: args_str

    type(inoutput) :: last_run_io
    type(cinoutput) :: force_cio
    character(len=FIELD_LENGTH) :: last_run_s
    integer :: this_run_i
    integer :: stat
    type(Atoms) :: for
    real(dp), pointer :: frc(:,:)

    call initialise(last_run_io, "cp2k_driver_fake_run", action=INPUT)
    last_run_s = read_line(last_run_io, status=stat)
    call finalise(last_run_io)
    if (stat /= 0) then
      this_run_i = 1
    else
      read (fmt=*,unit=last_run_s) this_run_i
      this_run_i = this_run_i + 1
    endif

    call print("do_cp2k_calc_fake run_i " // this_run_i, PRINT_ALWAYS)

    call initialise(force_cio, "cp2k_force_file_log")
    call read(force_cio, for, frame=this_run_i-1)
    !NB why does this crash now?
    ! call finalise(force_cio)
    if (.not. assign_pointer(for, 'frc', frc)) &
      call system_abort("do_cp2k_calc_fake couldn't find frc field in force log file")
    f = frc

    if (.not. get_value(for%params, "energy", e)) then
      if (.not. get_value(for%params, "Energy", e)) then
	if (.not. get_value(for%params, "E", e)) then
	  call system_abort("do_cp2k_calc_fake didn't find energy")
	endif
      endif
    endif

    e = e * HARTREE
    f  = f * HARTREE/BOHR 

    call initialise(last_run_io, "cp2k_driver_fake_run", action=OUTPUT)
    call print(""//this_run_i, file=last_run_io)
    call finalise(last_run_io)

  end subroutine do_cp2k_calc_fake


  subroutine do_cp2k_calc(at, f, e, args_str, error)
    type(Atoms), intent(inout) :: at
    real(dp), intent(out) :: f(:,:), e
    character(len=*), intent(in) :: args_str
    integer, intent(out), optional :: error

    type(Dictionary) :: cli
    character(len=FIELD_LENGTH) :: run_type, cp2k_template_file, psf_print, cp2k_program, link_template_file, topology_suffix
    logical :: clean_up_files, save_output_files, save_output_wfn_files, use_buffer
    integer :: clean_up_keep_n
    integer :: max_n_tries
    real(dp) :: max_force_warning
    real(dp) :: qm_vacuum
    real(dp) :: centre_pos(3), cp2k_box_centre_pos(3)
    logical :: auto_centre, has_centre_pos
    logical :: try_reuse_wfn
    character(len=FIELD_LENGTH) :: calc_qm_charges

    character(len=128) :: method

    type(Inoutput) :: template_io
    integer :: template_n_lines
    character(len=FIELD_LENGTH), allocatable :: cp2k_template_a(:)
    type(Inoutput) :: link_template_io
    integer :: link_template_n_lines
    character(len=FIELD_LENGTH), allocatable :: link_template_a(:)
    integer :: i_line

    character(len=FIELD_LENGTH) :: run_dir

    type(Table) :: qm_list
    type(Table) :: cut_bonds
    integer, pointer :: cut_bonds_p(:,:)
    integer, allocatable :: qm_list_a(:)
    integer, allocatable :: link_list_a(:)
    integer, allocatable :: qm_and_link_list_a(:)
    integer :: i_inner, i_outer
    logical :: inserted_atoms
    integer :: counter

    integer :: charge
    logical :: do_lsd

    logical :: can_reuse_wfn, qm_list_changed
    character(len=FIELD_LENGTH) :: qm_name_postfix

    logical :: use_QM, use_MM, use_QMMM
    logical :: cp2k_calc_fake

    integer, pointer :: isolated_atom(:)
    integer i, j, atno, insert_pos
    real(dp) :: cur_qmmm_qm_abc(3), old_qmmm_qm_abc(3)

    type(Atoms) :: at_cp2k

    character(len=TABLE_STRING_LENGTH), pointer :: from_str(:), to_str(:)
    character(len=TABLE_STRING_LENGTH) :: dummy_s
    real(dp), pointer :: from_dp(:), to_dp(:)
    integer, pointer :: from_i(:), to_i(:)

    integer, pointer :: old_cluster_mark_p(:), cluster_mark_p(:)
    logical :: dummy, have_silica_potential
    type(Table) :: intrares_impropers

    integer, pointer :: sort_index_p(:), saved_rev_sort_index_p(:)
    integer, allocatable :: rev_sort_index(:)
    integer :: at_i, iri_i
    type(Inoutput) :: rev_sort_index_io
    logical :: sorted

    integer :: run_dir_i, force_run_dir_i, delete_dir_i

    logical :: at_periodic
    integer :: form_bond(2), break_bond(2)
    integer :: form_bond_sorted(2), break_bond_sorted(2)

    character(len=FIELD_LENGTH) :: tmp_MM_param_filename, tmp_QM_pot_filename, tmp_QM_basis_filename
    character(len=FIELD_LENGTH) :: MM_param_filename, QM_pot_filename, QM_basis_filename
    logical :: truncate_parent_dir
    character(len=FIELD_LENGTH) :: dir, tmp_run_dir
    integer :: tmp_run_dir_i, stat
    logical :: exists

    INIT_ERROR(error)

    call system_timer('do_cp2k_calc')

    call initialise(cli)
      call param_register(cli, 'Run_Type', PARAM_MANDATORY, run_type, help_string="Type of run QS, MM, or QMMM")
      call param_register(cli, 'use_buffer', 'T', use_buffer, help_string="If true, use buffer as specified in relevant hybrid_mark")
      call param_register(cli, 'qm_name_postfix', '', qm_name_postfix, help_string="String to append to various marks and saved info to indicate distinct sets of calculations or QM/MM QM regions")
      call param_register(cli, 'cp2k_template_file', 'cp2k_input.template', cp2k_template_file, help_string="filename for cp2k input template")
      call param_register(cli, "qmmm_link_template_file", "", link_template_file, help_string="filename for cp2k link atoms template file")
      call param_register(cli, 'PSF_print', 'NO_PSF', psf_print, help_string="when to print PSF file: NO_PSF, DRIVER_PRINT_AND_SAVE, USE_EXISTING_PSF")
      call param_register(cli, "topology_suffix", "", topology_suffix, help_string="suffix to append to file containing topology info (for runs that do multiple topologies not to accidentally reuse PSF file")
      call param_register(cli, 'cp2k_program', PARAM_MANDATORY, cp2k_program, help_string="path to cp2k executable")
      call param_register(cli, 'clean_up_files', 'T', clean_up_files, help_string="if true, clean up run directory files")
      call param_register(cli, 'clean_up_keep_n', '1', clean_up_keep_n, help_string="number of old run directories to keep if cleaning up")
      call param_register(cli, 'save_output_files', 'T', save_output_files, help_string="if true, save the output files")
      call param_register(cli, 'save_output_wfn_files', 'F', save_output_wfn_files, help_string="if true, save output wavefunction files")
      call param_register(cli, 'max_n_tries', '2', max_n_tries, help_string="max number of times to run cp2k on failure")
      call param_register(cli, 'max_force_warning', '2.0', max_force_warning, help_string="generate warning if any force is larger than this")
      call param_register(cli, 'qm_vacuum', '6.0', qm_vacuum, help_string="amount of vacuum to add to size of qm region in hybrid (and nonperiodic?) runs")
      call param_register(cli, 'try_reuse_wfn', 'T', try_reuse_wfn, help_string="if true, try to reuse previous wavefunction file")
      call param_register(cli, 'have_silica_potential', 'F', have_silica_potential, help_string="if true, use 2.8A SILICA_CUTOFF for the connectivities")
      call param_register(cli, 'auto_centre', 'F', auto_centre, help_string="if true, automatically center configuration.  May cause energy/force fluctuations.  Mutually exclusive with centre_pos")
      call param_register(cli, 'centre_pos', '0.0 0.0 0.0', centre_pos, has_value_target=has_centre_pos, help_string="position to center around, mutually exclusive with auto_centre")
      call param_register(cli, 'cp2k_calc_fake', 'F', cp2k_calc_fake, help_string="if true, do fake cp2k runs that just read from old output files")
      call param_register(cli, 'form_bond', '0 0', form_bond, help_string="extra bond to form (for EVB)")
      call param_register(cli, 'break_bond', '0 0', break_bond, help_string="bond to break (for EVB)")
      call param_register(cli, 'qm_charges', '', calc_qm_charges, help_string="if not blank, name of property to put QM charges in")
      call param_register(cli, 'force_run_dir_i', '-1', force_run_dir_i, help_string="if > 0, force to run in this # run directory")
      call param_register(cli, 'tmp_run_dir_i', '-1', tmp_run_dir_i, help_string="if >0, the cp2k run directory will be /tmp/cp2k_run_$tmp_run_dir_i$, and all input files are also copied here when first called")
      call param_register(cli, 'MM_param_file', '', MM_param_filename, help_string="If tmp_run_dir>0, where to find MM parameter file to copy it to the cp2k run dir on /tmp.") !charmm.pot
      call param_register(cli, 'QM_potential_file', '', QM_pot_filename, help_string="If tmp_run_dir>0, where to find QM POTENTIAL file to copy it to the cp2k run dir on /tmp.") !POTENTIAL
      call param_register(cli, 'QM_basis_file', '', QM_basis_filename, help_string="If tmp_run_dir>0, where to find QM BASIS_SET file to copy it to the cp2k run dir on /tmp.") !BASIS_SET
      ! should really be ignore_unknown=false, but higher level things pass unneeded arguments down here
      if (.not.param_read_line(cli, args_str, ignore_unknown=.true.,task='cp2k_filepot_template args_str')) &
	call system_abort('cp2k_driver could not parse argument line')
    call finalise(cli)

    if (cp2k_calc_fake) then
      call print("do_fake cp2k calc calculation")
      call do_cp2k_calc_fake(at, f, e, args_str)
      return
    endif

    call print("do_cp2k_calc command line arguments")
    call print("  Run_Type " // Run_Type)
    call print("  use_buffer " // use_buffer)
    call print("  qm_name_postfix " // qm_name_postfix)
    call print("  cp2k_template_file " // cp2k_template_file)
    call print("  qmmm_link_template_file " // link_template_file)
    call print("  PSF_print " // PSF_print)
    call print("  clean_up_files " // clean_up_files)
    call print("  clean_up_keep_n " // clean_up_keep_n)
    call print("  save_output_files " // save_output_files)
    call print("  save_output_wfn_files " // save_output_wfn_files)
    call print("  max_n_tries " // max_n_tries)
    call print("  max_force_warning " // max_force_warning)
    call print("  qm_vacuum " // qm_vacuum)
    call print("  try_reuse_wfn " // try_reuse_wfn)
    call print('  have_silica_potential '//have_silica_potential)
    call print('  auto_centre '//auto_centre)
    call print('  centre_pos '//centre_pos)
    call print('  calc_qm_charges '//trim(calc_qm_charges))
    call print('  tmp_run_dir_i '//tmp_run_dir_i)
    call print('  MM_param_file '//trim(MM_param_filename))
    call print('  QM_potential_file '//trim(QM_pot_filename))
    call print('  QM_basis_file '//trim(QM_basis_filename))

    if (auto_centre .and. has_centre_pos) &
      call system_abort("do_cp2k_calc got both auto_centre and centre_pos, don't know which centre (automatic or specified) to shift to origin")

    if (tmp_run_dir_i>0 .and. clean_up_keep_n > 0) then
      call system_abort("do_cp2k_calc got both tmp_run_dir_i(only write on /tmp) and clean_up_keep_n (save in home).")
    endif

    !create run directory now, because it is needed if running on /tmp
    if (tmp_run_dir_i>0) then
      tmp_run_dir = "/tmp/cp2k_run_"//tmp_run_dir_i
      run_dir = link_run_directory(trim(tmp_run_dir), basename="cp2k_run", run_dir_i=run_dir_i)
      !and copy necessary files for access on /tmp if not yet present
      if (len_trim(MM_param_filename)>0) then
         tmp_MM_param_filename = trim(MM_param_filename)
         truncate_parent_dir=.true.
         do while(truncate_parent_dir)
            if (tmp_MM_param_filename(1:3)=="../") then
               tmp_MM_param_filename=trim(tmp_MM_param_filename(4:))
            else
               truncate_parent_dir=.false.
            endif
         enddo
         if (len_trim(tmp_MM_param_filename)==0) call system_abort("Empty tmp_MM_param_filename string")
         call print("if [ ! -s "//trim(tmp_run_dir)//"/"//trim(tmp_MM_param_filename)//" ] ; then echo 'copy charmm.pot' ; cp "//trim(MM_param_filename)//" "//trim(tmp_run_dir)//"/ ; else echo 'reuse charmm.pot' ; fi")
         call system_command("if [ ! -s "//trim(tmp_run_dir)//"/"//trim(tmp_MM_param_filename)//" ] ; then echo 'copy charmm.pot' ; cp "//trim(MM_param_filename)//" "//trim(tmp_run_dir)//"/ ; fi",status=stat)
         if ( stat /= 0 ) call system_abort("Something went wrong when tried to copy "//trim(MM_param_filename)//" into the tmp dir "//trim(tmp_run_dir))
      endif
      if (len_trim(QM_pot_filename)>0) then
         tmp_QM_pot_filename = trim(QM_pot_filename)
         truncate_parent_dir=.true.
         do while(truncate_parent_dir)
            if (tmp_QM_pot_filename(1:3)=="../") then
               tmp_QM_pot_filename=trim(tmp_QM_pot_filename(4:))
            else
               truncate_parent_dir=.false.
            endif
         enddo
         call print("if [ ! -s "//trim(tmp_run_dir)//"/"//trim(QM_pot_filename)//" ] ; then cp "//trim(QM_pot_filename)//" "//trim(tmp_run_dir)//"/ ; fi")
         call system_command("if [ ! -s "//trim(tmp_run_dir)//"/"//trim(tmp_QM_pot_filename)//" ] ; then echo 'copy QM potential' ; cp "//trim(QM_pot_filename)//" "//trim(tmp_run_dir)//"/ ; fi")
         if ( stat /= 0 ) call system_abort("Something went wrong when tried to copy "//trim(QM_pot_filename)//" into the tmp dir "//trim(tmp_run_dir))
      endif
      if (len_trim(QM_basis_filename)>0) then
         tmp_QM_basis_filename = trim(QM_basis_filename)
         truncate_parent_dir=.true.
         do while(truncate_parent_dir)
            if (tmp_QM_basis_filename(1:3)=="../") then
               tmp_QM_basis_filename=trim(tmp_QM_basis_filename(4:))
            else
               truncate_parent_dir=.false.
            endif
         enddo
         call print("if [ ! -s "//trim(tmp_run_dir)//"/"//trim(QM_basis_filename)//" ] ; then cp "//trim(QM_basis_filename)//" "//trim(tmp_run_dir)//"/ ; fi")
         call system_command("if [ ! -s "//trim(tmp_run_dir)//"/"//trim(tmp_QM_basis_filename)//" ] ; then echo 'copy QM basis' ; cp "//trim(QM_basis_filename)//" "//trim(tmp_run_dir)//"/ ; fi")
         if ( stat /= 0 ) call system_abort("Something went wrong when tried to copy "//trim(QM_basis_filename)//" into the tmp dir "//trim(tmp_run_dir))
      endif
    else
      run_dir = make_run_directory("cp2k_run", force_run_dir_i, run_dir_i)
    endif

    ! read template file
    if (tmp_run_dir_i>0) then
       call print("if [ ! -s "//trim(tmp_run_dir)//"/"//trim(cp2k_template_file)//" ] ; then cp "//trim(cp2k_template_file)//" "//trim(tmp_run_dir)//"/ ; fi")
       call system_command("if [ ! -s "//trim(tmp_run_dir)//"/"//trim(cp2k_template_file)//" ] ; then cp "//trim(cp2k_template_file)//" "//trim(tmp_run_dir)//"/ ; fi")
       if ( stat /= 0 ) call system_abort("Something went wrong when tried to copy "//trim(cp2k_template_file)//" into the tmp dir "//trim(tmp_run_dir))
       call initialise(template_io, trim(tmp_run_dir)//"/"//trim(cp2k_template_file), INPUT)
    else
       call initialise(template_io, trim(cp2k_template_file), INPUT)
    endif
    call read_file(template_io, cp2k_template_a, template_n_lines)
    call finalise(template_io)

    call prefix_cp2k_input_sections(cp2k_template_a(1:template_n_lines))

    if ( (trim(psf_print) /= 'NO_PSF') .and. &
         (trim(psf_print) /= 'DRIVER_PRINT_AND_SAVE') .and. &
         (trim(psf_print) /= 'USE_EXISTING_PSF')) &
      call system_abort("Unknown value for psf_print '"//trim(psf_print)//"'")

    ! parse run_type
    use_QM = .false.
    use_MM = .false.
    use_QMMM = .false.
    select case(trim(run_type))
      case("QS")
	use_QM=.true.
	method="QS"
      case("MM")
	use_MM=.true.
	method="Fist"
      case("QMMM")
	use_QM = .true.
	use_MM = .true.
	use_QMMM = .true.
	method="QMMM"
      case default
	call system_abort("Unknown run_type "//trim(run_type))
    end select

    ! prepare CHARMM params if necessary
    if (use_MM) then
      if (have_silica_potential) then
        call set_cutoff(at,SILICA_2body_CUTOFF)
      else
        call set_cutoff(at,0._dp)
      endif
      call calc_connect(at)

      call map_into_cell(at)
      call calc_dists(at)

    endif

    ! if writing PSF file, calculate residue labels, before sort
    if (run_type /= "QS") then
      if (trim(psf_print) == "DRIVER_PRINT_AND_SAVE") then
	call create_residue_labels_arb_pos(at,do_CHARMM=.true.,intrares_impropers=intrares_impropers,have_silica_potential=have_silica_potential,form_bond=form_bond,break_bond=break_bond)
      end if
    end if

    ! sort by molecule, residue ID
    call add_property(at, 'sort_index', 0, n_cols=1, ptr=sort_index_p, error=error)
    PASS_ERROR_WITH_INFO("Failed to add sort_index property", error)
    ! initialise sort index
    do at_i=1, at%N
      sort_index_p(at_i) = at_i
    end do

    ! do sort by read in order or labels
    sorted = .false.
    if (trim(psf_print) == 'USE_EXISTING_PSF') then ! read sort order
       ! add property for saved reverse sort indx
       call add_property(at, 'saved_rev_sort_index', 0, n_cols=1, ptr=saved_rev_sort_index_p, error=error)
       PASS_ERROR_WITH_INFO("Failed to add saved_rev_sort_index property", error)
       ! read it from file
       if (tmp_run_dir_i>0) then
         call system_command("if [ ! -s "//trim(tmp_run_dir)//"/quip_rev_sort_index"//trim(topology_suffix)// &
	                     " ] ; then cp quip_rev_sort_index"//trim(topology_suffix)//" /tmp/cp2k_run_"//tmp_run_dir_i//"/ ; fi",status=stat)
         if ( stat /= 0 ) then
	    call system_abort("Something went wrong when tried to copy quip_rev_sort_index"//trim(topology_suffix)// &
	                      " into the tmp dir "//trim(tmp_run_dir))
	 endif
         call initialise(rev_sort_index_io, trim(tmp_run_dir)//"/quip_rev_sort_index"//trim(topology_suffix), action=INPUT)
       else
         call initialise(rev_sort_index_io, "quip_rev_sort_index"//trim(topology_suffix), action=INPUT)
       endif
       call read_ascii(rev_sort_index_io, saved_rev_sort_index_p)
       call finalise(rev_sort_index_io)
       ! sort by it
       call atoms_sort(at, 'saved_rev_sort_index', error=error)
       PASS_ERROR_WITH_INFO ("do_cp2k_calc sorting atoms by read-in sort_index from quip_sort_order"//trim(topology_suffix), error)
       sorted = .true.
    endif
    if (trim(psf_print) == 'DRIVER_PRINT_AND_SAVE') then ! sort by labels
       if (has_property(at,'mol_id') .and. has_property(at,'atom_res_number')) then
	  if (has_property(at,'motif_atom_num')) then
	    call atoms_sort(at, 'mol_id', 'atom_res_number', 'motif_atom_num', error=error)
	  else
	    call atoms_sort(at, 'mol_id', 'atom_res_number', error=error)
	  endif
	  PASS_ERROR_WITH_INFO ("do_cp2k_calc sorting atoms by mol_id, atom_res_number, and motif_atom_num", error)
	  sorted = .true.
       endif
    endif

    allocate(rev_sort_index(at%N))
    do i=1, at%N
       rev_sort_index(sort_index_p(i)) = i
    end do

    if (sorted) then
      do at_i=1, at%N
	if (sort_index_p(at_i) /= at_i) then
	  call print("sort() of at%data reordered some atoms")
	  exit
	endif
      end do
      ! fix EVB bond forming/breaking indices for new sorted atom numbers
      call calc_connect(at)
      if ((all(form_bond > 0) .and. all(form_bond <= at%N)) .or. (all(break_bond > 0) .and. all(break_bond <= at%N))) then
	 if (all(form_bond > 0) .and. all(form_bond <= at%N)) form_bond_sorted(:) = rev_sort_index(form_bond(:))
	 if (all(break_bond > 0) .and. all(break_bond <= at%N)) break_bond_sorted(:) = rev_sort_index(break_bond(:))
      end if
      ! fix intrares impropers atom indices for new sorted atom numbers
      do iri_i=1, intrares_impropers%N
	intrares_impropers%int(1:4,iri_i) = rev_sort_index(intrares_impropers%int(1:4,iri_i))
      end do
    else
      call print("WARNING: didn't do sort_by_molecule - need saved sort_index or mol_id, atom_res_number, motif_atom_num.  CP2K may complain", PRINT_ALWAYS)
    end if

    ! write PSF file, if requested
    if (run_type /= "QS") then
      if (trim(psf_print) == "DRIVER_PRINT_AND_SAVE") then
	if (has_property(at, 'avgpos')) then
	  call write_psf_file_arb_pos(at, "quip_cp2k"//trim(topology_suffix)//".psf", run_type_string=trim(run_type),intrares_impropers=intrares_impropers, &
	    add_silica_23body=have_silica_potential,form_bond=form_bond_sorted,break_bond=break_bond_sorted)
	else if (has_property(at, 'pos')) then
	  call print("WARNING: do_cp2k_calc using pos for connectivity.  avgpos is preferred but not found.")
	  call write_psf_file_arb_pos(at, "quip_cp2k"//trim(topology_suffix)//".psf", run_type_string=trim(run_type),intrares_impropers=intrares_impropers, &
	    add_silica_23body=have_silica_potential,pos_field_for_connectivity='pos',form_bond=form_bond_sorted,break_bond=break_bond_sorted)
	else
	  call system_abort("do_cp2k_calc needs some pos field for connectivity (run_type='"//trim(run_type)//"' /= 'QS'), but found neither avgpos nor pos")
	endif
	! write sort order
	call initialise(rev_sort_index_io, "quip_rev_sort_index"//trim(topology_suffix), action=OUTPUT)
	call print(rev_sort_index, file=rev_sort_index_io)
	call finalise(rev_sort_index_io)
      endif
    endif

    ! set variables having to do with periodic configs
    if (.not. get_value(at%params, 'Periodic', at_periodic)) at_periodic = .true.
    insert_pos = 0
    if (at_periodic) then
      call insert_cp2k_input_line(cp2k_template_a, " @SET PERIODIC XYZ", after_line = insert_pos, n_l = template_n_lines); insert_pos = insert_pos + 1
    else
      call insert_cp2k_input_line(cp2k_template_a, " @SET PERIODIC NONE", after_line = insert_pos, n_l = template_n_lines); insert_pos = insert_pos + 1
    endif
    call insert_cp2k_input_line(cp2k_template_a, " @SET MAX_CELL_SIZE_INT "//int(max(norm(at%lattice(:,1)),norm(at%lattice(:,2)), norm(at%lattice(:,3)))), &
      after_line = insert_pos, n_l = template_n_lines); insert_pos = insert_pos + 1

    ! put in method
    insert_pos = find_make_cp2k_input_section(cp2k_template_a, template_n_lines, "", "&FORCE_EVAL")
    call insert_cp2k_input_line(cp2k_template_a, "&FORCE_EVAL METHOD "//trim(method), after_line = insert_pos, n_l = template_n_lines)

    ! get qm_list and link_list
    if (use_QMMM) then
      if (use_buffer) then
	call get_hybrid_list(at, qm_list, all_but_term=.true.,int_property="cluster_mark"//trim(qm_name_postfix))
      else
	call get_hybrid_list(at, qm_list, active_trans_only=.true.,int_property="cluster_mark"//trim(qm_name_postfix))
      endif
      allocate(qm_list_a(qm_list%N))
      if (qm_list%N > 0) qm_list_a = int_part(qm_list,1)
      !get link list

       if (assign_pointer(at,'cut_bonds'//trim(qm_name_postfix),cut_bonds_p)) then
          call initialise(cut_bonds,2,0,0,0,0)
          do i_inner=1,at%N
             do j=1,size(cut_bonds_p,1) !MAX_CUT_BONDS
		if (cut_bonds_p(j,i_inner) == 0) exit
		! correct for new atom indices resulting from sorting of atoms
                i_outer = rev_sort_index(cut_bonds_p(j,i_inner))
                call append(cut_bonds,(/i_inner,i_outer/))
             enddo
          enddo
          if (cut_bonds%N > 0) then
             call uniq(cut_bonds%int(2,1:cut_bonds%N),link_list_a)
             allocate(qm_and_link_list_a(size(qm_list_a)+size(link_list_a)))
             qm_and_link_list_a(1:size(qm_list_a)) = qm_list_a(1:size(qm_list_a))
             qm_and_link_list_a(size(qm_list_a)+1:size(qm_list_a)+size(link_list_a)) = link_list_a(1:size(link_list_a))
          else
             allocate(link_list_a(0))
             allocate(qm_and_link_list_a(size(qm_list_a)))
             if (size(qm_list_a) > 0) qm_and_link_list_a = qm_list_a
          endif
       else
          allocate(qm_and_link_list_a(size(qm_list_a)))
          if (size(qm_list_a) > 0) qm_and_link_list_a = qm_list_a
       endif

       !If needed, read QM/MM link_template_file
       if (size(link_list_a) > 0) then
          if (trim(link_template_file).eq."") call system_abort("There are QM/MM links, but qmmm_link_template is not defined.")
          call initialise(link_template_io, trim(link_template_file), INPUT)
          call read_file(link_template_io, link_template_a, link_template_n_lines)
          call finalise(link_template_io)
          call prefix_cp2k_input_sections(link_template_a)
       endif
    else
      allocate(qm_list_a(0))
      allocate(link_list_a(0))
      allocate(qm_and_link_list_a(0))
    endif

    if (auto_centre) then
      if (qm_list%N > 0) then
	centre_pos = pbc_aware_centre(at%pos(:,qm_list_a), at%lattice, at%g)
      else
	centre_pos = pbc_aware_centre(at%pos, at%lattice, at%g)
      endif
      call print("centering got automatic center " // centre_pos, PRINT_VERBOSE)
    endif
    ! move specified centre to origin (centre is already 0 if not specified)
    at%pos(1,:) = at%pos(1,:) - centre_pos(1)
    at%pos(2,:) = at%pos(2,:) - centre_pos(2)
    at%pos(3,:) = at%pos(3,:) - centre_pos(3)
    ! move origin into center of CP2K box (0.5 0.5 0.5 lattice coords)
    call map_into_cell(at)
    if (.not. at_periodic) then
      cp2k_box_centre_pos(1:3) = 0.5_dp*sum(at%lattice,2)
      at%pos(1,:) = at%pos(1,:) + cp2k_box_centre_pos(1)
      at%pos(2,:) = at%pos(2,:) + cp2k_box_centre_pos(2)
      at%pos(3,:) = at%pos(3,:) + cp2k_box_centre_pos(3)
    endif

    if (qm_list%N == at%N) then
      call print("WARNING: requested '"//trim(run_type)//"' but all atoms are in QM region, doing full QM run instead", PRINT_ALWAYS)
      run_type='QS'
      use_QM = .true.
      use_MM = .false.
      use_QMMM = .false.
      method = 'QS'
    endif

    can_reuse_wfn = .true.

    ! put in things needed for QMMM
    if (use_QMMM) then

      insert_pos = find_make_cp2k_input_section(cp2k_template_a, template_n_lines, "&FORCE_EVAL&QMMM", "&CELL")
      call print('INFO: The size of the QM cell is either the MM cell itself, or it will have at least '//(qm_vacuum/2.0_dp)// &
			' Angstrom around the QM atoms.')
      call print('WARNING! Please check if your cell is centreed around the QM region!',PRINT_ALWAYS)
      call print('WARNING! CP2K centreing algorithm fails if QM atoms are not all in the',PRINT_ALWAYS)
      call print('WARNING! 0,0,0 cell. If you have checked it, please ignore this message.',PRINT_ALWAYS)
      cur_qmmm_qm_abc = qmmm_qm_abc(at, qm_list_a, qm_vacuum)
      call insert_cp2k_input_line(cp2k_template_a, "&FORCE_EVAL&QMMM&CELL ABC " // cur_qmmm_qm_abc, after_line=insert_pos, n_l=template_n_lines); insert_pos = insert_pos + 1
      call insert_cp2k_input_line(cp2k_template_a, "&FORCE_EVAL&QMMM&CELL PERIODIC XYZ", after_line=insert_pos, n_l=template_n_lines); insert_pos = insert_pos + 1

      if (get_value(at%params, "QM_cell"//trim(qm_name_postfix), old_qmmm_qm_abc)) then
	if (cur_qmmm_qm_abc .fne. old_qmmm_qm_abc) can_reuse_wfn = .false.
      else
        can_reuse_wfn = .false.
      endif
      call set_value(at%params, "QM_cell"//trim(qm_name_postfix), cur_qmmm_qm_abc)
       call print('set_value QM_cell'//trim(qm_name_postfix)//' '//cur_qmmm_qm_abc)

      !check if QM list changed: compare cluster_mark and old_cluster_mark[_postfix]
!      if (get_value(at%params, "QM_list_changed", qm_list_changed)) then
!       if (qm_list_changed) can_reuse_wfn = .false.
!      endif
       if (.not.has_property(at, 'cluster_mark'//trim(qm_name_postfix))) call system_abort('no cluster_mark'//trim(qm_name_postfix)//' found in atoms object')
       if (.not.has_property(at, 'old_cluster_mark'//trim(qm_name_postfix))) call system_abort('no old_cluster_mark'//trim(qm_name_postfix)//' found in atoms object')
       dummy = assign_pointer(at, 'old_cluster_mark'//trim(qm_name_postfix), old_cluster_mark_p)
       dummy = assign_pointer(at, 'cluster_mark'//trim(qm_name_postfix), cluster_mark_p)

       qm_list_changed = .false.
       do i=1,at%N
          if (old_cluster_mark_p(i) /= cluster_mark_p(i)) then ! mark changed.  Does it matter?
	      if (use_buffer) then ! EXTENDED, check for transitions to/from HYBRID_NO_MARK
		if (any((/old_cluster_mark_p(i),cluster_mark_p(i)/) == HYBRID_NO_MARK)) qm_list_changed = .true.
	      else ! CORE, check for transitions between ACTIVE/TRANS and other
		if ( ( any(old_cluster_mark_p(i)  == (/ HYBRID_ACTIVE_MARK, HYBRID_TRANS_MARK /)) .and. &
		       all(cluster_mark_p(i) /= (/ HYBRID_ACTIVE_MARK, HYBRID_TRANS_MARK /)) ) .or. &
		     ( any(cluster_mark_p(i)  == (/ HYBRID_ACTIVE_MARK, HYBRID_TRANS_MARK /)) .and. &
		       all(old_cluster_mark_p(i) /= (/ HYBRID_ACTIVE_MARK, HYBRID_TRANS_MARK /)) ) ) qm_list_changed = .true.
              endif
	      if (qm_list_changed) exit
          endif
       enddo
       call set_value(at%params,'QM_list_changed',qm_list_changed)
       call print('set_value QM_list_changed '//qm_list_changed)

       if (qm_list_changed) can_reuse_wfn = .false.

      !Add QM atoms
      counter = 0
      do atno=minval(at%Z), maxval(at%Z)
	if (any(at%Z(qm_list_a) == atno)) then
	  insert_pos = find_make_cp2k_input_section(cp2k_template_a, template_n_lines, "&FORCE_EVAL&QMMM", "&QM_KIND "//ElementName(atno))
	  do i=1, size(qm_list_a)
	    if (at%Z(qm_list_a(i)) == atno) then
	      call insert_cp2k_input_line(cp2k_template_a, "&FORCE_EVAL&QMMM&QM_KIND-"//trim(ElementName(atno))// &
							" MM_INDEX "//qm_list_a(i), after_line = insert_pos, n_l = template_n_lines)
	      insert_pos = insert_pos + 1
	      counter = counter + 1
	    endif
	  end do
	end if
      end do
      if (size(qm_list_a) /= counter) &
	call system_abort("Number of QM list atoms " // size(qm_list_a) // " doesn't match number of QM_KIND atoms " // counter)

      !Add link sections from template file for each link
      if (size(link_list_a).gt.0) then
         do i=1,cut_bonds%N
            i_inner = cut_bonds%int(1,i)
            i_outer = cut_bonds%int(2,i)
            insert_pos = find_make_cp2k_input_section(cp2k_template_a, template_n_lines, "&FORCE_EVAL", "&QMMM")
            inserted_atoms = .false.
            do i_line=1,link_template_n_lines
               call insert_cp2k_input_line(cp2k_template_a, trim("&FORCE_EVAL&QMMM")//trim(link_template_a(i_line)), after_line = insert_pos, n_l = template_n_lines); insert_pos = insert_pos + 1
               if (.not.inserted_atoms) then
                  call insert_cp2k_input_line(cp2k_template_a, "&FORCE_EVAL&QMMM&LINK MM_INDEX "//i_outer, after_line = insert_pos, n_l = template_n_lines); insert_pos = insert_pos + 1
                  call insert_cp2k_input_line(cp2k_template_a, "&FORCE_EVAL&QMMM&LINK QM_INDEX "//i_inner, after_line = insert_pos, n_l = template_n_lines); insert_pos = insert_pos + 1
                  inserted_atoms = .true.
               endif
            enddo
         enddo
      endif
    endif

    ! put in things needed for QM
    if (use_QM) then
      if (try_reuse_wfn .and. can_reuse_wfn) then 
	insert_pos = find_make_cp2k_input_section(cp2k_template_a, template_n_lines, "&FORCE_EVAL", "&DFT")
        call insert_cp2k_input_line(cp2k_template_a, "&FORCE_EVAL&DFT WFN_RESTART_FILE_NAME ../wfn.restart.wfn"//trim(qm_name_postfix), after_line = insert_pos, n_l = template_n_lines); insert_pos = insert_pos + 1
	!insert_pos = find_make_cp2k_input_section(cp2k_template_a, template_n_lines, "&FORCE_EVAL&DFT", "&SCF")
	!call insert_cp2k_input_line(cp2k_template_a, "&FORCE_EVAL&DFT&SCF SCF_GUESS RESTART", after_line = insert_pos, n_l = template_n_lines); insert_pos = insert_pos + 1
      endif
      call calc_charge_lsd(at, qm_list_a, charge, do_lsd, error=error)
      PASS_ERROR(error)
      insert_pos = find_make_cp2k_input_section(cp2k_template_a, template_n_lines, "&FORCE_EVAL", "&DFT")
      call insert_cp2k_input_line(cp2k_template_a, "&FORCE_EVAL&DFT CHARGE "//charge, after_line = insert_pos, n_l = template_n_lines); insert_pos = insert_pos + 1
      if (do_lsd) call insert_cp2k_input_line(cp2k_template_a, "&FORCE_EVAL&DFT LSD ", after_line = insert_pos, n_l = template_n_lines); insert_pos = insert_pos + 1
      if (len_trim(calc_qm_charges) > 0) then
	 insert_pos = find_make_cp2k_input_section(cp2k_template_a, template_n_lines, "&FORCE_EVAL&DFT", "&PRINT")
	 insert_pos = find_make_cp2k_input_section(cp2k_template_a, template_n_lines, "&FORCE_EVAL&DFT&PRINT", "&MULLIKEN")
	 call insert_cp2k_input_line(cp2k_template_a, "&FORCE_EVAL&DFT&PRINT FILENAME qmcharges", after_line = insert_pos, n_l = template_n_lines); insert_pos = insert_pos + 1
      endif
    endif

    ! put in unit cell
    insert_pos = find_make_cp2k_input_section(cp2k_template_a, template_n_lines, "&FORCE_EVAL", "&SUBSYS")

    insert_pos = find_make_cp2k_input_section(cp2k_template_a, template_n_lines, "&FORCE_EVAL&SUBSYS", "&CELL")
    call insert_cp2k_input_line(cp2k_template_a, "&FORCE_EVAL&SUBSYS&CELL A " // at%lattice(:,1), after_line = insert_pos, n_l = template_n_lines); insert_pos = insert_pos + 1
    call insert_cp2k_input_line(cp2k_template_a, "&FORCE_EVAL&SUBSYS&CELL B " // at%lattice(:,2), after_line = insert_pos, n_l = template_n_lines); insert_pos = insert_pos + 1
    call insert_cp2k_input_line(cp2k_template_a, "&FORCE_EVAL&SUBSYS&CELL C " // at%lattice(:,3), after_line = insert_pos, n_l = template_n_lines); insert_pos = insert_pos + 1

    ! put in topology
    insert_pos = find_make_cp2k_input_section(cp2k_template_a, template_n_lines, "&FORCE_EVAL&SUBSYS", "&TOPOLOGY")
    insert_pos = find_make_cp2k_input_section(cp2k_template_a, template_n_lines, "&FORCE_EVAL&SUBSYS&TOPOLOGY", "&DUMP_PSF")
    insert_pos = find_make_cp2k_input_section(cp2k_template_a, template_n_lines, "&FORCE_EVAL&SUBSYS&TOPOLOGY", "&DUMP_PDB")
    insert_pos = find_make_cp2k_input_section(cp2k_template_a, template_n_lines, "&FORCE_EVAL&SUBSYS&TOPOLOGY", "&GENERATE")
    call insert_cp2k_input_line(cp2k_template_a, "&FORCE_EVAL&SUBSYS&TOPOLOGY&GENERATE REORDER F", after_line = insert_pos, n_l = template_n_lines); insert_pos = insert_pos + 1
    call insert_cp2k_input_line(cp2k_template_a, "&FORCE_EVAL&SUBSYS&TOPOLOGY&GENERATE CREATE_MOLECULES F", after_line = insert_pos, n_l = template_n_lines); insert_pos = insert_pos + 1
    insert_pos = find_make_cp2k_input_section(cp2k_template_a, template_n_lines, "&FORCE_EVAL&SUBSYS&TOPOLOGY&GENERATE", "&ISOLATED_ATOMS")
    if (use_QMMM) then
      do i=1, size(qm_list_a)
	call insert_cp2k_input_line(cp2k_template_a, "&FORCE_EVAL&SUBSYS&TOPOLOGY&GENERATE&ISOLATED_ATOMS LIST " // qm_list_a(i), after_line = insert_pos, n_l = template_n_lines); insert_pos = insert_pos + 1
      end do
    endif
    if (assign_pointer(at, "isolated_atom", isolated_atom)) then
      do i=1, at%N
	if (isolated_atom(i) /= 0) then
	  call insert_cp2k_input_line(cp2k_template_a, "&FORCE_EVAL&SUBSYS&TOPOLOGY&GENERATE&ISOLATED_ATOMS LIST " // i, after_line = insert_pos, n_l = template_n_lines); insert_pos = insert_pos + 1
	endif
      end do
    endif

    insert_pos = find_make_cp2k_input_section(cp2k_template_a, template_n_lines, "&FORCE_EVAL&SUBSYS", "&TOPOLOGY")
    call insert_cp2k_input_line(cp2k_template_a, "&FORCE_EVAL&SUBSYS&TOPOLOGY COORD_FILE_NAME quip_cp2k.xyz", after_line = insert_pos, n_l = template_n_lines); insert_pos = insert_pos + 1
    call insert_cp2k_input_line(cp2k_template_a, "&FORCE_EVAL&SUBSYS&TOPOLOGY COORDINATE XYZ", after_line = insert_pos, n_l = template_n_lines); insert_pos = insert_pos + 1
    if (trim(psf_print) == "DRIVER_PRINT_AND_SAVE" .or. trim(psf_print) == "USE_EXISTING_PSF") then
      if (tmp_run_dir_i>0) then
        call system_command("if [ ! -s "//trim(tmp_run_dir)//"/quip_cp2k"//trim(topology_suffix)//".psf ] ; then cp quip_cp2k"//trim(topology_suffix)//".psf /tmp/cp2k_run_"//tmp_run_dir_i//"/ ; fi",status=stat)
        if ( stat /= 0 ) call system_abort("Something went wrong when tried to copy quip_cp2k"//trim(topology_suffix)//".psf into the tmp dir "//trim(tmp_run_dir))
        call insert_cp2k_input_line(cp2k_template_a, "&FORCE_EVAL&SUBSYS&TOPOLOGY CONN_FILE_NAME quip_cp2k"//trim(topology_suffix)//".psf", after_line = insert_pos, n_l = template_n_lines); insert_pos = insert_pos + 1
      else
        call insert_cp2k_input_line(cp2k_template_a, "&FORCE_EVAL&SUBSYS&TOPOLOGY CONN_FILE_NAME ../quip_cp2k"//trim(topology_suffix)//".psf", after_line = insert_pos, n_l = template_n_lines); insert_pos = insert_pos + 1
      endif
      call insert_cp2k_input_line(cp2k_template_a, "&FORCE_EVAL&SUBSYS&TOPOLOGY CONN_FILE_FORMAT PSF", after_line = insert_pos, n_l = template_n_lines); insert_pos = insert_pos + 1
    endif

    ! put in global stuff to run a single force evalution, print out appropriate things
    insert_pos = find_make_cp2k_input_section(cp2k_template_a, template_n_lines, "", "&GLOBAL")
    call insert_cp2k_input_line(cp2k_template_a, "&GLOBAL   PROJECT quip", after_line = insert_pos, n_l = template_n_lines); insert_pos = insert_pos + 1
    call insert_cp2k_input_line(cp2k_template_a, "&GLOBAL   RUN_TYPE MD", after_line = insert_pos, n_l = template_n_lines); insert_pos = insert_pos + 1

    insert_pos = find_make_cp2k_input_section(cp2k_template_a, template_n_lines, "", "&MOTION")
    insert_pos = find_make_cp2k_input_section(cp2k_template_a, template_n_lines, "&MOTION", "&PRINT")
    insert_pos = find_make_cp2k_input_section(cp2k_template_a, template_n_lines, "&MOTION&PRINT", "&FORCES")
    call insert_cp2k_input_line(cp2k_template_a, "&MOTION&PRINT&FORCES       FORMAT XMOL", after_line = insert_pos, n_l = template_n_lines); insert_pos = insert_pos + 1
    insert_pos = find_make_cp2k_input_section(cp2k_template_a, template_n_lines, "&MOTION", "&MD")
    call insert_cp2k_input_line(cp2k_template_a, "&MOTION&MD     ENSEMBLE NVE", after_line = insert_pos, n_l = template_n_lines); insert_pos = insert_pos + 1
    call insert_cp2k_input_line(cp2k_template_a, "&MOTION&MD     STEPS 0", after_line = insert_pos, n_l = template_n_lines); insert_pos = insert_pos + 1

    call write_cp2k_input_file(cp2k_template_a(1:template_n_lines), trim(run_dir)//'/cp2k_input.inp')

    ! prepare xyz file for input to cp2k
    call write(at, trim(run_dir)//'/quip_cp2k.xyz', properties='species:pos')
    ! actually run cp2k
    call run_cp2k_program(trim(cp2k_program), trim(run_dir), max_n_tries)

    ! parse output
    call read_output(at, qm_and_link_list_a, cur_qmmm_qm_abc, trim(run_dir), "quip", e, f, trim(calc_qm_charges), error=error)

    at%pos(1,:) = at%pos(1,:) + centre_pos(1) - cp2k_box_centre_pos(1)
    at%pos(2,:) = at%pos(2,:) + centre_pos(2) - cp2k_box_centre_pos(2)
    at%pos(3,:) = at%pos(3,:) + centre_pos(3) - cp2k_box_centre_pos(3)
    call map_into_cell(at)

    ! unsort
    if (associated(sort_index_p)) then
      f(:,sort_index_p(:)) = f(:,:)
      call atoms_sort(at, 'sort_index', error=error)
      PASS_ERROR_WITH_INFO("do_cp2k_calc sorting atoms by sort_index",error)
    endif

    if (maxval(abs(f)) > max_force_warning) &
      call print('WARNING cp2k forces max component ' // maxval(abs(f)) // ' at ' // maxloc(abs(f)) // &
		 ' exceeds warning threshold ' // max_force_warning, PRINT_ALWAYS)

    ! save output

    if (use_QM) then
      call system_command('cp '//trim(run_dir)//'/quip-RESTART.wfn wfn.restart.wfn'//trim(qm_name_postfix))
      if (save_output_wfn_files) then
	call system_command('cp '//trim(run_dir)//'/quip-RESTART.wfn run_'//run_dir_i//'_end.wfn.restart.wfn'//trim(qm_name_postfix))
      endif
    endif

    if (save_output_files) then
      call system_command(&
        ' cat '//trim(run_dir)//'/cp2k_input.inp >> cp2k_input_log; echo "##############" >> cp2k_input_log;' // &
        ' cat '//trim(run_dir)//'/cp2k_output.out >> cp2k_output_log; echo "##############" >> cp2k_output_log;' // &
        ' cat filepot.0.xyz'//' >> cp2k_filepot_in_log.xyz;' // &
        ' cat '//trim(run_dir)//'/quip-frc-1.xyz'// ' >> cp2k_force_file_log')
    endif

    ! clean up

    if (tmp_run_dir_i>0) then
      if (clean_up_files) then
         !only delete files that need recreating, keep basis, potentials, psf
         call system_command('rm -f '//trim(tmp_run_dir)//"/quip-* "//trim(tmp_run_dir)//"/cp2k_input.inp "//trim(tmp_run_dir)//"/cp2k_output.out "//trim(run_dir))
      else !save dir
         exists = .true.
         i = 0
         do while (exists)
           i = i + 1
           dir = "cp2k_run_saved_"//i
           call system_command("bash -c '[ -e "//trim(dir)//" ]'", status=stat)
           exists = (stat == 0)
         end do
         call system_command("cp -r "//trim(run_dir)//" "//trim(dir), status=stat)
         if (stat /= 0) then
            RAISE_ERROR("Failed to copy "//trim(run_dir)//" to "//trim(dir)//" status " // stat, error)
         endif
         call system_command('rm -f '//trim(tmp_run_dir)//"/* "//trim(run_dir))
      endif
    else
       if (clean_up_files) then
	 if (clean_up_keep_n <= 0) then ! never keep any old directories around
	    call system_command('rm -rf '//trim(run_dir))
	 else ! keep some (>= 1) old directories around
	    delete_dir_i = mod(run_dir_i, clean_up_keep_n+1)+1
	    call system_command('rm -rf cp2k_run_'//delete_dir_i)
	 endif
       endif
    endif

    if (allocated(rev_sort_index)) deallocate(rev_sort_index)

    call system_timer('do_cp2k_calc')

  end subroutine do_cp2k_calc

  function find_make_cp2k_input_section(l_a, n_l, base_sec, new_sec) result(line_n)
    character(len=*), allocatable, intent(inout) :: l_a(:)
    integer, intent(inout) :: n_l
    character(len=*), intent(in) :: base_sec, new_sec
    integer :: line_n

    integer :: i, pamp, pspc
    character(len=FIELD_LENGTH) :: sec, word, arg, base_sec_root, base_sec_tail, new_sec_end

    line_n = 0

    do i=1, n_l
      call split_cp2k_input_line(trim(l_a(i)), sec, word, arg)

      if (trim(sec) == trim(base_sec) .and. trim(word) == trim(new_sec)) then
	line_n = i
	return
      endif
    end do

    if (len_trim(base_sec) == 0) then
      i = n_l
      call insert_cp2k_input_line(l_a, " "//trim(new_sec), after_line=i, n_l=n_l); i = i + 1
      pspc = index(trim(new_sec)," ")
      if (pspc == 0) then
	new_sec_end = "&END " // new_sec(2:len_trim(new_sec))
      else
	new_sec_end = "&END " // new_sec(2:pspc-1)
      endif
      call insert_cp2k_input_line(l_a, " " // trim(new_sec_end), after_line=i, n_l=n_l); i = i + 1
      line_n = n_l-1
      return
    endif

    pamp = index(base_sec,"&", back=.true.)
    if (pamp <= 1) then
      base_sec_root = ""
      base_sec_tail = trim(base_sec)
    else
      base_sec_root = base_sec(1:pamp-1)
      base_sec_tail = base_sec(pamp:len_trim(base_sec))
    endif


    do i=1, n_l
      call split_cp2k_input_line(trim(l_a(i)), sec, word, arg)
      if (trim(sec) == trim(base_sec_root) .and. trim(word) == trim(base_sec_tail)) then
	call insert_cp2k_input_line(l_a, trim(base_sec)//" "//trim(new_sec), after_line=i, n_l=n_l)
	pspc = index(trim(new_sec)," ")
	if (pspc == 0) then
	  new_sec_end = "&END " // new_sec(2:len_trim(new_sec))
	else
	  new_sec_end = "&END " // new_sec(2:pspc-1)
	endif
	call insert_cp2k_input_line(l_a, trim(base_sec)//" " // trim(new_sec_end), after_line=i+1, n_l=n_l)
	line_n = i+1
	return
      endif
    end do

    if (line_n == 0) &
      call system_abort("Could not find or make section '"//trim(new_sec)//" in base section '"//trim(base_sec))

  end function find_make_cp2k_input_section

  subroutine read_output(at, qm_list_a, cur_qmmm_qm_abc, run_dir, proj, e, f, calc_qm_charges, error)
    type(Atoms), intent(inout) :: at
    integer, intent(in) :: qm_list_a(:)
    real(dp), intent(in) :: cur_qmmm_qm_abc(3)
    character(len=*), intent(in) :: run_dir, proj
    real(dp), intent(out) :: e, f(:,:)
    real(dp), pointer :: force_p(:,:)
    character(len=*) :: calc_qm_charges
    integer, intent(out), optional :: ERROR

    real(dp), pointer :: qm_charges_p(:)
    type(Atoms) :: f_xyz, p_xyz
    integer :: m
    integer :: i, ti
    type(inoutput) :: qm_charges_io
    character(len=FIELD_LENGTH) :: species, qm_charges_l

    INIT_ERROR(error)

    call read(f_xyz, trim(run_dir)//'/'//trim(proj)//'-frc-1.xyz')
    call read(p_xyz, trim(run_dir)//'/'//trim(proj)//'-pos-1.xyz')
    nullify(qm_charges_p)
    if (len_trim(calc_qm_charges) > 0) then
      if (.not. assign_pointer(at, trim(calc_qm_charges), qm_charges_p)) then
	  call add_property(at, trim(calc_qm_charges), 0.0_dp, ptr=qm_charges_p)
      endif
      call initialise(qm_charges_io, trim(run_dir)//'/'//trim(proj)//'-qmcharges--1.mulliken',action=INPUT, error=error)
      PASS_ERROR_WITH_INFO("cp2k_driver read_output() failed to open qmcharges file", error)
      do i=1, 3
	qm_charges_l = read_line(qm_charges_io)
      end do
      do i=1, at%N
	qm_charges_l = read_line(qm_charges_io)
	read (unit=qm_charges_l,fmt=*) ti, species, qm_charges_p(i)
      end do
      call finalise(qm_charges_io)
    endif

    if (.not. get_value(f_xyz%params, "E", e)) then
      RAISE_ERROR('read_output failed to find E value in '//trim(run_dir)//'/quip-frc-1.xyz file', error)
    endif

    if (.not.(assign_pointer(f_xyz, 'frc', force_p))) then
      RAISE_ERROR("Did not find frc property in "//trim(run_dir)//'/quip-frc-1.xyz file', error)
    endif
    f = force_p

    e = e * HARTREE
    f  = f * HARTREE/BOHR 
    call reorder_if_necessary(at, qm_list_a, cur_qmmm_qm_abc, p_xyz%pos, f, qm_charges_p)

    call print('')
    call print('The energy of the system: '//e)
    call verbosity_push_decrement()
      call print('The forces acting on each atom (eV/A):')
      call print('atom     F(x)     F(y)     F(z)')
      do m=1,size(f,2)
        call print('  '//m//'    '//f(1,m)//'  '//f(2,m)//'  '//f(3,m))
      enddo
    call verbosity_pop()
    call print('Sum of the forces: '//sum(f,2))

  end subroutine read_output

  subroutine reorder_if_necessary(at, qm_list_a, qmmm_qm_abc, new_p, new_f, qm_charges_p)
    type(Atoms), intent(in) :: at
    integer, intent(in) :: qm_list_a(:)
    real(dp), intent(in) :: qmmm_qm_abc(3)
    real(dp), intent(in) :: new_p(:,:)
    real(dp), intent(inout) :: new_f(:,:)
    real(dp), intent(inout), pointer :: qm_charges_p(:)

    real(dp) :: shift(3)
    integer, allocatable :: reordering_index(:)
    integer :: i, j

    ! shifted cell in case of QMMM (cp2k/src/toplogy_coordinate_util.F)
    shift = 0.0_dp
    if (size(qm_list_a) > 0) then
      do i=1,3
	shift(i) = 0.5_dp * qmmm_qm_abc(i) - (minval(at%pos(i,qm_list_a)) + maxval(at%pos(i,qm_list_a)))*0.5_dp
      end do
    endif
    allocate(reordering_index(at%N))
    call check_reordering(at%pos, shift, new_p, at%g, reordering_index)
    if (any(reordering_index == 0)) then
      ! try again with shift of a/2 b/2 c/2 in case TOPOLOGY%CENTER_COORDINATES is set
      shift = sum(at%lattice(:,:),2)/2.0_dp - &
	      (minval(at%pos(:,:),2)+maxval(at%pos(:,:),2))/2.0_dp
      call check_reordering(at%pos, shift, new_p, at%g, reordering_index)
      if (any(reordering_index == 0)) then
	! try again with uniform shift (module periodic cell)
	shift = new_p(:,1) - at%pos(:,1)
	call check_reordering(at%pos, shift, new_p, at%g, reordering_index)
	if (any(reordering_index == 0)) &
	  call system_abort("Could not match original and read in atom objects")
      endif
    endif

    do i=1, at%N
      if (reordering_index(i) /= i) then
	 call print("WARNING: reorder_if_necessary indeed found reordered atoms", PRINT_ALWAYS)
	 exit
      endif
    end do

    new_f(1,reordering_index(:)) = new_f(1,:)
    new_f(2,reordering_index(:)) = new_f(2,:)
    new_f(3,reordering_index(:)) = new_f(3,:)
    if (associated(qm_charges_p)) then
      qm_charges_p(reordering_index(:)) = qm_charges_p(:)
    endif

    deallocate(reordering_index)
  end subroutine reorder_if_necessary

  subroutine check_reordering(old_p, shift, new_p, recip_lattice, reordering_index)
    real(dp), intent(in) :: old_p(:,:), shift(3), new_p(:,:), recip_lattice(3,3)
    integer, intent(out) :: reordering_index(:)

    integer :: N, i, j
    real(dp) :: dpos(3), dpos_i(3)

    N = size(old_p,2)

    reordering_index = 0
    do i=1, N
      ! check for same-order
      j = i
      dpos = matmul(recip_lattice(1:3,1:3), old_p(1:3,i) + shift(1:3) - new_p(1:3,j))
      dpos_i = nint(dpos)
      if (all(abs(dpos-dpos_i) <= 1.0e-4_dp)) then
        reordering_index(i) = j
      else ! not same order, search
	 do j=1, N
	   dpos = matmul(recip_lattice(1:3,1:3), old_p(1:3,i) + shift(1:3) - new_p(1:3,j))
	   dpos_i = nint(dpos)
	   if (all(abs(dpos-dpos_i) <= 1.0e-4_dp)) then
	     reordering_index(i) = j
	     exit
	   endif
	 end do
       end if
    end do
  end subroutine check_reordering

  subroutine run_cp2k_program(cp2k_program, run_dir, max_n_tries)
    character(len=*), intent(in) :: cp2k_program, run_dir
    integer, intent(in) :: max_n_tries

    integer :: n_tries
    logical :: converged
    character(len=FIELD_LENGTH) :: cp2k_run_command
    integer :: stat, error_stat

    n_tries = 0
    converged = .false.

    do while (.not. converged .and. (n_tries < max_n_tries))
      n_tries = n_tries + 1

      cp2k_run_command = 'cd ' // trim(run_dir)//'; '//trim(cp2k_program)//' cp2k_input.inp >> cp2k_output.out'
      call print("Doing '"//trim(cp2k_run_command)//"'")
      call system_timer('cp2k_run_command')
      call system_command(trim(cp2k_run_command), status=stat)
      call system_timer('cp2k_run_command')
      call print('grep -i warning '//trim(run_dir)//'/cp2k_output.out', PRINT_ALWAYS)
      call system_command("fgrep -i 'warning' "//trim(run_dir)//"/cp2k_output.out")
      call system_command("fgrep -i 'error' "//trim(run_dir)//"/cp2k_output.out", status=error_stat)
      if (stat /= 0) &
	call system_abort('cp2k_run_command has non zero return status ' // stat //'. check output file '//trim(run_dir)//'/cp2k_output.out')
      if (error_stat == 0) &
	call system_abort('cp2k_run_command generated ERROR message in output file '//trim(run_dir)//'/cp2k_output.out')

      call system_command('egrep "FORCE_EVAL.* QS " '//trim(run_dir)//'/cp2k_output.out',status=stat)
      if (stat == 0) then ! QS or QMMM run
	call system_command('grep "FAILED to converge" '//trim(run_dir)//'/cp2k_output.out',status=stat)
	if (stat == 0) then
	  call print("WARNING: cp2k_driver failed to converge, trying again",PRINT_ALWAYS)
	  converged = .false.
	else
	  call system_command('grep "outer SCF loop converged" '//trim(run_dir)//'/cp2k_output.out',status=stat)
	  if (stat == 0) then
	    converged = .true.
	  else
	    call print("WARNING: cp2k_driver couldn't find definitive sign of convergence or failure to converge in output file, trying again",PRINT_ALWAYS)
	    converged = .false.
	  endif
	end if
      else ! MM run
	converged = .true.
      endif
    end do

    if (.not. converged) &
      call system_abort('cp2k failed to converge after n_tries='//n_tries//'. see output file '//trim(run_dir)//'/cp2k_output.out')

  end subroutine run_cp2k_program

  subroutine write_cp2k_input_file(l_a, filename)
    character(len=*), intent(in) :: l_a(:)
    character(len=*), intent(in) :: filename

    integer :: i, pspc
    type(Inoutput) :: io

    call initialise(io, trim(filename), OUTPUT)
    do i=1, size(l_a)
      pspc = index(l_a(i), " ")
      call print(l_a(i)(pspc+1:len_trim(l_a(i))), file=io)
    end do
    call finalise(io)

  end subroutine write_cp2k_input_file

  subroutine split_cp2k_input_line(l, sec, word, arg)
    character(len=*) :: l, sec, word, arg

    character(len=len(l)) :: t_l
    integer :: pspc

    t_l = l

    pspc = index(trim(t_l), " ")
    if (pspc == 0 .or. pspc == 1) then ! top level
      sec = ""
    else
      sec = t_l(1:pspc-1)
    endif

    t_l = adjustl(t_l(pspc+1:len_trim(t_l)))

    pspc = index(trim(t_l), " ")
    if (pspc == 0) then ! no arg
      word = trim(t_l)
      arg = ""
    else
      word = t_l(1:pspc-1)
      arg = t_l(pspc+1:len_trim(t_l))
    endif

    sec = adjustl(sec)
    word = adjustl(word)
    arg = adjustl(arg)

  end subroutine split_cp2k_input_line

  subroutine substitute_cp2k_input(l_a, s_source, s_targ, n_l)
    character(len=*), intent(inout) :: l_a(:)
    character(len=*), intent(in) :: s_source, s_targ
    integer, intent(in) :: n_l

    character(len=len(l_a(1))) :: t
    integer :: i, p, l_targ, l_source, len_before, len_after

    l_targ = len(s_targ)
    l_source = len(s_source)

    do i=1, n_l
      p = index(l_a(i), s_source)
      if (p > 0) then
	len_before = p-1
	len_after = len_trim(l_a(i)) - l_source - p + 1
	t = ""
	if (len_before > 0) then
	  t(1:len_before) = l_a(i)(1:len_before)
	endif
	t(len_before+1:len_before+1+l_targ-1) = s_targ(1:l_targ)
	if (len_after > 0) then
	  t(len_before+1+l_targ:len_before+1+l_targ+len_after-1) = l_a(i)(len_before+l_source+1:len_before+l_source+1+len_after-1)
	endif
	l_a(i) = t
      end if
    end do
  end subroutine substitute_cp2k_input

  subroutine insert_cp2k_input_line(l_a, l, after_line, n_l)
    character(len=*), allocatable, intent(inout) :: l_a(:)
    character(len=*), intent(in) :: l
    integer, intent(in) :: after_line
    integer, intent(inout) :: n_l

    integer :: i

    if (n_l+1 > size(l_a)) call extend_char_array(l_a)

    do i=n_l, after_line+1, -1
      l_a(i+1) = l_a(i)
    end do
    l_a(after_line+1) = trim(l)
    n_l = n_l + 1
  end subroutine insert_cp2k_input_line

  subroutine prefix_cp2k_input_sections(l_a)
    character(len=*), intent(inout) :: l_a(:)

    integer :: pamp
    character(len=CP2K_LINE_LENGTH) :: section_str, new_section_str
    integer :: i, j, comment_i

    section_str = ""
    new_section_str = ""
    do i=1, size(l_a)
      comment_i = index(l_a(i), "#")
      if (comment_i /= 0) then
	l_a(i)(comment_i:len(l_a(i))) = ""
      endif
      if (index(l_a(i),"&END") /= 0) then
	pamp = index(section_str, "&", .true.)
	new_section_str = section_str(1:pamp-1)
	section_str = new_section_str
      else if (index(l_a(i),"&") /= 0) then
	pamp = index(l_a(i), "&")
	new_section_str = trim(section_str)
	do j=pamp, len_trim(l_a(i))
	  if (l_a(i)(j:j) == " ") then
	    new_section_str = trim(new_section_str) // "-"
	  else
	    new_section_str = trim(new_section_str) // l_a(i)(j:j)
	  endif
	end do
      endif
      l_a(i) = (trim(section_str) //" "//trim(l_a(i)))
      section_str = new_section_str
    end do
  end subroutine

  function qmmm_qm_abc(at, qm_list_a, qm_vacuum)
    type(Atoms), intent(in) :: at
    integer, intent(in) :: qm_list_a(:)
    real(dp), intent(in) :: qm_vacuum
    real(dp) :: qmmm_qm_abc(3)

    real(dp) :: qm_maxdist(3)
    integer i, j

    qm_maxdist = 0.0_dp
    do i=1, size(qm_list_a)
    do j=1, size(qm_list_a)
      qm_maxdist(1) = max(qm_maxdist(1), at%pos(1,qm_list_a(i))-at%pos(1,qm_list_a(j)))
      qm_maxdist(2) = max(qm_maxdist(2), at%pos(2,qm_list_a(i))-at%pos(2,qm_list_a(j)))
      qm_maxdist(3) = max(qm_maxdist(3), at%pos(3,qm_list_a(i))-at%pos(3,qm_list_a(j)))
    end do
    end do

    qmmm_qm_abc(1) = min(real(ceiling(qm_maxdist(1)))+qm_vacuum,at%lattice(1,1))
    qmmm_qm_abc(2) = min(real(ceiling(qm_maxdist(2)))+qm_vacuum,at%lattice(2,2))
    qmmm_qm_abc(3) = min(real(ceiling(qm_maxdist(3)))+qm_vacuum,at%lattice(3,3))

  end function

  subroutine calc_charge_lsd(at, qm_list_a, charge, do_lsd, error)
    type(Atoms), intent(in) :: at
    integer, intent(in) :: qm_list_a(:)
    integer, intent(out) :: charge
    logical, intent(out) :: do_lsd
    integer, intent(out), optional :: error

    real(dp), pointer :: atom_charge(:)
    integer, pointer  :: Z_p(:)
    integer           :: sum_Z
    integer           :: l_error

    INIT_ERROR(error)

    if (.not. assign_pointer(at, "Z", Z_p)) then
	RAISE_ERROR("calc_charge_lsd could not find Z property", error)
    endif

    if (size(qm_list_a) > 0) then
      if (.not. assign_pointer(at, "atom_charge", atom_charge)) then
	RAISE_ERROR("calc_charge_lsd could not find atom_charge", error)
      endif
      charge = nint(sum(atom_charge(qm_list_a)))
      !check if we have an odd number of electrons
      sum_Z = sum(Z_p(qm_list_a(1:size(qm_list_a))))
      do_lsd = (mod(sum_Z-charge,2) /= 0)
    else
      sum_Z = sum(Z_p)
      do_lsd = .false.
      charge = 0 
      call get_param_value(at, 'LSD', do_lsd, error=l_error) ! ignore error
      CLEAR_ERROR(error)
      !if charge is saved, also check if we have an odd number of electrons
      call get_param_value(at, 'Charge', charge, error=l_error)
      CLEAR_ERROR(error)
      if (l_error == 0) then
        call print("Using Charge " // charge)
        do_lsd = do_lsd .or. (mod(sum_Z-charge,2) /= 0)
      else !charge=0 is assumed by CP2K
        do_lsd = do_lsd .or. (mod(sum_Z,2) /= 0)
      endif
      if (do_lsd) call print("Using do_lsd " // do_lsd)
    endif

  end subroutine


end module cp2k_driver_template_module
