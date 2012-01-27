program convert

  use libatoms_module
  use gp_sparse_module
  use teach_sparse_mod

  implicit none

  type(teach_sparse) :: GAP

  integer :: i
  type(Dictionary) :: params, my_dictionary
  character(len=STRING_LENGTH) :: binary_file, xml_file

  logical :: has_xml_file
  integer :: quip_string_start, bracket_start, bracket_end, other_bracket_start
  character(len=10000) :: short_comment
  character(len=256) :: coordinates
  real(dp), dimension(:), allocatable :: z_eff, w_Z
  
  call system_initialise(verbosity=PRINT_SILENT)

  call initialise(params)
  call param_register(params, 'binary_file', 'gp.dat', binary_file, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'xml_file', 'gp.xml', xml_file, has_xml_file, help_string="No help yet.  This source file was $LastChangedBy$")

  if (.not. param_read_args(params) ) then
     call verbosity_push(PRINT_NORMAL)
     call print("Usage: covert_binary2xml binary_file=gp.dat xml_file=gp.xml")
     call system_abort('Exiting.')
  endif

  call finalise(params)

  call gp_read_binary(GAP%my_gp,trim(binary_file))

  quip_string_start = index(GAP%my_gp%comment,'quip_string')
  if( quip_string_start > 0 ) then
      bracket_start = index(GAP%my_gp%comment(quip_string_start:),'{')
      bracket_end = index(GAP%my_gp%comment(quip_string_start:),'}')
      other_bracket_start = bracket_start
      do
         other_bracket_start = index(GAP%my_gp%comment(quip_string_start+other_bracket_start:quip_string_start+bracket_end-2),'{')
         if(other_bracket_start > 0) then
            bracket_end = index(GAP%my_gp%comment(quip_string_start+other_bracket_start:),'}')
         else
            exit
         endif
      enddo

      call initialise(GAP%quip_string)
      call concat(GAP%quip_string,GAP%my_gp%comment(quip_string_start+bracket_start:quip_string_start+bracket_end-2))
      short_comment = GAP%my_gp%comment(:quip_string_start-1) // ' ' // GAP%my_gp%comment(quip_string_start+bracket_end:)
  else
     short_comment = GAP%my_gp%comment
  endif
                                         
  call read_string(my_dictionary,short_comment)

  coordinates = ''
  if (.not. get_value(my_dictionary, 'coordinates', coordinates)) &
     & call system_abort('convert: datafile_coordinates not found')

  if (trim(coordinates) == 'bispectrum') then
     if( .not. ( get_value(my_dictionary,'cutoff',GAP%r_cut) .and. &
               & get_value(my_dictionary,'j_max',GAP%j_max) .and. &
               & get_value(my_dictionary,'z0',GAP%z0) ) ) &
     & call system_abort('covert: did not find bispectrum parameters in gp data file, &
     & might be old version or format not correct')
!     GAP%do_qw_so3 = .false.
  elseif (trim(coordinates) == 'qw') then
     if (.not. get_value(my_dictionary, 'l_max', GAP%qw_l_max) .and. &
               get_value(my_dictionary, 'f_n', GAP%qw_f_n)) &
        call system_abort('convert: did not find qw parameters in gp data file')

     do i = 1, GAP%qw_f_n
        if (.not. (get_value(my_dictionary, 'cutoff_' // i, GAP%qw_cutoff(i)) .and. &
                   get_value(my_dictionary, 'cutoff_f_' // i, GAP%qw_cutoff_f(i)) .and. &
                   get_value(my_dictionary, 'cutoff_r1_' // i, GAP%qw_cutoff_r1(i)))) &
        call system_abort('convert: did not find qw parameters in gp data file')
     enddo

     if (.not. get_value(my_dictionary, 'do_q', GAP%qw_no_q)) GAP%qw_no_q = .true.
     if (.not. get_value(my_dictionary, 'do_w', GAP%qw_no_w)) GAP%qw_no_w = .true.

     GAP%qw_no_q = .not. GAP%qw_no_q
     GAP%qw_no_w = .not. GAP%qw_no_w

!     GAP%do_qw_so3 = .true.
     GAP%r_cut = maxval(GAP%qw_cutoff)
  else
     call system_abort('convert: datafile_coordinates '//trim(coordinates)//' unknown')
  endif

  if( get_value(my_dictionary,'n_species',GAP%n_species) ) then
     allocate( GAP%species_Z(GAP%n_species), w_Z(GAP%n_species), z_eff(GAP%n_species) )

     if( GAP%n_species == 1 ) then
        if( .not. ( get_value(my_dictionary,'Z',GAP%species_Z(1)) .and. get_value(my_dictionary,'w',w_Z(1))) ) &
        call system_abort('convert: no species information found')
       
        if( .not. get_value(my_dictionary,'z_eff',GAP%z_eff(GAP%species_Z(1))) ) &
        call system_abort('convert: no z_eff information found')
     else
        if( .not. ( get_value(my_dictionary,'Z',GAP%species_Z) .and. get_value(my_dictionary,'w',w_Z) ) ) &
        call system_abort('IPModel_GAP_Initialise_str: no species information found')


        if( .not. get_value(my_dictionary,'z_eff',z_eff) ) &
        call system_abort('IPModel_GAP_Initialise_str: no species information found')

     endif

     allocate( GAP%w_Z(maxval(GAP%species_Z)) )

     GAP%z_eff = 0.0_dp
     GAP%w_Z = 0.0_dp
     GAP%z_eff(GAP%species_Z) = z_eff
     GAP%w_Z(GAP%species_Z) = w_Z

  endif

  if( .not. get_value(my_dictionary,'do_ewald',GAP%do_ewald) ) GAP%do_ewald = .false.
  if( .not. get_value(my_dictionary,'do_ewald_corr',GAP%do_ewald_corr) ) GAP%do_ewald_corr = .false.
  if( .not. get_value(my_dictionary,'do_core',GAP%do_core) ) GAP%do_core = .false.
  if( .not. get_value(my_dictionary,'ip_args',GAP%ip_args) ) GAP%ip_args = ''
  if( .not. get_value(my_dictionary,'e0',GAP%e0) ) GAP%e0 = 0.0_dp
  if( .not. get_value(my_dictionary,'f0',GAP%f0) ) GAP%f0 = 0.0_dp

  call finalise(my_dictionary)

  if(.not.has_xml_file) xml_file = 'gp_'//GAP%my_gp%n//'.xml'
  call teach_sparse_print_xml(GAP,xml_file)

  call system_finalise()

end program convert
