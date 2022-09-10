program get_qw

  use libatoms_module
  use steinhardt_nelson_qw_module
  use system_module, only: isnan


  implicit none

  !local variables
  type(Dictionary)      :: params
  type(Atoms) :: at
  type(Cinoutput) :: infile, outfile       
  character(STRING_LENGTH)   :: atfile_in
  character(STRING_LENGTH)   :: central_mask, neighbour_mask
  integer  :: stat, n_string
  real(dp) :: r_cut, r_cut_min, summa_q, summa_w, q_global, w_glob
  real(dp), dimension(:), pointer :: q, w
  integer :: l, i , error, i_config, iat
  logical :: calc_QWave, print_QWxyz
  logical, allocatable, dimension(:) :: c_mask, n_mask

  call system_initialise(PRINT_NORMAL)
  call verbosity_push(PRINT_NORMAL)

  call initialise(params)
  call param_register(params, 'atfile_in', PARAM_MANDATORY, atfile_in,'') ! xyz input filename
  call param_register(params, 'r_cut', PARAM_MANDATORY, r_cut,'') ! cutoff for neighbours taken into account
  call param_register(params, 'l', PARAM_MANDATORY, l,'') ! l=2, 4, 6, 8, ....
  call param_register(params, 'r_cut_min', '0.0', r_cut_min,'') ! minimum cutoff within which atoms are not considered
  call param_register(params, 'calc_QWave', 'F', calc_QWave,'') ! average bond order params over the configuration
  call param_register(params, 'print_QWxyz', 'F', print_QWxyz,'') ! print an xyz with all QW as properties
  call param_register(params, 'central_mask', 'all', central_mask,'') ! which type of atom can be central
  call param_register(params, 'neighbour_mask', 'all', neighbour_mask,'') ! which type of atoms are neighbours

  if (.not. param_read_args(params, check_mandatory = .true.)) then
       call print("Usage: get_qw [atfile_in] [r_cut] [l] [r_cut_min=0.0] [calc_QWave=F] [print_QWxyz=F] &
                 &[central_mask=all] [neighbour_mask=all]" ,PRINT_SILENT)
       call system_abort('Exit: Mandatory argument(s) missing...')
  endif
  call finalise(params)

  call initialise(infile,atfile_in)

!     read(COMMAND_ARG(4),*) averageL
!     read(COMMAND_ARG(5),*) maximalQ 
!     read(COMMAND_ARG(6),*) atom_mask
!     read(COMMAND_ARG(7),*) print_qw

  if (print_QWxyz) then
     if (central_mask /= 'all' .or. neighbour_mask /= 'all') &
        & call system_abort("if qw to be printed to an xyz file, no mask should be used. Use 'all' option (default)")
     call initialise(outfile,"qw"//l//"_"//trim(atfile_in),action=output)
  endif


  i_config = 0
  do 
     call read(at, infile, error=error)
     if (error /= 0) then
        exit
     endif
     i_config = i_config + 1

     if (i_config == 1) then ! do the mask initialisation for the first config.

        allocate(c_mask(at%N), n_mask(at%N))
        c_mask = .true. ! default: all atoms taken into account
        n_mask = .true. ! default: all atoms taken into account
        if (central_mask /= 'all') then
           !call print('Selecting all '//trim(ElementName(at%Z(iat)))//' atoms (Z = '//at%Z(iat)//') as the central atom')
           do iat = 1,at%N
              if (central_mask /= trim(ElementName(at%Z(iat)))) c_mask(iat)=.false.
           enddo
        endif
        !write(*,*) c_mask
        if (neighbour_mask /= 'all') then
           !call print('Selecting all '//trim(ElementName(at%Z(iat)))//' atoms (Z = '//at%Z(iat)//') as the central atom')
           do iat = 1,at%N
              if (neighbour_mask /= trim(ElementName(at%Z(iat)))) n_mask(iat)=.false.
           enddo
        endif
        !write(*,*) n_mask

     endif

     call set_cutoff(at,r_cut)
     call calc_connect(at)
     call calc_qw(at,l,do_q=.true.,do_w=.true.,cutoff=r_cut,min_cutoff=r_cut_min,mask=c_mask, &
                  & mask_neighbour=n_mask,cutoff_transition_width=0.0_dp) ! l optional, default 4 (q4 es w4)
 
     if( .not. assign_pointer(at, 'q'//l, q) ) &
             & call system_abort('could not assign pointer to atoms object')
     if( .not. assign_pointer(at, 'w'//l, w) ) &
             & call system_abort('could not assign pointer to atoms object')

     do i = 1,size(w)
        if (isnan(w(i))) w(i) = 0.0_dp
     enddo
     do i = 1,size(q)
        if (isnan(q(i))) q(i) = 0.0_dp
     enddo

     if (calc_QWave) then ! print average QW parameters

        summa_q = 0.0_dp
        summa_w = 0.0_dp
        do iat = 1,at%N
           ! average only atoms with atom_mask true
           if (c_mask(iat)) then
              summa_q = summa_q + q(iat)
              summa_w = summa_w + w(iat)
           endif
        enddo
        write(*,*) summa_q/count(c_mask) , summa_w/count(c_mask) 

     else
        do iat = 1,at%N
           if (c_mask(iat) ) then
              write(*,*) i_config, q(iat), w(iat)
           endif
        enddo
        write(*,*)
        write(*,*)

     endif
     
     if (print_QWxyz) call write(at,outfile, properties="species:pos:q"//l//":w"//l)
     call finalise(at)
  enddo

  if (print_QWxyz) call finalise(outfile)
  call system_finalise()

       
end program get_qw
