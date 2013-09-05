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

program ts_main 
  use libAtoms_module
  use Potential_module 
  use ts_module
  use tsParams_module

  implicit none

  type(Dynamicalsystem) :: ds_in, ds_fin 
  type(Potential)     :: classicalpot, qmpot
  type(Potential) :: hybrid_pot
  type(Atoms)         :: at_in, at_fin, at_image
  type(inoutput)      :: xmlfile, file_res
  type(Cinoutput)     :: in_in, in_fin, in_image, outimage 
  type(TS)            :: tts 
  type(Dictionary)    :: pot_params
  type(MPI_Context)   :: mpi
  type(tsParams)      :: params

  real(dp), allocatable :: forces(:,:)
  integer               :: niter, steps, im
  character(len=STRING_LENGTH) :: xmlfilename
  real(dp), allocatable, dimension(:,:) :: conf

  call initialise(mpi)

  if (mpi%active) then
     call system_initialise (common_seed = .true., enable_timing=.true., mpi_all_inoutput=.false.)
     call print('MPI run with '//mpi%n_procs//' processes')
  else
     call system_initialise( enable_timing=.true.)
     call print('Serial run')
  end if

  xmlfilename = 'ts.xml'
  call print_title('Initialisation')
  call initialise(params)

  call print('Reading parameters from file '//trim(xmlfilename))
  call initialise(xmlfile,xmlfilename,INPUT)
  call read_xml(params,xmlfile)
  call verbosity_push(params%io_verbosity)   ! Set base verbosity
  call print(params)

  call print ("Initialising classical potential with args " // trim(params%classical_args) &
       // " from file " // trim(xmlfilename))
  call rewind(xmlfile)
  call initialise(classicalpot, params%classical_args, xmlfile, mpi_obj = mpi)
  call Print(classicalpot)

  if (params%simulation_hybrid) then
     call print ("Initialising QM potential with args " // trim(params%qm_args) &
          // " from file " // trim(xmlfilename))
     call rewind(xmlfile)
     call initialise(qmpot, params%qm_args, xmlfile, mpi_obj=mpi)
     call finalise(xmlfile)
     call Print(qmpot)

     call initialise(pot_params)
     call set_value(pot_params,'qm_args_str',params%qm_args_str)
     call set_value(pot_params,'method','force_mixing')
     call initialise(hybrid_pot, 'ForceMixing '//write_string(pot_params), &
             classicalpot, qmpot, mpi_obj=mpi)

     call print_title('Hybrid Potential')
     call print(hybrid_pot)

     call finalise(pot_params)
  end if

  call print_title('Initialising First and Last Image')
  call Initialise(in_in, trim(params%chain_first_conf), action=INPUT)
  call read(at_in, in_in)
  call Initialise(in_fin, trim(params%chain_last_conf), action=INPUT)
  call read(at_fin, in_fin)
  call Print('Setting neighbour cutoff to '//(cutoff(classicalpot))//' A.')
  call set_cutoff(at_in, cutoff(classicalpot))
  call set_cutoff(at_fin, cutoff(classicalpot))

  call initialise(ds_in,at_in)
  call initialise(ds_fin,at_fin)
! Fix atoms (if necessary) 
  if(params%chain_nfix.ne.0) then
    ds_in%atoms%move_mask(ds_in%atoms%N-params%chain_nfix+1:ds_in%atoms%N) = 0
  endif
  ds_in%Ndof = 3*count(ds_in%atoms%move_mask == 1)
  allocate(forces(3,at_in%N))

  tts%cos%N = params%chain_nimages 
! if you want to initialise interpolating between the first and last image
  if(.not.params%simulation_restart) then
     if(params%minim_end) then
        call Print_title('First Image Optimisation')
        if (.not. params%simulation_hybrid) then
          call calc_connect(ds_in%Atoms)
          steps = minim(classicalpot, ds_in%atoms, method=params%minim_end_method, convergence_tol=params%minim_end_tol, &
             max_steps=params%minim_end_max_steps, linminroutine=params%minim_end_linminroutine, &
             do_pos=.true., do_lat=.false., do_print=.false., &
             args_str=params%classical_args_str, eps_guess=params%minim_end_eps_guess)
        else
          steps = minim(hybrid_pot, ds_in%atoms, method=params%minim_end_method, convergence_tol=params%minim_end_tol, &
             max_steps=params%minim_end_max_steps, linminroutine=params%minim_end_linminroutine, &
             do_pos=.true., do_lat=.false., do_print=.false., &
             eps_guess=params%minim_end_eps_guess)
        end if

        call Print_title('Last Image Optimisation')
        if (.not. params%simulation_hybrid) then
          call calc_connect(ds_in%Atoms)
          steps = minim(classicalpot, ds_fin%atoms, method=params%minim_end_method, convergence_tol=params%minim_end_tol, &
             max_steps=params%minim_end_max_steps, linminroutine=params%minim_end_linminroutine, &
             do_pos=.true., do_lat=.false., do_print=.false., &
             args_str=params%classical_args_str, eps_guess=params%minim_end_eps_guess)
        else
          steps = minim(hybrid_pot, ds_fin%atoms, method=params%minim_end_method, convergence_tol=params%minim_end_tol, &
             max_steps=params%minim_end_max_steps, linminroutine=params%minim_end_linminroutine, &
             do_pos=.true., do_lat=.false., do_print=.false., &
             eps_guess=params%minim_end_eps_guess)
        end if
     endif
   
     call print_title('Initialisation of the chain of state interpolating between the first and last image')
     call initialise(tts,ds_in%atoms,ds_fin%atoms,params)

! if you want to start from previous configuration for the path
  else
     allocate(conf(tts%cos%N, 3 * at_in%N) )
     do im=1,tts%cos%N
       call Initialise(in_image, 'conf.'//im//'.xyz')
       call read(at_image, in_image)
       conf(im,:) = reshape(at_image%pos, (/3*at_image%N/) ) 
       call finalise(at_image)   
       call finalise(in_image)
     enddo
     call print_title('Initialisation of the chain of state using the guessed path')
     call initialise(tts,ds_in%atoms,conf, params)
  endif

  call print(tts)
  if (.not. mpi%active .or. (mpi%active .and.mpi%my_proc == 0)) then
     do im =1, tts%cos%N
       call initialise(outimage, 'image.'//im//'.xyz', action=OUTPUT)
       call write(outimage, tts%cos%image(im)%at)
       call finalise(outimage)
     enddo
  end if

  call initialise(file_res, "out.dat",OUTPUT)

  call print_title('Transition state calculation')
  if (.not. params%simulation_hybrid) then
     call calc(tts,classicalpot, niter, params, file_res, mpi)
  else
     call calc(tts,hybrid_pot, niter, params, file_res, mpi)
  endif
  call print('Number or Iterations :  ' // niter )

  if (.not. mpi%active .or. (mpi%active .and.mpi%my_proc == 0)) then
     do im =1, tts%cos%N 
       call initialise(outimage, 'final.'//im//'.xyz', action=OUTPUT)
       call write(outimage, tts%cos%image(im)%at)
       call finalise(outimage)
     enddo
  end if

  call finalise(ds_in)
  call finalise(ds_fin)
  call finalise(at_in)
  call finalise(at_fin)
  deallocate(forces)

  call system_finalise()

  end program ts_main 
