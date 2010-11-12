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

! RDFD
! Tool for calculating Radial Distribution Functions from xyz files
! The extra 'D' shows that the RDF is calculated as a function of distance from some centre and binned
! The centre is a position and not an atom
!
! rdfd <xyzfile> <symbol> <symbol> <centre x y z> <distance cutoff> <distance binwidth> <rdf cutoff> <rdf binwidth> <output file> <decimation> <from> <to> <IO_Rate>
program RDFD

    use libatoms_module

    implicit none

    type(Atoms)                           :: at
    type(InOutput)                        :: xyzfile, datafile
    integer                               :: Za, Zb
    integer                               :: frames_read, frames_processed
    integer                               :: status
    integer                               :: i, j

    !Binning & histograms
    integer                               :: rdf_bin, dist_bin
    integer                               :: num_rdf_bins, num_dist_bins
    integer                               :: num_A_atoms, num_B_atoms
    real(dp)                              :: r, B_density
    integer,  allocatable, dimension(:,:) :: rdf_data
    integer,  allocatable, dimension(:)   :: num_atoms

    !Input
    type(Dictionary)                      :: params_in
    character(FIELD_LENGTH)               :: xyzfilename, datafilename
    real(dp)                              :: d_cutoff, d_binwidth
    real(dp)                              :: rdf_cutoff, rdf_binwidth
    real(dp), dimension(3)                :: centre
    character(FIELD_LENGTH)               :: mask1, mask2
    integer                               :: IO_Rate
    integer                               :: decimation
    integer                               :: from, to
    logical                               :: Gaussian_smoothing
    real(dp)                              :: Gaussian_sigma

    !Start up LOTF, suppressing messages
    call system_initialise(PRINT_NORMAL)
    call verbosity_push(PRINT_NORMAL)

#ifdef DEBUG
    call print('********** DEBUG BUILD **********')
    call print('')
#endif

    call initialise(params_in)
    call param_register(params_in, 'xyzfile', param_mandatory, xyzfilename, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params_in, 'datafile', 'data.den1', datafilename, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params_in, 'AtomMask1', param_mandatory, mask1, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params_in, 'AtomMask2', param_mandatory, mask2, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params_in, 'centre', '0. 0. 0.', centre, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params_in, 'RDFCutoff', param_mandatory, rdf_cutoff, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params_in, 'RDFBinWidth', param_mandatory, rdf_binwidth, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params_in, 'DistCutoff', param_mandatory, d_cutoff, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params_in, 'DistBinWidth', param_mandatory, d_binwidth, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params_in, 'decimation', '1', decimation, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params_in, 'from', '0', from, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params_in, 'to', '0', to, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params_in, 'IO_Rate', '1', IO_Rate, help_string="No help yet.  This source file was $LastChangedBy$")
!    call param_register(params_in, 'Gaussian', 'F', Gaussian_smoothing, help_string="No help yet.  This source file was $LastChangedBy$")
!    call param_register(params_in, 'sigma', '0.0', Gaussian_sigma, help_string="No help yet.  This source file was $LastChangedBy$")
    if (.not. param_read_args(params_in, do_check = .true.)) then
       if (EXEC_NAME == '<UNKNOWN>') then
          call print_usage
       else
          call print_usage(EXEC_NAME)
       end if
      call system_abort('could not parse argument line')
    end if
    call finalise(params_in)

    call print('Run_parameters: ')
    call print('==================================')
    call print('    Input file: '//trim(xyzfilename))
    call print('   Output file: '//trim(datafilename))
    call print('     AtomMask1: '//trim(mask1))
    call print('     AtomMask2: '//trim(mask2))
    call print('        centre: '//round(centre(1),3)//' '//round(centre(2),3)//' '//round(centre(3),3))
    call print('     RDFCutoff: '//round(rdf_cutoff,3))
    call print('   RDFBinWidth: '//round(rdf_binwidth,3))
    call print('    DistCutoff: '//round(d_cutoff,3))
    call print('  DistBinWidth: '//round(d_binwidth,3))
    call print('    decimation: '//decimation)
    if (decimation == 1) then
       call print('             Processing every frame')
    else
       write(line,'(a,i0,a,a)')'             Processing every ',decimation,th(decimation),' frame'
       call print(line)
    end if
    call print('    from Frame: '//from)
    call print('      to Frame: '//to)
    call print('       IO_Rate: '//IO_Rate)
    call print('     Gaussians: '//Gaussian_smoothing)
    if (Gaussian_smoothing) call print('        sigma: '//round(Gaussian_sigma,3))
    call print('==================================')
    call print('')

    call initialise(xyzfile,xyzfilename,action=INPUT)

    call print('Mask 1:')
    call print('=======')
    Za = Atomic_Number(mask1)
    call print('')
    call print('Selecting all '//trim(ElementName(Za))//' atoms (Z = '//Za//')')
    call print('')

    call print('Mask 2:')
    call print('=======')
    Zb = Atomic_Number(mask2)
    call print('')
    call print('Selecting all '//trim(ElementName(Zb))//' atoms (Z = '//Zb//')')
    call print('')

    if (d_cutoff < 0.0_dp) call system_abort('Distance cutoff < 0.0 Angstroms')
    if (d_binwidth < 0.0_dp) call system_abort('Distance bin width < 0.0 Angstroms')
    if (rdf_cutoff < 0.0_dp) call system_abort('RDF cutoff < 0.0 Angstroms')
    if (rdf_binwidth < 0.0_dp) call system_abort('RDF bin width < 0.0 Angstroms')

    !Make the number of bins and bin widths consistent
    num_rdf_bins = ceiling(rdf_cutoff / rdf_binwidth)
    rdf_binwidth = rdf_cutoff / real(num_rdf_bins,dp)
    num_dist_bins = ceiling(d_cutoff / d_binwidth)
    d_binwidth = d_cutoff / real(num_dist_bins,dp)

    !Set up arrays
    allocate( rdf_data(num_rdf_bins,num_dist_bins), num_atoms(num_dist_bins) )

    rdf_data = 0
    num_atoms = 0

    call print('')
    call print('RDFs from ('//centre//') out to '//round(d_cutoff,3)//'A')
    call print('Each RDF has '//num_rdf_bins//' bins x '//round(rdf_binwidth,3)//'A per bin = '//round(rdf_cutoff,3)//'A cutoff')
    call print('')
  
    call print('Reading data...')
    frames_read = 0
    frames_processed = 0

    do
  
       call read_xyz(at,xyzfile,status=status)
       if (status/=0) then
          if (frames_processed > 1) call write_data
          exit
       end if
       frames_read = frames_read + 1
  
       ! Do the calculation between form and to:
       if (frames_read.ge.from) then
          if ((frames_read.gt.to) .and.(to.gt.0)) exit
          ! Loop over all atoms
          num_A_atoms = 0
          num_B_atoms = 0
          do i = 1, at%N
       
             !Count 'B' atoms
             if (at%Z(i) == Zb) num_B_atoms = num_B_atoms + 1
       
             !Do we have an 'A' atom?
             if (at%Z(i) /= Za) cycle
       
             num_A_atoms = num_A_atoms + 1
       
             !Is this 'A' atom within the cutoff of the centre?
             r = distance_min_image(at,i,centre)
             if (r >= d_cutoff) cycle
       
             !Increase the atom count at this distance
             dist_bin = bin(0.0_dp,d_cutoff,d_binwidth,r)
             num_atoms(dist_bin) = num_atoms(dist_bin) + 1
       
             !Loop over all other atoms
             do j = 1, at%N
       
                ! Skip self
                if (i==j) cycle
       
                !Do we have a 'B' atom?
                if (at%Z(j) /= Zb) cycle
       
                !Is this 'B' atom within the rdf cutoff of the 'A' atom?
                r = distance_min_image(at,i,j)
                if (r >= rdf_cutoff) cycle
       
                !Increase the count at this distance
                rdf_bin = bin(0.0_dp,rdf_cutoff,rdf_binwidth,r)
                rdf_data(rdf_bin,dist_bin) = rdf_data(rdf_bin,dist_bin) + 1
       
             end do
       
          end do  
       
          B_density = real(num_B_atoms,dp) / cell_volume(at)
          
          frames_processed = frames_processed + 1
          write(mainlog%unit,'(a,a,i0,$)') achar(13),'Frame ',frames_read
       
          if (mod(frames_processed,IO_Rate)==0) then
             !call print('writing data to datafile '//trim(datafile%filename))
             call write_data
          endif
       endif
       
       !Skip ahead (decimation-1) frames
       do i = 1, (decimation-1)
          call read_xyz(xyzfile,status)
          if (status/=0) exit
          frames_read = frames_read + 1
       end do
          if (status/=0) exit
  
    end do
  
    call print('Finished.')
    call print('Read '//frames_read//' frames, processed '//frames_processed//' frames.')
  
    call verbosity_pop
    call system_finalise()

contains

  subroutine write_data

    integer  :: i, j
    real(dp) :: dist, rdf_dist, rdf, rdf_int, dV

    call initialise(datafile,datafilename,action=OUTPUT)

    do i = 1, num_dist_bins

       dist = (real(i,dp)-0.5_dp)*d_binwidth

       do j = 1, num_rdf_bins

          rdf_dist = (real(j,dp)-0.5_dp)*rdf_binwidth
          dV = PI*rdf_binwidth * ( 4.0_dp*rdf_dist*rdf_dist + rdf_binwidth*rdf_binwidth/3.0_dp )
          
          if (num_atoms(i) > 0) then
             rdf_int = real(rdf_data(j,i),dp) / (real(num_atoms(i),dp))
          else
             rdf_int = 0
          end if
          rdf = rdf_int / ( dV * B_density )

          call print(round(dist,5)//' '//round(rdf_dist,5)//' '//round(rdf,5)//' '//round(rdf_int,5),file=datafile)

       end do

       call print('',file=datafile) !Insert blank line to separate data sets

    end do

    call finalise(datafile)

  end subroutine write_data


  subroutine print_usage(name)

    character(*), optional, intent(in) :: name

    if (present(name)) then
       write(line,'(3a)')'Usage: ',trim(name),' <xyzfile> <AtomMask1> <AtomMask2> <Centre> <DistCutoff> <DistBinwidth> <RDFCutoff> <RDFBinWidth> <datafile> [<decimation>] [<from>] [<to>] [<IO_Rate>]'
    else
       write(line,'(a)')'Usage: rdfd <xyzfile> <AtomMask1> <AtomMask2> <Centre> <DistCutoff> <DistBinwidth> <RDFCutoff> <RDFBinWidth> <datafile> [<decimation>] [<from>] [<to>] [<IO_Rate>]'
    end if
    call print(line)
    call print(' <xyzfile>       The input xyz file.')
    call print(' <AtomMask1>     An element symbol, e.g. H or Ca')
    call print(' <AtomMask2>     An element symbol, e.g. H or Ca')
    call print(' <Centre>        The central atom.')
    call print(' <DistBinWidth>  The width of bins into which rdfs are placed.') 
    call print(' <DistCutoff>    The cutoff radius for the distance from the centre.') 
    call print(' <RDFCutoff>    The cutoff radius for the RDF in Angstroms.')
    call print(' <RDFBinWidth>  The width of each bin in Angstroms used for individual rdf calculations.')
    call print(' <datafile>      The output data file.')
    call print(' <decimation>    Optional. Only process 1 out of every n frames')
    call print(' <fromn>         Optional. Only process frames from this frame')
    call print(' <to>            Optional. Only process frames until this frame')
    call print(' <IO_Rate>       Optional. Write data after every n processed frames')
    call print('')

    call verbosity_pop
    call system_finalise
        
  end subroutine print_usage

end program RDFD
