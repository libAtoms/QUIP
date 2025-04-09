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

!density <xyz> <numbins> <output> [<from>] [<to>] [<IO_Rate>]
program density

  use libatoms_module

  implicit none

    type(Atoms)                           :: at
    type(InOutput)                        :: xyzfile, datafile, datafile2
    integer                               :: frame_count, frames_processed
    integer                               :: status
    integer                               :: i

    !Input
    type(Dictionary)                      :: params_in
    character(STRING_LENGTH)               :: xyzfilename, datafilename
    integer                               :: numbins
    integer                               :: IO_Rate
    integer                               :: decimation
    integer                               :: from, to
    logical                               :: Gaussian_smoothing
    real(dp)                              :: Gaussian_sigma

    !Cell and subcell parameters
    real(dp)                              :: a, b, c
    real(dp)                              :: alpha, beta, gamma
    real(dp)                              :: cellvol

    !Density binning
    real(dp)                              :: t(3), r(3)
    integer                               :: box_a, box_b, box_c
    real(dp), allocatable                 :: mass_cuml(:,:,:)
    !Gaussian density binning
    real(dp), dimension(3)                :: increment
    real(dp)                              :: x_min_bin, y_min_bin, z_min_bin
    real(dp)                              :: x_min, y_min, z_min
    real(dp)                              :: x_max, y_max, z_max
    real(dp), dimension(3)                :: Gaussian_centre
    integer                               :: j,k,l
    integer                               :: bin_j,bin_k,bin_l

    call system_initialise(PRINT_SILENT)
    call verbosity_push(PRINT_NORMAL)

    call initialise(params_in)
    call param_register(params_in, 'xyzfile', param_mandatory, xyzfilename, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params_in, 'datafile', 'data.den1', datafilename, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params_in, 'NumBins', param_mandatory, numbins, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params_in, 'decimation', '1', decimation, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params_in, 'from', '0', from, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params_in, 'to', '0', to, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params_in, 'IO_Rate', '1', IO_Rate, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params_in, 'Gaussian', 'F', Gaussian_smoothing, help_string="No help yet.  This source file was $LastChangedBy$")
    call param_register(params_in, 'sigma', '0.0', Gaussian_sigma, help_string="No help yet.  This source file was $LastChangedBy$")
    if (.not. param_read_args(params_in)) then
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
    call print('       NumBins: '//numbins)
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
 
    call initialise(xyzfile,trim(xyzfilename),action=INPUT)
    if (Gaussian_smoothing.and.(Gaussian_sigma.le.0._dp)) call system_abort('sigma must be > 0._dp')

    call read_xyz(at,xyzfile,status=status)
    call get_lattice_params(at%lattice,a,b,c,alpha,beta,gamma)
  
    allocate( mass_cuml(numbins,numbins,numbins) )

    mass_cuml = 0.0_dp

    call print('grid size '//numbins//'x'//numbins//'x'//numbins)

    cellvol = abs(matrix3x3_det(at%lattice)) / real(numbins*numbins*numbins,dp)

    call print('Subcell volume = '//round(cellvol,6)//'A^3')

    call print('Reading data...')

    frame_count = 0
    frames_processed = 0

    do
     
       if (status/=0) exit
  
       !Skip ahead (decimation-1) frames in the xyz file
       frame_count = frame_count + 1
       do i = 1, (decimation-1)
          call read_xyz(xyzfile,status)
          if (status/=0) exit
          frame_count = frame_count + 1
       end do

       write(mainlog%unit,'(a,a,i0,$)') achar(13),'Frame ',frame_count
  
       if (frame_count.gt.from) then
          if ((frame_count.gt.to).and.(to.ne.0)) exit
          call calc_connect(at)     
          call map_into_cell(at)

          !spatial binning of the cell
          if (Gaussian_smoothing) then
             increment = at%lattice .mult. (/1.0_dp/real(numbins,dp),1.0_dp/real(numbins,dp),1.0_dp/real(numbins,dp)/)
             x_min_bin = -0.5*at%lattice(1,1)
             y_min_bin = -0.5*at%lattice(2,2)
             z_min_bin = -0.5*at%lattice(3,3)
          endif

          do i = 1, at%N
          
             !calculate bin / central bin of the Gaussian
             t = at%g .mult. at%pos(:,i)
             box_a = floor(real(numbins,dp) * (t(1)+0.5_dp)) + 1
             box_b = floor(real(numbins,dp) * (t(2)+0.5_dp)) + 1
             box_c = floor(real(numbins,dp) * (t(3)+0.5_dp)) + 1
            
             if (.not.Gaussian_smoothing) then

                mass_cuml(box_a,box_b,box_c) = mass_cuml(box_a,box_b,box_c) + ElementMass(at%Z(i))

             else

                !binning is centred around the Gaussian
                Gaussian_centre = at%pos(1:3,i)

                do j=int(real(box_a,dp)+real(numbins,dp)/2._dp)-numbins+1,int(real(box_a,dp)+real(numbins,dp)/2._dp)

                   !calculate bin in the range of [1:numbins]
                   bin_j = mod(j+numbins,numbins)
                   if (bin_j.eq.0) bin_j = numbins

                   do k=int(real(box_b,dp)+real(numbins,dp)/2._dp)-numbins+1,int(real(box_b,dp)+real(numbins,dp)/2._dp)

                      !calculate bin in the range of [1:numbins]
                      bin_k = mod(k+numbins,numbins)
                      if (bin_k.eq.0) bin_k = numbins

                      do l=int(real(box_c,dp)+real(numbins,dp)/2._dp)-numbins+1,int(real(box_c,dp)+real(numbins,dp)/2._dp)

                         !calculate bin in the range of [1:numbins]
                         bin_l = mod(l+numbins,numbins)
                         if (bin_l.eq.0) bin_l = numbins

                         !spatial coordinates of this bin
                         x_min = x_min_bin + real(j-1,dp)*increment(1)
                         y_min = y_min_bin + real(k-1,dp)*increment(2)
                         z_min = z_min_bin + real(l-1,dp)*increment(3)
                         x_max = x_min + increment(1)
                         y_max = y_min + increment(2)
                         z_max = z_min + increment(3)

                         !add the integral of the Gaussian on the range [x_min,y_min,z_min:x_max,y_max,z_max]
                         mass_cuml(bin_j,bin_k,bin_l) = mass_cuml(bin_j,bin_k,bin_l) + ElementMass(at%Z(i)) * &
                           0.5_dp * ( - erf((x_min-Gaussian_centre(1))/(Gaussian_sigma*sqrt(2.0_dp))) + erf((x_max-Gaussian_centre(1))/(Gaussian_sigma*sqrt(2.0_dp))) ) * &
                           0.5_dp * ( - erf((y_min-Gaussian_centre(2))/(Gaussian_sigma*sqrt(2.0_dp))) + erf((y_max-Gaussian_centre(2))/(Gaussian_sigma*sqrt(2.0_dp))) ) * &
                           0.5_dp * ( - erf((z_min-Gaussian_centre(3))/(Gaussian_sigma*sqrt(2.0_dp))) + erf((z_max-Gaussian_centre(3))/(Gaussian_sigma*sqrt(2.0_dp))) )

                      enddo
                   enddo
                enddo

             endif
          end do
  
          if (mod(frame_count-from,IO_RATE)==0) call write_data
       endif
  
       !Try to read another frame
       call read_xyz(at,xyzfile,status=status)
       
       if (status.ne.0) then
          exit
       endif

    end do
    
    call print('')
    call write_data
    call print('Read '//frame_count//' frames')
  
    !Free up memory
    call finalise(at)
    call finalise(xyzfile)
  
    deallocate(mass_cuml)
  
    call print('Finished.')
  
!    call verbosity_pop
    call system_finalise

contains

  subroutine print_usage(name)

    character(*), optional, intent(in) :: name

    if (present(name)) then
       write(line,'(3a)')'Usage: ',trim(name),' xyzfile datafile NumBins [decimation] [from] [to] [IO_Rate] [Gaussian] [sigma]'
    else
       write(line,'(a)')'Usage: density xyzfile datafile NumBins [decimation] [from] [to] [IO_Rate] [Gaussian] [sigma]'
    end if
    call print(line)
    call print(' <xyzfile>       The input xyz file.')
    call print(' <datafile>      The output data file.')
    call print(' <NumBins>       The number of bins in the cell.')
    call print(' <decimation>    Optional. Only process 1 out of every n frames.')
    call print(' <from>          Optional. Only process frames from this frame.')
    call print(' <to>            Optional. Only process frames until this frame.')
    call print(' <IO_Rate>       Optional. Write data after every n processed frames.')
    call print(' <Gaussian>      Optional. Use Gaussians instead of delta functions.')
    call print(' <sigma>         Optional. The sigma is the sqrt(variance) of the Gaussian function.')
    call print('')
    call print('Pressing Ctrl-C during execution will leave the output file with the rdf averaged over the frames read so far')
    call print('')

    !call verbosity_pop
    call system_finalise
    stop
        
  end subroutine print_usage

  subroutine write_data

    integer :: i, j, k
    real(dp) :: density

    call initialise(datafile,trim(datafilename),action=OUTPUT)

    do k = 1, numbins
       t(3) = (real(k,dp) - 0.5_dp) / real(numbins,dp) - 0.5_dp
       do j = 1, numbins
          t(2) = (real(j,dp) - 0.5_dp) / real(numbins,dp) - 0.5_dp
          do i = 1, numbins
             t(1) = (real(i,dp) - 0.5_dp) / real(numbins,dp) - 0.5_dp
             r = at%lattice .mult. t
             density = mass_cuml(i,j,k) / (cellvol * real(frame_count-from,dp))
             density = density * abs(matrix3x3_det(at%lattice)) / sum(ElementMass(at%Z(1:at%N)))
             call print(r//' '//density,file=datafile)
          end do
          call print('',file=datafile)
       end do
       call print('',file=datafile)
    end do

    call finalise(datafile)

    call initialise(datafile2,trim(datafilename)//'_smooth',action=OUTPUT)
    do k = 1, numbins
       t(3) = (real(k,dp) - 0.5_dp) / real(numbins,dp) - 0.5_dp
       do j = 1, numbins
          t(2) = (real(j,dp) - 0.5_dp) / real(numbins,dp) - 0.5_dp
          do i = 1, numbins
             t(1) = (real(i,dp) - 0.5_dp) / real(numbins,dp) - 0.5_dp
             r = at%lattice .mult. t
             density = 0.125_dp * ( mass_cuml(i,j,k) + mass_cuml(numbins+1-i,j,k) + &
                                    mass_cuml(i,numbins+1-j,k) + mass_cuml(numbins+1-i,numbins+1-j,k) + &
                                    mass_cuml(i,j,numbins+1-k) + mass_cuml(numbins+1-i,j,numbins+1-k) + &
                                    mass_cuml(i,numbins+1-j,numbins+1-k) + mass_cuml(numbins+1-i,numbins+1-j,numbins+1-k) ) &
                                  / (cellvol * real(frame_count-from,dp))
             density = density * abs(matrix3x3_det(at%lattice)) / sum(ElementMass(at%Z(1:at%N)))
             call print(r//' '//density,file=datafile2)
          end do
          call print('',file=datafile2)
       end do
       call print('',file=datafile2)
    end do
    call finalise(datafile2)

  end subroutine write_data

end program density
