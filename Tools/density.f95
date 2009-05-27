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
    character(FIELD_LENGTH)               :: xyzfilename, datafilename
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

    call system_initialise(SILENT)
    call verbosity_push(NORMAL)

    call initialise(params_in)
    call param_register(params_in, 'xyzfile', param_mandatory, xyzfilename)
    call param_register(params_in, 'datafile', 'data.den1', datafilename)
    call param_register(params_in, 'NumBins', param_mandatory, numbins)
    call param_register(params_in, 'decimation', '1', decimation)
    call param_register(params_in, 'from', '0', from)
    call param_register(params_in, 'to', '0', to)
    call param_register(params_in, 'IO_Rate', '1', IO_Rate)
    call param_register(params_in, 'Gaussian', 'F', Gaussian_smoothing)
    call param_register(params_in, 'sigma', '0.0', Gaussian_sigma)
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

       if (status.ne.0) then
          call print('double exit')
          exit
       endif

       write(mainlog%unit,'(a,a,i0,$)') achar(13),'Frame ',frame_count
  
       if (frame_count.gt.from) then
          if ((frame_count.gt.to).and.(to.ne.0)) exit
          call calc_connect(at)     
          call map_into_cell(at)
          do i = 1, at%N
          
             if (.not.Gaussian_smoothing) then
                t = at%g .mult. at%pos(:,i)
                box_a = floor(real(numbins,dp) * (t(1)+0.5_dp)) + 1
                box_b = floor(real(numbins,dp) * (t(2)+0.5_dp)) + 1
                box_c = floor(real(numbins,dp) * (t(3)+0.5_dp)) + 1
             
                mass_cuml(box_a,box_b,box_c) = mass_cuml(box_a,box_b,box_c) + ElementMass(at%Z(i))
             else
!call print('hi!'//i)
                Gaussian_centre = at%pos(1:3,i)
                increment = at%lattice .mult. (/1.0_dp/real(numbins,dp),1.0_dp/real(numbins,dp),1.0_dp/real(numbins,dp)/)
                x_min_bin = -0.5*at%lattice(1,1)
                y_min_bin = -0.5*at%lattice(2,2)
                z_min_bin = -0.5*at%lattice(3,3)
!call print('hi!'//increment)
                do j=1,numbins
                   do k=1,numbins
                      do l=1,numbins
!call print('hi! '//j//' '//k//' '//l)
                         x_min = x_min_bin + real(j-1,dp)*increment(1)
!call print('x_min '//x_min)
                         y_min = y_min_bin + real(k-1,dp)*increment(2)
!call print('y_min '//y_min)
                         z_min = z_min_bin + real(l-1,dp)*increment(3)
!call print('z_min '//z_min)
                         x_max = x_min + increment(1)
!call print('x_max '//x_max)
                         y_max = y_min + increment(2)
!call print('y_max '//y_max)
                         z_max = z_min + increment(3)
!call print('z_max '//z_max)
                         !call print('min_bin: '//x_min_bin//', max_bin: '//x_max)
                         !call print('ERF(min_bin) = '//erf((Gaussian_centre(1)-min_bin)/(Gaussian_sigma*sqrt(2._dp))))
                         !call print('ERF(max_bin) = '//erf((max_bin-Gaussian_centre(1))/(Gaussian_sigma*sqrt(2.0_dp))))
                         !call print('Adding to bin '//j//' '//erf((Gaussian_centre(1)-min_bin)/(Gaussian_sigma*sqrt(2._dp))) + erf((max_bin-Gaussian_centre(1))/(Gaussian_sigma*sqrt(2.0_dp))))
!call print('hi! '//size(mass_cuml,1)//' '//size(mass_cuml,2)//' '//size(mass_cuml,3))
!call print('hi! '//mass_cuml(j,k,l))
!call print('Gaussian centre'//Gaussian_centre)
!call print('ERF '//( (x_min-Gaussian_centre(1))/(Gaussian_sigma*sqrt(2.0_dp))) )
!call print('ERF term '//( - erf((x_min-Gaussian_centre(1))/(Gaussian_sigma*sqrt(2.0_dp))) ))
!call print('ERF '//( (x_max-Gaussian_centre(1))/(Gaussian_sigma*sqrt(2.0_dp))) )
!call print('ERF term '//( + erf((x_max-Gaussian_centre(1))/(Gaussian_sigma*sqrt(2.0_dp))) ))
!call print('ERF '//( (y_min-Gaussian_centre(2))/(Gaussian_sigma*sqrt(2.0_dp))) )
!call print('ERF term '//( - erf((y_min-Gaussian_centre(2))/(Gaussian_sigma*sqrt(2.0_dp))) ))
!call print('ERF '//( (y_max-Gaussian_centre(2))/(Gaussian_sigma*sqrt(2.0_dp))) )
!call print('ERF term '//( + erf((y_max-Gaussian_centre(2))/(Gaussian_sigma*sqrt(2.0_dp))) ))
!call print('ERF '//( (z_min-Gaussian_centre(3))/(Gaussian_sigma*sqrt(2.0_dp))) )
!call print('ERF term '//( - erf((z_min-Gaussian_centre(3))/(Gaussian_sigma*sqrt(2.0_dp))) ))
!call print('ERF '//( (z_max-Gaussian_centre(3))/(Gaussian_sigma*sqrt(2.0_dp))) )
!call print('ERF term '//( + erf((z_max-Gaussian_centre(3))/(Gaussian_sigma*sqrt(2.0_dp))) ))
                         mass_cuml(j,k,l) = mass_cuml(j,k,l) + ElementMass(at%Z(i)) * &
                           0.5_dp * ( - erf((x_min-Gaussian_centre(1))/(Gaussian_sigma*sqrt(2.0_dp))) + erf((x_max-Gaussian_centre(1))/(Gaussian_sigma*sqrt(2.0_dp))) ) * &
                           0.5_dp * ( - erf((y_min-Gaussian_centre(2))/(Gaussian_sigma*sqrt(2.0_dp))) + erf((y_max-Gaussian_centre(2))/(Gaussian_sigma*sqrt(2.0_dp))) ) * &
                           0.5_dp * ( - erf((z_min-Gaussian_centre(3))/(Gaussian_sigma*sqrt(2.0_dp))) + erf((z_max-Gaussian_centre(3))/(Gaussian_sigma*sqrt(2.0_dp))) )
!call print('hi! '//mass_cuml(j,k,l))
                      enddo
                   enddo
                enddo
!call print('hi!'//i)
             endif
          end do
  
!call print('hi!')
          if (mod(frame_count-from,IO_RATE)==0) call write_data
       endif
  
       !Try to read another frame
       call read_xyz(at,xyzfile,status=status)
       
    end do
    
    call print('')
    call write_data
    call print('Read '//frame_count//' frames')
  
    !Free up memory
    call finalise(at)
    call finalise(xyzfile)
  
    deallocate(mass_cuml)
  
    call print('Finished.')
  
    call verbosity_pop
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
