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

! compute Gaussian fit for a distribution(s)
! compute Kullback--Leibler divergence of two distributions
! compute the derivative of the free energy at UI_xi (using UI_use_xi, which is mostly 0)

program histogram_process

use libatoms_module

implicit none

  !input
  character(len=STRING_LENGTH) :: infile1_name, infile2_name, outfile_name, listfilename
  logical                     :: file_is_list
  integer                     :: n_data
  logical                     :: KL_divergence, gaussian_fit, UI_unbias
  type(InOutput)              :: outfile

  !reading in infiles
  type(Dictionary)              :: cli_params
  type(InOutput)                :: list_file
  character(STRING_LENGTH)       :: inline
  character(1024), dimension(2) :: fields
  integer                       :: num_fields,status
  logical                       :: read_more

  real(dp), allocatable :: histogram1(:,:), histogram2(:,:)
  real(dp)              :: sum_hist1, sum_hist2
  real(dp)              :: amplitude, mean, variance, chisq
  real(dp)              :: KL_1_2, KL_2_1, KL
  integer               :: i
  real(dp)              :: UI_xi, UI_use_xi
  real(dp)              :: dA_dxi_unbiased, temp

  call system_initialise()

  call initialise(cli_params)
  call param_register(cli_params, "infile1"      , "stdin"        , infile1_name , help_string="infile1 - containing histogram with no blank lines or comments")
  call param_register(cli_params, "infile2"      , ""             , infile2_name , help_string="infile2 (needed for KL divergence)")
  call param_register(cli_params, "outfile"      , "stdout"       , outfile_name , help_string="outfile")
  call param_register(cli_params, "n_data"       , PARAM_MANDATORY, n_data       , help_string="number of data points in the histogram")
  call param_register(cli_params, "KL_divergence", "F"            , KL_divergence, help_string="calculate Kullback-Leibler divergence")
  call param_register(cli_params, "file_is_list" , "F"            , file_is_list , help_string="if infile1 is just a list of infile names (comments starting with # are allowed")
  call param_register(cli_params, "gaussian_fit" , "F"            , gaussian_fit , help_string="calculate Gaussian fit for the histogram")
  call param_register(cli_params, "UI_unbias"    , "F"            , UI_unbias    , help_string="calculate unbiased derivative of the free energy using Umbrella Integration")
  call param_register(cli_params, "UI_xi"        , "0.0"          , UI_xi        , help_string="the set value of the reaction coordinate (used as a label), for Umbrella Integration, default=0.0")
  call param_register(cli_params, "UI_use_xi"    , "0.0"          , UI_use_xi    , help_string="value of the reaction coordinate, around which the histogram is meant to be centred, for Umbrella Integration (often the histogram is shifted to 0), default=0.0")
  call param_register(cli_params, "temperature"  , "300.0"        , temp  , help_string="simulation temperature, needed for Umbrella Integration, default=300.0")
  if (.not.param_read_args(cli_params)) then
    call param_print_help(cli_params)
    call print("Usage: "//trim(EXEC_NAME)//" infile1=stdin infile2 outfile=stdout n_data KL_divergence=F gaussian_fit=F", PRINT_ALWAYS)
    call system_abort("Unable to parse command line")
  endif
  call finalise(cli_params)

  !check arguments
  if (len_trim(infile2_name)==0 .and. KL_divergence) &
    call system_abort("KL_divergence needs both infile1 and infile2")
  if (n_data<1) &
    call system_abort("At least 1 data point should be in the input files, specified in n_data. n_data="//n_data)
  if (.not.KL_divergence .and. .not.gaussian_fit .and. .not. UI_unbias) &
    call system_abort("At least 1 task of KL_divergence and gaussian_fit and UI_unbias should be chosen.")
  if (UI_unbias .and. temp<=0._dp) &
    call system_abort("UI_unbias needs temperature > 0, got "//temp)


  allocate(histogram1(n_data,1:2))
  allocate(histogram2(n_data,1:2))

  if (file_is_list) then
    listfilename=trim(infile1_name)
    call initialise(list_file,trim(listfilename),action=input)
    status=0
    read_more=.true.
    do while(status==0.and.read_more)
       inline=read_line(list_file,status)
       if (status==0) then
            if (inline(1:1)/="#") then
               read_more=.false.
            else
               call print("",file=outfile)
               call print("",file=outfile)
            endif
       endif
    enddo
    if (num_fields==1 .and. KL_divergence) &
    call system_abort("KL_divergence needs 2 infiles in the list file.")
  endif

  call initialise(outfile,trim(outfile_name),action=output,append=.false.)
  if (gaussian_fit) call print("GAUSS| #amplitude mean variance chisq",file=outfile)
  if (KL_divergence) call print("DKL| #D_KL(1|2) D_KL(2|1) D_KL(1,2)",file=outfile)

  do while (status==0)

    if (file_is_list) then
      call parse_string(inline,' ',fields,num_fields)
      infile1_name=fields(1)
      infile2_name=fields(2)
    endif

    call print("Hist1: "//trim(infile1_name))
    call print("Hist2: "//trim(infile2_name))
    histogram1 = 0._dp
    histogram2 = 0._dp
    !read in histogram(s)
    histogram1 = read_histogram(trim(infile1_name),n_data)
    if (len_trim(infile2_name)/=0) histogram2 = read_histogram(trim(infile2_name),n_data)

    !Gaussian fit if needed
    if (gaussian_fit .or. UI_unbias) then
       call gaussian(histogram1(1:n_data,1:2),amplitude,mean,variance,chisq)
       call print("Gaussian fit of histogram(1):")
       call print("  amplitude: "//amplitude)
       call print("  mean:      "//mean)
       call print("  variance:  "//variance)
       call print("")
       if (UI_unbias) then
          dA_dxi_unbiased = 1.0_dp * BOLTZMANN_K * temp * ( UI_use_xi - mean ) / variance
          call print("UI| "//trim(infile1_name)//" "//UI_xi//" "//dA_dxi_unbiased//" "//mean//" "//variance//" "//chisq,file=outfile)
       else
          call print("GAUSS| "//trim(infile1_name)//" "//amplitude//" "//mean//" "//variance//" "//chisq,file=outfile)
       endif

       if (len_trim(infile2_name)/=0) then
          call gaussian(histogram2(1:n_data,1:2),amplitude,mean,variance,chisq)
       call print("Gaussian fit of histogram(2):")
       call print("  amplitude: "//amplitude)
       call print("  mean:      "//mean)
       call print("  variance:  "//variance)
       if (UI_unbias) then
          dA_dxi_unbiased = 1.0_dp * BOLTZMANN_K * temp * ( UI_use_xi - mean ) / variance
          call print("UI| "//trim(infile2_name)//" "//UI_xi//" "//dA_dxi_unbiased//" "//mean//" "//variance//" "//chisq,file=outfile)
       else
          call print("GAUSS| "//trim(infile2_name)//" "//amplitude//" "//mean//" "//variance//" "//chisq,file=outfile)
       endif
       call print("")
       endif

    endif

    !Kullback--Leibler divergence if needed
    if (KL_divergence) then
      !distributions must add up to 1 (not the same as being normalised)
      sum_hist1 = sum(histogram1(1:n_data,2))
      sum_hist2 = sum(histogram2(1:n_data,2))
      histogram1(1:n_data,2) = histogram1(1:n_data,2) / sum_hist1
      histogram2(1:n_data,2) = histogram2(1:n_data,2) / sum_hist2
      sum_hist1=0._dp ; sum_hist2=0._dp 
      KL_1_2 = 0._dp
      KL_2_1 = 0._dp
      do i=1,n_data
        sum_hist1 = sum_hist1 + histogram1(i,2) ; sum_hist2 = sum_hist2 + histogram2(i,2)
        if (histogram1(i,2)>epsilon(0._dp) .and. histogram2(i,2)>epsilon(0._dp)) then
          !call print("Adding "//histogram1(i,2)//"*log("//histogram1(i,2)//"/"//histogram2(i,2)//") to D_KL(1|2)")
          KL_1_2 = KL_1_2 + histogram1(i,2) * log ( histogram1(i,2)/histogram2(i,2) )
          !call print("KL_1_2 so far "//KL_1_2)
          !call print("Adding "//histogram2(i,2)//"*log("//histogram2(i,2)//"/"//histogram1(i,2)//") to D_KL(2|1)")
          KL_2_1 = KL_2_1 + histogram2(i,2) * log ( histogram2(i,2)/histogram1(i,2) )
          !call print("KL_2_1 so far "//KL_2_1)
        else
          !call print("Skipping "//histogram1(i,2)//" and "//histogram2(i,2))
        endif
      enddo
      KL = 0.5_dp * (KL_1_2 + KL_2_1)
      call print("Kullback--Leibler divergence:")
      call print("  asym. D_KL(hist1|hist2): "//KL_1_2)
      call print("  asym. D_KL(hist2|hist1): "//KL_2_1)
      call print("  symm. D_KL(hist2,hist1): "//KL)
      call print("  Sum(Hist1): "//sum_hist1)
      call print("  Sum(Hist2): "//sum_hist2)
      call print("DKL| "//trim(infile1_name)//" "//trim(infile2_name)//" "//KL_1_2//" "//KL_2_1//" "//KL,file=outfile)
      call print("")
    endif

    if (file_is_list) then
      status=0
      read_more=.true.
      do while(status==0.and.read_more)
         inline=read_line(list_file,status)
         if (status==0) then
            if (inline(1:1)/="#") then
               read_more=.false.
            else
               call print("",file=outfile)
               call print("",file=outfile)
            endif
         endif
      enddo
    else
      exit
    endif

  enddo

  deallocate(histogram1)
  deallocate(histogram2)
  if (file_is_list) call finalise(list_file)
  call finalise(outfile)

  call system_finalise()

contains

  function read_histogram(filename,n_data) result(histo)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: n_data
    real(dp) :: histo(1:n_data,2)
    type(InOutput) :: infile
    integer :: i

    call initialise(infile,trim(filename),action=INPUT)
    do i=1,n_data
      call read_ascii(infile, histo(i,1:2))
    enddo
    call finalise(infile)

  end function read_histogram
    

  subroutine gaussian(histo,amplitude,mean,variance,chisq)
    real(dp), intent(in)  :: histo(:,:)
    real(dp), intent(out) :: amplitude, mean, variance
    real(dp), intent(out) :: chisq
    real(dp) :: params(3) !(/A,B,C/) of fit: Ax^2+Bx+C
    real(dp) :: x,x2,x3,x4,y,yx,yx2,one
    integer :: n_data, i
    real(dp) :: inversematrix(3,3)


    n_data = size(histo,1)

    !fit a second order polynomial to (ln y)(x)
    !least squares fit of log(y) (x) = log(1/2variance) - (x-mean)**2.0/(2*variance)
    x = 0._dp
    x2 = 0._dp
    x3 = 0._dp
    x4 = 0._dp
    y = 0._dp
    yx = 0._dp
    yx2 = 0._dp
    one = 0._dp
    do i=1,n_data
!    call print(histo(i,1)//" "//histo(i,2))
      if (histo(i,2)>epsilon(1._dp)) then
        x = x + histo(i,1)
        x2 = x2 + histo(i,1)*histo(i,1)
        x3 = x3 + histo(i,1)*histo(i,1)*histo(i,1)
        x4 = x4 + histo(i,1)*histo(i,1)*histo(i,1)*histo(i,1)
!call print(histo(i,1)//" "//log(histo(i,2)))
        y = y + log(histo(i,2))
        yx = yx + log(histo(i,2))*histo(i,1)
        yx2 = yx2 + log(histo(i,2))*histo(i,1)*histo(i,1)
        one = one + 1._dp
      !weight=1/histo(i,2)
      endif
    enddo


    !call print((/x4,x3,x2,x3,x2,x,x2,x,one/))
    !call print( reshape((/x4,x3,x2,x3,x2,x,x2,x,one/),(/3,3/)) )
    call matrix3x3_inverse( reshape((/x4,x3,x2,x3,x2,x,x2,x,one/),(/3,3/)) , inversematrix )
    !call print(inversematrix)
    !call print( inversematrix .mult. reshape((/x4,x3,x2,x3,x2,x,x2,x,one/),(/3,3/)) )
    !call print((/yx2,yx,y/))
    params = inversematrix .mult. (/yx2,yx,y/)
    !call print(params)

    chisq = 0._dp
    do i=1,n_data
      if (histo(i,2)>epsilon(1._dp)) then
        chisq = chisq + ( log(histo(i,2)) - (params(1)*histo(i,1)*histo(i,1) + params(2)*histo(i,1) + params(3) ) )**2.0_dp
      endif
    enddo

    call print("Chi^2 of fitting: "//chisq)
    call print("est. variance: "//(chisq/(one-3)))

!    least_squares(histo(x=1:size(histo,1),1), &
!                        y=1:size(histo,1),2), &
!                        sig=0.01_dp, &
!                        a=(/A,B,C/), &
!                        chisq=chisq, &
!                        funcs=func)

    variance= -1._dp/2._dp*params(1)
    mean=params(2)/variance
    amplitude=exp(params(3)+mean**2._dp/(2._dp*variance)) * 2._dp*variance

  end subroutine gaussian

!  function func(x,afunc)
!    real(dp) :: x,afunc
!    afunc = A*x**2.0_dp + B*x + C
!  end function func

end program histogram_process
