program test

  use libatoms_module
  use bispectrum_module

  implicit none

  type(cinoutput) :: xyzfile
  type(atoms) :: at
  type(qw_so3) :: qw
  type(grad_qw_so3) :: dqw
  type(fourier_so3) :: f
  type(grad_fourier_so3) :: df
  logical :: status

  call system_initialise()

  call initialise(xyzfile,'test.xyz')

  call initialise(f, 10, (/5.0_dp/), (/3/), (/1.0_dp/))
  call initialise(df, 10, (/5.0_dp/), (/3/), (/1.0_dp/))
  call initialise(qw, 10, 1, do_q = .true., do_w = .true.)
  call initialise(dqw, 10, 1, do_q = .true., do_w = .true.)

  call read(xyzfile,at,frame=0)

  call set_cutoff(at, 5.0_dp)
  call calc_connect(at)

  status = test_qw_gradient(f, df, qw, dqw, at, (/6.0_dp, 6.0_dp, 6.0_dp/))

endprogram test
