program xyz_to_SPECSYM_start
use libatoms_module
implicit none
  type(Atoms) :: at
  integer i
  real(dp) :: dipole_mom(3)

  call system_initialise()
  call read_xyz(at, "stdin")


  mainlog%prefix = "SPECSYM"

  call print(at%N//" 0.02 0.005 0.0 0.0")
  call print("300.0 0.0 1 6.0")
  do i=1, at%N
    call print(ElementMass(at%Z(i))/MASSCONVERT)
  end do
  call print("0.0 0")
  do i=1, at%N
    call print (at%pos(:,i)/BOHR // " " // at%Z(i) // " " // at%Z(i))
    call print ("0.0 0.0 0.0")
  end do

  if (.not. get_value(at%params, 'Dipole_Moment', dipole_mom)) dipole_mom = 0.0_dp

  call print (dipole_mom/BOHR)
  call print ("0.0 0.0 0.0")

  mainlog%prefix = ""

  call system_finalise()
end program
