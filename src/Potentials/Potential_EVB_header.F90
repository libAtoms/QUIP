
  public :: Potential_EVB
  type Potential_EVB
     type(MPI_context) :: mpi

     type(Potential), pointer :: pot1 => null()    !% The underlying MM potential, pot1 and pot2

     character(STRING_LENGTH) :: mm_args_str        !% Args string to be passed to 'calc' method of pot1 and pot2

     character(STRING_LENGTH) :: topology_suffix1   !% Suffix in topology filename, to be added to mm_args_str of pot1
     character(STRING_LENGTH) :: topology_suffix2   !% Suffix in topology filename, to be added to mm_args_str of pot2

     integer :: form_bond(2)                       !% Atom pair that is bonded in EVB1 only
     integer :: break_bond(2)                      !% Atom pair that is bonded in EVB2 only

     real(dp) :: diagonal_dE2                      !% The energy offset to E2
     real(dp) :: offdiagonal_A12                   !% The offdiagonal pre-exponent factor of EVB Hamiltonian
     real(dp) :: offdiagonal_mu12                  !% The offdiagonal exponent factor of the EVB Hamiltonian
     real(dp) :: offdiagonal_mu12_square           !% The offdiagonal exponent factor of the EVB Hamiltonian:  A exp(-mu12(r-r0)-mu12_square(r-r0)^2)
     real(dp) :: offdiagonal_r0                    !% The offdiagonal exponent of the EVB Hamiltonian

     logical :: save_forces                        !% Whether to save forces from the 2 MM calculations in the Atoms object (needed for force on E_GAP)
     logical :: save_energies                      !% Whether to save energies from the 2 MM calculations in the Atoms object (needed for E_GAP)

  end type Potential_EVB

  interface Initialise
     module procedure Potential_EVB_Initialise
  end interface

  interface Finalise
     module procedure Potential_EVB_Finalise
  end interface

  interface Print
     module procedure Potential_EVB_Print
  end interface

  interface Cutoff
     module procedure Potential_EVB_Cutoff
  end interface

  interface Calc
     module procedure Potential_EVB_Calc
  end interface

