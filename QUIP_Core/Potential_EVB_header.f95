
  public :: Potential_EVB
  type Potential_EVB
     type(MPI_context) :: mpi

     type(Potential), pointer :: pot1 => null() 
     type(Potential), pointer :: pot2 => null() 

     character(FIELD_LENGTH) :: evb_step !% step number containing string to print energies
     character(FIELD_LENGTH) :: mm_args_str !% Args string to be passed to 'calc' method of EVB%pot1 and EVB%pot2
     character(FIELD_LENGTH) :: topology_suffix1 !% suffix in topology filename, to be added to mm_args_str of EVB%pot1
     character(FIELD_LENGTH) :: topology_suffix2 !% suffix in topology filename, to be added to mm_args_str of EVB%pot2
     integer :: form_bond(2)    !atom pair that is bonded in EVB1 only
     integer :: break_bond(2)   !atom pair that is bonded in EVB2 only
!     real(dp) :: energy_offset  !the energy offset of the 2 resonance states
!     type(Constraint) :: off_diagonal

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

