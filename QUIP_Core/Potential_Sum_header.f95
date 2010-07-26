
  public :: Potential_Sum
  type Potential_Sum
     type(MPI_context) :: mpi

     type(Potential), pointer :: pot1 => null() 
     type(Potential), pointer :: pot2 => null() 

     logical  :: subtract_pot1
     logical  :: subtract_pot2
  end type Potential_Sum

  interface Initialise
     module procedure Potential_Sum_Initialise
  end interface

  interface Finalise
     module procedure Potential_Sum_Finalise
  end interface

  interface Print
     module procedure Potential_Sum_Print
  end interface

  interface Cutoff
     module procedure Potential_Sum_Cutoff
  end interface

  interface Calc
     module procedure Potential_Sum_Calc
  end interface

