  public :: MetaPotential_ONIOM
  type MetaPotential_ONIOM
     type(Potential), pointer :: pot_region1
     type(Potential), pointer :: pot_region2
     logical :: terminate = .true.
     real(dp):: r_scale_pot1 = 1.0_dp, E_scale_pot1 = 1.0_dp

     type(MetaPotential) :: relax_metapot

     logical       :: minimise_mm 
     character(FIELD_LENGTH) :: minim_mm_method 
     real(dp)      :: minim_mm_tol, minim_mm_eps_guess
     integer       :: minim_mm_max_steps
     character(FIELD_LENGTH) :: minim_mm_linminroutine 
     logical       :: minim_mm_do_pos, minim_mm_do_lat
     logical       :: minim_mm_do_print, minim_mm_use_n_minim
     character(FIELD_LENGTH) :: minim_mm_args_str

      type(Inoutput), pointer :: minim_inoutput_movie
      type(CInoutput), pointer :: minim_cinoutput_movie

  end type MetaPotential_ONIOM

  interface Print
     module procedure metapotential_oniom_print
  end interface Print

  interface Cutoff
     module procedure metapotential_oniom_cutoff
  end interface Cutoff

  interface Calc
     module procedure metapotential_oniom_calc
  end interface Calc

  interface Initialise
     module procedure metapotential_oniom_initialise
  end interface Initialise

  interface Finalise
     module procedure metapotential_oniom_finalise
  end interface Finalise
