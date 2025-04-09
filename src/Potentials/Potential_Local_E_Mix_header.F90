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

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X Local Energy Mixing header stuff to be included in Potential.f95
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


  public :: Potential_Local_E_Mix
  type Potential_Local_E_Mix
     type(Potential), pointer :: pot_region1
     type(Potential), pointer :: pot_region2
     logical :: terminate = .true.
     real(dp):: r_scale_pot1 = 1.0_dp, E_scale_pot1 = 1.0_dp

     type(Potential) :: relax_pot

     character(STRING_LENGTH) :: run_suffix
     logical       :: minimise_mm 
     character(STRING_LENGTH) :: minim_mm_method 
     real(dp)      :: minim_mm_tol, minim_mm_eps_guess
     integer       :: minim_mm_max_steps
     character(STRING_LENGTH) :: minim_mm_linminroutine 
     logical       :: minim_mm_do_pos, minim_mm_do_lat
     logical       :: minim_mm_do_print
     character(STRING_LENGTH) :: minim_mm_args_str

     type(Dictionary) :: create_hybrid_weights_params

     type(Inoutput), pointer :: minim_inoutput_movie
     type(CInoutput), pointer :: minim_cinoutput_movie

  end type Potential_Local_E_Mix

  interface Print
     module procedure potential_local_e_mix_print
  end interface Print

  interface Cutoff
     module procedure potential_local_e_mix_cutoff
  end interface Cutoff

  interface Calc
     module procedure potential_local_e_mix_calc
  end interface Calc

  interface Initialise
     module procedure potential_local_e_mix_initialise
  end interface Initialise

  interface Finalise
     module procedure potential_local_e_mix_finalise
  end interface Finalise
