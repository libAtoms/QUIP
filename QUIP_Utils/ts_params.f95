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

module tsParams_module

use libAtoms_module
use QUIP_Common_module
use potential_module

implicit none

private

public :: tsParams

type tsParams

   character(STRING_LENGTH) :: classical_args !% Arguments used to initialise classical potential
   character(STRING_LENGTH) :: classical_args_str !% Arguments used by Calc Potential
   real(dp)                 :: classical_force_reweight = 1.0_dp   ! Fraction

   character(STRING_LENGTH) :: qm_args  !% Arguments used to initialise QM potential
   character(STRING_LENGTH) :: qm_args_str  !% Arguments used by QM potential

   logical  :: minim_end 
   character(STRING_LENGTH) :: minim_end_method !% Minimisation method: use 'cg' for conjugate gradients or 'sd' for steepest descent. 
                                            !% See 'minim()' in 'libAtoms/minimisation.f95' for details.
   real(dp) :: minim_end_tol                !% Target force tolerance - geometry optimisation is considered to be 
                                        !% converged when $|\mathbf{f}|^2 <$ 'tol'
   real(dp) :: minim_end_eps_guess          !% Initial guess for line search step size $\epsilon$.
   integer  :: minim_end_max_steps          !% Maximum number of minimisation steps.
   character(STRING_LENGTH) :: minim_end_linminroutine !% Linmin routine, e.g. 'FAST_LINMIN' for classical potentials with total energy, or 
                                                   !% 'LINMIN_DERIV' when doing a LOTF hybrid simulation and only forces are available.
   
   real(dp) :: minim_gfac 
   real(dp) :: minim_force_tol 
   real(dp) :: minim_energy_tol 
   integer  :: minim_max_steps

   integer :: io_verbosity              !% Output verbosity. In XML file, this should be specified as one of
   integer :: io_print_interval 
                                       
   logical  :: simulation_hybrid
   logical  :: simulation_restart 
   character(STRING_LENGTH) :: simulation_method 
   logical  :: simulation_climbing 
   logical  :: simulation_newtangent
   integer  :: simulation_climbing_steps 
   integer  :: simulation_freq_rep   !% Frequency of reparametrization for string method 
   real(dp) :: simulation_spring_constant

   integer  :: chain_nfix 
   integer  :: chain_nimages 
   logical  :: chain_mobile_last 
   logical  :: chain_mobile_first
   character(STRING_LENGTH) :: chain_first_conf 
   character(STRING_LENGTH) :: chain_last_conf 
  
end type tsparams

type(tsParams), pointer :: parse_ts
logical :: parse_in_ts

public :: initialise
interface initialise
   module procedure tsParams_initialise 
end interface

public :: print
interface print
   module procedure tsParams_print
end interface

public :: read_xml
interface read_xml
   module procedure tsParams_read_xml
end interface


contains

subroutine tsParams_initialise(this)
  type(tsParams), intent(inout) :: this

   this%classical_args     = 'IP SW' 
   this%classical_args_str = ' ' 
   this%classical_force_reweight = 1.0_dp   ! Fraction

   ! QM parameters
   this%qm_args                  = 'FilePot command=./castep_driver.py property_list=pos:embed'
!'
   this%qm_args_str              = ' '

   this%minim_end               = .false.
   this%minim_end_method            = 'cg'
   this%minim_end_tol               = 1e-3_dp  ! normsq(force) eV/A
   this%minim_end_eps_guess         = 0.01_dp  ! Angstrom
   this%minim_end_max_steps         = 1000     ! number

   this%io_verbosity            = PRINT_NORMAL
   this%io_print_interval          = 10

   this%simulation_hybrid     = .false. 
   this%simulation_restart    = .false.
   this%simulation_method     = 'neb'
   this%simulation_spring_constant = 0.1 
   this%simulation_newtangent      = .true. 
   this%simulation_climbing        = .false. 
   this%simulation_climbing_steps  = 10 
   this%simulation_freq_rep        = 5

   this%minim_gfac       = 0.1
   this%minim_force_tol  = 0.001_dp
   this%minim_energy_tol = 0.01_dp 
   this%minim_max_steps  = 10 

   this%chain_first_conf  = 'first.xyz'
   this%chain_last_conf   = 'last.xyz'
   this%chain_nfix       = 0
   this%chain_nimages    = 10
   this%chain_mobile_last  = .false.
   this%chain_mobile_first = .false.
end subroutine tsParams_initialise

subroutine tsParams_read_xml(this, xmlfile)
  type(tsParams), intent(inout), target :: this
  type(Inoutput) :: xmlfile

  type (xml_t) :: fxml
  type(extendable_str) :: ss

  call initialise(this) 

  call Initialise(ss)
  call read(ss, xmlfile%unit, convert_to_string=.true.)

  if (len(trim(string(ss))) <= 0) return

  call open_xml_string(fxml, string(ss))

  parse_ts => this
  parse_in_ts = .false.

  call parse(fxml, &
       startElement_handler = tsParams_startElement_handler, &
       endElement_handler = tsParams_endElement_handler)

  call close_xml_t(fxml)

  call Finalise(ss)

end subroutine tsParams_read_xml

subroutine tsParams_startElement_handler(URI, localname, name, attributes)
  character(len=*), intent(in)   :: URI  
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name 
  type(dictionary_t), intent(in) :: attributes

  integer status
  character(len=1024) :: value

  if (name == 'ts_params') then ! new crack_params stanza

     parse_in_ts = .true.

 elseif (parse_in_ts .and. name == 'chain') then
     call QUIP_FoX_get_value(attributes, "first_conf", value, status)
     if (status == 0) then
        parse_ts%chain_first_conf = value
     end if

     call QUIP_FoX_get_value(attributes, "last_conf", value, status)
     if (status == 0) then
        parse_ts%chain_last_conf = value
     end if

     call QUIP_FoX_get_value(attributes, "nfix", value, status)
     if (status == 0) then
        read (value, *) parse_ts%chain_nfix
     end if

     call QUIP_FoX_get_value(attributes, "nimages", value, status)
     if (status == 0) then
        read (value, *) parse_ts%chain_nimages
     end if

     call QUIP_FoX_get_value(attributes, "mobile_last", value, status)
     if (status == 0) then
        read (value, *)  parse_ts%chain_mobile_last 
     end if

     call QUIP_FoX_get_value(attributes, "mobile_first", value, status)
     if (status == 0) then
        read (value, *)  parse_ts%chain_mobile_first
     end if

 elseif (parse_in_ts .and. name == 'simulation') then

     call QUIP_FoX_get_value(attributes, "hybrid", value, status)
     if (status == 0) then
       read (value, *) parse_ts%simulation_hybrid
     end if
  
     call QUIP_FoX_get_value(attributes, "restart", value, status)
     if (status == 0) then
        read (value,*) parse_ts%simulation_restart 
     end if

     call QUIP_FoX_get_value(attributes, "method", value, status)
     if (status == 0) then
        parse_ts%simulation_method = value
     end if

     call QUIP_FoX_get_value(attributes, "newtangent", value, status)
     if (status == 0) then
        read(value,*) parse_ts%simulation_newtangent
     end if

     call QUIP_FoX_get_value(attributes, "climbing", value, status)
     if (status == 0) then
        read(value,*) parse_ts%simulation_climbing 
     end if

     call QUIP_FoX_get_value(attributes, "climbing_steps", value, status)
     if (status == 0) then
        read(value,*) parse_ts%simulation_climbing_steps 
     end if

     call QUIP_FoX_get_value(attributes, "freq_rep", value, status)
     if (status == 0) then
        read (value, *) parse_ts%simulation_freq_rep
     end if

     call QUIP_FoX_get_value(attributes, "spring_constant", value, status)
     if (status == 0) then
        read (value, *) parse_ts%simulation_spring_constant
     end if

 elseif (parse_in_ts .and. name == 'io') then

     call QUIP_FoX_get_value(attributes, "verbosity", value, status)
     if (status == 0) then
       parse_ts%io_verbosity = verbosity_of_str(value)
     end if

     call QUIP_FoX_get_value(attributes, "print_interval", value, status)
     if (status == 0) then
        read (value, *) parse_ts%io_print_interval
     end if

 elseif (parse_in_ts .and. name == 'qm') then

    call QUIP_FoX_get_value(attributes, "args", value, status)
    if (status == 0) then
       parse_ts%qm_args = value
    end if

    call QUIP_FoX_get_value(attributes, "args_str", value, status)
    if (status == 0) then
       parse_ts%qm_args_str = value
    end if

 elseif (parse_in_ts .and. name == 'classical') then

    call QUIP_FoX_get_value(attributes, "args", value, status)
    if (status == 0) then
       parse_ts%classical_args = value
    end if

    call QUIP_FoX_get_value(attributes, "args_str", value, status)
    if (status == 0) then
       parse_ts%classical_args_str = value
    end if

    call QUIP_FoX_get_value(attributes, "force_reweight", value, status)
    if (status == 0) then
       read (value, *) parse_ts%classical_force_reweight
    end if

 elseif (parse_in_ts .and. name == 'minim') then
      print *, 'AAA'
   

     call QUIP_FoX_get_value(attributes, "gfac", value, status)
     if (status == 0) then
        read (value, *) parse_ts%minim_gfac
     end if

     call QUIP_FoX_get_value(attributes, "force_tol", value, status)
     if (status == 0) then
        read (value, *) parse_ts%minim_force_tol
     end if

     call QUIP_FoX_get_value(attributes, "energy_tol", value, status)
     if (status == 0) then
        read (value, *) parse_ts%minim_energy_tol
     end if

     call QUIP_FoX_get_value(attributes, "max_steps", value, status)
     if (status == 0) then
        read (value, *) parse_ts%minim_max_steps
     end if

 elseif (parse_in_ts .and. name == 'minim_end') then

      call QUIP_FoX_get_value(attributes, "optimise_end", value, status)
      if (status == 0) then
         read(value,*) parse_ts%minim_end
      end if

      call QUIP_FoX_get_value(attributes, "method", value, status)
      if (status == 0) then
         parse_ts%minim_end_method = value
      end if

      call QUIP_FoX_get_value(attributes, "tol", value, status)
      if (status == 0) then     
         read (value, *) parse_ts%minim_end_tol
      end if

      call QUIP_FoX_get_value(attributes, "eps_guess", value, status)
      if (status == 0) then     
         read (value, *) parse_ts%minim_end_eps_guess
      end if

      call QUIP_FoX_get_value(attributes, "max_steps", value, status)
      if (status == 0) then
         read (value, *) parse_ts%minim_end_max_steps
      end if

      call QUIP_FoX_get_value(attributes, "linminroutine", value, status)
      if (status == 0) then
         parse_ts%minim_end_linminroutine = value 
      end if

  endif

end subroutine tsParams_startElement_handler

subroutine tsParams_endElement_handler(URI, localname, name)
  character(len=*), intent(in)   :: URI  
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name 

  if (parse_in_ts) then
     if (name == 'ts_params') then
        parse_in_ts = .false.
     endif
  endif

end subroutine tsParams_endElement_handler


subroutine tsParams_print(this,file)
  type(tsParams), intent(in) :: this
  type(Inoutput), optional, intent(in) :: file

end subroutine tsParams_print

end module tsParams_module
