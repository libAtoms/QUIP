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

!X
!X IP module 
!X
!% General object which manages all the possible interatomic potentials (IP), addressing
!% the calls to the right modules. The available IPs are the following:
!%   \begin{itemize}
!%    \item    Gaussian Approximation Potential           ({\bf IPModel_GAP}) 
!%    \item    Lennard-Jones Potential                    ({\bf IPModel_LJ}) 
!%    \item    Morse Potential                            ({\bf IPModel_Morse}) 
!%    \item    Force-Constant Potential                   ({\bf IPModel_FC}) 
!%    \item    Stillinger and Weber Potential for Silicon ({\bf IPModel_SW}) 
!%    \item    Tersoff potential for C/Si/Ge              ({\bf IPModel_Tersoff})  
!%    \item    Embedded Atom Potential of Ercolessi Adam  ({\bf IPModel_EAM_ErcolAd})
!%    \item    Brenner Potential                          ({\bf IPModel_Brenner})
!%    \item    FB Potential for SiO2                      ({\bf IPModel_FB})
!%    \item    Silicon Modified Embedded Atom Potentia l  ({\bf IPModel_Si_MEAM})
!%    \item    Finnis-Sinclair                            ({\bf IPModel_FS})
!%    \item    Bond Order Potential                       ({\bf IPModel_BOP})
!%    \item    Screened Brenner Potential (interface)     ({\bf IPModel_Brenner_Screened})
!%    \item    2nd generation Brenner Potential (interface) ({\bf IPModel_Brenner_2002})
!%    \item    Template potential ({\bf IPModel_Template})
!%   \end{itemize}
!%  The IP_type object contains details regarding the selected IP.
!%  When a type Potential is defined
!%>   type(Potential) :: pot
!%  it is then necessary to initialise the IP parameters (readable from an external input file or 
!%  from an internal string) and to define the IP type:
!%>   call Initialise(pot, IP_type, params) 
!%  where IP_type can be 
!%   \begin{itemize}
!%    \item    'IP GAP'
!%    \item    'IP LJ'
!%    \item    'IP Morse'
!%    \item    'IP FC'
!%    \item    'IP SW' 
!%    \item    'IP Tersoff'
!%    \item    'IP EAM_ErcolAd'
!%    \item    'IP Brenner'
!%    \item    'IP FB'
!%    \item    'IP Si_MEAM'
!%    \item    'IP FS'
!%    \item    'IP BOP'
!%    \item    'IP Brenner_Screened'
!%    \item    'IP_Brenner_2002'
!%    \item    'IP_Template'
!%   \end{itemize}
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#include "error.inc"
module IP_module

use libatoms_module

use MPI_context_module
use IPModel_GAP_module
use IPModel_LJ_module
use IPModel_Morse_module
use IPModel_FC_module
use IPModel_SW_module
use IPModel_Tersoff_module
use IPModel_EAM_ErcolAd_module
use IPModel_Brenner_module
use IPModel_FB_module
use IPModel_Si_MEAM_module
use IPModel_FS_module
use IPModel_BOP_module
use IPModel_Brenner_Screened_module
use IPModel_Brenner_2002_module
use IPModel_ASAP_module
use IPModel_ASAP2_module
use IPModel_Glue_module
! Add new IPs here
use IPModel_Template_module
use QUIP_Common_module

implicit none

private

integer, parameter :: FF_LJ = 1, FF_SW = 2, FF_Tersoff = 3, FF_EAM_ErcolAd = 4, &
     FF_Brenner = 5, FF_GAP = 6, FF_FS = 7, FF_BOP = 8, FF_FB = 9, FF_Si_MEAM = 10, FF_Brenner_Screened = 11, &
     FF_Brenner_2002 = 12, FF_ASAP = 13, FF_ASAP2 = 14, FF_FC = 15, FF_Morse = 16, FF_GLUE = 17, & ! Add new IPs here
     FF_Template = 99

public :: IP_type
type IP_type
  integer :: functional_form = 0

  type(IPModel_GAP) ip_gap
  type(IPModel_LJ) ip_lj
  type(IPModel_Morse) ip_morse
  type(IPModel_FC) ip_fc
  type(IPModel_SW) ip_sw
  type(IPModel_Tersoff) ip_Tersoff
  type(IPModel_EAM_ErcolAd) ip_EAM_ErcolAd
  type(IPModel_Brenner) ip_Brenner
  type(IPModel_FB) ip_FB
  type(IPModel_Si_MEAM) ip_Si_MEAM
  type(IPModel_FS) ip_fs
  type(IPModel_BOP) ip_BOP
  type(IPModel_Brenner_Screened) ip_Brenner_Screened
  type(IPModel_Brenner_2002) ip_Brenner_2002
  type(IPModel_ASAP) ip_ASAP
  type(IPModel_ASAP2) ip_ASAP2
  type(IPModel_Glue) ip_Glue
     ! Add new IPs here  
  type(IPModel_Template) ip_Template

  type(mpi_context) :: mpi_glob

end type IP_type

!% Initialise IP_type object defining the IP model and the corresponding parameters. If necessary
!% it initialises the input file for the potential parameters. 
public :: Initialise
public :: IP_Initialise_filename
interface Initialise
  module procedure IP_Initialise_inoutput, IP_Initialise_str
end interface Initialise

public :: Finalise
interface Finalise
  module procedure IP_Finalise
end interface Finalise

!% Return the cutoff of this interatomic potential
public :: cutoff
interface cutoff
   module procedure IP_cutoff
end interface

!% Print the parameters of the selected IP model
public :: Print
interface Print
  module procedure IP_Print
end interface Print

!% Hook routine called before calc() is invoked. Can be used to add required properties.
public :: setup_atoms
interface setup_atoms
   module procedure ip_setup_atoms
end interface

!% Call the potential calculator for the selected IP model
public :: Calc
private :: IP_Calc
interface Calc
  module procedure IP_Calc
end interface Calc

!% set up optimal parameters for parallel computation for a given system
public :: setup_parallel
private :: IP_setup_parallel
interface setup_parallel
  module procedure IP_setup_parallel
end interface setup_parallel

contains

!% OMIT
subroutine IP_Initialise_filename(this, args_str, filename, mpi_obj, error)
  type(IP_type), intent(inout) :: this
  character(len=*), intent(in) :: args_str
  character(len=*), intent(in) :: filename  !% File name containing the IP parameters
  type(MPI_context), intent(in), optional :: mpi_obj
  integer, intent(out), optional :: error

  type(inoutput) io

  INIT_ERROR(error)

  call Initialise(io, filename, INPUT)

  call Initialise(this, args_str, io, mpi_obj)

  call Finalise(io)

end subroutine


subroutine IP_Initialise_inoutput(this, args_str, io_obj, mpi_obj, error)
  type(IP_type), intent(inout) :: this
  character(len=*), intent(in) :: args_str
  type(inoutput), intent(inout), optional :: io_obj
  type(MPI_context), intent(in), optional :: mpi_obj
  integer, intent(out), optional :: error

  type(extendable_str) :: ss

  INIT_ERROR(error)

  call Initialise(ss)
  if (present(io_obj)) then
    if (present(mpi_obj)) then
      call read(ss, io_obj%unit, convert_to_string=.true., mpi_comm=mpi_obj%communicator)
    else
      call read(ss, io_obj%unit, convert_to_string=.true.)
    endif
  endif
  call Initialise(this, args_str, string(ss), mpi_obj)
  call Finalise(ss)

end subroutine IP_Initialise_inoutput

subroutine IP_Initialise_str(this, args_str, param_str, mpi_obj, error)
  type(IP_type), intent(inout) :: this
  character(len=*), intent(in) :: args_str, param_str
  type(MPI_context), intent(in), optional :: mpi_obj
  integer, intent(out), optional :: error

  type(Dictionary) :: params
  logical is_GAP, is_LJ, is_FC, is_Morse, is_SW, is_Tersoff, is_EAM_ErcolAd, is_Brenner, is_FS, is_BOP, is_FB, is_Si_MEAM, &
       is_Brenner_Screened, is_Brenner_2002, is_ASAP, is_ASAP2, is_Glue, is_template

  INIT_ERROR(error)

  call Finalise(this)

  call initialise(params)
  call param_register(params, 'GAP', 'false', is_GAP)
  call param_register(params, 'LJ', 'false', is_LJ)
  call param_register(params, 'Morse', 'false', is_Morse)
  call param_register(params, 'FC', 'false', is_FC)
  call param_register(params, 'SW', 'false', is_SW)
  call param_register(params, 'Tersoff', 'false', is_Tersoff)
  call param_register(params, 'EAM_ErcolAd', 'false', is_EAM_ErcolAd)
  call param_register(params, 'Brenner', 'false', is_Brenner)
  call param_register(params, 'FB', 'false', is_FB)
  call param_register(params, 'Si_MEAM', 'false', is_Si_MEAM)
  call param_register(params, 'FS', 'false', is_FS)
  call param_register(params, 'BOP', 'false', is_BOP)
  call param_register(params, 'Brenner_Screened', 'false', is_Brenner_Screened)
  call param_register(params, 'Brenner_2002', 'false', is_Brenner_2002)
  call param_register(params, 'ASAP', 'false', is_ASAP)
  call param_register(params, 'ASAP2', 'false', is_ASAP2)
  call param_register(params, 'Glue', 'false', is_Glue)
  ! Add new IPs here
  call param_register(params, 'Template', 'false', is_template)

  if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='IP_Initialise_str args_str')) then
    RAISE_ERROR("IP_Initialise_str failed to parse args_str='"//trim(args_str)//"'", error)
  endif
  call finalise(params)

  if (count((/is_GAP, is_LJ, is_FC, is_Morse, is_SW, is_Tersoff, is_EAM_ErcolAd, is_Brenner, is_FS, is_BOP, is_FB, is_Si_MEAM, &
       is_Brenner_Screened, is_Brenner_2002, is_ASAP, is_ASAP2, is_Glue,  &        ! add new IPs here
       is_Template /)) /= 1) then
    RAISE_ERROR("IP_Initialise_str found too few or too many IP Model types args_str='"//trim(args_str)//"'", error)
  endif

  if (is_GAP) then
    this%functional_form = FF_GAP
    call Initialise(this%ip_gap, args_str, param_str)
  else if (is_LJ) then
    this%functional_form = FF_LJ
    call Initialise(this%ip_lj, args_str, param_str)
  else if (is_Morse) then
    this%functional_form = FF_Morse
    call Initialise(this%ip_morse, args_str, param_str)
  else if (is_FC) then
    this%functional_form = FF_FC
    call Initialise(this%ip_fc, args_str, param_str)
  else if (is_SW) then
    this%functional_form = FF_SW
    call Initialise(this%ip_sw, args_str, param_str)
  else if (is_Tersoff) then
    this%functional_form = FF_Tersoff
    call Initialise(this%ip_tersoff, args_str, param_str)
  else if (is_EAM_ErcolAd) then
    this%functional_form = FF_EAM_ErcolAd
    call Initialise(this%ip_EAM_ErcolAd, args_str, param_str)
  else if (is_Brenner) then
    this%functional_form = FF_Brenner
    call Initialise(this%ip_Brenner, args_str, param_str)
  else if (is_FB) then
    this%functional_form = FF_FB
    call Initialise(this%ip_FB, args_str, param_str)
  else if (is_Si_MEAM) then
    this%functional_form = FF_Si_MEAM
    call Initialise(this%ip_Si_MEAM, args_str, param_str)
  else if (is_FS) then
    this%functional_form = FF_FS
    call Initialise(this%ip_fs, args_str, param_str)
  else if (is_BOP) then
    this%functional_form = FF_BOP
    call Initialise(this%ip_bop, args_str, param_str)
  else if (is_Brenner_Screened) then
    this%functional_form = FF_Brenner_Screened
    call Initialise(this%ip_brenner_screened, args_str, param_str)
  else if (is_Brenner_2002) then
    this%functional_form = FF_Brenner_2002
    call Initialise(this%ip_brenner_2002, args_str, param_str)
  else if (is_ASAP) then
    this%functional_form = FF_ASAP
    call Initialise(this%ip_ASAP, args_str, param_str) 
  else if (is_ASAP2) then
    this%functional_form = FF_ASAP2
    call Initialise(this%ip_ASAP2, args_str, param_str) 
  else if (is_Glue) then
    this%functional_form = FF_GLUE
    call Initialise(this%ip_Glue, args_str, param_str) 
    ! add new IPs here
  else if (is_Template) then
    this%functional_form = FF_Template
    call Initialise(this%ip_Template, args_str, param_str) 
  end if

  if (present(mpi_obj)) this%mpi_glob = mpi_obj

end subroutine IP_Initialise_str

subroutine IP_Finalise(this)
  type(IP_type), intent(inout) :: this

  select case (this%functional_form)
    case (FF_GAP)
      if (this%ip_gap%mpi%active) call free_context(this%ip_gap%mpi)
      call Finalise(this%ip_gap)
    case (FF_LJ)
      if (this%ip_lj%mpi%active) call free_context(this%ip_lj%mpi)
      call Finalise(this%ip_lj)
    case (FF_Morse)
      if (this%ip_morse%mpi%active) call free_context(this%ip_morse%mpi)
      call Finalise(this%ip_morse)
    case (FF_FC)
      if (this%ip_fc%mpi%active) call free_context(this%ip_fc%mpi)
      call Finalise(this%ip_fc)
    case (FF_SW)
      if (this%ip_sw%mpi%active) call free_context(this%ip_sw%mpi)
      call Finalise(this%ip_sw)
    case (FF_Tersoff)
      if (this%ip_tersoff%mpi%active) call free_context(this%ip_tersoff%mpi)
      call Finalise(this%ip_tersoff)
    case (FF_EAM_ErcolAd)
      if (this%ip_EAM_ErcolAd%mpi%active) call free_context(this%ip_EAM_ErcolAd%mpi)
      call Finalise(this%ip_EAM_ErcolAd)
    case(FF_Brenner)
      if (this%ip_brenner%mpi%active) call free_context(this%ip_brenner%mpi)
      call Finalise(this%ip_Brenner)
    case(FF_FB)
      if (this%ip_fb%mpi%active) call free_context(this%ip_fb%mpi)
      call Finalise(this%ip_fb)
    case(FF_Si_MEAM)
      if (this%ip_Si_MEAM%mpi%active) call free_context(this%ip_Si_MEAM%mpi)
      call Finalise(this%ip_Si_MEAM)
    case (FF_FS)
      if (this%ip_fs%mpi%active) call free_context(this%ip_fs%mpi)
      call Finalise(this%ip_fs)
    case (FF_BOP)
      if (this%ip_bop%mpi%active) call free_context(this%ip_bop%mpi)
      call Finalise(this%ip_bop)
    case (FF_Brenner_Screened)
      if(this%ip_brenner_screened%mpi%active) call free_context(this%ip_brenner_screened%mpi)
      call Finalise(this%ip_brenner_screened)
    case (FF_Brenner_2002)
      if(this%ip_brenner_2002%mpi%active) call free_context(this%ip_brenner_2002%mpi)
      call Finalise(this%ip_brenner_2002)
    case (FF_ASAP)
      if(this%ip_ASAP%mpi%active) call free_context(this%ip_ASAP%mpi)
      call Finalise(this%ip_ASAP)
    case (FF_ASAP2)
      if(this%ip_ASAP2%mpi%active) call free_context(this%ip_ASAP2%mpi)
      call Finalise(this%ip_ASAP2)
    case (FF_GLUE)
      if(this%ip_Glue%mpi%active) call free_context(this%ip_Glue%mpi)
      call Finalise(this%ip_Glue)
      ! add new IP here
    case (FF_Template)
      if(this%ip_template%mpi%active) call free_context(this%ip_template%mpi)
      call Finalise(this%ip_template)
  end select

end subroutine IP_Finalise

function IP_cutoff(this)
  type(IP_type), intent(in) :: this
  real(dp) :: IP_cutoff

  select case (this%functional_form)
  case (FF_GAP)
     IP_cutoff = this%ip_gap%cutoff
  case (FF_LJ)
     IP_cutoff = this%ip_lj%cutoff
  case (FF_Morse)
     IP_cutoff = this%ip_morse%cutoff
  case (FF_FC)
     IP_cutoff = this%ip_fc%cutoff
  case (FF_SW)
     IP_cutoff = this%ip_sw%cutoff
  case (FF_Tersoff)
     IP_cutoff = this%ip_tersoff%cutoff
  case (FF_EAM_ErcolAd)
     IP_cutoff = this%ip_EAM_ErcolAd%cutoff
  case(FF_Brenner)
     IP_cutoff = this%ip_Brenner%cutoff
  case(FF_FB)
     IP_cutoff = this%ip_FB%cutoff
  case(FF_Si_MEAM)
     IP_cutoff = this%ip_Si_MEAM%cutoff
  case (FF_FS)
     IP_cutoff = this%ip_fs%cutoff
  case (FF_BOP)
     IP_cutoff = this%ip_bop%cutoff
  case (FF_Brenner_Screened)
     IP_cutoff = this%ip_brenner_screened%cutoff
  case (FF_Brenner_2002)
     IP_cutoff = this%ip_brenner_2002%cutoff
  case (FF_ASAP)
     IP_cutoff = maxval(this%ip_asap%cutoff)*BOHR
  case (FF_ASAP2)
     IP_cutoff = max(this%ip_asap2%cutoff_ms, this%ip_asap2%cutoff_coulomb)*BOHR
  case (FF_GLUE)
     IP_cutoff = this%ip_Glue%cutoff
  ! Add new IP here
  case (FF_Template)
     IP_cutoff = this%ip_template%cutoff
  case default
     IP_cutoff = 0.0_dp
  end select
end function IP_cutoff

subroutine IP_setup_atoms(this, at)
  type(IP_type), intent(in) :: this
  type(Atoms), intent(inout) :: at

  select case (this%functional_form)
  case (FF_ASAP2)
     call setup_atoms(this%ip_asap2, at)
  end select

end subroutine IP_setup_atoms

subroutine IP_Calc(this, at, energy, local_e, f, virial, args_str, error)
  type(IP_type), intent(inout) :: this
  type(Atoms), intent(inout) :: at                
  real(dp), intent(out), optional :: energy, local_e(:) !% \texttt{energy} = System total energy, \texttt{local_e} = energy of each atom, vector dimensioned as \texttt{at%N}.
  real(dp), intent(out), optional :: f(:,:)
  real(dp), intent(out), optional :: virial(3,3)
  character(len=*), intent(in), optional      :: args_str 
  integer, intent(out), optional :: error

  logical mpi_active

  INIT_ERROR(error)

  call system_timer("IP_Calc")

  select case (this%functional_form)
    case (FF_GAP)
      mpi_active = this%ip_gap%mpi%active
    case (FF_LJ)
      mpi_active = this%ip_lj%mpi%active
    case (FF_Morse)
      mpi_active = this%ip_morse%mpi%active
    case (FF_FC)
      mpi_active = this%ip_fc%mpi%active
    case (FF_SW)
      mpi_active = this%ip_sw%mpi%active
    case (FF_Tersoff)
      mpi_active = this%ip_tersoff%mpi%active
    case (FF_EAM_ErcolAd)
      mpi_active = this%ip_EAM_ErcolAd%mpi%active
    case(FF_Brenner)
      mpi_active = this%ip_Brenner%mpi%active
    case(FF_FB)
      mpi_active = this%ip_FB%mpi%active
    case(FF_Si_MEAM)
      mpi_active = this%ip_Si_MEAM%mpi%active
    case(FF_FS)
      mpi_active = this%ip_fs%mpi%active
    case(FF_BOP)
      mpi_active = this%ip_bop%mpi%active
    case(FF_Brenner_Screened)
      mpi_active = this%ip_brenner_screened%mpi%active
    case(FF_Brenner_2002)
      mpi_active = this%ip_brenner_2002%mpi%active
    case(FF_ASAP)
      mpi_active = this%ip_asap%mpi%active
    case(FF_ASAP2)
      mpi_active = this%ip_asap2%mpi%active
    case(FF_GLUE)
      mpi_active = this%ip_Glue%mpi%active
    ! add new IP here
    case(FF_Template)
      mpi_active = this%ip_template%mpi%active
    case default
      RAISE_ERROR("IP_Calc confused by functional_form " // this%functional_form, error)
  end select

  if (this%mpi_glob%active .and. .not. mpi_active) then
    call setup_parallel(this, at)
  endif

  select case (this%functional_form)
    case (FF_GAP)
      call calc(this%ip_gap, at, energy, local_e, f, virial,args_str)
    case (FF_LJ)
      call calc(this%ip_lj, at, energy, local_e, f, virial, args_str, error=error)
      PASS_ERROR(error)
    case (FF_Morse)
      call calc(this%ip_morse, at, energy, local_e, f, virial, args_str)
    case (FF_FC)
      call calc(this%ip_fc, at, energy, local_e, f, virial, args_str)
    case (FF_SW)
      call calc(this%ip_sw, at, energy, local_e, f, virial)
    case (FF_Tersoff)
      call calc(this%ip_tersoff, at, energy, local_e, f, virial, args_str)
    case (FF_EAM_ErcolAd)
      call calc(this%ip_EAM_ErcolAd, at, energy, local_e, f, virial)
    case(FF_Brenner)
      call calc(this%ip_Brenner, at, energy, local_e, f, virial)
    case(FF_FB)
      call calc(this%ip_FB, at, energy, local_e, f, virial)
    case(FF_Si_MEAM)
      call calc(this%ip_Si_MEAM, at, energy, local_e, f, virial)
    case (FF_FS)
      call calc(this%ip_fs, at, energy, local_e, f, virial)
    case (FF_BOP)
      call calc(this%ip_bop, at, energy, local_e, f, virial, args_str)
    case (FF_Brenner_Screened)
      call calc(this%ip_brenner_screened, at, energy, local_e, f, virial, args_str)
    case (FF_Brenner_2002)
      call calc(this%ip_brenner_2002, at, energy, local_e, f, virial, args_str)
    case (FF_ASAP)
      call calc(this%ip_asap, at, energy, local_e, f, virial, args_str)
    case (FF_ASAP2)
      call calc(this%ip_asap2, at, energy, local_e, f, virial, args_str)
    case (FF_GLUE)
      call calc(this%ip_Glue, at, energy, local_e, f, virial, args_str)
    ! add new IP here
    case (FF_Template)
      call calc(this%ip_template, at, energy, local_e, f, virial, args_str)
    case default
      RAISE_ERROR("IP_Calc confused by functional_form " // this%functional_form, error)
  end select

  call system_timer("IP_Calc")
end subroutine IP_Calc


subroutine IP_Print(this, file, error)
  type(IP_type), intent(inout) :: this
  type(Inoutput), intent(inout),optional :: file
  integer, intent(out), optional :: error

  INIT_ERROR(error)

  call Print ("IP : " // this%functional_form, file=file)

  select case (this%functional_form)
    case (FF_GAP)
      call Print(this%ip_gap, file=file)
    case (FF_LJ)
      call Print(this%ip_lj, file=file)
    case (FF_Morse)
      call Print(this%ip_morse, file=file)
    case (FF_FC)
      call Print(this%ip_fc, file=file)
    case (FF_SW)
      call Print(this%ip_sw, file=file)
    case (FF_Tersoff)
      call Print(this%ip_tersoff, file=file)
    case (FF_EAM_ErcolAd)
      call Print(this%ip_EAM_ErcolAd, file=file)
    case (FF_Brenner)
      call Print(this%ip_Brenner, file=file)
    case (FF_FB)
      call Print(this%ip_FB, file=file)
    case (FF_Si_MEAM)
      call Print(this%ip_Si_MEAM, file=file)
    case (FF_FS)
      call Print(this%ip_fs, file=file)
    case (FF_BOP)
      call Print(this%ip_bop, file=file)
    case (FF_Brenner_Screened)
      call Print(this%ip_brenner_screened, file=file)
    case (FF_Brenner_2002)
      call Print(this%ip_brenner_2002, file=file)
    case (FF_ASAP)
      call Print(this%ip_asap, file=file)
    case (FF_ASAP2)
      call Print(this%ip_asap2, file=file)
    case (FF_GLUE)
      call Print(this%ip_Glue, file=file)
    ! add new IP here
    case (FF_Template)
      call Print(this%ip_template, file=file)
    case default
      RAISE_ERROR("IP_Print confused by functional_form " // this%functional_form, error)
  end select

end subroutine IP_Print

subroutine IP_setup_parallel(this, at, energy, local_e, f, virial, args_str)
  type(IP_type), intent(inout) :: this
  type(Atoms), intent(inout) :: at
  real(dp), intent(out), optional :: energy, local_e(:) !% \texttt{energy} = System total energy, \texttt{local_e} = energy of each atom, vector dimensioned as \texttt{at%N}.  
  real(dp), intent(out), optional :: f(:,:)
  real(dp), intent(out), optional :: virial(3,3)
  character(len=STRING_LENGTH), intent(in), optional      :: args_str 


  integer :: pgroup_size, prev_pgroup_size
  integer :: n_groups
  real(dp) :: prev_time, this_time

  prev_time = 1.0e38_dp
  prev_pgroup_size = 0

  call print('setup_parallel timings', PRINT_VERBOSE)
  call print('group_size  time/sec', PRINT_VERBOSE)
  
  call setup_atoms(this, at)
  do pgroup_size=1, this%mpi_glob%n_procs
    n_groups = this%mpi_glob%n_procs / pgroup_size
    if (n_groups*pgroup_size == this%mpi_glob%n_procs) then
      call setup_parallel_groups(this, this%mpi_glob, pgroup_size)
      call system_timer("IP_parallel", do_always = .true.)
      call calc(this, at, energy, local_e, f, virial, args_str)
      call system_timer("IP_parallel", do_always = .true., time_elapsed = this_time)
      this_time = max(this%mpi_glob, this_time)
      call print(pgroup_size//' '//this_time, PRINT_VERBOSE)
      if (this_time > prev_time) then
	call setup_parallel_groups(this, this%mpi_glob, prev_pgroup_size)
	exit
      else
	prev_time = this_time
	prev_pgroup_size = pgroup_size
      endif
    endif
  end do

  call print("Parallelizing IP using group_size " // prev_pgroup_size, PRINT_ALWAYS)
end subroutine IP_setup_parallel

subroutine setup_parallel_groups(this, mpi, pgroup_size, error)
  type(IP_type), intent(inout) :: this
  type(mpi_context), intent(in) :: mpi
  integer, intent(in) :: pgroup_size
  integer, intent(out), optional :: error

  type(mpi_context) :: mpi_local
  integer :: split_index

  INIT_ERROR(error)

  if (mpi%active) then
    split_index = mpi%my_proc/pgroup_size
    call split_context(mpi, split_index, mpi_local)
  endif

  select case (this%functional_form)
    case (FF_GAP)
      if (this%ip_gap%mpi%active) call free_context(this%ip_gap%mpi)
      this%ip_gap%mpi = mpi_local
    case (FF_LJ)
      if (this%ip_lj%mpi%active) call free_context(this%ip_lj%mpi)
      this%ip_lj%mpi = mpi_local
    case (FF_Morse)
      if (this%ip_morse%mpi%active) call free_context(this%ip_morse%mpi)
      this%ip_morse%mpi = mpi_local
    case (FF_FC)
      if (this%ip_fc%mpi%active) call free_context(this%ip_fc%mpi)
      this%ip_fc%mpi = mpi_local
    case (FF_SW)
      if (this%ip_sw%mpi%active) call free_context(this%ip_sw%mpi)
      this%ip_sw%mpi = mpi_local
    case (FF_Tersoff)
      if (this%ip_tersoff%mpi%active) call free_context(this%ip_tersoff%mpi)
      this%ip_tersoff%mpi = mpi_local
    case (FF_EAM_ErcolAd)
      if (this%ip_EAM_ErcolAd%mpi%active) call free_context(this%ip_EAM_ErcolAd%mpi)
      this%ip_EAM_ErcolAd%mpi = mpi_local
    case(FF_Brenner)
      if (this%ip_brenner%mpi%active) call free_context(this%ip_brenner%mpi)
      this%ip_Brenner%mpi = mpi_local
    case(FF_FB)
      if (this%ip_FB%mpi%active) call free_context(this%ip_FB%mpi)
      this%ip_FB%mpi = mpi_local
    case(FF_Si_MEAM)
      if (this%ip_Si_MEAM%mpi%active) call free_context(this%ip_Si_MEAM%mpi)
      this%ip_Si_MEAM%mpi = mpi_local
    case(FF_FS)
      if (this%ip_fs%mpi%active) call free_context(this%ip_fs%mpi)
      this%ip_fs%mpi = mpi_local
    case(FF_BOP)
      if (this%ip_bop%mpi%active) call free_context(this%ip_bop%mpi)
      this%ip_bop%mpi = mpi_local
    case(FF_Brenner_Screened)
      if (this%ip_brenner_screened%mpi%active) call free_context(this%ip_brenner_screened%mpi)
      this%ip_brenner_screened%mpi = mpi_local
    case(FF_Brenner_2002)
      if (this%ip_brenner_2002%mpi%active) call free_context(this%ip_brenner_2002%mpi)
      this%ip_brenner_2002%mpi = mpi_local
    case(FF_ASAP)
      if (this%ip_asap%mpi%active) call free_context(this%ip_asap%mpi)
      this%ip_asap%mpi = mpi_local
    case(FF_ASAP2)
      if (this%ip_asap2%mpi%active) call free_context(this%ip_asap2%mpi)
      this%ip_asap2%mpi = mpi_local
    case(FF_GLUE)
      if (this%ip_Glue%mpi%active) call free_context(this%ip_Glue%mpi)
      this%ip_Glue%mpi = mpi_local
    ! add new IP here
    case(FF_Template)
      if (this%ip_template%mpi%active) call free_context(this%ip_template%mpi)
      this%ip_template%mpi = mpi_local
    case default
      RAISE_ERROR("setup_parallel_groups confused by functional_form " // this%functional_form, error)
  end select

end subroutine setup_parallel_groups

end module IP_module
