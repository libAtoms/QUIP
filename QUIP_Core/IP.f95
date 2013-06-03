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
!%    \item    Partridge-Schwenke potential               ({\bf IPModel_PartridgeSchwenke})
!%    \item    Einstein crystal potential ({\bf IPModel_Einstein})
!%    \item    Coulomb potential ({\bf IPModel_Coulomb})
!%    \item    Sutton-Chen potential ({\bf IPModel_Sutton_Chen})
!%    \item    Simple potential for an HF dimer ({\bf IPModel_HFdimer})
!%    \item    2-body potential for water dimer ({\bf IPModel_WaterDimer_Gillan})
!%    \item    Born-Mayer potential ({\bf IPModel_BornMayer})
!%    \item    Customised (hardwired) potential (\{\bf IPModel_Custom})
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
!%    \item    'IP Brenner_2002'
!%    \item    'IP PartridgeSchwenke'
!%    \item    'IP Einstein'
!%    \item    'IP Coulomb'
!%    \item    'IP Sutton_Chen'
!%    \item    'IP HFdimer'
!%    \item    'IP BornMayer'
!%    \item    'IP_WaterDimer_Gillan'
!%    \item    'IP_WaterTrimer_Gillan'
!%    \item    'IP Custom'
!%    \item    'IP Template'
!%   \end{itemize}
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#include "error.inc"
module IP_module

! use system_module, only : dp
use error_module
use system_module, only : dp, inoutput, initialise, INPUT, print, PRINT_VERBOSE, PRINT_ALWAYS, system_timer
use units_module
use mpi_context_module
use extendable_str_module
use dictionary_module
use paramreader_module
use atoms_module

use QUIP_Common_module
use IPModel_GAP_module, only : ipmodel_gap, initialise, finalise, calc, print
use IPModel_LJ_module, only : ipmodel_lj, initialise, finalise, calc, print
use IPModel_Morse_module, only : ipmodel_morse, initialise, finalise, calc, print
use IPModel_FC_module, only : ipmodel_fc, initialise, finalise, calc, print
use IPModel_SW_module, only : ipmodel_sw, initialise, finalise, calc, print
use IPModel_Tersoff_module, only : ipmodel_tersoff, initialise, finalise, calc, print
use IPModel_EAM_ErcolAd_module, only : ipmodel_eam_ercolad, initialise, finalise, calc, print
use IPModel_Brenner_module, only : ipmodel_brenner, initialise, finalise, calc, print
use IPModel_FB_module, only : ipmodel_fb, initialise, finalise, calc, print
use IPModel_Si_MEAM_module, only : ipmodel_si_meam, initialise, finalise, calc, print
use IPModel_FS_module, only : ipmodel_fs, initialise, finalise, calc, print
use IPModel_BOP_module, only : ipmodel_bop, initialise, finalise, calc, print
use IPModel_Brenner_Screened_module, only : ipmodel_brenner_screened, initialise, finalise, calc, print
use IPModel_Brenner_2002_module, only : ipmodel_brenner_2002, initialise, finalise, calc, print
#ifdef HAVE_ASAP
use IPModel_ASAP_module, only : ipmodel_asap, initialise, finalise, calc, print
#endif
use IPModel_TS_module, only : ipmodel_ts, initialise, finalise, calc, print, setup_atoms
use IPModel_Glue_module, only : ipmodel_glue, initialise, finalise, calc, print
use IPModel_PartridgeSchwenke_module, only : ipmodel_partridgeschwenke, initialise, finalise, calc, print
use IPModel_Einstein_module, only : ipmodel_einstein, initialise, finalise, calc, print
use IPModel_Coulomb_module, only : ipmodel_coulomb, initialise, finalise, calc, print
use IPModel_Sutton_Chen_module, only : ipmodel_sutton_chen, initialise, finalise, calc, print
use IPModel_Template_module, only : ipmodel_template, initialise, finalise, calc, print
#ifdef HAVE_KIM
use IPModel_KIM_module, only : ipmodel_kim, initialise, finalise, calc, print
#endif
use IPModel_FX_module, only : ipmodel_fx, initialise, finalise, calc, print
use IPModel_HFdimer_module, only : ipmodel_hfdimer, initialise, finalise, calc, print
use IPModel_WaterDimer_Gillan_module, only : ipmodel_waterdimer_gillan, initialise, finalise, calc, print
use IPModel_WaterTrimer_Gillan_module, only : ipmodel_watertrimer_gillan, initialise, finalise, calc, print
use IPModel_BornMayer_module, only : ipmodel_bornmayer, initialise, finalise, calc, print
use IPModel_Custom_module, only : ipmodel_custom, initialise, finalise, calc, print
use IPModel_SW_VP_module, only : ipmodel_sw_vp, initialise, finalise, calc, print
! Add new IP here

implicit none

private

integer, parameter :: FF_LJ = 1, FF_SW = 2, FF_Tersoff = 3, FF_EAM_ErcolAd = 4, &
     FF_Brenner = 5, FF_GAP = 6, FF_FS = 7, FF_BOP = 8, FF_FB = 9, FF_Si_MEAM = 10, FF_Brenner_Screened = 11, &
     FF_Brenner_2002 = 12, FF_ASAP = 13, FF_TS = 14, FF_FC = 15, FF_Morse = 16, FF_GLUE = 17, FF_PartridgeSchwenke = 18, &
     FF_Einstein = 19, FF_Coulomb = 20, FF_Sutton_Chen = 21, FF_KIM = 22, FF_FX = 23, FF_HFdimer = 24, FF_Custom = 25, FF_SW_VP=26, &
     FF_BornMayer = 27, FF_WaterDimer_Gillan=28, FF_WaterTrimer_Gillan=29, &! Add new IPs here
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
#ifdef HAVE_ASAP
  type(IPModel_ASAP) ip_ASAP
#endif
  type(IPModel_TS) ip_TS
  type(IPModel_Glue) ip_Glue
  type(IPModel_PartridgeSchwenke) ip_PartridgeSchwenke
  type(IPModel_Einstein) ip_Einstein
  type(IPModel_Coulomb) ip_Coulomb
  type(IPModel_Sutton_Chen) ip_Sutton_Chen
#ifdef HAVE_KIM
  type(IPModel_KIM) ip_KIM
#endif
  type(IPModel_FX) ip_FX
  type(IPModel_BornMayer) ip_BornMayer
  type(IPModel_Custom) ip_Custom
  type(IPModel_HFdimer) ip_HFdimer 
  type(IPModel_SW_VP) ip_sw_vp
  type(IPModel_WaterDimer_Gillan) ip_waterdimer_gillan
  type(IPModel_WaterTrimer_Gillan) ip_watertrimer_gillan
     ! Add new IP here
  type(IPModel_Template) ip_Template
  type(mpi_context) :: mpi_glob, mpi_local

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
       is_Brenner_Screened, is_Brenner_2002, is_ASAP, is_TS, is_Glue, is_PartridgeSchwenke, is_Einstein, is_Coulomb, &
       is_Sutton_Chen, is_KIM, is_FX, is_HFdimer, is_BornMayer, is_Custom, is_Template, is_SW_VP, is_WaterDimer_Gillan , is_WaterTrimer_Gillan
  ! Add new IPs here

  INIT_ERROR(error)

  call Finalise(this)

  call initialise(params)
  call param_register(params, 'GAP', 'false', is_GAP, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'LJ', 'false', is_LJ, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'Morse', 'false', is_Morse, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'FC', 'false', is_FC, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'SW', 'false', is_SW, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'Tersoff', 'false', is_Tersoff, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'EAM_ErcolAd', 'false', is_EAM_ErcolAd, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'Brenner', 'false', is_Brenner, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'FB', 'false', is_FB, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'Si_MEAM', 'false', is_Si_MEAM, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'FS', 'false', is_FS, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'BOP', 'false', is_BOP, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'Brenner_Screened', 'false', is_Brenner_Screened, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'Brenner_2002', 'false', is_Brenner_2002, help_string="No help yet.  This source file was $LastChangedBy$")
#ifdef HAVE_ASAP
  call param_register(params, 'ASAP', 'false', is_ASAP, help_string="No help yet.  This source file was $LastChangedBy$")
#else
  is_ASAP = .false.
#endif
  call param_register(params, 'TS', 'false', is_TS, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'Glue', 'false', is_Glue, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'PartridgeSchwenke', 'false', is_PartridgeSchwenke, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'Einstein', 'false', is_Einstein, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'Coulomb', 'false', is_Coulomb, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'Sutton_Chen', 'false', is_Sutton_Chen, help_string="No help yet.  This source file was $LastChangedBy$")
#ifdef HAVE_KIM
  call param_register(params, 'KIM', 'false', is_KIM, help_string="No help yet.  This source file was $LastChangedBy$")
#else
  is_KIM = .false.
#endif
  call param_register(params, 'FX', 'false', is_FX, help_string="Fanourgakis-Xantheas model of water")
  call param_register(params, 'BornMayer', 'false', is_BornMayer, help_string="Born-Mayer potential")
  call param_register(params, 'Custom', 'false', is_Custom, help_string="Customised (hard wired) potential")
  call param_register(params, 'HFdimer', 'false', is_HFdimer, help_string="HF dimer potential")
  call param_register(params, 'SW_VP', 'false', is_SW_VP, help_string="No help yet.  This source file was $LastChangedBy$")
  call param_register(params, 'WaterDimer_Gillan', 'false', is_WaterDimer_Gillan, help_string='Water Dimer 2-body potential by Mike Gillan')
  call param_register(params, 'WaterTrimer_Gillan', 'false', is_WaterTrimer_Gillan, help_string='Water Trimer 3-body potential by Mike Gillan')
  
 ! Add new IP here
  call param_register(params, 'Template', 'false', is_Template, help_string="No help yet.  This source file was $LastChangedBy$")

  if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='IP_Initialise_str args_str')) then
    RAISE_ERROR("IP_Initialise_str failed to parse args_str='"//trim(args_str)//"'", error)
  endif
  call finalise(params)

  if (count((/is_GAP, is_LJ, is_FC, is_Morse, is_SW, is_Tersoff, is_EAM_ErcolAd, is_Brenner, is_FS, is_BOP, is_FB, is_Si_MEAM, &
       is_Brenner_Screened, is_Brenner_2002, is_ASAP, is_TS, is_Glue, is_PartridgeSchwenke, is_Einstein, is_Coulomb, &
       is_Sutton_Chen, is_KIM, is_FX, is_HFdimer, is_BornMayer, is_Custom, is_SW_VP, is_WaterDimer_Gillan,is_WaterTrimer_Gillan,&       ! add new IPs here
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
#ifdef HAVE_ASAP
  else if (is_ASAP) then
    this%functional_form = FF_ASAP
    call Initialise(this%ip_ASAP, args_str, param_str) 
#endif
  else if (is_TS) then
    this%functional_form = FF_TS
    call Initialise(this%ip_TS, args_str, param_str) 
  else if (is_Glue) then
    this%functional_form = FF_GLUE
    call Initialise(this%ip_Glue, args_str, param_str) 
 else if (is_PartridgeSchwenke) then
    this%functional_form = FF_PartridgeSchwenke
    call Initialise(this%ip_PartridgeSchwenke, args_str, param_str)
  else if (is_Einstein) then
    this%functional_form = FF_Einstein
    call Initialise(this%ip_Einstein, args_str, param_str) 
  else if (is_Coulomb) then
    this%functional_form = FF_Coulomb
    call Initialise(this%ip_Coulomb, args_str, param_str) 
  else if (is_Sutton_Chen) then
    this%functional_form = FF_Sutton_Chen
    call Initialise(this%ip_Sutton_Chen, args_str, param_str) 
#ifdef HAVE_KIM
  else if (is_KIM) then
    this%functional_form = FF_KIM
    call Initialise(this%ip_KIM, args_str, param_str) 
#endif
  else if (is_FX) then
     this%functional_form = FF_FX
     call Initialise(this%ip_fx, args_str, param_str)
  else if (is_BornMayer) then
     this%functional_form = FF_BornMayer
     call Initialise(this%ip_bornmayer, args_str, param_str)
  else if (is_Custom) then
     this%functional_form = FF_Custom
     call Initialise(this%ip_custom, args_str, param_str)
  else if (is_HFdimer) then
     this%functional_form = FF_HFdimer
     call Initialise(this%ip_HFdimer, args_str, param_str)
  else if (is_SW_VP) then
    this%functional_form = FF_SW_VP
    call Initialise(this%ip_sw_vp, args_str, param_str)
 else if (is_WaterDimer_Gillan) then
    this%functional_form = FF_WaterDimer_Gillan
    call Initialise(this%ip_waterdimer_gillan, args_str, param_str)
 else if (is_WaterTrimer_Gillan) then
    this%functional_form = FF_WaterTrimer_Gillan
    call Initialise(this%ip_watertrimer_gillan, args_str, param_str)
    ! Add new IP here
  else if (is_Template) then
    this%functional_form = FF_Template
    call Initialise(this%ip_Template, args_str, param_str) 
  end if

  if (present(mpi_obj)) this%mpi_glob = mpi_obj

end subroutine IP_Initialise_str

subroutine IP_Finalise(this)
  type(IP_type), intent(inout) :: this

  if (this%mpi_local%active) call free_context(this%mpi_local)
  select case (this%functional_form)
    case (FF_GAP)
      call Finalise(this%ip_gap)
    case (FF_LJ)
      call Finalise(this%ip_lj)
    case (FF_Morse)
      call Finalise(this%ip_morse)
    case (FF_FC)
      call Finalise(this%ip_fc)
    case (FF_SW)
      call Finalise(this%ip_sw)
    case (FF_Tersoff)
      call Finalise(this%ip_tersoff)
    case (FF_EAM_ErcolAd)
      call Finalise(this%ip_EAM_ErcolAd)
    case(FF_Brenner)
      call Finalise(this%ip_Brenner)
    case(FF_FB)
      call Finalise(this%ip_fb)
    case(FF_Si_MEAM)
      call Finalise(this%ip_Si_MEAM)
    case (FF_FS)
      call Finalise(this%ip_fs)
    case (FF_BOP)
      call Finalise(this%ip_bop)
    case (FF_Brenner_Screened)
      call Finalise(this%ip_brenner_screened)
    case (FF_Brenner_2002)
      call Finalise(this%ip_brenner_2002)
#ifdef HAVE_ASAP
    case (FF_ASAP)
      call Finalise(this%ip_ASAP)
#endif
    case (FF_TS)
      call Finalise(this%ip_TS)
    case (FF_GLUE)
      call Finalise(this%ip_Glue)
   case (FF_PartridgeSchwenke)
      call Finalise(this%ip_PartridgeSchwenke)
    case (FF_Einstein)
      call Finalise(this%ip_Einstein)
    case (FF_Coulomb)
      call Finalise(this%ip_Coulomb)
    case (FF_Sutton_Chen)
      call Finalise(this%ip_Sutton_Chen)
#ifdef HAVE_KIM
    case (FF_KIM)
      call Finalise(this%ip_KIM)
#endif
   case (FF_FX)
      call Finalise(this%ip_FX)
   case(FF_BornMayer)
      call Finalise(this%ip_BornMayer)
   case (FF_Custom)
      call Finalise(this%ip_Custom)
   case (FF_HFdimer)
      call Finalise(this%ip_HFdimer)
   case (FF_SW_VP)
      call Finalise(this%ip_sw_vp)
   case (FF_WaterDimer_Gillan)
      call Finalise(this%ip_waterdimer_gillan)
   case (FF_WaterTrimer_Gillan)
      call Finalise(this%ip_watertrimer_gillan)
      ! add new IP here
   case (FF_Template)
      call Finalise(this%ip_Template)
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
#ifdef HAVE_ASAP
  case (FF_ASAP)
     IP_cutoff = max(this%ip_asap%cutoff_ms, this%ip_asap%cutoff_coulomb)*BOHR
#endif
  case (FF_TS)
     IP_cutoff = max(this%ip_TS%cutoff_ms, this%ip_TS%cutoff_coulomb)*BOHR
  case (FF_GLUE)
     IP_cutoff = this%ip_Glue%cutoff
  case (FF_PartridgeSchwenke)
     IP_cutoff = this%ip_PartridgeSchwenke%cutoff
  case (FF_Einstein)
     IP_cutoff = this%ip_Einstein%cutoff
  case (FF_Coulomb)
     IP_cutoff = this%ip_Coulomb%cutoff
  case (FF_Sutton_Chen)
     IP_cutoff = this%ip_Sutton_Chen%cutoff
#ifdef HAVE_KIM
  case (FF_KIM)
     IP_cutoff = this%ip_kim%cutoff
#endif
  case (FF_FX)
     IP_cutoff = this%ip_fx%cutoff
  case (FF_BornMayer)
     IP_cutoff = this%ip_BornMayer%cutoff
  case (FF_Custom)
     IP_cutoff = this%ip_custom%cutoff
  case (FF_HFdimer)
     IP_cutoff = this%ip_HFdimer%cutoff
  case (FF_SW_VP)
     IP_cutoff = this%ip_sw_vp%cutoff 
  case (FF_WaterDimer_Gillan)
     IP_cutoff = this%ip_waterdimer_gillan%cutoff
  case (FF_WaterTrimer_Gillan)
     IP_cutoff = this%ip_watertrimer_gillan%cutoff
  ! Add new IP here
  case (FF_Template)
     IP_cutoff = this%ip_Template%cutoff
  case default
     IP_cutoff = 0.0_dp
  end select
end function IP_cutoff

subroutine IP_setup_atoms(this, at)
  type(IP_type), intent(in) :: this
  type(Atoms), intent(inout) :: at

  select case (this%functional_form)
  case (FF_TS)
     call setup_atoms(this%ip_TS, at)
  end select

end subroutine IP_setup_atoms

subroutine IP_Calc(this, at, energy, local_e, f, virial, local_virial, args_str, error)
  type(IP_type), intent(inout) :: this
  type(Atoms), intent(inout) :: at                
  real(dp), intent(out), optional :: energy, local_e(:) !% \texttt{energy} = System total energy, \texttt{local_e} = energy of each atom, vector dimensioned as \texttt{at%N}.
  real(dp), intent(out), optional :: f(:,:), local_virial(:,:)
  real(dp), intent(out), optional :: virial(3,3)
  character(len=*), intent(in), optional      :: args_str 
  integer, intent(out), optional :: error

  INIT_ERROR(error)

  call system_timer("IP_Calc")

  if (this%mpi_glob%active .and. .not. this%mpi_local%active) then
    call setup_parallel(this, at, energy, local_e, f, virial, local_virial)
  endif

  select case (this%functional_form)
    case (FF_GAP)
      call calc(this%ip_gap, at, energy, local_e, f, virial, local_virial, args_str, mpi=this%mpi_local, error=error)
    case (FF_LJ)
      call calc(this%ip_lj, at, energy, local_e, f, virial, local_virial, args_str, mpi=this%mpi_local, error=error)
      PASS_ERROR(error)
    case (FF_Morse)
      call calc(this%ip_morse, at, energy, local_e, f, virial, local_virial, args_str, mpi=this%mpi_local, error=error)
    case (FF_FC)
      call calc(this%ip_fc, at, energy, local_e, f, virial, local_virial, args_str, mpi=this%mpi_local, error=error)
    case (FF_SW)
      call calc(this%ip_sw, at, energy, local_e, f, virial, local_virial, args_str, mpi=this%mpi_local, error=error)
    case (FF_Tersoff)
      call calc(this%ip_tersoff, at, energy, local_e, f, virial, local_virial, args_str, mpi=this%mpi_local, error=error)
    case (FF_EAM_ErcolAd)
      call calc(this%ip_EAM_ErcolAd, at, energy, local_e, f, virial, local_virial, args_str, mpi=this%mpi_local, error=error)
    case(FF_Brenner)
      call calc(this%ip_Brenner, at, energy, local_e, f, virial, local_virial, args_str, mpi=this%mpi_local, error=error)
    case(FF_FB)
      call calc(this%ip_FB, at, energy, local_e, f, virial, local_virial, args_str, mpi=this%mpi_local, error=error)
    case(FF_Si_MEAM)
      call calc(this%ip_Si_MEAM, at, energy, local_e, f, virial, local_virial, args_str, mpi=this%mpi_local, error=error)
    case (FF_FS)
      call calc(this%ip_fs, at, energy, local_e, f, virial, local_virial, args_str, mpi=this%mpi_local, error=error)
    case (FF_BOP)
      call calc(this%ip_bop, at, energy, local_e, f, virial, local_virial, args_str, mpi=this%mpi_local, error=error)
    case (FF_Brenner_Screened)
      call calc(this%ip_brenner_screened, at, energy, local_e, f, virial, local_virial, args_str, mpi=this%mpi_local, error=error)
    case (FF_Brenner_2002)
      call calc(this%ip_brenner_2002, at, energy, local_e, f, virial, local_virial, args_str, mpi=this%mpi_local, error=error)
#ifdef HAVE_ASAP
    case (FF_ASAP)
      call calc(this%ip_asap, at, energy, local_e, f, virial, local_virial, args_str, mpi=this%mpi_local, error=error)
#endif
    case (FF_TS)
      call calc(this%ip_TS, at, energy, local_e, f, virial, local_virial, args_str, mpi=this%mpi_local, error=error)
    case (FF_GLUE)
      call calc(this%ip_Glue, at, energy, local_e, f, virial, local_virial, args_str, mpi=this%mpi_local, error=error)
   case (FF_PartridgeSchwenke)
      call calc(this%ip_PartridgeSchwenke, at, energy, local_e, f, virial, local_virial, args_str, mpi=this%mpi_local, error=error)
    case (FF_Einstein)
      call calc(this%ip_Einstein, at, energy, local_e, f, virial, local_virial, args_str, mpi=this%mpi_local, error=error)
    case (FF_Coulomb)
      call calc(this%ip_Coulomb, at, energy, local_e, f, virial, local_virial, args_str, mpi=this%mpi_local, error=error)
    case (FF_Sutton_Chen)
      call calc(this%ip_Sutton_Chen, at, energy, local_e, f, virial, local_virial, args_str, mpi=this%mpi_local, error=error)
#ifdef HAVE_KIM
    case (FF_KIM)
      call calc(this%ip_KIM, at, energy, local_e, f, virial, local_virial, args_str, mpi=this%mpi_local, error=error)
      PASS_ERROR(error)
#endif
   case (FF_FX)
      call calc(this%ip_FX, at, energy, local_e, f, virial, local_virial, args_str, mpi=this%mpi_local, error=error)
   case (FF_BornMayer)
      call calc(this%ip_BornMayer, at, energy, local_e, f, virial, local_virial, args_str, mpi=this%mpi_local, error=error)
   case (FF_Custom)
      call calc(this%ip_Custom, at, energy, local_e, f, virial, local_virial, args_str, mpi=this%mpi_local, error=error)
   case (FF_HFdimer)
      call calc(this%ip_HFdimer, at, energy, local_e, f, virial, local_virial, args_str, mpi=this%mpi_local, error=error)
   case (FF_SW_VP)
      call calc(this%ip_sw_vp, at, energy, local_e, f, virial, local_virial, args_str, mpi=this%mpi_local, error=error)
   case (FF_WaterDimer_Gillan)
      call calc(this%ip_waterdimer_gillan, at, energy, local_e, f, virial, local_virial, args_str, mpi=this%mpi_local, error=error)
   case (FF_WaterTrimer_Gillan)
      call calc(this%ip_watertrimer_gillan, at, energy, local_e, f, virial, local_virial, args_str, mpi=this%mpi_local, error=error)
   ! Add new IP here   
   case (FF_Template)
      call calc(this%ip_Template, at, energy, local_e, f, virial, local_virial, args_str, mpi=this%mpi_local, error=error)
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
#ifdef HAVE_ASAP
    case (FF_ASAP)
      call Print(this%ip_asap, file=file)
#endif
    case (FF_TS)
      call Print(this%ip_TS, file=file)
    case (FF_GLUE)
      call Print(this%ip_Glue, file=file)
   case (FF_PartridgeSchwenke)
      call Print(this%ip_PartridgeSchwenke, file=file)
    case (FF_Einstein)
      call Print(this%ip_Einstein, file=file)
    case (FF_Coulomb)
      call Print(this%ip_Coulomb, file=file)
    case (FF_Sutton_Chen)
      call Print(this%ip_Sutton_Chen, file=file)
#ifdef HAVE_KIM
    case (FF_KIM)
      call Print(this%ip_kim, file=file)
#endif
   case (FF_FX)
      call Print(this%ip_fx, file=file)
   case (FF_BornMayer)
      call Print(this%ip_BornMayer, file=file)
   case (FF_Custom)
      call Print(this%ip_custom, file=file)
   case (FF_HFdimer)
      call Print(this%ip_HFdimer, file=file)
   case (FF_SW_VP)
      call Print(this%ip_sw_vp, file=file)
   case (FF_WaterDimer_Gillan)
      call Print(this%ip_waterdimer_gillan, file=file)
   case (FF_WaterTrimer_Gillan)
      call Print(this%ip_watertrimer_gillan, file=file)
    ! add new IP here
   case (FF_Template)
      call Print(this%ip_Template, file=file)
   case default
      RAISE_ERROR("IP_Print confused by functional_form " // this%functional_form, error)
  end select

end subroutine IP_Print

subroutine IP_setup_parallel(this, at, energy, local_e, f, virial, local_virial, args_str)
  type(IP_type), intent(inout) :: this
  type(Atoms), intent(inout) :: at
  real(dp), intent(out), optional :: energy, local_e(:) !% \texttt{energy} = System total energy, \texttt{local_e} = energy of each atom, vector dimensioned as \texttt{at%N}.  
  real(dp), intent(out), optional :: f(:,:), local_virial(:,:)   !% Forces, dimensioned as \texttt{f(3,at%N)}, local virials, dimensioned as \texttt{local_virial(9,at%N)} 
  real(dp), intent(out), optional :: virial(3,3)
  character(len=*), intent(in), optional      :: args_str 


  integer :: pgroup_size, prev_pgroup_size
  integer :: n_groups
  real(dp) :: prev_time, this_time

  prev_time = 1.0e38_dp
  prev_pgroup_size = 0

  call print('IP_Setup_Parallel timings', PRINT_VERBOSE)
  call print('IP_Setup_Parallel group_size  time/sec', PRINT_VERBOSE)

  call print("IP_Setup_Parallel at%N "//at%N// " present(energy) "//present(energy)// " present(local_e) "//present(local_e)// " present(f) "//present(f)// " present(virial) "//present(virial)// " present(local_virial) "//present(local_virial), PRINT_VERBOSE)

  call setup_atoms(this, at)
  do pgroup_size=this%mpi_glob%n_procs, 1, -1
    n_groups = this%mpi_glob%n_procs / pgroup_size
    if (n_groups*pgroup_size == this%mpi_glob%n_procs) then
      call setup_parallel_groups(this, this%mpi_glob, pgroup_size)
      call system_timer("IP_parallel", do_always = .true.)
      call calc(this, at, energy, local_e, f, virial, local_virial, args_str)
      call system_timer("IP_parallel", do_always = .true., time_elapsed = this_time)
      this_time = max(this%mpi_glob, this_time)
      call print("IP_Setup_Parallel "//pgroup_size//' '//this_time, PRINT_VERBOSE)
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

  integer :: split_index

  INIT_ERROR(error)

  if (this%mpi_local%active) call free_context(this%mpi_local)

  if (mpi%active) then
    split_index = mpi%my_proc/pgroup_size
    call split_context(mpi, split_index, this%mpi_local)
  endif

end subroutine setup_parallel_groups

end module IP_module
