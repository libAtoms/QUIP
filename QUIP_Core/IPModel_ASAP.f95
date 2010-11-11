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
!X IPModel_ASAP
!X
!% Interface to ASAP potential.
!% P. Tangney and S. Scandolo,
!% An ab initio parametrized interatomic force field for silica
!% J. Chem. Phys, 117, 8898 (2002). 
!%
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#include "error.inc"

module IPModel_ASAP_module

use error_module
use libatoms_module, only: dp, FIELD_LENGTH
use mpi_context_module
use QUIP_Common_module

implicit none
private

include 'IPModel_interface.h'

logical, private :: asap_initialised = .false.

public :: IPModel_ASAP
type IPModel_ASAP
  integer :: n_types = 0, n_atoms = 0
  real(dp) :: betapol, tolpol, yukalpha, yuksmoothlength
  integer :: maxipol, pred_order
  integer, allocatable :: atomic_num(:), type_of_atomic_num(:)
  real(dp), allocatable, dimension(:) :: pol, z
  real(dp), allocatable, dimension(:,:) :: D_ms, gamma_ms, R_ms, B_pol, C_pol
  logical :: tewald, tdip_sr
  real :: raggio, a_ew, gcut
  integer :: iesr(3)

  real(dp) :: cutoff(4)

  character(len=FIELD_LENGTH) :: label
  logical :: initialised

end type IPModel_ASAP

logical, private :: parse_in_ip, parse_matched_label
type(IPModel_ASAP), private, pointer :: parse_ip

interface Initialise
  module procedure IPModel_ASAP_Initialise_str
end interface Initialise

interface Finalise
  module procedure IPModel_ASAP_Finalise
end interface Finalise

interface Print
  module procedure IPModel_ASAP_Print
end interface Print

interface Calc
  module procedure IPModel_ASAP_Calc
end interface Calc

contains

#ifdef HAVE_ASAP
  subroutine asap_singlepoint_finalise()

    use atoms
    use verlet
    use neighbour 
    use pot_parameters
    use stuff 
    use energy_i
    use print65
    use electric_field
    use nndim
    use neighbour3
    use fixpar
    use distortion, only : taimsp,rmin_aim
    use thermal_cond
    use shakes
    use compress
    use quench
    use forcetest
    use iteration
    use changepress
    use metric
    use netquant
    use universal
    use initstuff
    use plot_efield
    use polar
    use structopt
    use parameters
    use presstemp
    use testf
    use minimiser
    implicit none

    deallocate(mass,spind,objind,numobj,numobjsp)
    deallocate(nnlist,nnlist2,nnlist3)
    deallocate(nnat,nnat2,nnat3,nn_imag)
    deallocate(nn_imag2,nn_imag3)
    deallocate(s,sm)
    deallocate(r,rm)
    deallocate(Z,bij,Cij)
    deallocate(Dij,alphaij)
    deallocate(Eij,Nij)
    deallocate(pol,bpol,cpol)
    deallocate(theta0jik,bjik)
    deallocate(b_tt,r_ms)
    deallocate(gamma_ms,d_ms)
    deallocate(c_harm,rmin,rmin_aim)
    deallocate(adist,bdist,cdist,ddist)
    deallocate(bu1,alphau1)
    deallocate(sigma1,sigma2,sigma3,sigma4)
    deallocate(bu2,alphau2)
    deallocate(taimsp,minr,mindistl)
    deallocate(elements)
    deallocate(aelements)
    deallocate(ielements)
    deallocate(dt2bym, aumass)
    deallocate(parf, tparf)
    deallocate(posobj,velobj,dipobj)
    deallocate(timposepos,timposevel,timposedip)
    deallocate(zerovec,falseimp)
    if (allocated(dip_perm)) deallocate(dip_perm)
    if (t_therm) deallocate(nzave,ztempave)

    if (allocated(efield_old)) deallocate(efield_old)
    if (allocated(efield)) deallocate(efield)
    if (allocated(fqdip)) deallocate(fqdip)
    if (allocated(fdipdip)) deallocate(fdipdip)
    if (allocated(force_ew)) deallocate(force_ew)
    if (allocated(efield_ew)) deallocate(efield_ew)
    if (allocated(dip)) deallocate(dip)
    if (allocated(stress_ew)) deallocate(stress_ew)
    if (allocated(dip_sr)) deallocate(dip_sr)
    if (allocated(efield_old)) deallocate(efield_old)
    if (allocated(force_pol)) deallocate(force_pol)
    if (allocated(stress_pol)) deallocate(stress_pol)
    if (allocated(testt)) deallocate(testt)

  end subroutine asap_singlepoint_finalise

#endif


subroutine IPModel_ASAP_Initialise_str(this, args_str, param_str)
  use libatoms_module
  type(IPModel_ASAP), intent(inout) :: this
  character(len=*), intent(in) :: args_str, param_str

  type(Dictionary) :: params

  this%initialised = .false.

  call Finalise(this)

  call initialise(params)
  this%label=''
  call param_register(params, 'label', '', this%label)
  if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='IPModel_ASAP_Initialise_str args_str')) then
    call system_abort("IPModel_ASAP_Initialise_str failed to parse label from args_str="//trim(args_str))
  endif
  call finalise(params)

  call IPModel_ASAP_read_params_xml(this, param_str)
  this%n_atoms = 0
  this%initialised = .true.
  
end subroutine IPModel_ASAP_Initialise_str

subroutine IPModel_ASAP_Finalise(this)
  use libatoms_module
  type(IPModel_ASAP), intent(inout) :: this

#ifdef HAVE_ASAP
  if (asap_initialised) then
     call asap_singlepoint_finalise()
     asap_initialised = .false.
  end if
  this%initialised = .false.
#else
  call system_abort('ASAP potential is not compiled in. Recompile with HAVE_ASAP=1')
#endif

  if (allocated(this%atomic_num)) deallocate(this%atomic_num)
  if (allocated(this%type_of_atomic_num)) deallocate(this%type_of_atomic_num)
  if (allocated(this%pol)) deallocate(this%pol)
  if (allocated(this%z)) deallocate(this%z)
  if (allocated(this%D_ms)) deallocate(this%D_ms)
  if (allocated(this%gamma_ms)) deallocate(this%gamma_ms)
  if (allocated(this%R_ms)) deallocate(this%R_ms)
  if (allocated(this%B_pol)) deallocate(this%B_pol)
  if (allocated(this%C_pol)) deallocate(this%C_pol)
  this%n_types = 0
  this%label = ''
end subroutine IPModel_ASAP_Finalise


subroutine IPModel_ASAP_Calc(this, at, e, local_e, f, virial, local_virial, args_str, mpi, error)
#ifdef HAVE_ASAP
  use neighbour!, only : nsp,nat,iesr,spind,rcut
  use neighbour3
  use pot_parameters
  use logical_stuff
  use polar
  use stuff
  use distortion, only :  xgmin,xgmax,taimsp
  use parameters
  use testf
  use initstuff
  use atoms
  use verlet
  use iteration
  use nndim
  use metric
  use fixpar, only: ntype, tpbc, treadpar
  use energy_i
  use print65
  use electric_field
  use nndim
  use distortion, only : taimsp,rmin_aim
  use thermal_cond
  use shakes
  use compress
  use quench
  use forcetest
  use iteration
  use changepress
  use metric
  use netquant
  use universal, only: au_to_kbar, small
  use initstuff
  use plot_efield
  use polar
  use structopt
  use parameters
  use presstemp
  use testf
  use minimiser

#endif
   use libAtoms_module, only : Atoms, Dictionary, BOHR, HARTREE, &
        ElementMass, massconvert, initialise, param_read_line, finalise, &
        ElementName, print, operator(//), cell_volume, has_property, &
        assign_pointer, add_property, system_abort, fit_box_in_cell, &
        param_register
	
   type(IPModel_ASAP), intent(inout):: this
   type(Atoms), intent(inout)      :: at
   real(dp), intent(out), optional :: e, local_e(:)
   real(dp), intent(out), optional :: f(:,:), local_virial(:,:)   !% Forces, dimensioned as \texttt{f(3,at%N)}, local virials, dimensioned as \texttt{local_virial(9,at%N)} 
   real(dp), intent(out), optional :: virial(3,3)
   character(len=*), optional, intent(in) :: args_str
   type(MPI_Context), intent(in), optional :: mpi

   integer, intent(out), optional :: error
   type(Dictionary) :: params	 
   real(dp) :: asap_e, asap_stress(3,3), a_ew, f_ew	 
   logical smlra,smlrc,smlgc,smlrc2,smlrc3,smlrc4	 
   logical nsmlra,nsmlrc,nsmlgc	 
   real*8 rarec,gcrec,rcrec,root2	 
   real*8 raggio_in,rcut_in,gcut_in	 
   real(dp), allocatable :: asap_f(:,:)	 
   real(dp), pointer :: dipoles_ptr(:,:)	 
   integer :: i, ti, tj	 
   logical :: do_restart, calc_dipoles	 
   integer at0, atf	 
   integer idebug	 
   real(dp) dtold,dtnew	 
   logical  tzeroc	 
   logical sumewald	 
   real(dp) mass_cel	 
   logical tscaled	 
   logical readnat	 
   logical texist	 
   logical tpow,tgmin	 
   integer nesr

   logical :: has_atom_mask_name
   character(FIELD_LENGTH) :: atom_mask_name

#ifdef HAVE_ASAP

   INIT_ERROR(error)

   if (present(e)) e = 0.0_dp
   if (present(local_e)) then
      call check_size('Local_E',local_e,(/at%N/),'IPModel_ASAP_Calc', error)
      local_e = 0.0_dp
   endif
   if (present(f)) then
      call check_size('Force',f,(/3,at%Nbuffer/),'IPModel_ASAP_Calc', error)
      f = 0.0_dp
   end if
   if (present(virial)) virial = 0.0_dp
   if (present(local_virial)) then
      call check_size('Local_virial',local_virial,(/9,at%Nbuffer/),'IPModel_ASAP_Calc', error)
      local_virial = 0.0_dp
      RAISE_ERROR("IPModel_ASAP_Calc: local_virial calculation requested but not supported yet.", error)
   endif

   allocate(asap_f(3,at%N))

   call initialise(params)
   this%label=''
   call param_register(params, 'restart', 'F', do_restart)
   call param_register(params, 'calc_dipoles', 'F', calc_dipoles)
   call param_register(params, 'atom_mask_name', 'NONE', atom_mask_name, has_value_target=has_atom_mask_name)
   if (.not. param_read_line(params, args_str, ignore_unknown=.true.,task='IPModel_ASAP_Calc args_str')) then
      call system_abort("IPModel_ASAP_Initialise_str failed to parse args_str="//trim(args_str))
   endif
   call finalise(params)

   if(has_atom_mask_name) then
      RAISE_ERROR('IPModel_ASAP_Calc: atom_mask_name found, but not supported', error)
   endif

   if (.not. asap_initialised .or. this%n_atoms /= at%n) then
      
      if (asap_initialised) then
         call asap_singlepoint_finalise()
         asap_initialised = .false.
      end if

      this%n_atoms = at%n

      idebug = 0 
      nat = this%n_atoms
      nsp = this%n_types
      at0 = 1
      atf = nat
      nobj = 1

      nthermo = min(atf-at0+1,nat-1)
      allocate(mass(this%n_types))

      tsinglepoint = .true.
      tangstrom = .false.
      irestart = -1
      itsave = 0
      readunit = 90
      saveunit = 90
      readnat = .false.
      dtold = 0.0d0
      dtnew = 0.0d0
      ntstep = 0
      iprint = 0
      iprintxyz = 0
      deltat = dtnew
      tcenter = .false.
      timposepbc = .false.
      tremovetrans = .false.
      treadxyz = .false.
      tfirst_xyzprint = .true.
      treadpar = .false.
      tappendfiles = .false.
      tstructopt = .false.
      tsdp = .false.
      tpbc = .true.
      tpow = .false.
      tgmin = .false.
      testewald = .false.
      ttime = .true.
      tforcetest = .false. 
      tscaled = .false.
      tangstrom = .false.
      texist = .false.
      tprint65 = .false.
      tenergy_i = .false.
      tbin65 = .false.
      iprint65 = 1
      tzeroc = .false.
      tcel = .false.
      tsdc = .false.
      cell_dir = 3
      press = 0.0d0
      mass_cel = 2.0d6
      press = press / au_to_kbar
      wc = mass_cel
      trescale = .false.
      tnosep = .false.
      tboltz = .false.
      tshake = .false.
      tcompress = .false.
      tquench = .false.
      tchangeP= .false.
      tpolwrite = .false.
      t_therm = .false.
      tshowforce = .false.
      tplot_efield = .false.
      tfirst_efield_plot = .true.
      tewald = .false.
      write_sr = .false.
      read_sr= .false.
      calc_sr=.false.
      calc_pol = .false.
      write_pol=.false.
      read_pol=.false.
      calc_gvec=.false.
      write_gvec=.false.
      read_gvec=.false.
      write_aim = .false.
      read_aim= .false.
      calc_aim=.false.
      tzvar=.false.
      tpolvar=.false.
      tsrvar=.false.
      tpol=.false.
      tsr=.false.
      taim = .false.
      tz=.false.
      hafta=.false.
      tbegin = .true.
      tyukawa = .true.
      
      sumewald = .true.
      allocate(spind(nat),objind(nat),numobj(nobj),numobjsp(nobj,nsp))
      allocate(z(nsp),bij(nsp,nsp),Cij(nsp,nsp))
      allocate(Dij(nsp,nsp),alphaij(nsp,nsp))
      allocate(Eij(nsp,nsp),Nij(nsp,nsp))
      allocate(pol(nsp),bpol(nsp,nsp),cpol(nsp,nsp))
      allocate(theta0jik(nsp,nsp,nsp),bjik(nsp,nsp,nsp))
      allocate(b_tt(3,nsp,nsp),r_ms(nsp,nsp))
      allocate(gamma_ms(nsp,nsp),d_ms(nsp,nsp))
      allocate(c_harm(nsp,nsp),rmin(nsp,nsp))
      allocate(rmin_aim(nsp,nsp))
      rmin_aim = 0.d0
      allocate(adist(nsp,nsp),bdist(nsp,nsp))
      allocate(cdist(nsp),ddist(nsp,nsp))
      allocate(bu1(nsp,nsp),alphau1(nsp,nsp))
      allocate(sigma1(nsp),sigma2(nsp),sigma3(nsp),sigma4(nsp))
      allocate(bu2(nsp,nsp),alphau2(nsp,nsp))
      allocate(taimsp(nsp))
      allocate(minr(nsp,nsp),mindistl(nsp,nsp))
      mindistl = 0.d0

      allocate(s(3,nat),sm(3,nat))
      allocate(r(3,nat),rm(3,nat))
      allocate(dt2bym(nat),aumass(nat))
      allocate(elements(nsp))
      allocate(aelements(nat))
      allocate(ielements(nat))
      nparm = (nsp*nsp*nsp+21*nsp*nsp+28*nsp)/2 + 4*nsp
      allocate (parf(nparm),tparf(nparm))
      tgen = .false.
      testforce = .false.
      npar = 1

      ntnlist = (/1,1/)
      nnatmax = 100000
      nnatmax_sr = 30000
      allocate(posobj(3,nobj),velobj(3,nobj),dipobj(3,nobj))
      allocate(timposepos(nobj),timposevel(nobj),timposedip(nobj))
      allocate(zerovec(3,nobj),falseimp(nobj))
      zerovec = 0.0d0!

      falseimp = .false.
      XNOS2M = 0.D0
      XNOSM  = 0.D0
      XNOS0  = 0.D0
      VRNOS  = 0.D0 
      
      it = 1

      c_harm = 0.0_dp
      rmin = 0.0_dp

      ! Per type parameters
      do ti=1,this%n_types
         mass(ti)     = ElementMass(this%atomic_num(ti))/MASSCONVERT
         elements(ti) = ElementName(this%atomic_num(ti))
         z(ti) = this%z(ti)
         pol(ti) = this%pol(ti)
         ntype(ti) = count(at%Z == this%atomic_num(ti))
      end do

      write(6,'(/," Nos of each species:",10i6)') (ntype(i),i=1,nsp)
      write(6,'(/," Masses :",2(1x,f10.5))') mass

      ! Ewald parameters
      tewald = this%tewald
      raggio = this%raggio
      a_ew = this%a_ew
      gcut = this%gcut
      rcut = this%cutoff
      iesr = this%iesr

      call print('tewald :'//tewald)
      write(6,'(/," Raggio :",1x,f10.5)') raggio
      write(6,'(" A_ew   :",1x,e10.5)') a_ew
      write(6,'(" iesr   :",3i3)') iesr
      write(6,'(" Gcut   :",1x,f10.5,/)') gcut

      smlgc = gcut.lt.small
      smlrc = rcut(1).lt.small
      smlra = raggio.lt.small
      smlrc2 = rcut(2).lt.small
      smlrc3 = rcut(3).lt.small
      smlrc4 = rcut(4).lt.small

      nsmlgc = .not.smlgc      ! iesr will be computed later from rcut and cell

      nsmlrc = .not.smlrc
      nsmlra = .not.smlra

      f_ew = dsqrt(-1.0d0*log(a_ew))
      root2 = dsqrt(2.0d0)

      rarec = 0.0d0
      rcrec = 0.0d0
      gcrec = 0.0d0

      raggio_in = raggio
      gcut_in = gcut
      rcut_in = rcut(1)

      if (smlrc.and.smlgc.and.smlra) then
       rarec = 3.0d0
       gcrec = f_ew*root2/rarec
       rcrec = gcut*rarec*rarec
      else if (nsmlrc.and.smlgc.and.smlra) then
       rarec = rcut(1)/f_ew/root2
       gcrec = f_ew*root2/rarec
       raggio = rarec
       gcut = gcrec
      else if (smlrc.and.nsmlgc.and.smlra) then
       rarec  = f_ew*root2/gcut
       rcrec  = gcut*rarec*rarec
       raggio = rarec
       rcut(1) = rcrec
      else if (smlrc.and.smlgc.and.nsmlra) then
       gcrec = f_ew*dsqrt(2.0d0)/raggio
       rcrec = gcrec*raggio*raggio
       gcut = gcrec
       rcut(1) = rcrec
      else if (nsmlrc.and.nsmlgc.and.smlra) then
       rarec = rcut(1)/f_ew/root2
       gcrec = f_ew*root2/raggio
       raggio = rarec
      else if (nsmlrc.and.smlgc.and.nsmlra) then
       rarec = rcut(1)/f_ew/root2
       gcrec = f_ew*root2/raggio
       gcut = gcrec
      else if (smlrc.and.nsmlgc.and.nsmlra) then
       gcrec = f_ew*root2/raggio
       rcrec = gcrec*raggio*raggio
       rcut(1) = rcrec
      else if (nsmlrc.and.nsmlgc.and.nsmlra) then
       write(6,*) 'raggio, rcut, gcut, are all defined!!'
      endif

      if (smlrc2) rcut(2) = rcut(1)
      if (smlrc3) rcut(3) = rcut(1)
      if (smlrc4) rcut(4) = rcut(1)

      write(6,'(/," Input Raggio                : ",f10.5)')raggio_in
      write(6,'(" Input r-space cutoff        : ",f10.5)')rcut_in
      write(6,'(" Input g-space cutoff        : ",f10.5)')gcut_in
      write(6,'(/," Recommended Raggio          : ",f10.5)')rarec
      write(6,'(" Recommended r-space cutoff  : ",f10.5)')rcrec
      write(6,'(" Recommended g-space cutoff  : ",f10.5)')gcrec
      write(6,'(/," Using Raggio                : ",f10.5)')raggio
      write(6,'(" Using r-space cutoff        : ",f10.5)')rcut(1)
      write(6,'(" Using g-space cutoff        : ",f10.5)')gcut  
      write(6,'(/," Secondary r-space cutoff    : ",f10.5)')rcut(2)
      write(6,'(" Tertiary  r-space cutoff    : ",f10.5)')rcut(3)
      write(6,'(" Quaternary r-space cutoff   : ",f10.5)')rcut(4)
      
      netcharge = 0.0d0
      do i=1,this%n_types
         netcharge = netcharge + z(i)*dfloat(ntype(i))
         if (dabs(z(i)).gt.1.d-10) tz = .true.
      end do
      
      write(6,'(/," Charges : ",100(e13.5))') z
      write(6,'(/," Net charge ",e13.5)') netcharge
      write(6,*)
      
      ! Short range parameters
      tsr = .true.
      alphaij = 0.0_dp
      bij = 0.0_dp
      cij = 0.0_dp
      dij = 0.0_dp
      eij = 0.0_dp
      nij = 0.0_dp
      b_tt = 0.0_dp
      gamma_ms = 0.d0
      r_ms = 0.0d0
      d_ms= 0.0d0
      do ti=1,this%n_types
         do tj=1,this%n_types
            d_ms(ti,tj) = this%d_ms(ti,tj)
            gamma_ms(ti,tj) = this%gamma_ms(ti,tj)
            R_ms(ti,tj) = this%R_ms(ti,tj)
         end do
      end do

      ! Polarisation
      betapol = this%betapol
      maxipol = this%maxipol
      tolpol = this%tolpol
      pred_order = this%pred_order

      tpol = any(dabs(pol) > 1.0e-6_dp)
      call Print('Polarisation: '//pol)

      bpol = 0.0_dp
      cpol = 0.0_dp
      do ti=1,this%n_types
         do tj=1,this%n_types
            bpol(ti,tj) = this%B_pol(ti,tj)
            cpol(ti,tj) = this%C_pol(ti, tj)
         end do
      end do

      ! Smoothing
      smooth = .false.
      
      ! Yukawa parameters
      yukalpha = this%yukalpha
      yuksmoothlength = this%yuksmoothlength
      tdip_sr = this%tdip_sr
      
      asap_initialised = .true.
   end if

   spind = 0
   do i=1,this%n_types
      where(at%Z == this%atomic_num(i)) spind = i
   end do

   ! ASAP uses atomic units - lengths are in Bohr, energies in Hartree,
   ! forces in Hartree/Bohr and stress in Hartree/Bohr**3

   r(:,:) = at%pos(:,:)/BOHR         ! positions
   rm(:,:) = at%pos(:,:)/BOHR
   htm = transpose(at%lattice/BOHR)  ! lattice
   ht = transpose(at%lattice/BOHR) 
   restart = do_restart   ! reset electric field
   
   objind = 1
   numobj(1) = nat
   ntype = 1
   numobjsp = 0
   do i=1,nat
      ntype(spind(i)) = ntype(spind(i)) + 1
      numobjsp(1,spind(i)) = numobjsp(1,spind(i)) + 1
   end do

   call inv3(ht,htm1,omega)

   if (all(iesr == -1)) then
      ! given rcut(1) and lattice, compute iesr
      call fit_box_in_cell(rcut(1),rcut(1),rcut(1), at%lattice, iesr(1), iesr(2), iesr(3))
      iesr = iesr/2
   end if

   nesr = (2*iesr(1)+1)*(2*iesr(2)+1)*(2*iesr(3)+1)
   nnatmax = min(nnatmax,nat*nesr)
   nnatmax = nnatmax + 10
  
   if (allocated(nnat)) deallocate(nnat)
   if (allocated(nnat2)) deallocate(nnat2)
   if (allocated(nnat3)) deallocate(nnat3)
   if (allocated(nnlist)) deallocate(nnlist)
   if (allocated(nnlist2)) deallocate(nnlist2)
   if (allocated(nnlist3)) deallocate(nnlist3)
   if (allocated(nn_imag)) deallocate(nn_imag)
   if (allocated(nn_imag2)) deallocate(nn_imag2)
   if (allocated(nn_imag3)) deallocate(nn_imag3)

   call label_elements
   call init
   call prepare_traj

   if((mod(it,ntnlist(1)).eq.1).or.(ntnlist(1).eq.1)) then
      call nbrlist(r,nnatmax)
      call nbrlist3(r)
   endif

   if(mod(it,ntnlist(2)).eq.1) then
      call nbrlist3(r)
   endif

   call force_ft(asap_e, asap_f, asap_stress, 0)
   it = it + 1

   ! Convert to {eV, A, fs} units
   if (present(e)) e = asap_e*HARTREE
   if (present(f)) f = asap_f*(HARTREE/BOHR)
   if (present(virial)) virial = asap_stress*(HARTREE/(BOHR**3))*cell_volume(at)
   if (present(local_e)) local_e = 0.0_dp

   if (calc_dipoles .and. tpol) then
      if (.not. has_property(at, 'dipoles')) call add_property(at, 'dipoles', 0.0_dp, n_cols=3)
      if (.not. assign_pointer(at, 'dipoles', dipoles_ptr)) call system_abort('IPModel_ASAP_calc: assign_pointer dipoles failed')
      dipoles_ptr = dip
   end if
   
   deallocate(asap_f)
#else
  call system_abort('ASAP potential is not compiled in. Recompile with HAVE_ASAP=1')
#endif

end subroutine IPModel_ASAP_Calc


subroutine IPModel_ASAP_Print(this, file)
  use libatoms_module
  type(IPModel_ASAP), intent(in) :: this
  type(Inoutput), intent(inout),optional :: file

  integer :: ti, tj

  call Print("IPModel_ASAP : ASAP Potential", file=file)
  call Print("IPModel_ASAP : n_types = " // this%n_types //" n_atoms = "//this%n_atoms, file=file)
  call Print("IPModel_ASAP : betapol = "//this%betapol//" maxipol = "//this%maxipol//" tolpol = "//this%tolpol//" pred_order = "//this%pred_order, file=file)
  call Print("IPModel_ASAP : yukalpha = "//this%yukalpha//" yuksmoothlength = "//this%yuksmoothlength, file=file)

  do ti=1, this%n_types
    call Print ("IPModel_ASAP : type " // ti // " atomic_num " // this%atomic_num(ti), file=file)
    call Print ("IPModel_ASAP : pol = "//this%pol(ti), file=file)
    call Print ("IPModel_ASAP : z   = "//this%z(ti), file=file)
   call verbosity_push_decrement()
    do tj =1,this%n_types
       call Print ("IPModel_ASAP : pair interaction ti tj " // ti // " " // tj // " Zi Zj " // this%atomic_num(ti) //&
            " " // this%atomic_num(tj), file=file)
       call Print ("IPModel_ASAP : pair " // this%D_ms(ti,tj) // " " // this%gamma_ms(ti,tj) // " " &
            // this%R_ms(ti,tj) // " " // this%B_pol(ti,tj) // " " // this%C_pol(ti, tj), file=file)
    end do
   call verbosity_pop()
  end do

end subroutine IPModel_ASAP_Print

subroutine IPModel_ASAP_read_params_xml(this, param_str)
  use libatoms_module
  type(IPModel_ASAP), intent(inout), target :: this
  character(len=*), intent(in) :: param_str

  type(xml_t) :: fxml

  if (len(trim(param_str)) <= 0) return

  parse_in_ip = .false. 
  parse_matched_label = .false.
  parse_ip => this

  call open_xml_string(fxml, param_str)
  call parse(fxml,  &
    startElement_handler = IPModel_startElement_handler, &
    endElement_handler = IPModel_endElement_handler)
  call close_xml_t(fxml)

  if (this%n_types == 0) then
    call system_abort("IPModel_ASAP_read_params_xml parsed file, but n_types = 0")
  endif

end subroutine IPModel_ASAP_read_params_xml

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X 
!% XML param reader functions
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
subroutine IPModel_startElement_handler(URI, localname, name, attributes)
  use libatoms_module
  character(len=*), intent(in)   :: URI
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name
  type(dictionary_t), intent(in) :: attributes

  integer :: status
  character(len=FIELD_LENGTH) :: val

  integer ti, tj, Zi, Zj

  if (name == 'ASAP_params') then ! new ASAP stanza

    if (parse_matched_label) return ! we already found an exact match for this label

    call Quip_fox_get_value(attributes, 'label', val, status)
    if (status /= 0) val = ''

    if (len(trim(parse_ip%label)) > 0) then ! we were passed in a label
      if (val == parse_ip%label) then ! exact match
        parse_matched_label = .true.
        parse_in_ip = .true.
      else ! no match
        parse_in_ip = .false.
      endif
    else ! no label passed in
      parse_in_ip = .true.
    endif

    if (parse_in_ip) then
      if (parse_ip%n_types /= 0) then
        call finalise(parse_ip)
      endif

      call Quip_fox_get_value(attributes, 'n_types', val, status)
      if (status == 0) then
        read (val, *), parse_ip%n_types
      else
        call system_abort("Can't find n_types in ASAP_params")
      endif

      allocate(parse_ip%atomic_num(parse_ip%n_types))
      parse_ip%atomic_num = 0

      allocate(parse_ip%pol(parse_ip%n_types))
      parse_ip%pol = 0.0_dp
      allocate(parse_ip%z(parse_ip%n_types))

      allocate(parse_ip%D_ms(parse_ip%n_types,parse_ip%n_types))
      parse_ip%D_ms = 0.0_dp
      allocate(parse_ip%gamma_ms(parse_ip%n_types,parse_ip%n_types))
      parse_ip%gamma_ms = 0.0_dp
      allocate(parse_ip%R_ms(parse_ip%n_types,parse_ip%n_types))
      parse_ip%R_ms = 0.0_dp
      allocate(parse_ip%B_pol(parse_ip%n_types,parse_ip%n_types))
      parse_ip%B_pol = 0.0_dp
      allocate(parse_ip%C_pol(parse_ip%n_types,parse_ip%n_types))
      parse_ip%C_pol = 0.0_dp

      call Quip_fox_get_value(attributes, "cutoff", val, status)
      if (status /= 0) call system_abort ("IPModel_ASAP_read_params_xml cannot find cutoff")
      read (val, *) parse_ip%cutoff

      call Quip_fox_get_value(attributes, "betapol", val, status)
      if (status == 0) read (val, *) parse_ip%betapol

      call Quip_fox_get_value(attributes, "maxipol", val, status)
      if (status == 0) read (val, *) parse_ip%maxipol

      call Quip_fox_get_value(attributes, "tolpol", val, status)
      if (status == 0) read (val, *) parse_ip%tolpol

      call Quip_fox_get_value(attributes, "pred_order", val, status)
      if (status == 0) read (val, *) parse_ip%pred_order

      call Quip_fox_get_value(attributes, "yukalpha", val, status)
      if (status == 0) read (val, *) parse_ip%yukalpha

      call Quip_fox_get_value(attributes, "yuksmoothlength", val, status)
      if (status == 0) read (val, *) parse_ip%yuksmoothlength

      parse_ip%tewald = .false.
      call Quip_fox_get_value(attributes, "tewald", val, status)
      if (status == 0) read (val, *), parse_ip%tewald

      parse_ip%tdip_sr = .true.
      call Quip_fox_get_value(attributes, "tdip_sr", val, status)
      if (status == 0) read (val, *), parse_ip%tdip_sr

      parse_ip%raggio = 0.0_dp
      call Quip_fox_get_value(attributes, "raggio", val, status)
      if (status == 0) read (val, *), parse_ip%raggio

      parse_ip%a_ew = 0.0_dp
      call Quip_fox_get_value(attributes, "a_ew", val, status)
      if (status == 0) read (val, *), parse_ip%a_ew

      parse_ip%gcut = 0.0_dp
      call Quip_fox_get_value(attributes, "gcut", val, status)
      if (status == 0) read (val, *), parse_ip%gcut

      parse_ip%iesr = 0
      call Quip_fox_get_value(attributes, "iesr", val, status)
      if (status == 0) read (val, *), parse_ip%iesr

    endif

  elseif (parse_in_ip .and. name == 'per_type_data') then

    call Quip_fox_get_value(attributes, "type", val, status)
    if (status /= 0) call system_abort ("IPModel_ASAP_read_params_xml cannot find type")
    read (val, *) ti

    call Quip_fox_get_value(attributes, "atomic_num", val, status)
    if (status /= 0) call system_abort ("IPModel_ASAP_read_params_xml cannot find atomic_num")
    read (val, *) parse_ip%atomic_num(ti)

    call Quip_fox_get_value(attributes, "pol", val, status)
    if (status /= 0) call system_abort ("IPModel_ASAP_read_params_xml cannot find pol")
    read (val, *) parse_ip%pol(ti)

    call Quip_fox_get_value(attributes, "z", val, status)
    if (status /= 0) call system_abort ("IPModel_ASAP_read_params_xml cannot find z")
    read (val, *) parse_ip%z(ti)

    if (allocated(parse_ip%type_of_atomic_num)) deallocate(parse_ip%type_of_atomic_num)
    allocate(parse_ip%type_of_atomic_num(maxval(parse_ip%atomic_num)))
    parse_ip%type_of_atomic_num = 0
    do ti=1, parse_ip%n_types
      if (parse_ip%atomic_num(ti) > 0) &
        parse_ip%type_of_atomic_num(parse_ip%atomic_num(ti)) = ti
    end do

  elseif (parse_in_ip .and. name == 'per_pair_data') then

    call Quip_fox_get_value(attributes, "atnum_i", val, status)
    if (status /= 0) call system_abort ("IPModel_SW_read_params_xml cannot find atnum_i")
    read (val, *) Zi
    call Quip_fox_get_value(attributes, "atnum_j", val, status)
    if (status /= 0) call system_abort ("IPModel_SW_read_params_xml cannot find atnum_j")
    read (val, *) Zj

    ti = get_type(parse_ip%type_of_atomic_num,Zi)
    tj = get_type(parse_ip%type_of_atomic_num,Zj)

    call Quip_fox_get_value(attributes, "D_ms", val, status)
    if (status /= 0) call system_abort ("IPModel_ASAP_read_params_xml cannot find D_ms")
    read (val, *) parse_ip%D_ms(ti,tj)
    call Quip_fox_get_value(attributes, "gamma_ms", val, status)
    if (status /= 0) call system_abort ("IPModel_ASAP_read_params_xml cannot find gamma_ms")
    read (val, *) parse_ip%gamma_ms(ti,tj)
    call Quip_fox_get_value(attributes, "R_ms", val, status)
    if (status /= 0) call system_abort ("IPModel_ASAP_read_params_xml cannot find R_ms")
    read (val, *) parse_ip%R_ms(ti,tj)
    call Quip_fox_get_value(attributes, "B_pol", val, status)
    if (status /= 0) call system_abort ("IPModel_ASAP_read_params_xml cannot find B_pol")
    read (val, *) parse_ip%B_pol(ti,tj)
    call Quip_fox_get_value(attributes, "C_pol", val, status)
    if (status /= 0) call system_abort ("IPModel_ASAP_read_params_xml cannot find C_pol")
    read (val, *) parse_ip%C_pol(ti,tj)

    if (ti /= tj) then
      parse_ip%D_ms(tj,ti) = parse_ip%D_ms(ti,tj)
      parse_ip%gamma_ms(tj,ti) = parse_ip%gamma_ms(ti,tj)
      parse_ip%R_ms(tj,ti) = parse_ip%R_ms(ti,tj)
      parse_ip%B_pol(tj,ti) = parse_ip%B_pol(ti,tj)
      parse_ip%C_pol(tj,ti) = parse_ip%C_pol(ti,tj)
    endif

  endif

end subroutine IPModel_startElement_handler

subroutine IPModel_endElement_handler(URI, localname, name)
  use libatoms_module
  character(len=*), intent(in)   :: URI
  character(len=*), intent(in)   :: localname
  character(len=*), intent(in)   :: name

  if (parse_in_ip) then
    if (name == 'ASAP_params') then
      parse_in_ip = .false.
    end if
  endif

end subroutine IPModel_endElement_handler

end module IPModel_ASAP_module
