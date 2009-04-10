      module ts_module
      
      use libatoms_module
      use QUIP_module
      use MetaPotential_module

      implicit none

      type Nudged_Elastic_Band
        real(dp)        :: spring_constant
        logical         :: climbing 
      end type Nudged_Elastic_Band 

      type String_method 
        integer         :: rep
        integer         :: freq_rep
      end type String_method 

      type Images 
        type(Atoms) :: at
        logical     :: mobile
        real(dp)    :: energy
        real(dp)    :: frac 
      end type Images 
      
      type Chain_of_states 
        integer         :: N                          !% Number of images belonging to the elastic band
        type(Images), allocatable, dimension(:)    :: image   !% images 
      end type Chain_of_states 
      
      type TS
        type (Chain_of_states)     :: cos 
        logical                    :: lnewtangent 
        logical                    :: lneb
        logical                    :: lsm
        type (Nudged_Elastic_Band) :: neb
        type (String_method)       :: sm 
        real(dp)                   :: gfac
      end type TS 
 
      interface initialise 
        module procedure TS_Initialise, TS_Initialise_interp,  TS_Initialise_nointerp
      end interface initialise 

      interface Finalise
        module procedure TS_Finalise
      end interface Finalise
      
      interface Calc
        module procedure TS_Calc
      end interface Calc

      interface Print 
        module procedure TS_print_energy, TS_print
      end interface Print 

      interface fix 
        module procedure TS_fix
      end interface fix 

      contains

!--------------------------------------------------------------------------      
      subroutine TS_Initialise(this,method,freq_rep,lclimbing) 
       type(TS), intent(inout) :: this
       character(len=*)        :: method
       integer, optional       :: freq_rep
       logical, optional       :: lclimbing

       select case(method)

       case("neb")
          this%lneb = .true.
          this%lnewtangent = .true.
          this%neb%spring_constant = 0.1d0
          if(present(lclimbing)) then
            this%neb%climbing        = lclimbing 
          else
            this%neb%climbing        = .false. 
          endif
       case("sm")
          this%lsm = .true.
          this%lnewtangent = .true.
          if(present(freq_rep)) then
            this%sm%freq_rep = freq_rep 
          else
            this%sm%freq_rep = 5
          endif
       case("default")
          stop "Warning:: I do not recognise the method!"   

       end select 

       this%gfac = 0.1d0

      end subroutine TS_Initialise
!--------------------------------------------------------------------------
      subroutine TS_Initialise_interp(this,method,at_in,at_fin,interpolating_method,lmobile_first,lmobile_last,freq_rep,gfac,lclimbing)
        type(TS), intent(inout)               :: this 
        character(len=*)                      :: method
        type(Atoms), intent(inout)            :: at_in 
        type(Atoms), intent(inout)            :: at_fin
        character(*), optional, intent(inout) :: interpolating_method 
        logical, optional, intent(inout)      :: lmobile_first,lmobile_last, lclimbing 
        integer, optional, intent(inout)      :: freq_rep
        real(dp), optional, intent(inout)     :: gfac 
        integer                               :: im
        logical, allocatable                  :: lmobile(:)
    
        call finalise(this) 

        call Initialise(this,method,freq_rep,lclimbing)
     
        allocate(this%cos%image(this%cos%N))
        this%cos%image(1:this%cos%N)%mobile = .true.

        call initialise(this%cos%image(1)%at,at_in%N,at_in%lattice)
        call initialise(this%cos%image(this%cos%N)%at,at_fin%N,at_fin%lattice)
        
        this%cos%image(1)%at          = at_in
        this%cos%image(this%cos%N)%at = at_fin
        
        do im = 2, this%cos%N-1
          call initialise(this%cos%image(im)%at,at_in%N,at_in%lattice)
          this%cos%image(im)%at = at_in
        enddo 

        if(this%cos%image(1)%at%N.ne. this%cos%image(this%cos%N)%at%N) then
          call system_abort('NEB Initialise: # of Atoms in first and last images are different')
        endif
     
        if(present(interpolating_method)) then
          call Print ('Interpolating method --> ' // adjustl(trim(interpolating_method)) ) 
          call Interpolate_images(this,interpolating_method)
        else
          call Print ('Interpolating method --> ' // 'Linear interpolation' ) 
          call Interpolate_images(this)
        endif

        allocate(lmobile(this%cos%N))
        lmobile          = .true.
        if(present(lmobile_first)) lmobile(1)          = lmobile_first
        if(present(lmobile_last))  lmobile(this%cos%N) = lmobile_last
        call fix(this,lmobile)
        deallocate(lmobile)

        if(present(gfac)) then 
            this%gfac        = gfac
        else
            this%gfac        = 0.1_dp
        endif
        
        call TS_initilise_fractionate_chain(this)
 
      end subroutine TS_Initialise_interp
!-----------------------------------------------------------------------------
      subroutine TS_Initialise_nointerp(this,method,at_ref,conf,lmobile_first,lmobile_last,freq_rep,gfac,lclimbing)
        type(TS), intent(inout)               :: this 
        type(Atoms), intent(in)               :: at_ref
        real(dp), dimension(this%cos%N,3*at_ref%N) :: conf 
        integer                               :: im
        character(len=120)                    :: method
        logical, optional, intent(inout)      :: lmobile_first,lmobile_last,lclimbing 
        integer, optional, intent(inout)      :: freq_rep
        real(dp), optional, intent(inout)     :: gfac 
        logical, allocatable                  :: lmobile(:)

        call finalise(this)

        call Initialise(this,method,freq_rep,lclimbing)

        allocate(this%cos%image(this%cos%N))
        this%cos%image(:)%mobile = .true.

        do im = 1, this%cos%N
          call initialise(this%cos%image(im)%at,at_ref%N,at_ref%lattice)
          this%cos%image(im)%at = at_ref
          this%cos%image(im)%at%pos = reshape(conf(im,:), (/3, at_ref%N/))   
          this%cos%image(im)%at%Z(:)   = at_ref%Z(:) 
          this%cos%image(im)%at%cutoff = at_ref%cutoff
          this%cos%image(im)%at%use_uniform_cutoff  = at_ref%use_uniform_cutoff
          if(associated(at_ref%move_mask)) &
                  this%cos%image(im)%at%move_mask = at_ref%move_mask 
        enddo

        allocate(lmobile(this%cos%N))
        lmobile          = .true.
        if(present(lmobile_first)) lmobile(1)          = lmobile_first
        if(present(lmobile_last))  lmobile(this%cos%N) = lmobile_last
        call fix(this,lmobile)
        deallocate(lmobile)

        if(present(gfac)) then 
            this%gfac        = gfac
        else
            this%gfac        = 0.1_dp
        endif
        

        call TS_initilise_fractionate_chain(this)
 
      end subroutine TS_Initialise_nointerp
!-----------------------------------------------------------------------------
      subroutine TS_initilise_fractionate_chain(this)
        type(TS), intent(inout) :: this 
        real(dp)                :: dx
       integer                  :: i

        dx = 1.d0 / dble(this%cos%N - 1._dp)
        do i = 1, this%cos%N
          this%cos%image(i)%frac = dble(i - 1._dp) * dx 
        enddo

      end subroutine TS_initilise_fractionate_chain
!-----------------------------------------------------------------------------
      subroutine Interpolate_images(this,interpolating_method)
        type(TS), intent(inout)                :: this
        character(*), optional, intent(inout)  :: interpolating_method
        real(dp), allocatable, dimension(:,:)  :: increm 
        integer                 :: im,i,j 
        real(dp), dimension(3)  :: pos_in, pos_fin
       
        if(allocated(increm)) deallocate(increm)  
        allocate(increm(3,this%cos%image(1)%at%N) )
       
        if(.not.present(interpolating_method)) then
           interpolating_method = ''
        endif
       
        select case(interpolating_method)
       
        case('to do')
       
        case default 

          do i = 1, this%cos%image(1)%at%N
             pos_in = this%cos%image(1)%at%pos(:,i)
             pos_fin = this%cos%image(this%cos%N)%at%pos(:,i)
             increm(:,i) = diff_min_image(this%cos%image(1)%at,pos_in,pos_fin) 
          enddo
          increm = increm / real(this%cos%N-1,dp)

          do im = 2, this%cos%N -1
             this%cos%image(im)%at%pos(:,:)     = this%cos%image(1)%at%pos(:,:) + real(im-1) * increm(:,:)  
             this%cos%image(im)%at%Z(:)         = this%cos%image(1)%at%Z(:) 
             this%cos%image(im)%at%cutoff       = this%cos%image(1)%at%cutoff
             this%cos%image(im)%at%use_uniform_cutoff  = this%cos%image(1)%at%use_uniform_cutoff
             if(associated(this%cos%image(1)%at%move_mask)) &
                  this%cos%image(im)%at%move_mask = this%cos%image(1)%at%move_mask 
          enddo 
        end select

      end subroutine Interpolate_images
!-----------------------------------------------------------------------------
      subroutine TS_Finalise(this)
       type(TS), intent(inout) :: this
       integer :: im

       if(allocated(this%cos%image)) then
          do im = 1, size(this%cos%image(:))
             call finalise(this%cos%image(im)%at)
          enddo
          deallocate(this%cos%image)  
       end if

       this%lneb = .false.
       this%lsm  = .false.
       this%lnewtangent  = .false.
       this%neb%climbing = .false.
        
      end subroutine TS_Finalise 

!------------------------------------------------------------------------------
      subroutine TS_Calc(this,pot,niter,convergence_tol,force_tol,max_steps,nprint,file,args_str)
        type(TS), intent(inout)               :: this 
        type(Potential), intent(inout)        :: pot
        real(dp),     intent(in)              :: convergence_tol, force_tol
        integer,      intent(in)              :: max_steps  !% Maximum number of steps
        integer,   intent(inout)              :: niter 
        integer, optional, intent(inout)      :: nprint
        type(Inoutput), intent(inout),optional:: file
        character(len=STRING_LENGTH), optional:: args_str
        real(dp), dimension(this%cos%N,3,this%cos%image(1)%at%N) :: forces
        real(dp), dimension(this%cos%N,3,this%cos%image(1)%at%N) :: tau 
        real(dp), dimension(this%cos%N,3,this%cos%image(1)%at%N) :: forces_pot, forces_spring 
        real(dp), dimension(this%cos%N,3,this%cos%image(1)%at%N) :: forces2
        real(dp), dimension(this%cos%N,3,this%cos%image(1)%at%N) :: forces_4c
        real(dp), dimension(this%cos%N)                          :: energy, gfac, error
        real(dp)                                                 :: fnorma, gmaxstepfact, rmaxstep, energy_save
        logical                                                  :: lcheck(this%cos%N) 
        integer                                                  :: i, iat, isaddle, is
        real(dp)                                                 :: emax, emin
        character(len=STRING_LENGTH)                             :: string_arg,args_str_calc

        rmaxstep = 0.1d0
        gfac     = this%gfac 
        do i = 1, this%cos%N
          call calc_connect(this%cos%image(i)%at)
          if(present(args_str)) then
           ! write (string_arg,'(a,i0)') 'project=image_',i  
           ! args_str_calc=trim(args_str)//" "//trim(string_arg)//" "//  'index='//i
           ! args_str_calc='index='//i
            call calc(pot,this%cos%image(i)%at,e=this%cos%image(i)%energy,f=forces_pot(i,:,:),args_str=args_str_calc)
          else
            call calc(pot,this%cos%image(i)%at,e=this%cos%image(i)%energy,f=forces_pot(i,:,:))
          endif
          energy(i) = this%cos%image(i)%energy
          if(present(file)) call print(this,i, file)
        enddo

        lcheck = .false.

        steepest_loop : DO niter = 1, max_steps         

          tau = 0._dp
          call calc_tangent(this,this%lnewtangent,energy,tau)

          call calc_force_perp(this, forces_pot, tau, forces)

          if(this%lneb) then
            call calc_spring(this, tau, forces_spring) 
            if(this%neb%climbing) then
              is = isaddle(this) 
              do i = 1, this%cos%N
                if(i.eq.is) then
                    call calc_force_4c(this,is,forces_pot,tau,forces_4c)
                    forces(i,:,:) =   forces_4c(i,:,:)
                else
                    forces(i,:,:) = forces(i,:,:) + forces_spring(i,:,:)
                endif
              enddo
            else
              forces = forces + forces_spring
            endif
          endif

          forces2 = forces * forces

          do i = 1, this%cos%N

             if(.not.this%cos%image(i)%mobile) then 
               if(mod(niter,nprint).eq.0.and.present(file)) call print(this,i,file)
               cycle
             endif
             fnorma  = 0.0_dp
             do iat = 1, this%cos%image(i)%at%N
               fnorma = fnorma + sum(forces2(i,:,iat))
             enddo
             fnorma = dsqrt(fnorma)
             gmaxstepfact = 1.d0
             if (gfac(i) * fnorma > rmaxstep) then
               gmaxstepfact = 1.d0/fnorma/gfac(i) * rmaxstep
             end if
             
             energy_save = energy(i)
             do iat =  1, this%cos%image(i)%at%N
                if(associated(this%cos%image(i)%at%move_mask)) then
                   if(this%cos%image(i)%at%move_mask(iat).eq.0)  cycle 
                end if
               this%cos%image(i)%at%pos(:,iat) =  this%cos%image(i)%at%pos(:,iat) + gfac(i) * forces(i,:,iat) * gmaxstepfact
             enddo 

             if(.not.this%cos%image(i)%mobile) cycle
             call calc_connect(this%cos%image(i)%at)
             
             if(present(args_str)) then
               !write (string_arg,'(a,i0)') 'project=image_',i 
               !args_str_calc=trim(args_str)//" "//trim(string_arg)//" "//  'index='//i
               args_str_calc='index='//i
               call calc(pot,this%cos%image(i)%at,e=this%cos%image(i)%energy,f=forces_pot(i,:,:),args_str=args_str_calc)
             else
               call calc(pot,this%cos%image(i)%at,e=this%cos%image(i)%energy,f=forces_pot(i,:,:))
             endif
             energy(i) = this%cos%image(i)%energy
          
             if(niter.gt.1) call check(this, energy_save, energy(i), convergence_tol, force_tol, fnorma, error(i), lcheck(i))

             if(present(nprint).and.mod(niter,nprint).eq.0) then
               if(present(file)) call print(this,i,file)
             endif
          enddo

          if(all(lcheck(:))) then
             call print ('convergence reached at step ' // niter)
             exit
          endif
     
          if(this%lsm.and. mod(niter,this%sm%freq_rep).eq.0) then
            call calc_reparameterisation(this)
          endif

          if(present(nprint).and.mod(niter,nprint).eq.0) then
             call TS_print_xyz(this, niter)
             is = isaddle(this)
             call print("Iter " // niter // "; Saddle in = " // is // "; Energy error for the saddle = " // error(is))
          endif
!"
        end do steepest_loop
        niter = niter - 1 
                
      end subroutine TS_Calc 
!----------------------------------------------------------------------------
      subroutine calc_tangent(this,lnewtangent,energy, tau)
!----------------------------------------------------------------------------
! Tangent computed according to G. Henkelman, H Jonsson, JCP 113, 9978 (2000)
! Notice: on the first and last image tau = 0.
! It should be = tauP (tauM) but it is necessary to handle when E_0 > E_i
! and E_N > E_N-1
!----------------------------------------------------------------------------
        type(TS), intent(in)     :: this 
        logical,  intent(in)     :: lnewtangent
        real(dp), intent(inout), dimension(this%cos%N,3,this%cos%image(1)%at%N) :: tau
        real(dp), intent(in), dimension(this%cos%N)                             :: energy
        real(dp), dimension(this%cos%image(1)%at%N)                             :: taunorm
        real(dp), dimension(this%cos%N,3,this%cos%image(1)%at%N)                :: tauP, tauM
        real(dp), dimension(3)                                                  :: pos_1, pos_2
        integer                                                                :: im, iat
        real(dp)                                                               :: E1, E2, Emin, Emax

        if(.not.lnewtangent) return

        tau = 0._dp
        do im = 1, this%cos%N 
          do iat = 1, this%cos%image(1)%at%N
            if(im.ne.this%cos%N) then 
              pos_1 = this%cos%image(im+1)%at%pos(:,iat)
              pos_2 = this%cos%image(im)%at%pos(:,iat)
              tauP(im,:,iat) = diff_min_image(this%cos%image(1)%at,pos_2,pos_1) 
            endif
            if(im.ne.1) then
              pos_1 = this%cos%image(im)%at%pos(:,iat)
              pos_2 = this%cos%image(im-1)%at%pos(:,iat)
              tauM(im,:,iat) = diff_min_image(this%cos%image(1)%at,pos_2,pos_1) 
            endif
           enddo
        enddo
        do im = 2, this%cos%N-1
          if(energy(im+1).gt.energy(im).and.energy(im).gt.energy(im-1)) then
             tau(im,1:3,:) = tauP(im,1:3,:)
          elseif(energy(im+1).lt.energy(im).and.energy(im).lt.energy(im-1)) then
             tau(im,1:3,:) = tauM(im,1:3,:)
          else
             E1 = abs(energy(im+1) - energy(im))
             E2 = abs(energy(im-1) - energy(im))
             if(E1.gt.E2) then
               Emax = E1
               Emin = E2
             else
               Emin = E1
               Emax = E2
             endif
             if(energy(im+1).gt.energy(im-1)) then
               tau(im,1:3,:) = tauP(im,1:3,:) * Emax + tauM(im,1:3,:) * Emin
             elseif(energy(im+1).lt.energy(im-1)) then
               tau(im,1:3,:) = tauP(im,1:3,:) * Emin + tauM(im,1:3,:) * Emax
             endif
          endif
          do iat = 1,  this%cos%image(1)%at%N
            taunorm(iat) = sqrt(sum(tau(im,:,iat) * tau(im,:,iat) ) )
            if(taunorm(iat).gt.0.00001_dp) tau(im,:,iat) = tau(im,:,iat) / taunorm(iat)
            if(associated(this%cos%image(1)%at%move_mask)) then
               if(this%cos%image(1)%at%move_mask(iat).eq.0) tau(im,:,iat) = 0._dp
            end if
          enddo
        enddo

      end subroutine calc_tangent

!--------------------------------------------------------------------
!%  Compute the force perpendicular to the path for each image
      subroutine calc_force_perp(this, forces, tau, forces_perp)
        type(TS), intent(in)     :: this 
        real(dp), intent(in),    dimension(this%cos%N,3,this%cos%image(1)%at%N) :: tau
        real(dp), intent(inout), dimension(this%cos%N,3,this%cos%image(1)%at%N) :: forces, forces_perp 
        integer   :: im, iat
        real(dp)  :: ftemp

        do im = 1, this%cos%N
         do iat = 1, this%cos%image(1)%at%N
            ftemp = sum(forces(im,:,iat) * tau(im,:,iat))
            forces_perp(im,:,iat) =  forces(im,:,iat) - ftemp * tau(im,:,iat) 
         enddo
        enddo

      end  subroutine calc_force_perp
!--------------------------------------------------------------------
!%   Compute the effective force changing the sign of the parallel component 
!%   Used for the climbing image method 
      subroutine calc_force_4c(this, im, forces, tau, forces_4c)
        type(TS), intent(in)     :: this
        real(dp), intent(in),    dimension(this%cos%N,3,this%cos%image(1)%at%N) :: tau
        real(dp), intent(inout), dimension(this%cos%N,3,this%cos%image(1)%at%N) :: forces, forces_4c
        integer   :: im, iat
        real(dp)  :: ftemp

        forces_4c = 0._dp
        do iat = 1, this%cos%image(1)%at%N
            ftemp = sum(forces(im,:,iat) * tau(im,:,iat))
            forces_4c(im,:,iat) =  forces(im,:,iat) - 2._dp * ftemp * tau(im,:,iat) 
        enddo

      end  subroutine calc_force_4c
!--------------------------------------------------------------------
!%  Compute the force due to the spring for each image
      subroutine calc_spring(this, tau, forces)
        type(TS), intent(in)     :: this 
        real(dp), intent(in),    dimension(this%cos%N,3,this%cos%image(1)%at%N) :: tau
        real(dp), intent(inout), dimension(this%cos%N,3,this%cos%image(1)%at%N) :: forces
        real(dp), dimension(3)   :: pos_im,pos_imP, pos_imM, diff_p, diff_m
        integer   :: im, iat
        real(dp)  :: ftemp

       forces = 0._dp
       do im = 2, this%cos%N-1 
         do iat = 1, this%cos%image(1)%at%N
           pos_imP = this%cos%image(im+1)%at%pos(:,iat)
           pos_im  = this%cos%image(im)%at%pos(:,iat)
           pos_imM = this%cos%image(im-1)%at%pos(:,iat)
           diff_p = diff_min_image(this%cos%image(1)%at,pos_im,pos_imP) 
           diff_m = diff_min_image(this%cos%image(1)%at,pos_imM,pos_im) 
           ftemp = sqrt(sum(diff_p(:)*diff_p(:))) -  sqrt(sum(diff_m(:)*diff_m(:)))
           forces(im,:,iat) = this%neb%spring_constant * ftemp * tau(im,:,iat)
         enddo
       enddo

      end subroutine calc_spring
!--------------------------------------------------------------------
      subroutine calc_reparameterisation(this)
        type(TS), intent(inout)           :: this 
        real(dp), dimension(this%cos%N,3,this%cos%image(1)%at%N) :: dis, pos_temp_im
        real(dp), dimension(this%cos%N,3,this%cos%image(1)%at%N) :: dis2
        real(dp), dimension(this%cos%N,this%cos%image(1)%at%N)   :: ldist
        real(dp), dimension(3)                                   :: pos_im, pos_imM
        real(dp)  :: psum, psum_temp
        integer   :: nim, nat, i, iat

        nim = this%cos%N 
        nat = this%cos%image(1)%at%N
        dis = 0._dp
        do i=2, nim 
         do iat = 1, nat 
      !     pos_im  = realpos(this%cos%image(i)%at,iat)
      !     pos_imM = realpos(this%cos%image(i-1)%at,iat)
      !     dis(i,:,iat) = pos_im - pos_imM 
             
           pos_im  = this%cos%image(i)%at%pos(:,iat)
           pos_imM = this%cos%image(i-1)%at%pos(:,iat)
           dis(i,:,iat) = diff_min_image(this%cos%image(1)%at,pos_imM,pos_im) 
         enddo
        enddo
      
        psum =0._dp
        ldist = 0._dp
        do iat = 1, nat 
          psum = 0_dp
          do i = 1, nim
            psum = psum + sqrt(sum (dis(i,:,iat) * dis(i,:,iat) ) )
            ldist(i,iat) = psum
            pos_temp_im(i,:,iat) = realpos(this%cos%image(i)%at,iat)
          enddo
          ldist(:,iat) = ldist(:,iat)/ldist(nim,iat) 
        enddo

        do i = 1, nim
          do iat = 1, nat  
             !if(associated(this%cos%image(i)%at%move_mask)) then
             !   if(this%cos%image(i)%at%move_mask(iat).eq.0) cycle
             !end if
            pos_im(1) = interp1(nim,ldist(:,iat),pos_temp_im(:,1,iat),this%cos%image(i)%frac)
            pos_im(2) = interp1(nim,ldist(:,iat),pos_temp_im(:,2,iat),this%cos%image(i)%frac)
            pos_im(3) = interp1(nim,ldist(:,iat),pos_temp_im(:,3,iat),this%cos%image(i)%frac)
            this%cos%image(i)%at%pos(:,iat) = pos_im(:) - (this%cos%image(i)%at%lattice .mult. this%cos%image(i)%at%travel(:,iat)) 
          enddo
        enddo 

      end subroutine calc_reparameterisation
!---------------------------------------------------------------------
      function interp1(n,xa,ya,x)
      integer  :: n
      real(dp) :: xa(n), ya(n)
      real(dp) :: x, y
      real(dp) :: a, b
      integer  :: k, klo, khi
      real(dp) :: interp1

      klo=1
      khi=n
1     if (khi-klo.gt.1) then
       k=(khi+klo)/2
       if(xa(k).gt.x)then
         khi=k
       else
         klo=k
       endif
       goto 1
      endif      

      if(abs(xa(klo)-xa(khi)).lt.0.0000001_dp) then
       interp1 =  (ya(klo) - ya(khi))/2._dp
      else
       b = (ya(klo) - ya(khi))/(xa(klo) - xa(khi))
       a = ya(klo) - b * xa(klo)
       interp1 = a + b * x
      endif

      return
      end function interp1
!--------------------------------------------------------------------
      subroutine check(this, e_old, e_new, toll, force_tol, fnorma, error, lcheck)
       type(TS) :: this
       real(dp) :: e_old, e_new
       real(dp) :: fnorma, toll
       real(dp) :: force_tol, error
       logical  :: lback, lcheck, lcheck_ene, lcheck_force

       lcheck       = .false.
       lcheck_ene   = .false.
       lcheck_force = .false.
       lback        = .false.

       if((e_old-e_new).lt.toll) then 
         lcheck_ene = .true.
         error = e_old-e_new
       else
         error = e_old-e_new
       endif
       if(fnorma/real(this%cos%image(1)%at%N,dp) .lt. force_tol) lcheck_force = .true.
       if(lcheck_ene.and.lcheck_force) lcheck = .true.

      end subroutine check
!------------------------------------------------------------------
      subroutine TS_fix(this,lfix) 
        type(TS) :: this
        logical  :: lfix(this%cos%N)
  
        this%cos%image(:)%mobile = lfix(:)
         
      end subroutine TS_fix
!------------------------------------------------------------------
      subroutine TS_print_xyz(this, iter) 
        type(TS) :: this
        integer  :: iter
        integer  :: im

        do im =1, this%cos%N
            call print_xyz(this%cos%image(im)%at, &
                        'neb.'//im//'.'//iter//'.xyz')
        enddo

      end subroutine TS_print_xyz
!-----------------------------------------------------------------
      subroutine TS_print_energy(this, im, file) 
        type(TS) :: this
        type(Inoutput), intent(inout),optional:: file
        integer  :: im

         call Print ("image= "// im // " Energy= "// this%cos%image(im)%energy, file=file)
         if(im.eq.this%cos%N) call Print ("    ", file=file)

      end subroutine TS_print_energy
!-----------------------------------------------------------------
      subroutine TS_print(this) 
        type(TS) :: this
        
       call print( '       ')
        call Print_title('Parameters')
        call Print('Chain of state with ' // this%cos%N // ' images')
!'
        if(this%cos%image(1)%mobile) then 
             call print ('First image: mobile')
        else
             call print ('First image: frozen')
        endif
        if(this%cos%image(this%cos%N)%mobile) then 
             call print ('Last image : mobile')
        else
             call print ('Last image : frozen')
        endif

       call print( '       ')
       if(this%lneb) then
           call print ('Transition state compute with NEB')
           if(this%neb%climbing) call print ('Climbing image active')
           call print ('Spring Constant ' // this%neb%spring_constant)
           call print ('Climbing Image method ' // this%neb%climbing)
       else
           call print ('Transition state compute with String method')
           call print ('Reparameterisation performed every ' // this%sm%freq_rep // ' steps')
       endif 
       if(this%lnewtangent) call print ('Tangent computed with the newtangent method')

       call print( '       ')
       call print( 'Minimization :')
       call print( 'gfac = ' // this%gfac)
 
      end subroutine TS_print
!------------------------------------------------------------------
      function isaddle(this)
        type(TS) :: this
        real(dp) :: Emax
        integer  :: i
        integer  :: isaddle  

        Emax = -1.d13
        do i  = 1, this%cos%N
          if(this%cos%image(i)%energy.gt.Emax) then
             isaddle = i
             Emax = this%cos%image(i)%energy 
          endif
        enddo
        if(isaddle.gt. this%cos%N) stop 'isaddle function error'
      end function

!-----------------------------------------------------------------
      end module ts_module
