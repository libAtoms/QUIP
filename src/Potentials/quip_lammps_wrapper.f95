module QUIP_LAMMPS_wrapper_module

   use system_module, only : dp, system_initialise, PRINT_SILENT, system_abort, a2s
   use linearalgebra_module
   use connection_module
   use atoms_types_module
   use atoms_module

   use potential_module
   implicit none

   private

   public :: quip_lammps_wrapper, quip_lammps_potential_initialise

   type quip_lammps_potential
      type(Potential), pointer :: pot
   endtype quip_lammps_potential

   contains

   subroutine quip_lammps_wrapper(nlocal, nghost, atomic_numbers, &
      inum, sum_num_neigh, ilist, &
      quip_num_neigh, quip_neigh, lattice, &
      quip_potential, n_quip_potential, quip_x, &
      quip_e, quip_local_e, quip_virial, quip_local_virial, quip_force) bind(c)

      use iso_c_binding, only: c_double, c_int

      integer(kind=c_int), intent(in) :: nlocal                                      ! number of local atoms in process
      integer(kind=c_int), intent(in) :: nghost                                      ! number of ghost atoms in process
      integer(kind=c_int), intent(in) :: inum                                        ! size of ilist, the list of local atoms (should be nlocal in most cases)
      integer(kind=c_int), intent(in) :: sum_num_neigh                               ! number of all neighbours

      integer(kind=c_int), intent(in), dimension(nlocal+nghost) :: atomic_numbers    ! list of atomic numbers of all atoms
      integer(kind=c_int), intent(in), dimension(inum) :: ilist                      ! list of local atoms
      integer(kind=c_int), intent(in), dimension(inum) :: quip_num_neigh             ! number of neighbours of each local atom
      integer(kind=c_int), intent(in), dimension(sum_num_neigh) :: quip_neigh        ! list of neighbours of local atoms, packaged
      real(kind=c_double), intent(in), dimension(3,3) :: lattice                    ! lattice parameters
      integer(kind=c_int), intent(in), dimension(n_quip_potential) :: quip_potential ! integer array transfering the location of the QUIP potential in the memory
      integer(kind=c_int), intent(in) :: n_quip_potential                            ! size of quip_potential
      real(kind=c_double), intent(in), dimension(3,(nlocal+nghost)) :: quip_x       ! atomic positions
      real(kind=c_double), intent(out) :: quip_e                                    ! energy
      real(kind=c_double), intent(out), dimension(nlocal+nghost) :: quip_local_e    ! atomic energies
      real(kind=c_double), intent(out), dimension(3,3) :: quip_virial               ! virial
      real(kind=c_double), intent(out), dimension(9,(nlocal+nghost)) :: quip_local_virial ! atomic virials
      real(kind=c_double), intent(out), dimension(3,(nlocal+nghost)) :: quip_force  ! force

      integer :: i, j, n, nn, ni, i_n1n, j_n2n, error
      integer, save :: nn_guess = 0
      type(atoms), save     :: at
      type(quip_lammps_potential) :: am
      type(Potential), pointer :: pot
      logical, dimension(:), pointer :: local

      if( n_quip_potential == 0 ) then
         call system_abort('quip_lammps_wrapper: quip_potential not initialised')
      else
         am = transfer(quip_potential,am)
         pot => am%pot
      endif

      ! if this region is vacuum, don't do anything.
      if(nlocal > 0) then

         ! Initialise atoms object. If number of atoms does not change, keep the
         ! object. Allocate connection table, be generous with the number of
         ! neighbours.
         if( .not. is_initialised(at) .or. at%N /= (nlocal+nghost) .or. nn_guess < maxval(quip_num_neigh) ) then
            call finalise(at)
            call initialise(at,(nlocal+nghost),lattice)
            call set_cutoff(at,cutoff(pot))
            nn_guess = 2*maxval(quip_num_neigh)
            call connection_fill(at%connect,at%N,at%N,nn_guess=nn_guess) 
         endif

         ! Transport atomic numbers and positions.
         at%Z = atomic_numbers
         at%pos = quip_x

         ! Add atom mask to atoms object. Potentials will only include these atoms
         ! in the calculation.
         call add_property(at,'local',.false.,ptr = local, overwrite=.true., error=error) 

         ! Zero all connections.
         do i = 1, at%N
            at%connect%neighbour1(i)%t%N = 0
            at%connect%neighbour2(i)%t%N = 0
         enddo

         ! Fill connection table.
         nn = 0
         do ni = 1, inum
            ! Look up ni-th atom as atom i.
            i = ilist(ni) + 1
            ! Set local flag - include in potential calculation.
            local(i) = .true.

            i_n1n = 0
            do n = 1, quip_num_neigh(ni)
               nn = nn + 1;
               j = quip_neigh(nn)

               ! j is atom i's n-th neighbour. Fill the connection table for both i
               ! and j at the same time. As LAMMPS repeats periodic images
               ! explicitly, shift is set to 0. For the same reason, the lattice
               ! variable is not really important, except in calculating the
               ! virial.
               if( i <= j ) then
                  i_n1n = i_n1n + 1
                  at%connect%neighbour1(i)%t%N = at%connect%neighbour1(i)%t%N + 1    ! Increase size of neighbour list by one
                  at%connect%neighbour1(i)%t%int(1,i_n1n) = j                        ! The last neighbour is j
                  at%connect%neighbour1(i)%t%int(2:4,i_n1n) = 0                      ! Set the shift to zero
                  at%connect%neighbour1(i)%t%real(1,i_n1n) = norm(at%pos(:,i) - at%pos(:,j))

                  at%connect%neighbour2(j)%t%N = at%connect%neighbour2(j)%t%N + 1 ! Fill the connection for the other atom in the pair.
                  j_n2n = at%connect%neighbour2(j)%t%N
                  at%connect%neighbour2(j)%t%int(1,j_n2n) = i
                  at%connect%neighbour2(j)%t%int(2,j_n2n) = i_n1n
               endif
            enddo
         enddo

         ! Call the QUIP potential.
         call calc(pot,at,energy=quip_e,local_energy=quip_local_e,force=quip_force,virial=quip_virial,local_virial=quip_local_virial,args_str="atom_mask_name=local lammps do_calc_connect=F")
      else
         quip_e = 0.0_dp
         quip_local_e = 0.0_dp
         quip_virial = 0.0_dp
         quip_local_virial = 0.0_dp
         quip_force = 0.0_dp
      endif

   endsubroutine quip_lammps_wrapper

   subroutine quip_lammps_potential_initialise(quip_potential,n_quip_potential,quip_cutoff,quip_file,n_quip_file,quip_string,n_quip_string) bind(c)
      
      use iso_c_binding, only: c_double, c_int, c_char

      integer(kind=c_int), intent(out), dimension(n_quip_potential) :: quip_potential
      integer(kind=c_int), intent(inout) :: n_quip_potential
      real(kind=c_double), intent(out) :: quip_cutoff
      character(kind=c_char), intent(in) :: quip_file(n_quip_file)
      integer(kind=c_int), intent(in) :: n_quip_file
      character(kind=c_char), intent(in) :: quip_string(n_quip_string)
      integer(kind=c_int), intent(in) :: n_quip_string

      type(Potential), target, save :: pot
      type(quip_lammps_potential) :: am

      integer, dimension(1) :: ammold

      ! In this routine we initialise the QUIP potential and set a pointer to
      ! it. The pointer will then be transferred as an integer array to the
      ! wrapper.
      ! Routine has to be called twice. The integer array used to transfer the
      ! pointer will be allocated by LAMMPS. The first call sends the size of
      ! the array, and the second call transfers the pointer via the array.
      if( n_quip_potential == 0 ) then
         call system_initialise(verbosity=PRINT_SILENT)
         call Potential_Filename_Initialise(pot, trim(a2s(quip_string)), trim(a2s(quip_file)))
         am%pot => pot
         n_quip_potential = size(transfer(am,ammold))
      else
         am%pot => pot
         quip_potential = transfer(am,ammold)
      endif

      ! Cutoff returned.
      quip_cutoff = cutoff(pot)

   endsubroutine quip_lammps_potential_initialise

endmodule QUIP_LAMMPS_wrapper_module
