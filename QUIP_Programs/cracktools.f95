module CrackTools_module

  use libAtoms_module
  use QUIP_module

  use CrackParams_module

  implicit none

  !% Print crack slab to XYZ file, using properties defined in 'params%io_print_properties'
  !% or all properties if 'params%io_print_all_properties' is true.
  interface crack_print
     module procedure crack_print_cio
     module procedure crack_print_filename
  end interface

contains

  subroutine crack_fix_pointers(crack_slab, nn, changed_nn, load, move_mask, edge_mask, md_old_changed_nn, &
       old_nn, hybrid, hybrid_mark, u_disp, k_disp)

    type(Atoms), intent(in) :: crack_slab
    real(dp), pointer, dimension(:,:) :: load, k_disp, u_disp
    integer, pointer, dimension(:) :: move_mask, nn, changed_nn, edge_mask, md_old_changed_nn, &
         old_nn, hybrid, hybrid_mark

    if (has_property(crack_slab, 'nn')) then
       if (.not. assign_pointer(crack_slab, 'nn', nn)) &
            call system_abort('nn pointer assignment failed')
    end if

    if (has_property(crack_slab, 'changed_nn')) then
       if (.not. assign_pointer(crack_slab, 'changed_nn', changed_nn)) &
            call system_abort('changed_nn pointer assignment failed')
    end if
  
    if (has_property(crack_slab, 'load')) then
       if (.not. assign_pointer(crack_slab, 'load', load)) &
            call system_abort('load pointer assignment failed')
    end if

    if (has_property(crack_slab, 'move_mask')) then
       if (.not. assign_pointer(crack_slab, 'move_mask', move_mask)) &
            call system_abort('move_mask pointer assignment failed')
    end if

    if (has_property(crack_slab, 'edge_mask')) then
       if (.not. assign_pointer(crack_slab, 'edge_mask', edge_mask)) &
            call system_abort('edge_mask pointer assignment failed')
    end if

    if (has_property(crack_slab, 'md_old_changed_nn')) then
       if (.not. assign_pointer(crack_slab, 'md_old_changed_nn', md_old_changed_nn)) &
            call system_abort('md_old_changed_nn pointer assignment failed')
    end if

    if (has_property(crack_slab, 'old_nn')) then
       if (.not. assign_pointer(crack_slab, 'old_nn', old_nn)) &
            call system_abort('old_nn pointer assignment failed')
    end if
    
    if (has_property(crack_slab, 'hybrid')) then
       if (.not. assign_pointer(crack_slab, 'hybrid', hybrid)) &
            call system_abort('hybrid pointer assignment failed')
    end if

    if (has_property(crack_slab, 'hybrid_mark')) then
       if (.not. assign_pointer(crack_slab, 'hybrid_mark', hybrid_mark)) &
            call system_abort('hybrid_mark pointer assignment failed')
    end if

    if (has_property(crack_slab, 'uniform_disp')) then
       if (.not. assign_pointer(crack_slab, 'uniform_disp', u_disp)) &
            call system_abort('u_disp pointer assignment failed')
    end if
    
    if (has_property(crack_slab, 'k_disp')) then
       if (.not. assign_pointer(crack_slab, 'k_disp', k_disp)) &
            call system_abort('k_disp pointer assignment failed')
    end if

  end subroutine crack_fix_pointers




  !% Parse crack name in the format is (ijk)[lmn], with negative numbers
  !% denoted by a trailing 'b' (short for 'bar'), e.g. '(111)[11b0]'
  !% Axes of crack slab returned as $3\times3$ matrix with columns
  !% $\mathbf{x}$,$\mathbf{y}$,$\mathbf{z}$.
  subroutine crack_parse_name(crackname, axes)
    character(len=*), intent(in) :: crackname
    real(dp), intent(out), dimension(3,3) :: axes

    real(dp), dimension(3) :: x,y,z
    integer :: crack(6), i, j

    j = 1
    do i=1, len(crackname)
       if (crackname(i:i) == 'b') crack(j-1) = -crack(j-1)
       if (verify(crackname(i:i), '0123456789') == 0) then
          crack(j) = ichar(crackname(i:i)) - ichar('0')
          j = j + 1
       end if
    end do

    y = (/ (real(crack(i),dp), i=1,3) /)
    z = (/ (real(crack(i),dp), i=4,6) /)
    x = y .cross. z

    write (line,'(a,3f8.3)') 'x =', x
    call print(line)
    write (line,'(a,3f8.3)') 'y =', y
    call print(line)
    write (line,'(a,3f8.3)') 'z =', z
    call print(line)

    axes(:,1) = x;  axes(:,2) = y;  axes(:,3) = z

  end subroutine crack_parse_name


  !% Convert a string representation of the material to
  !% an atomic number array suitable for passing to 'diamond()', \emph{e.g.}
  !% 'Si' $\to$ '14', 'C' $\to$ '6' and 'SiC' $\to$ '(/14,6/)'.
  subroutine crack_parse_atomic_numbers(str, Z)
    character(len=*), intent(in) :: str
    integer, intent(out) :: Z(:)

    integer i, p

    i = 0
    p=1
    do while (p <= len(trim(str)))
       i = i + 1
       if (i > size(Z)) call system_abort("Too many elements in str='"//trim(str)//"' for size(Z)="//size(Z))
       if (len(trim(str)) > p .and. scan(str(p+1:p+1),'abcdefghijklmnopqrstuvwxzy') == 1)  then ! 2 letter element
          Z(i) = Atomic_Number(str(p:p+1))
          p = p + 1
       else
          Z(i) = Atomic_Number(str(p:p))
       endif
       p = p + 1
    end do

    if (i == 0) then
       call system_abort("parse_atomic_numbers failed to parse anything, str='"//trim(str)//"'")
    else if (i == 1) then
       Z(2:) = Z(1)
    else
       if (i /= size(Z)) then
          call system_abort("parse_atomic_numbers found "//i//" elements, but size(Z)="//size(Z))
       endif
    endif
  end subroutine crack_parse_atomic_numbers


  !% Find the distinct elements in the array 'at_Z' and
  !% return a list of them in the array 'Z'.
  subroutine crack_find_elements(at_Z, Z)
    integer, intent(in) :: at_Z(:)
    integer, intent(out), allocatable :: Z(:)

    integer :: i, ii
    integer :: Zs(size(ElementName))

    Zs = 0
    Zs(at_Z) = 1

    allocate(Z(sum(Zs)))
    ii = 1
    do i=1, size(Zs)
       if (Zs(i) /= 0) then
          Z(ii) = i
          ii = ii + 1
       end if
    end do
  end subroutine crack_find_elements

  !% Calculate energy release rate $G$ from strain using 
  !% $$G = \frac{1}{2} \frac{E}{1-\nu^2} \epsilon^2 h$$
  !% from thin strip result. Quantities are:
  !% 'strain',$\epsilon$, dimensionless ratio $\frac{\Delta y}{y}$;
  !% 'E', $E$,  Young's modulus, GPa;
  !% 'v', $\nu$, Poisson ratio, dimensionless;
  !% 'height', $h$ \AA{}, 10$^{-10}$~m;
  !% 'G', Energy release rate, J/m$^2$.
  function crack_strain_to_g(strain, E, v, height)
    real(dp), intent(in) :: strain, E, v, height
    real(dp) :: crack_strain_to_g

    crack_strain_to_g = 1.0_dp/2.0_dp*E/(1-v*v)*strain*strain*height*0.1_dp

  end function crack_strain_to_g


  !% Calculate $epsilon$ from $G$, inverse of above formula.
  !% Units are as the same as 'crack_strain_to_g'
  function crack_g_to_strain(G, E, v, height)
    real(dp), intent(in) :: G, E, v, height
    real(dp) :: crack_g_to_strain

    crack_g_to_strain = sqrt(2.0_dp*G*(1.0_dp-v*v)/(E*height*0.1_dp))

  end function crack_g_to_strain

  !% Convert from energy release rate $G$ to stress intensity factor $K$
  !% Units: G (J/m$^2$), E (GPa), K (Pa sqrt(m))
  function crack_g_to_k(G, E, v, mode) result(K)
    real(dp), intent(in) :: G, E, v
    character(*), optional, intent(in) :: mode
    real(dp) :: K

    real(dp) :: Ep
    character(20) :: use_mode

    use_mode = optional_default("plane strain", mode)

    if (trim(use_mode) == 'plane stress') then
       Ep = E
    else if (trim(use_mode) == 'plane strain') then
       Ep = E/(1.0_dp-v*v)
    else
       call system_abort('crack_k_to_g: bad mode '//trim(use_mode))
    end if

    K = sqrt(G*Ep*1e9_dp)

  end function crack_g_to_k

  !% Convert from stress intensity factor $K$ to energy release rate $G$
  !% Units: G (J/m$^2$), E (GPa), K (Pa sqrt(m))
  function crack_k_to_g(K, E, v, mode) result(G)
    real(dp), intent(in) :: K, E, v
    character(*), optional, intent(in) :: mode
    real(dp) :: G

    real(dp) :: Ep
    character(20) :: use_mode

    use_mode = optional_default("plane strain", mode)

    if (trim(use_mode) == 'plane stress') then
       Ep = E
    else if (trim(use_mode) == 'plane strain') then
       Ep = E/(1.0_dp-v*v)
    else
       call system_abort('crack_k_to_g: bad mode '//trim(use_mode))
    end if

    G = K*K/(Ep*1e9_dp)

  end function crack_k_to_g


  !% Measure the current height of slab and calculate
  !% energy release rate $G$ from current and original
  !% heights and elastic constants $E$ and $\nu$, using the equation
  !% $$ G = \frac{1}{2} \frac{E}{1-\nu^2} \frac{{h - h_0}^2}{h_0} $$
  !% where $h_0$ is the original height and $h$ the new height.
  !% Otherwise, symbols and units are the same as in 'crack_strain_to_g'.
  function crack_measure_g(at, E, v, orig_height) result(G)
    type(Atoms), intent(in) :: at
    real(dp), intent(in) :: E, v, orig_height
    real(dp) :: G

    real(dp) :: new_height

    new_height = maxval(at%pos(2,:))-minval(at%pos(2,:))
    G = 0.1_dp*0.5_dp*E/(1.0_dp-v*v)*(new_height - orig_height)**2.0_dp/orig_height

  end function crack_measure_g


  !% Rescale atoms in slab, with atoms in front of either crack tip
  !% strained in y direction by 'strain' and atoms behind crack tip
  !% rigidly shifted to keep top and bottom edges flat. A transition
  !% zone is created in between with linearly varying strain to 
  !% avoid creation of defects.
  !%
  !%>  --------------------------------------
  !%>  |           |   |   |   |            |
  !%>  |           |   |___|   |            |
  !%>  |           |   |   |   |            |
  !%>  |           |   |___|   |            |
  !%>  |           |   |   |   |            |
  !%>  |    1      | 2 | 3 | 4 |     5      |
  !%>  --------------------------------------
  !% \begin{center}\begin{tabular}{clc}
  !% \hline \hline
  !% Region &  Positon & Load \\
  !% \hline
  !%   1    &  $x <$ 'l_crack_pos'-'zone_width' & $G$ \\
  !%   2    &  'l_crack_pos'-'zone_width' $\le x <$ 'l_crack_pos' & $G \to 0$ \\
  !%   3    &  'l_crack_pos' $< x <$ 'r_crack_pos'& $0$ \\
  !%   4    &  'r_crack_pos' $< x \le$ 'r_crack_pos'+'zone_width' & $0 \to G$ \\
  !%   5    &  $x >$ 'r_crack_pos'+'zone_width' & $G$ \\
  !% \hline 
  !% \hline
  !% \end{tabular}\end{center}
  subroutine crack_uniform_load(at, l_crack_pos, r_crack_pos, zone_width, G, apply_load)
    type(Atoms), intent(inout) :: at
    real(dp), intent(in) :: l_crack_pos, r_crack_pos, zone_width
    real(dp), intent(inout) :: G
    logical, optional :: apply_load

    integer ::  j
    real(dp) top_old, top_new, bottom_old, bottom_new, x, y, q, strain, E, v, &
         orig_height, new_height
    real(dp), pointer, dimension(:,:) :: disp
    logical :: do_apply_load

    if (.not. get_value(at%params, 'OrigHeight', orig_height)) &
         call system_abort('crack_uniform_load: "OrigHeight" parameter missing')

    if (.not. get_value(at%params, 'YoungsModulus', E)) &
         call system_abort('crack_uniform_load: "YoungsModulus" missing')

    if (.not. get_value(at%params, 'PoissonRatio_yx', v)) &
         call system_abort('crack_uniform_load: "PoissonRatio_yx" missing')

    call add_property(at, 'uniform_disp', 0.0_dp, n_cols=3)
    if (.not. assign_pointer(at, 'uniform_disp', disp)) &
         call system_abort('crack_uniform_load: failed to assign pointer to uniform_disp')

    ! Calculate strain corresponding to given G
    strain = crack_g_to_strain(G, E, v, orig_height)

    call print('Requested G = '//G)
    call Print('Applying strain '//strain)

    top_old = 0.0; top_new = 0.0; bottom_old = 0.0; bottom_new = 0.0

    ! Find top and bottom
    do j=1, at%N
       if (at%pos(2,j) > top_old)       top_old = at%pos(2,j)
       if (at%pos(2,j) < bottom_old) bottom_old = at%pos(2,j);
    end do

    disp = 0.0_dp

    do j=1, at%N
       x = at%pos(1,j);
       y = at%pos(2,j);

       ! Strain regions 1 and 5
       if (x < l_crack_pos-zone_width .or. x > r_crack_pos+zone_width) then
          disp(2,j) = strain*y
       end if

       ! Strain region 2
       if (x >= l_crack_pos-zone_width .and. x < l_crack_pos) then
          q = (x - (l_crack_pos - zone_width))/zone_width
          if (y >= 0.0) then
             disp(2,j) = strain*(y*(1.0-q) + top_old*q)
          else
             disp(2,j) = strain*(y*(1.0-q) + bottom_old*q)
          end if
       end if

       ! Strain region 4
       if (x > r_crack_pos .and. x <= r_crack_pos+zone_width) then
          q = (x - r_crack_pos)/zone_width;
          if (y >= 0.0) then
             disp(2,j) = strain*(y*q + top_old*(1.0-q))
          else
             disp(2,j) = strain*(y*q + bottom_old*(1.0-q))
          end if
       end if

       ! Rigidly shift region 3
       if (x >= l_crack_pos .and. x <= r_crack_pos) then
          if(y >= 0.0) then
             disp(2,j) = strain*top_old
          else
             disp(2,j) = strain*bottom_old
          end if
       end if
    end do

    do_apply_load = optional_default(.false.,apply_load)
    if (do_apply_load) then

       at%pos(2,:) = at%pos(2,:) + disp(2,:)

       new_height = maxval(at%pos(2,:))-minval(at%pos(2,:))
       call print('Measured slab height after strain = '//new_height)
       G = 0.1_dp*0.5_dp*E/(1.0_dp-v*v)*(new_height - orig_height)**2.0_dp/orig_height
       call print('Measured G = '//G)

       call set_value(at%params, 'G', G)
    end if

  end subroutine crack_uniform_load

  !% Calculate Irwin K-field stresses and/or displacements for all atoms in 'at'.
  !% Atomic positions should be the original undistorted bulk crystal positions. 
  !% 'YoungsModulus' and 'PoissonRatio_yx' parameters are extracted from 'at', along
  !% with 'CrackPos' to specify the location of the crack tip. If neither 'sig' nor 'disp'
  !% are present thenn properties are added to at if do_disp or do_sig are true.
  !% Stress is in 6 component Voigt notation: $1=xx, 2=yy, 3=zz, 4=yz, 5=zx$ and $6=xy$, and 
  !% displacement is a Cartesian vector $(u_x,u_y,u_z)$.
  subroutine crack_k_field(at, K, mode, sig, disp, do_sig, do_disp)
    type(Atoms), intent(inout) :: at
    real(dp), intent(in) :: K
    character(*), intent(in), optional :: mode
    real(dp), dimension(:,:), optional, target :: sig, disp
    logical, optional, intent(in) :: do_sig, do_disp

    real(dp), pointer, dimension(:,:) :: psig, pdisp
    real(dp), allocatable, dimension(:,:) :: mysig
    real(dp) :: E, v, r, theta, kappa, vp, vpp, crack_pos(3), pos(3)
    character(20) :: use_mode
    logical :: dummy, my_do_sig, my_do_disp
    integer :: i

    use_mode = optional_default("plane strain", mode)
    my_do_sig  = optional_default(.false., do_sig)
    my_do_disp = optional_default(.false., do_disp)

    if (.not. get_value(at%params, 'YoungsModulus', E)) &
         call system_abort('crack_k_field: "YoungsModulus" missing')

    if (.not. get_value(at%params, 'PoissonRatio_yx', v)) &
         call system_abort('crack_k_field: "PoissonRatio_yx" missing')

    crack_pos = 0.0_dp
    if (.not. get_value(at%params, 'CrackPos', crack_pos(1))) &
         call system_abort('crack_update_selection: CrackPos parameter missing from atoms')

    if (trim(use_mode) == 'plane stress') then
       kappa = (3.0_dp - v)/(1.0_dp + v)
       vp = 0.0_dp
       vpp = v
    else if (trim(use_mode) == 'plane strain') then
       kappa = 3.0_dp - 4.0_dp*v
       vp = v
       vpp = 0.0_dp
    else
       call system_abort('crack_k_field: bad mode '//trim(use_mode))
    end if

    nullify(psig,pdisp)
    allocate(mysig(6,at%N))

    if (present(sig) .or. present(disp)) then
       if (present(sig)) psig => sig
       if (present(disp)) pdisp => disp
    else
       if (my_do_sig) then
          call add_property(at, 'k_sig', 0.0_dp, n_cols=6)
          dummy = assign_pointer(at, 'k_sig', psig)
       end if

       if (my_do_disp) then
          call add_property(at, 'k_disp', 0.0_dp, n_cols=3)    
          dummy = assign_pointer(at, 'k_disp', pdisp)
       end if
    end if

    do i=1,at%N

       pos = at%pos(:,i) - crack_pos
       r = sqrt(pos(1)*pos(1) + pos(2)*pos(2))*1e-10_dp
       theta = atan2(pos(2),pos(1))

       mysig(1,i) = K/sqrt(2.0_dp*pi*r)*cos(theta/2.0_dp)*(1.0_dp - sin(theta/2.0_dp)* &
            sin(3.0_dp*theta/2.0_dp))

       mysig(2,i) = K/sqrt(2.0_dp*pi*r)*cos(theta/2.0_dp)*(1.0_dp + sin(theta/2.0_dp)* &
            sin(3.0_dp*theta/2.0_dp))

       mysig(3,i) = vp*(mysig(1,i) + mysig(2,i))
       mysig(4,i) = 0.0_dp
       mysig(5,i) = 0.0_dp

       mysig(6,i) = K/sqrt(2.0_dp*pi*r)*sin(theta/2.0_dp)*cos(theta/2.0_dp)* &
            cos(3.0_dp*theta/2.0_dp)

       if (associated(pdisp)) then
          pdisp(1,i) = K/(2.0_dp*E*1e9_dp)*sqrt(r/(2.0_dp*pi))*((1.0_dp+v)*(2.0_dp*kappa-1.0_dp)*cos(theta/2.0_dp) - &
               cos(3.0_dp*theta/2.0_dp))/1e-10
          pdisp(2,i) = K/(2.0_dp*E*1e9_dp)*sqrt(r/(2.0_dp*pi))*((1.0_dp+v)*(2.0_dp*kappa-1.0_dp)*sin(theta/2.0_dp) - &
               sin(3.0_dp*theta/2.0_dp))/1e-10
          pdisp(3,i) = -vpp*pos(3)*1e-10_dp/(E*1e9_dp)*(mysig(1,i) + mysig(2,i))/1e-10
       end if

    end do

    if (associated(psig)) psig = mysig/1e9_dp ! Convert to GPa
    deallocate(mysig)

  end subroutine crack_k_field



  subroutine crack_apply_strain_ramp(at, G1, G2, d1, d2, d3)
    type(Atoms), intent(inout) :: at
    real(dp), intent(in) :: G1, G2, d1, d2, d3

    integer ::  j
    real(dp) x, y, q, strain1, strain2, &
         E, v, orig_height, orig_width, new_height, G

    if (.not. get_value(at%params, 'OrigHeight', orig_height)) &
         call system_abort('crack_apply_strain_ramp: "OrigHeight" parameter missing')

    if (.not. get_value(at%params, 'OrigWidth', orig_width)) &
         call system_abort('crack_apply_strain_ramp: "OrigWidth" parameter missing')

    if (.not. get_value(at%params, 'YoungsModulus', E)) &
         call system_abort('crack_apply_strain_ramp: "YoungsModulus" missing')

    if (.not. get_value(at%params, 'PoissonRatio_yx', v)) &
         call system_abort('crack_apply_strain_ramp: "PoissonRatio_yx" missing')

    ! Calculate strains corresponding to given Gs
    strain1 = crack_g_to_strain(G1, E, v, orig_height)
    strain2 = crack_g_to_strain(G2, E, v, orig_height)

    call print('Starting load G1 = '//G1//' --> strain '//strain1)
    call print('Final load    G2 = '//G2//' --> strain '//strain2)


    do j=1, at%N
       x = at%pos(1,j)
       y = at%pos(2,j)

       ! Rigidly shift first region
       if (x <= d1) then
          at%pos(2,j) = at%pos(2,j) + strain1*orig_height/2.0_dp*sign(1.0_dp,y)
       end if

       ! Linearly increasing strain in region 2
       ! Interpolate between 0 at x=d1 and strain1 at x=d2
       ! Add shift to align with first region at x=d1
       if (x > d1 .and. x <= d2) then
          q = (x - d1)/(d2-d1)
          at%pos(2,j) = at%pos(2,j) * (1.0_dp + strain1*q) + &
               strain1*orig_height/2.0_dp*sign(1.0_dp,y)*(1.0_dp-q)
       end if

       ! Decreasing strain in region 3
       ! Interpolate between strain1 at x=d2 and strain2 at x=d3
       if (x > d2 .and. x <= d3) then
          q = (x - d2)/(d3-d2)
          at%pos(2,j) = at%pos(2,j) * (1.0_dp +  (strain1*(1.0_dp-q) + strain2*q))
       end if

       ! Constant strain2 in last region
       if (x > d3) then
          at%pos(2,j) = at%pos(2,j) * (1.0_dp + strain2)
       end if

    end do

    new_height = maxval(at%pos(2,:))-minval(at%pos(2,:))
    call print('Measured slab height after strain = '//new_height)
    G = 0.1_dp*0.5_dp*E/(1.0_dp-v*v)*(new_height - orig_height)**2.0_dp/orig_height
    call print('Measured G = '//G)

    call set_value(at%params, 'G', G)

  end subroutine crack_apply_strain_ramp


  subroutine crack_setup_marks(crack_slab, params)
    type(Atoms) :: crack_slab
    type(CrackParams) :: params

    integer :: i
    real(dp) :: crack_pos

    ! Pointers into Atoms data structure
    real(dp), pointer, dimension(:,:) :: load, k_disp, u_disp
    integer, pointer, dimension(:) :: move_mask, nn, changed_nn, edge_mask, md_old_changed_nn, &
         old_nn, hybrid, hybrid_mark

    call crack_fix_pointers(crack_slab, nn, changed_nn, load, move_mask, edge_mask, md_old_changed_nn, &
         old_nn, hybrid, hybrid_mark, u_disp, k_disp)

    ! Setup edge_mask to allow easy exclusion of edge atoms
    do i=1,crack_slab%N
       if (crack_is_edge_atom(crack_slab, i, params%selection_edge_tol)) &
	    edge_mask(i) = 1
    end do

    ! Calculate connectivity and numbers of nearest neighbours
    call crack_update_connect(crack_slab, params)

    ! Find rightmost undercoordinated atoms in bulk - this is initial crack tip position
    crack_pos = crack_find_crack_pos(crack_slab, params)

    ! Artificially set changed_nn to 1 for atoms near to crack tip
    do i = 1, crack_slab%N
       if (distance_min_image(crack_slab, i, (/crack_pos,0.0_dp,0.0_dp/)) < params%crack_seed_embed_tol) &
	    changed_nn(i) = 1
    end do

    call Print('Seeded embed region with '//count(changed_nn /= 0)//' atoms.')

    call crack_update_selection(crack_slab, params)

  end subroutine crack_setup_marks


  subroutine crack_calc_load_field(crack_slab, params, metapot, load_method, overwrite_pos, mpi)
    type(Atoms), intent(inout) :: crack_slab
    type(CrackParams), intent(in) :: params
    type(MetaPotential), intent(inout) :: metapot
    character(*), intent(in) :: load_method
    logical, intent(in) :: overwrite_pos
    type(MPI_Context), intent(in) :: mpi

    real(dp), pointer, dimension(:,:) :: load, k_disp, u_disp
    integer, pointer, dimension(:) :: move_mask, nn, changed_nn, edge_mask, md_old_changed_nn, &
         old_nn, hybrid, hybrid_mark
    type(CInOutput) :: movie
    real(dp), allocatable :: pos1(:,:), pos2(:,:)
    integer :: i, k, steps
    real(dp) :: G, G1, E, v, v2, Orig_Width, Orig_Height, K1, r, l_crack_pos, r_crack_pos, crack_pos

    if (.not. get_value(crack_slab%params, 'OrigHeight', orig_height)) &
         call system_abort('crack_calc_load_field: "OrigHeight" parameter missing')

    if (.not. get_value(crack_slab%params, 'OrigWidth', orig_width)) &
         call system_abort('crack_calc_load_field: "OrigWidth" parameter missing')

    if (.not. get_value(crack_slab%params, 'YoungsModulus', E)) &
         call system_abort('crack_calc_load_field: "YoungsModulus" missing')

    if (.not. get_value(crack_slab%params, 'PoissonRatio_yx', v)) &
         call system_abort('crack_calc_load_field: "PoissonRatio_yx" missing')

    if (.not. get_value(crack_slab%params, 'PoissonRatio_yz', v2)) &
         call system_abort('crack_calc_load_field: "PoissonRatio_yz" missing')

    if (.not. get_value(crack_slab%params, 'CrackPos', crack_pos)) &
         call system_abort('crack_calc_load_field: "CrackPos" missing')

    if (.not. get_value(crack_slab%params, 'G', G)) &
         call system_abort('crack_calc_load_field: "G" missing')

    call add_property(crack_slab, 'load', 0.0_dp, n_cols=3)

    call crack_fix_pointers(crack_slab, nn, changed_nn, load, move_mask, edge_mask, md_old_changed_nn, &
         old_nn, hybrid, hybrid_mark, u_disp, k_disp)

    if (.not. mpi%active .or. (mpi%active .and.mpi%my_proc == 0)) then
       if (params%io_netcdf) then
          call initialise(movie, 'crack_relax_loading_movie.nc', action=OUTPUT)
       else
          call initialise(movie, 'crack_relax_loading_movie.xyz', action=OUTPUT)
       end if
    end if

    allocate(pos1(3,crack_slab%N), pos2(3,crack_slab%N))

    l_crack_pos = -orig_width
    r_crack_pos = crack_pos
    
!!$    if (trim(params%crack_structure) == 'graphene') then
!!$       r_crack_pos = -orig_width
!!$    end if

    ! Save original positions
    pos1 = crack_slab%pos
    
    if (params%crack_relax_loading_field) then
       ! Geometry optimise
       steps = minim(metapot, crack_slab, method=params%minim_mm_method, convergence_tol=params%minim_mm_tol, &
	     max_steps=params%minim_mm_max_steps, linminroutine=params%minim_mm_linminroutine, do_print=.true., &
	     do_pos=.true.,do_lat=.false., use_fire=(trim(params%minim_mm_method)=='fire'), &
             print_cinoutput=movie)
    end if
    
    ! Save relaxed positions
    pos2 = crack_slab%pos
        
    ! Apply loading field
    
    if (trim(load_method) == 'saved') then

       call print_title('Applying saved load increment')
       call crack_apply_load_increment(crack_slab)

    else if (trim(load_method) == 'uniform') then

       call print_title('Applying uniform load increment')
       ! strain it a bit
       do i=1,crack_slab%N
          crack_slab%pos(2,i) = crack_slab%pos(2,i)*(1.0_dp+params%crack_initial_loading_strain)
       end do

    else if (trim(load_method) == 'kfield') then

       call print_title('Applying K-field load increment')

       if (.not. has_property(crack_slab, 'k_disp')) &
            call system_abort('crack_calc_load_field: kfield mode and k_disp property missing.')

       ! remove old load k_disp
       do i=1,crack_slab%N
          crack_slab%pos(:,i) = crack_slab%pos(:,i) - k_disp(:,i)
       end do

       ! apply load increment
       K1 = crack_g_to_k(crack_strain_to_g( &
            crack_g_to_strain(G, E, v, orig_height) + &
            params%crack_initial_loading_strain, E, v, orig_height),E,v)

       call print('Stress Intensity Factor K_1 = '//(K1/1e6_dp)//' MPa.sqrt(m)')
       call crack_k_field(crack_slab, K1, do_disp=.true.)
       call crack_fix_pointers(crack_slab, nn, changed_nn, load, move_mask, edge_mask, md_old_changed_nn, &
            old_nn, hybrid, hybrid_mark, u_disp, k_disp)

       do i=1,crack_slab%N
          crack_slab%pos(:,i) = crack_slab%pos(:,i) + k_disp(:,i)
       end do
       
    else if (trim(load_method) == 'interp_kfield_uniform') then

       call print_title('Applying interp. kfield-uniform load increment')

       if (.not. has_property(crack_slab, 'uniform_disp')) &
            call system_abort('crack_calc_load_field: interp_kfield_uniform mode and uniform_disp property missing.')
       
       if (.not. has_property(crack_slab, 'k_disp')) &
            call system_abort('crack_calc_load_field: interp_kfield_uniform mode and k_disp property missing.')

       ! remove old load u_disp + k_disp
       do i=1,crack_slab%N
          r = sqrt((crack_slab%pos(1,i) - crack_pos)**2.0_dp + &
               crack_slab%pos(2,i)**2.0_dp)
          if (r > params%crack_load_interp_length) then
             crack_slab%pos(:,i) = crack_slab%pos(:,i) - u_disp(:,i)
          else
             do k=1,3
                crack_slab%pos(k,i) = crack_slab%pos(k,i) -  &
                     linear_interpolate(0.0_dp, k_disp(k,i), params%crack_load_interp_length, u_disp(k,i), r)
             end do
          end if
       end do

       ! apply load increment to u_disp
       G1 = crack_strain_to_g( &
            crack_g_to_strain(G, E, v, orig_height) + &
            params%crack_initial_loading_strain, E, v, orig_height)

       call crack_uniform_load(crack_slab, l_crack_pos, r_crack_pos, &
            params%crack_strain_zone_width, G1, apply_load=.false.)

       ! apply load increment to K_disp
       K1 = crack_g_to_k(G1,E,v)

       call print('Energy release rate     G_1 = '//G1//' J/m^2')
       call print('Stress Intensity Factor K_1 = '//(K1/1e6_dp)//' MPa.sqrt(m)')
       call crack_k_field(crack_slab, K1, do_sig=.false., do_disp=.true.)
       call crack_fix_pointers(crack_slab, nn, changed_nn, load, move_mask, edge_mask, md_old_changed_nn, &
            old_nn, hybrid, hybrid_mark, u_disp, k_disp)

       do i=1,crack_slab%N
          r = sqrt((crack_slab%pos(1,i) - crack_pos)**2.0_dp + &
               crack_slab%pos(2,i)**2.0_dp)
          if (r > params%crack_load_interp_length) then
             crack_slab%pos(:,i) = crack_slab%pos(:,i) + u_disp(:,i)
          else
             do k=1,3
                crack_slab%pos(k,i) = crack_slab%pos(k,i) +  &
                     linear_interpolate(0.0_dp, k_disp(k,i), params%crack_load_interp_length, u_disp(k,i), r)
             end do
          end if
       end do
       
    end if
    
    if (params%crack_relax_loading_field) then
       
       ! now re-relax
       steps = minim(metapot, crack_slab, method=params%minim_mm_method, convergence_tol=params%minim_mm_tol, &
            max_steps=params%minim_mm_max_steps, linminroutine=params%minim_mm_linminroutine, do_print=.true., &
            do_pos=.true.,do_lat=.false., use_fire=(trim(params%minim_mm_method)=='fire'), &
            print_cinoutput=movie)
    end if
    
    ! work out displacement field, using relaxed positions
    do i=1,crack_slab%N
       load(:,i) = crack_slab%pos(:,i) - pos2(:,i)
    end do
    
    call print('Displacement field generated. Max disp: '//maxval(load))
    call print('                              RMS disp: '//sqrt(norm2(reshape(load,(/3*crack_slab%N/)))/(3.0_dp*crack_slab%N)))
     
    if (overwrite_pos) then
       crack_slab%pos = pos2
    else
       crack_slab%pos = pos1
    end if
    call calc_connect(crack_slab)

    deallocate(pos1,pos2)
    call finalise(movie)
    
  end subroutine crack_calc_load_field


  subroutine crack_make_seed(crack_slab, params) 
    type(Atoms), intent(inout) :: crack_slab
    type(CrackParams), intent(in) :: params

    ! Pointers into Atoms data structure
    real(dp), pointer, dimension(:,:) :: load, k_disp, u_disp
    integer, pointer, dimension(:) :: move_mask, nn, changed_nn, edge_mask, md_old_changed_nn, &
         old_nn, hybrid, hybrid_mark
    real(dp) :: G, E, v, v2, Orig_Width, Orig_Height,  r, l_crack_pos, r_crack_pos, strain
    integer :: i, k

    call crack_fix_pointers(crack_slab, nn, changed_nn, load, move_mask, edge_mask, md_old_changed_nn, &
         old_nn, hybrid, hybrid_mark, u_disp, k_disp)

    if (.not. get_value(crack_slab%params, 'OrigHeight', orig_height)) &
         call system_abort('crack_make_seed: "OrigHeight" parameter missing')

    if (.not. get_value(crack_slab%params, 'OrigWidth', orig_width)) &
         call system_abort('crack_make_seed: "OrigWidth" parameter missing')

    if (.not. get_value(crack_slab%params, 'YoungsModulus', E)) &
         call system_abort('crack_make_seed: "YoungsModulus" missing')

    if (.not. get_value(crack_slab%params, 'PoissonRatio_yx', v)) &
         call system_abort('crack_make_seed: "PoissonRatio_yx" missing')

    if (.not. get_value(crack_slab%params, 'PoissonRatio_yz', v2)) &
         call system_abort('crack_make_seed: "PoissonRatio_yz" missing')

    ! Determine position of seed crack
    l_crack_pos = -orig_width ! single ended crack
    r_crack_pos = -orig_width/2.0_dp + params%crack_seed_length
    
!!$    if (trim(params%crack_structure) == 'graphene') then
!!$       r_crack_pos = -orig_width
!!$    end if

    if (trim(params%crack_loading) == 'uniform') then
       if (params%crack_G > 0.0_dp) then

          call print_title('Seed crack - Uniform Load')

          call add_property(crack_slab, 'uniform_disp', 0.0_dp, n_cols=3)
          call crack_fix_pointers(crack_slab, nn, changed_nn, load, move_mask, edge_mask, md_old_changed_nn, &
               old_nn, hybrid, hybrid_mark, u_disp, k_disp)

          G = params%crack_G
          call crack_uniform_load(crack_slab, l_crack_pos, r_crack_pos, &
               params%crack_strain_zone_width, G, apply_load=.true.)
          call set_value(crack_slab%params, 'G', G)
          call set_value(crack_slab%params, 'CrackPos', r_crack_pos + 0.85_dp*params%crack_strain_zone_width)
          call set_value(crack_slab%params, 'OrigCrackPos', r_crack_pos + 0.85_dp*params%crack_strain_zone_width)

          if(params%crack_rescale_x_z) then
             !  Rescale in x direction by v and in z direction by v2
             if (.not. get_value(crack_slab%params,'OrigHeight',orig_height)) orig_height = 0.0_dp
             strain = crack_g_to_strain(params%crack_G, E, v, orig_height)
             crack_slab%pos(1,:) = crack_slab%pos(1,:)*(1.0_dp-v*strain)
             crack_slab%pos(3,:) = crack_slab%pos(3,:)*(1.0_dp-v2*strain)
             crack_slab%lattice(3,3) = crack_slab%lattice(3,3)*(1.0_dp-v2*strain)
             call atoms_set_lattice(crack_slab, crack_slab%lattice)
          elseif(params%crack_rescale_x) then 
             !  Rescale in x direction by v 
             if (.not. get_value(crack_slab%params,'OrigHeight',orig_height)) orig_height = 0.0_dp
             strain = crack_g_to_strain(params%crack_G, E, v, orig_height)
             crack_slab%pos(1,:) = crack_slab%pos(1,:)*(1.0_dp-v*strain)
          endif
       end if

    else if (trim(params%crack_loading) == 'ramp') then

       call print_title('Seed crack - Loading Ramp')

       call crack_apply_strain_ramp(crack_slab, params%crack_G, params%crack_ramp_end_G, r_crack_pos, &
            r_crack_pos+params%crack_strain_zone_width, &
            r_crack_pos+params%crack_strain_zone_width+params%crack_ramp_length)

    else if (trim(params%crack_loading) == 'kfield') then

       call print_title('Seed crack - Irwin K-field Loading')

       call add_property(crack_slab, 'k_disp', 0.0_dp, n_cols=3)
       call crack_fix_pointers(crack_slab, nn, changed_nn, load, move_mask, edge_mask, md_old_changed_nn, &
            old_nn, hybrid, hybrid_mark, u_disp, k_disp)

       if (.not. get_value(crack_slab%params,'OrigHeight',orig_height)) orig_height = 0.0_dp

       call print('Initial stress intesity factor K_0 = '//crack_g_to_k(params%crack_G, E, v)/1e6_dp//' MPa.sqrt(m)')

       call set_value(crack_slab%params, 'CrackPos', r_crack_pos + 0.85_dp*params%crack_strain_zone_width)
       call set_value(crack_slab%params, 'OrigCrackPos', r_crack_pos + 0.85_dp*params%crack_strain_zone_width)
       call set_value(crack_slab%params, 'G', params%crack_G)
       call crack_k_field(crack_slab, crack_g_to_k(params%crack_G, E, v), do_disp=.true.)

       call crack_fix_pointers(crack_slab, nn, changed_nn, load, move_mask, edge_mask, md_old_changed_nn, &
            old_nn, hybrid, hybrid_mark, u_disp, k_disp)

       do i=1,crack_slab%N
          crack_slab%pos(:,i) = crack_slab%pos(:,i) + k_disp(:,i)
       end do

    else if (trim(params%crack_loading) == 'interp_kfield_uniform') then

       ! interpolate linearly between K field (near tip) and uniform loading (near edge)

       call print_title('Seed crack - K-field and Uniform Loading')

       call add_property(crack_slab, 'uniform_disp', 0.0_dp, n_cols=3)
       call add_property(crack_slab, 'k_disp', 0.0_dp, n_cols=3)
       call crack_fix_pointers(crack_slab, nn, changed_nn, load, move_mask, edge_mask, md_old_changed_nn, &
            old_nn, hybrid, hybrid_mark, u_disp, k_disp)

       call print('Interpolation length '//params%crack_load_interp_length//' A')

       if (.not. get_value(crack_slab%params,'OrigHeight',orig_height)) orig_height = 0.0_dp

       call print('Initial energy release rate    G_0 = '//params%crack_G//' J/m^2')
       call print('Initial stress intesity factor K_0 = '//crack_g_to_k(params%crack_G, E, v)/1e6_dp//' MPa.sqrt(m)')

       call set_value(crack_slab%params, 'CrackPos', r_crack_pos + 0.85_dp*params%crack_strain_zone_width)
       call set_value(crack_slab%params, 'OrigCrackPos', r_crack_pos + 0.85_dp*params%crack_strain_zone_width)
       call crack_k_field(crack_slab, crack_g_to_k(params%crack_G, E, v), do_disp=.true.)

       G = params%crack_G
       call crack_uniform_load(crack_slab, l_crack_pos, r_crack_pos, &
            params%crack_strain_zone_width, G, apply_load=.false.)
       call set_value(crack_slab%params, 'G', G)

       call crack_fix_pointers(crack_slab, nn, changed_nn, load, move_mask, edge_mask, md_old_changed_nn, &
            old_nn, hybrid, hybrid_mark, u_disp, k_disp)

       do i=1,crack_slab%N
          r = sqrt((crack_slab%pos(1,i) - (r_crack_pos + 0.85*params%crack_strain_zone_width))**2.0_dp + &
               crack_slab%pos(2,i)**2.0_dp)
          if (r > params%crack_load_interp_length) then
             crack_slab%pos(:,i) = crack_slab%pos(:,i) + u_disp(:,i)
          else
             do k=1,3
                crack_slab%pos(k,i) = crack_slab%pos(k,i) +  &
                     linear_interpolate(0.0_dp, k_disp(k,i), params%crack_load_interp_length, u_disp(k,i), r)
             end do
          end if
       end do
    else
       call system_abort('Unknown loading type '//trim(params%crack_loading))
    end if

  end subroutine crack_make_seed


  ! displacement field subtracted must be identical to that added - save interp field




  !% Increase the load by adding the the load displacement field
  !% to the atomic positions. The routine recalculates the loading 
  !% G and stores it in the atom parameter dictionary.
  subroutine crack_apply_load_increment(at, G_increment)
    type(Atoms) :: at
    real(dp), optional :: G_increment

    integer :: i
    real(dp), pointer, dimension(:,:) :: load
    real(dp) :: orig_height, new_height, G, E, v
    real(dp) :: load_scale, my_G_increment
    real(dp) :: cur_height, cur_G, test_height, load_dheight, prefactor, target_G

    if (.not. get_value(at%params, 'OrigHeight', orig_height)) &
         call system_abort('crack_apply_load_increment: "OrigHeight" parameter missing from atoms structure')

    if (.not. get_value(at%params, 'YoungsModulus', E)) &
         call system_abort('crack_apply_load_increment: "YoungsModulus" missing')

    if (.not. get_value(at%params, 'PoissonRatio_yx', v)) &
         call system_abort('crack_apply_load_increment: "PoissonRatio_yx" missing')

    if (.not. assign_pointer(at, 'load', load)) &
         call system_abort('crack_apply_load_increment: load field is missing')

    my_G_increment = optional_default(0.0_dp, G_increment)
    if (my_G_increment .fne. 0.0_dp) then
       ! calc cur properties
       cur_height = maxval(at%pos(2,:))-minval(at%pos(2,:))
       cur_G = 0.1_dp*0.5_dp*E/(1.0_dp-v*v)*(cur_height - orig_height)**2.0_dp/orig_height
       target_G = cur_G + G_increment

       call print("crack_apply_load: doing G_increment, cur_G="//cur_G)

       ! calc properties after full load_step
       do i=1,at%N
          at%pos(:,i) = at%pos(:,i) + load(:,i)
       end do
       test_height = maxval(at%pos(2,:))-minval(at%pos(2,:))
       load_dheight = test_height-cur_height
       ! step back to orig positions
       do i=1,at%N
          at%pos(:,i) = at%pos(:,i) - load(:,i)
       end do

       ! cur_G = A (cur_h-orig_h)^2/orig_h
       ! load_dheight = test_height - cur_height
       ! target_G = A (cur_h + x load_dh - orig_h)^2/orig_h
       !
       ! A = cur_G*orig_h / (cur_h-orig_h)^2
       ! x = (sqrt(target_G*orig_h/A)-cur_h+orig_h)/load_dh
       prefactor = cur_G*orig_height / (cur_height-orig_height)**2
       load_scale = (sqrt(target_G*orig_height/prefactor) - cur_height + orig_height) / load_dheight
    else
       load_scale = 1.0_dp
    endif

    do i=1,at%N
       at%pos(:,i) = at%pos(:,i) + load_scale*load(:,i)
    end do
    call calc_dists(at)

    new_height = maxval(at%pos(2,:))-minval(at%pos(2,:))
    call print('crack_apply_load: new height = '//new_height)
    G = 0.1_dp*0.5_dp*E/(1.0_dp-v*v)*(new_height - orig_height)**2.0_dp/orig_height
    call print('crack_apply_load: new loading G = '//G)

    call set_value(at%params, 'G', G)

  end subroutine crack_apply_load_increment



  !% Return true if the point 'd' is within an ellipse centred at the origin
  !% with the $x$, $y$, and $z$ radii specifeid in the vector 'ellipse'.
  function in_ellipse(d, ellipse)
    real(dp), dimension(3), intent(in) :: d, ellipse
    logical :: In_Ellipse

    In_Ellipse = (d(1)/ellipse(1))*(d(1)/ellipse(1))+ &
         (d(2)/ellipse(2))*(d(2)/ellipse(2))+ &
         (d(3)/ellipse(3))*(d(3)/ellipse(3)) < 1.0_dp
  end function in_ellipse


  !% Select atoms in ellipse around atom 'c'.
  !% Principal radii of ellipse in $x$,$y$ and $z$ directions
  !% are given by the components of the vector 'ellipse'
  !% 'ellipse_bias' shifts ellipse, positive values forward
  !% On exit 'list' contains indexes of selected atoms.
  subroutine select_ellipse(at, ellipse, ellipse_bias, list, c)
    type(Atoms), intent(in) :: at
    real(dp), intent(in) :: ellipse(3), ellipse_bias(3)
    type(Table), intent(inout) :: list
    integer, intent(in) :: c

    integer :: i
    real(dp) :: p(3), cutoff

    call table_allocate(list, 4, 0, 0, 0)
    call append(list, (/c,0,0,0/))

    if (at%use_uniform_cutoff) then
       cutoff = at%cutoff*bond_length(at%Z(c),at%Z(c))
    else
       cutoff = at%cutoff
    endif

    ! Grow in all directions by enough steps to select sphere of 2*principal radius.
    ! With nneigh_only set to false, each hop will increase radius by about the
    ! neighbour cutoff distance
    call bfs_grow(at, list, max(nint(2.5_dp*maxval(ellipse)/cutoff),1), &
         nneighb_only=.false., min_images_only = .true.)

    ! Remove things we've added that are outside ellipse
    i = 1
    do while (i <= list%N)
       p = diff(at,c,list%int(1,i),list%int(2:4,i))
       p = p - ellipse_bias ! bias ellipse (arbitrary direction)
       if (.not. In_Ellipse(p, ellipse)) then
          call delete(list, i)
       else
          i = i + 1
       end if
    end do

  end subroutine select_ellipse


  !% Returns true if atom 'i' is near to an open surface of slab.
  !% Open surfaces are planes at $x = \pm$'OrigWidth/2' and $y = \pm$'OrigHeight/2'
  !% Near to means within 'edge_gap' of the surface.
  function crack_is_edge_atom(slab, i, edge_gap)
    type(Atoms), intent(in) :: slab
    integer, intent(in) :: i
    real(dp), intent(in) :: edge_gap
    logical :: crack_is_edge_atom
    real(dp) :: width, height

    if (.not. get_value(slab%params, 'OrigWidth', width)) &
         call system_abort('crack_is_edge_atom: "OrigWidth" parameter missing')

    if (.not. get_value(slab%params, 'OrigHeight', height)) &
         call system_abort('crack_is_edge_atom: "OrigWidth" parameter missing')

    crack_is_edge_atom =  abs(slab%pos(2,i)) > height/2.0_dp - edge_gap .or. &
         abs(slab%pos(1,i)) > width/2.0_dp - edge_gap

  end function crack_is_edge_atom


  !% Update the connectivity of a crack slab. calc_connect is only called if 
  !% necessary (i.e. if the maximal atomic displacement is bigger than
  !% 'params%md_recalc_connect_factor*params%md_crust'
  !% The 'nn' and 'changed_nn' properties are updated each call, with
  !% the (cheaper) nearest neighbour calc_connect always being perforemd.
  subroutine crack_update_connect(at, params)
    type(Atoms), intent(inout) :: at
    type(CrackParams), intent(in) :: params

    type(Atoms), save :: nn_atoms
    real(dp) :: max_disp
    real(dp), allocatable, save, dimension(:,:) :: stored_pos
    integer, pointer, dimension(:) :: old_nn, nn, changed_nn, edge_mask
    integer :: i
    logical :: first_time

    call system_timer('connectivity update')

    if (.not. assign_pointer(at, 'nn', nn)) &
         call system_abort('crack_update_connect: nn property missing from atoms')

    if (.not. assign_pointer(at, 'old_nn', old_nn)) &
         call system_abort('crack_update_connect: old_nn property missing from atoms')

    if (.not. assign_pointer(at, 'changed_nn', changed_nn)) &
         call system_abort('crack_update_connect: changed_nn property missing from atoms')

    if (.not. assign_pointer(at, 'edge_mask', edge_mask)) &
         call system_abort('crack_update_connect: edge property missing from atoms')

    ! First time we need to allocate stored_pos
    first_time = .false.
    if (.not. allocated(stored_pos)) then
       allocate(stored_pos(3,at%N))
       stored_pos = 0.0_dp
       first_time = .true.
       nn_atoms = at  ! First time we copy entire atoms structure
    end if

    if (.not. first_time) then
       max_disp = maxval(norm2(stored_pos - at%pos, 1))
       call print('Maximum atomic displacement since last calc_connect is '//max_disp)
    end if

    if (first_time .or. (max_disp > params%md_recalc_connect_factor*params%md_crust)) then
       call print('Recalculating connectivity')
       call calc_connect(at)
       stored_pos = at%pos ! Save positions for next time
    end if

    call print('Recalculating nearest neighbour table')
    call atoms_set_cutoff_factor(nn_atoms, params%md_nneigh_tol)
    if (trim(params%simulation_task) == 'md' .and. associated(at%avgpos)) then
       nn_atoms%pos = at%avgpos
    else
       nn_atoms%pos = at%pos
    end if
    call calc_connect(nn_atoms)

    nn = 0
    do i = 1,nn_atoms%N
       nn(i) = atoms_n_neighbours(nn_atoms, i)
    end do

    if (all(old_nn == 0)) old_nn = nn ! Special case for first time

    ! Age all the NN changes by 1
    where (changed_nn /= 0) changed_nn = changed_nn + 1

    ! Update changed_nn flag (excludes edge atoms)
    where (edge_mask == 0 .and. nn /= old_nn) changed_nn = 1

    if (count(changed_nn == 1) /= 0) then
       call print('CONNECT '//count(changed_nn == 1)//' atoms changed their neighbour count.')
    end if
    old_nn = nn

    call system_timer('connectivity update')

  end subroutine crack_update_connect


  !% Update QM selection region for a crack configuration using the 'nn' and 'changed_nn'
  !% properties and the 'CrackPos' parameter from the atoms structure, as well as the
  !% selection parameters in 'params'. If 'update_embed' is true then the embed region is
  !% updated, otherwise we simply recompute the fit region from the embed region.
  !% The value of 'num_directionality' returned can be passed to adjustable_potential_init.
  subroutine crack_update_selection(at, params)
    type(Atoms), intent(inout) :: at
    type(CrackParams), intent(in) :: params

    integer :: p, i, j, k, surface, age, n
    type(Table) :: old_embed, selectlist(2), tmp_select, embedlist
    type(Table), dimension(2) :: new_embed
    integer, allocatable, dimension(:) :: sorted, sindex
    real(dp), dimension(2,3) :: selection_ellipse
    real(dp) :: ellipse_bias(3), crack_pos, real_pos(3), lattice_coord(3)

    integer, pointer, dimension(:) :: nn, changed_nn, hybrid

    call system_timer('selection update')

    if (.not. assign_pointer(at, 'nn', nn)) &
         call system_abort('crack_update_selection: nn property missing from atoms')

    if (.not. assign_pointer(at, 'changed_nn', changed_nn)) &
         call system_abort('crack_update_selection: changed_nn property missing from atoms')

    if (.not. assign_pointer(at, 'hybrid', hybrid)) &
         call system_abort('crack_update_selection: atoms structure is missing hybrid property')


    if (.not. get_value(at%params, 'CrackPos', crack_pos)) &
         call system_abort('crack_update_selection: CrackPos parameter missing from atoms')

    call print('Building QM selection zone...')

    call allocate(embedlist, 1,0,0,0)
    call allocate(old_embed, 1,0,0,0)
    call print('count(changed_nn /= 0) = '//count(changed_nn /= 0))
    call print('Got '//count(hybrid == HYBRID_ACTIVE_MARK)//' old embed atoms')
    if (count(hybrid == HYBRID_ACTIVE_MARK) /= 0) &
         call append(old_embed, find(hybrid == HYBRID_ACTIVE_MARK))

    selection_ellipse(1,:) = params%selection_ellipse
    selection_ellipse(2,:) = params%selection_ellipse
    do i=1,3
       selection_ellipse(2,i) = selection_ellipse(2,i) + params%selection_ellipse_buffer
    end do
    ellipse_bias = 0.0_dp
    ! bias ellipse forward by arbitrary fraction of a radius (0.5 => 1/4 back, 3/4 ahead)
    ellipse_bias(1) = params%selection_ellipse_bias*selection_ellipse(1,1)

    ! Do selection twice, once to get inner and once to get outer surface
    do surface=1,2

       call table_allocate(selectlist(surface), 5, 0, 0, 0)

       ! Mark ellipsoid around each real atom with changed_nn /= 0 with its age
       !  - If central (active) atom already marked, keep the newer mark
       !  - If embedded atoms already marked, also keep newer mark

       do i=1,at%N
          if (changed_nn(i) == 0) cycle
          
          if (abs(at%pos(1,i)-crack_pos) < params%selection_cutoff_plane .and. &
               abs(at%pos(2,i)) < params%selection_cutoff_plane) then

             p = Find_in_array(selectlist(surface)%int(1,1:selectlist(surface)%N), i)
             if (p == 0) then
                call append(selectlist(surface), (/i,0,0,0,changed_nn(i)/))
             else
                selectlist(surface)%int(5,p) = min(selectlist(surface)%int(5,p), changed_nn(i))
             end if

             call select_ellipse(at, selection_ellipse(surface,:), ellipse_bias, tmp_select, i)

             do j = 1, tmp_select%N
                p = Find_in_array(int_part(selectlist(surface),(/1,2,3,4/)), tmp_select%int(:,j))
                if (p == 0) then
                   ! Marking for first time
                   call append(selectlist(surface), (/tmp_select%int(:,j), changed_nn(i)/))
                else
                   ! Keep newer mark
                   selectlist(surface)%int(5,p) = min(selectlist(surface)%int(5,p), changed_nn(i))
                end if
             end do
          end if
       end do

       ! Sort by age of NN changes, most recent are smallest values
       allocate(sorted(selectlist(surface)%N))
       allocate(sindex(selectlist(surface)%N))
       sorted = selectlist(surface)%int(5,1:selectlist(surface)%N)

       call insertion_sort(sorted, sindex)

       i = 1
       do while (i <= selectlist(surface)%N .and. new_embed(surface)%N < params%selection_max_qm_atoms)
          age = sorted(i)
          write (line, '(a,i0)') '  Selecting changed_nn age ', age
          call print(line)

          do while(i <= selectlist(surface)%N)
             if (sorted(i) /= age) exit
             call append(new_embed(surface), selectlist(surface)%int(1:4,sindex(i)))
             i = i + 1
          end do

          write (line,'(a,i0,a,i0,a)') 'Surface ',surface,' Now embedding ', new_embed(surface)%N, ' atoms'
          call print(line)
       end do

       deallocate(sorted)
       deallocate(sindex)
    end do


    call wipe(embedlist)

    ! Keep old embed atoms unless they're now outside outer surface
    do i=1,old_embed%N
       if (is_in_array(new_embed(2)%int(1,1:new_embed(2)%N), old_embed%int(1,i))) &
            call append(embedlist, old_embed%int(:,i))
    end do

    ! Add atoms inside inner surface
    do i=1,new_embed(1)%N
       if (.not. is_in_array(embedlist%int(1,1:embedlist%N), new_embed(1)%int(1,i))) &
            call append(embedlist, new_embed(1)%int(1,i))
    end do

    call Print('Embedding '//embedlist%N//' atoms.')

    call finalise(old_embed)
    call finalise(new_embed(1))
    call finalise(new_embed(2))
    call finalise(selectlist(1))
    call finalise(selectlist(2))
    call finalise(tmp_select)


!!$    ! Grow QM region to form fit region
!!$    call print('Building fit zone...')
!!$
!!$    if (params%hack_fit_on_eqm_coordination_only) then
!!$       fitlist = embedlist
!!$
!!$       ! Only add atoms which aren't undercoordinated to fit list
!!$       do i=1, params%fit_hops
!!$          call wipe(tmplist)
!!$          call bfs_step(at, fitlist, tmplist, nneighb_only = .true., min_images_only = .true.)
!!$          do j=1,tmplist%N
!!$             if (nn(tmplist%int(1,j)) == params%md_eqm_coordination) call append(fitlist, tmplist%int(:,j)) 
!!$          end do
!!$       end do
!!$    else
!!$       fitlist = embedlist
!!$       call bfs_grow(at, fitlist, params%fit_hops, min_images_only = .true.)
!!$    end if
!!$    call print('Fitting on '//fitlist%N//' atoms')
!!$
!!$    ! How many atoms should we require good directionality on, i.e. good spring
!!$    ! spanning of 3D space?
!!$    if (params%selection_directionality) then
!!$       num_directionality = embedlist%N
!!$    else
!!$       num_directionality = 0
!!$    end if

    ! Update crack position: set to average of embed atom time-averaged x coordinate
    ! (or just normal positions if we're not doing MD)
    if (embedlist%N /= 0) then
       crack_pos = 0.0_dp
       do i=1,embedlist%N
          if (trim(params%simulation_task) == 'md' .and. associated(at%avgpos)) then

             lattice_coord = at%g .mult. at%avgpos(:,embedlist%int(1,i))
             do n=1,3
                if ((lattice_coord(n) < -0.5_dp) .or. (lattice_coord(n) >= 0.5_dp)) then
                   k = floor(lattice_coord(n)+0.5_dp)
                   lattice_coord(n) = lattice_coord(n) - k
                end if
             end do

             real_pos = at%lattice .mult. lattice_coord
             crack_pos = crack_pos + real_pos(1)
          else
             crack_pos = crack_pos + at%pos(1,embedlist%int(1,i))
          end if
       end do
       crack_pos = crack_pos/embedlist%N
    end if
    call Print('Crack position = '//crack_pos)
    call set_value(at%params, 'CrackPos', crack_pos)

    ! Copy embedlist to 'hybrid' property
    hybrid = 0
    hybrid(int_part(embedlist,1)) = 1

    call system_timer('selection update')

    call finalise(embedlist)

  end subroutine crack_update_selection


  !% Return $x$ coordinate of rightmost undercoordinated atom
  function crack_find_crack_pos(at, params) result(crack_pos)
    type(Atoms), intent(inout) :: at
    type(CrackParams) :: params
    real(dp) :: crack_pos

    integer :: i, crack_tip_atom
    integer, pointer, dimension(:) :: nn, edge_mask
    real(dp) :: orig_width

    if (.not. assign_pointer(at, 'nn', nn)) &
         call system_abort('crack_find_crack_pos: nn property missing from atoms')

    if (.not. assign_pointer(at, 'edge_mask', edge_mask)) &
         call system_abort('crack_find_crack_pos: edge property missing from atoms')

    if (.not. get_value(at%params, 'OrigWidth', orig_width)) &
         call system_abort('crack_find_crack_pos: "OrigWidth" parameter missing from atoms')

    crack_pos = -orig_width
    crack_tip_atom = 0
    do i=1,at%N
       if (nn(i) == params%md_eqm_coordination) cycle
       if (edge_mask(i) == 1) cycle
       if (at%pos(1,i) > crack_pos) then
          crack_pos = at%pos(1,i)
          crack_tip_atom = i
       end if
    end do

    call Print('Crack position = '//crack_pos//' near atom '//crack_tip_atom)
    call set_value(at%params, 'CrackPos', crack_pos)

  end function crack_find_crack_pos

!!$  subroutine crack_print_xyz_name(at, xyzfilename, params)
!!$    type(Atoms), intent(inout) :: at
!!$    character(len=*), intent(in) :: xyzfilename
!!$    type(CrackParams), intent(in) :: params
!!$
!!$    if (params%io_print_all_properties) then
!!$       call print_xyz(at, xyzfilename, all_properties=.true.)
!!$    else
!!$       call print_xyz(at, xyzfilename, properties=params%io_print_properties)
!!$    end if
!!$  end subroutine crack_print_xyz_name
!!$
!!$  subroutine crack_print_xyz_file(at, xyzfile, params)
!!$    type(Atoms), intent(inout) :: at
!!$    type(inoutput), intent(inout) :: xyzfile
!!$    type(CrackParams), intent(in) :: params
!!$
!!$    if (params%io_print_all_properties) then
!!$       call print_xyz(at, xyzfile, all_properties=.true.)
!!$    else
!!$       call print_xyz(at, xyzfile, properties=params%io_print_properties)
!!$    end if
!!$  end subroutine crack_print_xyz_file

  subroutine crack_print_cio(at, cio, params, mpi)
    type(Atoms), intent(inout) :: at
    type(CInoutput), intent(inout) :: cio
    type(CrackParams), intent(in) :: params
    type(MPI_context) :: mpi

    if (.not. mpi%active .or. (mpi%active .and.mpi%my_proc == 0)) then
       if (params%io_print_all_properties) then
          call write(cio, at)
       else
          call write(cio, at, properties=params%io_print_properties)
       end if
    end if
  end subroutine crack_print_cio

  subroutine crack_print_filename(at, filename, params, mpi)
    type(Atoms), intent(inout) :: at
    character(*), intent(in) :: filename
    type(CrackParams), intent(in) :: params
    type(MPI_context) :: mpi

    type(CInOutput) :: cio

    if (.not. mpi%active .or. (mpi%active .and.mpi%my_proc == 0)) then
       call initialise(cio, filename, action=OUTPUT)
       if (params%io_print_all_properties) then
          call write(cio, at)
       else
          call write(cio, at, properties=params%io_print_properties)
       end if
    end if
  end subroutine crack_print_filename


end module CrackTools_module
