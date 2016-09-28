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

! do various structural analysis to a trajectory, save an intermediate file
! to be postprocessed (mean, variance, correlation, etc)
! #/vol density on a radial mesh
! #/vol density on a grid
! RDF, possibly as a function of distance of center from a fixed point
! to be added: ADF

module structure_analysis_routines_module
use libatoms_module
implicit none
private

public :: density_sample_radial_mesh_Gaussians, rdfd_calc, xrd_calc, propdf_radial_calc, adfd_calc, &
   density_sample_rectilinear_mesh_Gaussians, density_bin_rectilinear_mesh, &
   geometry_calc, shift_silica_to_edges, density_axial_calc, num_hbond_calc, &
   water_orientation_calc

real(dp), parameter, dimension(0:116) :: ASF_c = (/ &
    0.0, &
    0.001305, &
    0.0064, &
    0.0377, &
    0.0385, &
    -0.1932, &
    0.2156, &
    -11.529, &
    0.2508, &
    0.2776, &
    0.3515, &
    0.676, &
    0.8584, &
    1.1151, &
    1.1407, &
    1.1149, &
    0.8669, &
    -9.5574, &
    1.4445, &
    1.4228, &
    1.3751, &
    1.3329, &
    1.2807, &
    1.2199, &
    1.1832, &
    1.0896, &
    1.0369, &
    1.0118, &
    1.0341, &
    1.191, &
    1.3041, &
    1.7189, &
    2.1313, &
    2.531, &
    2.8409, &
    2.9557, &
    2.825, &
    3.4873, &
    2.5064, &
    1.91213, &
    2.06929, &
    3.75591, &
    4.3875, &
    5.40428, &
    5.37814, &
    5.328, &
    5.26593, &
    5.179, &
    5.0694, &
    4.9391, &
    4.7821, &
    4.5909, &
    4.352, &
    4.0712, &
    3.7118, &
    3.3352, &
    2.7731, &
    2.14678, &
    1.86264, &
    2.0583, &
    1.98486, &
    2.02876, &
    2.20963, &
    2.5745, &
    2.4196, &
    3.58324, &
    4.29728, &
    4.56796, &
    5.92046, &
    1.63929, &
    7.56672, &
    7.97628, &
    8.58154, &
    9.24354, &
    9.8875, &
    10.472, &
    11.0005, &
    11.4722, &
    11.6883, &
    12.0658, &
    12.6089, &
    13.1746, &
    13.4118, &
    13.5782, &
    13.677, &
    13.7108, &
    13.6905, &
    13.7247, &
    13.6211, &
    13.5266, &
    13.4314, &
    13.4287, &
    13.3966, &
    13.3573, &
    13.3812, &
    13.3592, &
    13.2887, &
    13.2754, &
    13.2674, &
    1.0, &
    1.0, &
    1.0, &
    1.0, &
    1.0, &
    1.0, &
    1.0, &
    1.0, &
    1.0, &
    1.0, &
    1.0, &
    1.0, &
    1.0, &
    1.0, &
    1.0, &
    1.0, &
    1.0, &
    1.0 /)

real(dp), parameter, dimension(4,0:116) :: ASF_a = reshape( (/ &
    0.0, 0.0, 0.0, 0.0, &
    0.489918, 0.262003, 0.196767, 0.049879, &
    0.8734, 0.6309, 0.3112, 0.178, &
    1.1282, 0.7508, 0.6175, 0.4653, &
    1.5919, 1.1278, 0.5391, 0.7029, &
    2.0545, 1.3326, 1.0979, 0.7068, &
    2.31, 1.02, 1.5886, 0.865, &
    12.2126, 3.1322, 2.0125, 1.1663, &
    3.0485, 2.2868, 1.5463, 0.867, &
    3.5392, 2.6412, 1.517, 1.0243, &
    3.9553, 3.1125, 1.4546, 1.1251, &
    4.7626, 3.1736, 1.2674, 1.1128, &
    5.4204, 2.1735, 1.2269, 2.3073, &
    6.4202, 1.9002, 1.5936, 1.9646, &
    6.2915, 3.0353, 1.9891, 1.541, &
    6.4345, 4.1791, 1.78, 1.4908, &
    6.9053, 5.2034, 1.4379, 1.5863, &
    11.4604, 7.1964, 6.2556, 1.6455, &
    7.4845, 6.7723, 0.6539, 1.6442, &
    8.2186, 7.4398, 1.0519, 0.8659, &
    8.6266, 7.3873, 1.5899, 1.0211, &
    9.189, 7.3679, 1.6409, 1.468, &
    9.7595, 7.3558, 1.6991, 1.9021, &
    10.2971, 7.3511, 2.0703, 2.0571, &
    10.6406, 7.3537, 3.324, 1.4922, &
    11.2819, 7.3573, 3.0193, 2.2441, &
    11.7695, 7.3573, 3.5222, 2.3045, &
    12.2841, 7.3409, 4.0034, 2.3488, &
    12.8376, 7.292, 4.4438, 2.38, &
    13.338, 7.1676, 5.6158, 1.6735, &
    14.0743, 7.0318, 5.1652, 2.41, &
    15.2354, 6.7006, 4.3591, 2.9623, &
    16.0816, 6.3747, 3.7068, 3.683, &
    16.6723, 6.0701, 3.4313, 4.2779, &
    17.0006, 5.8196, 3.9731, 4.3543, &
    17.1789, 5.2358, 5.6377, 3.9851, &
    17.3555, 6.7286, 5.5493, 3.5375, &
    17.1784, 9.6435, 5.1399, 1.5292, &
    17.5663, 9.8184, 5.422, 2.6694, &
    17.776, 10.2946, 5.72629, 3.26588, &
    17.8765, 10.948, 5.41732, 3.65721, &
    17.6142, 12.0144, 4.04183, 3.53346, &
    3.7025, 17.2356, 12.8876, 3.7429, &
    19.1301, 11.0948, 4.64901, 2.71263, &
    19.2674, 12.9182, 4.86337, 1.56756, &
    19.2957, 14.3501, 4.73425, 1.28918, &
    19.3319, 15.5017, 5.29537, 0.605844, &
    19.2808, 16.6885, 4.8045, 1.0463, &
    19.2214, 17.6444, 4.461, 1.6029, &
    19.1624, 18.5596, 4.2948, 2.0396, &
    19.1889, 19.1005, 4.4585, 2.4663, &
    19.6418, 19.0455, 5.0371, 2.6827, &
    19.9644, 19.0138, 6.14487, 2.5239, &
    20.1472, 18.9949, 7.5138, 2.2735, &
    20.2933, 19.0298, 8.9767, 1.99, &
    20.3892, 19.1062, 10.662, 1.4953, &
    20.3361, 19.297, 10.888, 2.6959, &
    20.578, 19.599, 11.3727, 3.28719, &
    21.1671, 19.7695, 11.8513, 3.33049, &
    22.044, 19.6697, 12.3856, 2.82428, &
    22.6845, 19.6847, 12.774, 2.85137, &
    23.3405, 19.6095, 13.1235, 2.87516, &
    24.0042, 19.4258, 13.4396, 2.89604, &
    24.6274, 19.0886, 13.7603, 2.9227, &
    25.0709, 19.0798, 13.8518, 3.54545, &
    25.8976, 18.2185, 14.3167, 2.95354, &
    26.507, 17.6383, 14.5596, 2.96577, &
    26.9049, 17.294, 14.5583, 3.63837, &
    27.6563, 16.4285, 14.9779, 2.98233, &
    28.1819, 15.8851, 15.1542, 2.98706, &
    28.6641, 15.4345, 15.3087, 2.98963, &
    28.9476, 15.2208, 15.1, 3.71601, &
    29.144, 15.1726, 14.7586, 4.30013, &
    29.2024, 15.2293, 14.5135, 4.76492, &
    29.0818, 15.43, 14.4327, 5.11982, &
    28.7621, 15.7189, 14.5564, 5.44174, &
    28.1894, 16.155, 14.9305, 5.67589, &
    27.3049, 16.7296, 15.6115, 5.83377, &
    27.0059, 17.7639, 15.7131, 5.7837, &
    16.8819, 18.5913, 25.5582, 5.86, &
    20.6809, 19.0417, 21.6575, 5.9676, &
    27.5446, 19.1584, 15.538, 5.52593, &
    31.0617, 13.0637, 18.442, 5.9696, &
    33.3689, 12.951, 16.5877, 6.4692, &
    34.6726, 15.4733, 13.1138, 7.02588, &
    35.3163, 19.0211, 9.49887, 7.42518, &
    35.5631, 21.2816, 8.0037, 7.4433, &
    35.9299, 23.0547, 12.1439, 2.11253, &
    35.763, 22.9064, 12.4739, 3.21097, &
    35.6597, 23.1032, 12.5977, 4.08655, &
    35.5645, 23.4219, 12.7473, 4.80703, &
    35.8847, 23.2948, 14.1891, 4.17287, &
    36.0228, 23.4128, 14.9491, 4.188, &
    36.1874, 23.5964, 15.6402, 4.1855, &
    36.5254, 23.8083, 16.7707, 3.47947, &
    36.6706, 24.0992, 17.3415, 3.49331, &
    36.6488, 24.4096, 17.399, 4.21665, &
    36.7881, 24.7736, 17.8919, 4.23284, &
    36.9185, 25.1995, 18.3317, 4.24391, &
    0.0, 0.0, 0.0, 0.0, &
    0.0, 0.0, 0.0, 0.0, &
    0.0, 0.0, 0.0, 0.0, &
    0.0, 0.0, 0.0, 0.0, &
    0.0, 0.0, 0.0, 0.0, &
    0.0, 0.0, 0.0, 0.0, &
    0.0, 0.0, 0.0, 0.0, &
    0.0, 0.0, 0.0, 0.0, &
    0.0, 0.0, 0.0, 0.0, &
    0.0, 0.0, 0.0, 0.0, &
    0.0, 0.0, 0.0, 0.0, &
    0.0, 0.0, 0.0, 0.0, &
    0.0, 0.0, 0.0, 0.0, &
    0.0, 0.0, 0.0, 0.0, &
    0.0, 0.0, 0.0, 0.0, &
    0.0, 0.0, 0.0, 0.0, &
    0.0, 0.0, 0.0, 0.0, &
    0.0, 0.0, 0.0, 0.0 /), (/ 4, 117 /) )

real(dp), parameter, dimension(4,0:116) :: ASF_b = reshape( (/ &
    0.0, 0.0, 0.0, 0.0, &
    20.6593, 7.74039, 49.5519, 2.20159, &
    9.1037, 3.3568, 22.9276, 0.9821, &
    3.9546, 1.0524, 85.3905, 68.261, &
    43.6427, 1.8623, 103.483, 0.542, &
    23.2185, 1.021, 60.3498, 0.1403, &
    20.8439, 10.2075, 0.5687, 51.6512, &
    0.0057, 9.8933, 28.9975, 0.5826, &
    13.2771, 5.7011, 0.3239, 32.9089, &
    10.2825, 4.2944, 0.2615, 26.1476, &
    8.4042, 3.4262, 0.2306, 21.7184, &
    3.285, 8.8422, 0.3136, 129.424, &
    2.8275, 79.2611, 0.3808, 7.1937, &
    3.0387, 0.7426, 31.5472, 85.0886, &
    2.4386, 32.3337, 0.6785, 81.6937, &
    1.9067, 27.157, 0.526, 68.1645, &
    1.4679, 22.2151, 0.2536, 56.172, &
    0.0104, 1.1662, 18.5194, 47.7784, &
    0.9072, 14.8407, 43.8983, 33.3929, &
    12.7949, 0.7748, 213.187, 41.6841, &
    10.4421, 0.6599, 85.7484, 178.437, &
    9.0213, 0.5729, 136.108, 51.3531, &
    7.8508, 0.5, 35.6338, 116.105, &
    6.8657, 0.4385, 26.8938, 102.478, &
    6.1038, 0.392, 20.2626, 98.7399, &
    5.3409, 0.3432, 17.8674, 83.7543, &
    4.7611, 0.3072, 15.3535, 76.8805, &
    4.2791, 0.2784, 13.5359, 71.1692, &
    3.8785, 0.2565, 12.1763, 66.3421, &
    3.5828, 0.247, 11.3966, 64.8126, &
    3.2655, 0.2333, 10.3163, 58.7097, &
    3.0669, 0.2412, 10.7805, 61.4135, &
    2.8509, 0.2516, 11.4468, 54.7625, &
    2.6345, 0.2647, 12.9479, 47.7972, &
    2.4098, 0.2726, 15.2372, 43.8163, &
    2.1723, 16.5796, 0.2609, 41.4328, &
    1.9384, 16.5623, 0.2261, 39.3972, &
    1.7888, 17.3151, 0.2748, 164.934, &
    1.5564, 14.0988, 0.1664, 132.376, &
    1.4029, 12.8006, 0.125599, 104.354, &
    1.27618, 11.916, 0.117622, 87.6627, &
    1.18865, 11.766, 0.204785, 69.7957, &
    0.2772, 1.0958, 11.004, 61.6584, &
    0.864132, 8.14487, 21.5707, 86.8472, &
    0.80852, 8.43467, 24.7997, 94.2928, &
    0.751536, 8.21758, 25.8749, 98.6062, &
    0.698655, 7.98929, 25.2052, 76.8986, &
    0.6446, 7.4726, 24.6605, 99.8156, &
    0.5946, 6.9089, 24.7008, 87.4825, &
    0.5476, 6.3776, 25.8499, 92.8029, &
    5.8303, 0.5031, 26.8909, 83.9571, &
    5.3034, 0.4607, 27.9074, 75.2825, &
    4.81742, 0.420885, 28.5284, 70.8403, &
    4.347, 0.3814, 27.766, 66.8776, &
    3.9282, 0.344, 26.4659, 64.2658, &
    3.569, 0.3107, 24.3879, 213.904, &
    3.216, 0.2756, 20.2073, 167.202, &
    2.94817, 0.244475, 18.7726, 133.124, &
    2.81219, 0.226836, 17.6083, 127.113, &
    2.77393, 0.222087, 16.7669, 143.644, &
    2.66248, 0.210628, 15.885, 137.903, &
    2.5627, 0.202088, 15.1009, 132.721, &
    2.47274, 0.196451, 14.3996, 128.007, &
    2.3879, 0.1942, 13.7546, 123.174, &
    2.25341, 0.181951, 12.9331, 101.398, &
    2.24256, 0.196143, 12.6648, 115.362, &
    2.1802, 0.202172, 12.1899, 111.874, &
    2.07051, 0.19794, 11.4407, 92.6566, &
    2.07356, 0.223545, 11.3604, 105.703, &
    2.02859, 0.238849, 10.9975, 102.961, &
    1.9889, 0.257119, 10.6647, 100.417, &
    1.90182, 9.98519, 0.261033, 84.3298, &
    1.83262, 9.5999, 0.275116, 72.029, &
    1.77333, 9.37046, 0.295977, 63.3644, &
    1.72029, 9.2259, 0.321703, 57.056, &
    1.67191, 9.09227, 0.3505, 52.0861, &
    1.62903, 8.97948, 0.382661, 48.1647, &
    1.59279, 8.86553, 0.417916, 45.0011, &
    1.51293, 8.81174, 0.424593, 38.6103, &
    0.4611, 8.6216, 1.4826, 36.3956, &
    0.545, 8.4484, 1.5729, 38.3246, &
    0.65515, 8.70751, 1.96347, 45.8149, &
    0.6902, 2.3576, 8.618, 47.2579, &
    0.704, 2.9238, 8.7937, 48.0093, &
    0.700999, 3.55078, 9.55642, 47.0045, &
    0.68587, 3.97458, 11.3824, 45.4715, &
    0.6631, 4.0691, 14.0422, 44.2473, &
    0.646453, 4.17619, 23.1052, 150.645, &
    0.616341, 3.87135, 19.9887, 142.325, &
    0.589092, 3.65155, 18.599, 117.02, &
    0.563359, 3.46204, 17.8309, 99.1722, &
    0.547751, 3.41519, 16.9235, 105.251, &
    0.5293, 3.3253, 16.0927, 100.613, &
    0.511929, 3.25396, 15.3622, 97.4908, &
    0.499384, 3.26371, 14.9455, 105.98, &
    0.483629, 3.20647, 14.3136, 102.273, &
    0.465154, 3.08997, 13.4346, 88.4834, &
    0.451018, 3.04619, 12.8946, 86.003, &
    0.437533, 3.00775, 12.4044, 83.7881, &
    0.0, 0.0, 0.0, 0.0, &
    0.0, 0.0, 0.0, 0.0, &
    0.0, 0.0, 0.0, 0.0, &
    0.0, 0.0, 0.0, 0.0, &
    0.0, 0.0, 0.0, 0.0, &
    0.0, 0.0, 0.0, 0.0, &
    0.0, 0.0, 0.0, 0.0, &
    0.0, 0.0, 0.0, 0.0, &
    0.0, 0.0, 0.0, 0.0, &
    0.0, 0.0, 0.0, 0.0, &
    0.0, 0.0, 0.0, 0.0, &
    0.0, 0.0, 0.0, 0.0, &
    0.0, 0.0, 0.0, 0.0, &
    0.0, 0.0, 0.0, 0.0, &
    0.0, 0.0, 0.0, 0.0, &
    0.0, 0.0, 0.0, 0.0, &
    0.0, 0.0, 0.0, 0.0, &
    0.0, 0.0, 0.0, 0.0 /), (/ 4, 117 /) )


contains

subroutine is_in_mask(mask_a, at, mask_str)
  type(Atoms), intent(in) :: at
  logical, intent(out) :: mask_a(at%N)
  character(len=*), intent(in) :: mask_str

  integer :: i_at, i_Z
  integer :: Zmask
  type(Table) :: atom_indices
  character(len=4) :: species(128)
  integer :: n_species
  character(len=len(mask_str)) :: fields(2)
  integer :: n_fields

  integer :: i_val
  logical :: l_val

  integer, pointer :: p_i(:)
  logical, pointer :: p_l(:)

  if (len_trim(mask_str) == 0) then
    mask_a = .true.
    return
  endif

  mask_a = .false.
  if (mask_str(1:1)=='@') then ! list of indices
    call parse_atom_mask(mask_str,atom_indices)
    do i_at=1, atom_indices%N
      mask_a(atom_indices%int(1,i_at)) = .true.
    end do
  else if (scan(mask_str,'=')/=0) then ! arbitrary property
    call split_string(mask_str,'=','""', fields,n_fields)
    if (assign_pointer(at, trim(fields(1)), p_i)) then 
       ! integer match
       read (unit=fields(2), fmt=*) i_val
       mask_a = (p_i == i_val)
    else if (assign_pointer(at, trim(fields(1)), p_l)) then 
       ! integer match
       read (unit=fields(2), fmt=*) l_val
       mask_a = (p_l .eqv. l_val)
    else
       call system_abort("mask is arbitrary property match, but apparently not integer or logical, so unsupported")
    endif
  else ! species (via atomic number)
    call split_string(mask_str, ' ,', '""', species, n_species)
    do i_Z=1, n_species
      Zmask = Atomic_Number(species(i_Z))
      do i_at=1, at%N
        if (at%Z(i_at) == Zmask) mask_a(i_at) = .true.
      end do
    end do
  end if
end subroutine is_in_mask

subroutine density_sample_radial_mesh_Gaussians(histogram, at, center_pos, center_i, rad_bin_width, n_rad_bins, gaussian_sigma, mask_str, radial_pos, accumulate, quantity)
  real(dp), intent(inout) :: histogram(:)
  type(Atoms), intent(inout) :: at
  real(dp), intent(in), optional :: center_pos(3)
  integer, intent(in), optional :: center_i
  real(dp), intent(in) :: rad_bin_width
  integer, intent(in) :: n_rad_bins
  real(dp), intent(in) :: gaussian_sigma
  character(len=*), optional, intent(in) :: mask_str
  real(dp), intent(out), optional :: radial_pos(:)
  logical, optional, intent(in) :: accumulate
  character(len=*), optional, intent(in) :: quantity

  logical :: my_accumulate

  real(dp) :: use_center_pos(3), d, r0, s_sq, s_cu, h_val
  real(dp), parameter :: SQROOT_PI = sqrt(PI), PI_THREE_HALVES = PI**(3.0_dp/2.0_dp)
  logical, allocatable :: mask_a(:)
  integer :: at_i, rad_sample_i
  real(dp) :: rad_sample_r, exp_arg, ep, em
  logical :: quantity_1, quantity_KE
  real(dp) :: sample_weight

  if (present(center_pos) .and. present(center_i)) then
    call system_abort("density_sample_radial_mesh_Gaussians received both center_pos and center_i")
  else if (.not. present(center_pos) .and. .not. present(center_i)) then
    call system_abort("density_sample_radial_mesh_Gaussians received neither center_pos nor center_i")
  endif

  quantity_1 = .false.
  quantity_KE = .false.
  if (present(quantity)) then
    select case (quantity)
      case("1")
	quantity_1 = .true.
      case("KE")
	quantity_KE = .true.
      case default
	call system_abort("density_sample_radial_mesh_Gaussians called with unknown quantity='"//trim(quantity)//"'")
    end select
  else
    quantity_1 = .true.
  endif
  if (quantity_1) sample_weight = 1.0_dp

  my_accumulate = optional_default(.false., accumulate)
  if (.not. my_accumulate) histogram = 0.0_dp

  allocate(mask_a(at%N))
  call is_in_mask(mask_a, at, mask_str)

  if (present(radial_pos)) then
    do rad_sample_i=1, n_rad_bins
      rad_sample_r = (rad_sample_i-1)*rad_bin_width
      radial_pos(rad_sample_i) = rad_sample_r
    end do
  endif

  ! ratio of 20/4=5 is bad
  ! ratio of 20/3=6.66 is bad
  ! ratio of 20/2.5=8 is borderline (2e-4)
  ! ratio of 20/2.22=9 is fine  (error 2e-5)
  ! ratio of 20/2=10 is definitely fine
  if (min(norm(at%lattice(:,1)), norm(at%lattice(:,2)), norm(at%lattice(:,3))) < 9.0_dp*gaussian_sigma) &
    call print("WARNING: at%lattice may be too small for sigma, errors (noticeably too low a density) may result", PRINT_ALWAYS)

  if (present(center_i)) then
    use_center_pos = at%pos(:,center_i)
  else
    use_center_pos = center_pos
  end if

  s_sq = gaussian_sigma**2
  s_cu = s_sq*gaussian_sigma
  do at_i=1, at%N
    if (.not. mask_a(at_i)) cycle
    if (present(center_i)) then
      if (at_i == center_i) cycle
    endif
    if (quantity_KE) then
      if (associated(at%velo)) then
	sample_weight = kinetic_energy(ElementMass(at%Z(at_i)), at%velo(1:3,at_i))
      else
	sample_weight = 0.0_dp
      endif
    endif
    d = distance_min_image(at,use_center_pos,at%pos(:,at_i))
    ! skip if atom is outside range
    if (d > (n_rad_bins-1)*rad_bin_width+6.0_dp*gaussian_sigma) cycle
    do rad_sample_i=1, n_rad_bins
      rad_sample_r = (rad_sample_i-1)*rad_bin_width
      r0 = rad_sample_r
      ep = 0.0_dp
      exp_arg = -(r0+d)**2/s_sq
      if (exp_arg > -20.0_dp) ep = exp(exp_arg)
      em = 0.0_dp
      exp_arg = -(r0-d)**2/s_sq
      if (exp_arg > -20.0_dp) em = exp(exp_arg)
      ! should really fix d->0 limit
      if (d .feq. 0.0_dp) then
	h_val = 0.0_dp
	exp_arg = -r0**2/s_sq
	if (exp_arg > -20.0_dp) &
	  h_val = exp(exp_arg) / (PI_THREE_HALVES * s_cu)
      else if (r0 .feq. 0.0_dp) then
	h_val = 0.0_dp
	exp_arg = -d**2/s_sq
	if (exp_arg > -20.0_dp) &
	  h_val = exp(exp_arg) / (PI_THREE_HALVES * s_cu)
      else
	h_val = (r0/(SQROOT_PI * gaussian_sigma * d) * (em - ep)) /  (4.0_dp * PI * r0**2)
      endif
      histogram(rad_sample_i) = histogram(rad_sample_i) +  sample_weight*h_val
    end do ! rad_sample_i
  end do ! at_i

end subroutine density_sample_radial_mesh_Gaussians

! from LAMMPS "compute xrd" documentation
subroutine xrd_calc(xrd, at, lambda, range_2theta, n_2theta, pos_2theta)
  real(dp), intent(inout) :: xrd(:)
  type(Atoms), intent(inout) :: at
  real(dp), intent(in) :: lambda
  real(dp), intent(in) :: range_2theta(2)
  integer, intent(in) :: n_2theta
  real(dp) , intent(inout), optional :: pos_2theta(:)

  integer :: i_at, i_2theta
  real(dp) :: min_norm_k, max_norm_k
  integer :: k1_min, k1_max, k2_min, k2_max, k3_min, k3_max, n_k, k1, k2, k3
  real(dp) :: k(3), theta
  integer :: theta_i
  integer, allocatable :: theta_i_of_k(:), Z_of_type(:)
  real(dp), allocatable :: ks(:,:), sin_sq_theta_over_lambda_sq(:), ASF(:), k_dot_r(:,:)
  complex(dp), allocatable :: exp_2_pi_i_k_dot_r(:,:), F_k(:)
  integer :: i_type, i_Z, n_types, i_k
  complex(dp) :: two_pi_i

  ! print *, "xrd_calc for lambda ", lambda, "range_2theta ", range_2theta, "n_2theta ", n_2theta
  if (present(pos_2theta)) then
    do i_2theta=1, n_2theta
       pos_2theta(i_2theta) = range_2theta(1) + &
         (i_2theta-1)*(range_2theta(2)-range_2theta(1)) / real(n_2theta-1.0_dp,dp)
    end do
    ! print *, "pos_2theta ", pos_2theta
  end if

  min_norm_k = 2.0*sin(range_2theta(1)/2.0)/lambda
  max_norm_k = 2.0*sin(range_2theta(2)/2.0)/lambda
  ! print *, "norm k range ", min_norm_k, max_norm_k
  k1_min = min_norm_k/norm(at%g(:,1))-1
  k1_max = max_norm_k/norm(at%g(:,1))+1
  k2_min = min_norm_k/norm(at%g(:,2))-1
  k2_max = max_norm_k/norm(at%g(:,2))+1
  k3_min = min_norm_k/norm(at%g(:,3))-1
  k3_max = max_norm_k/norm(at%g(:,3))+1
  k1_min = max(k1_min, 0)
  k2_min = max(k2_min, 0)
  k3_min = max(k3_min, 0)
  ! print *, "k range ", k1_min, k1_max, k2_min, k2_max, k3_min, k3_max
  n_k = 0
  do k1 = -k1_max, k1_max
  do k2 = -k2_max, k2_max
  do k3 = -k3_max, k3_max
    k = k1*at%g(:,1) + k2*at%g(:,2) + k3*at%g(:,3)
    if (norm(k) >= min_norm_k .and. norm(k) <= max_norm_k) n_k = n_k + 1
  end do
  end do
  end do

  ! print *, "n_k ", n_k

  allocate(ks(n_k,3))
  allocate(theta_i_of_k(n_k), sin_sq_theta_over_lambda_sq(n_k))

  n_k = 0
  do k1 = -k1_max, k1_max
  do k2 = -k2_max, k2_max
  do k3 = -k3_max, k3_max
    k = k1*at%g(:,1) + k2*at%g(:,2) + k3*at%g(:,3)
    if (norm(k) >= min_norm_k .and. norm(k) <= max_norm_k) then
        n_k = n_k + 1
        ks(n_k,:) = k
        theta = asin(norm(k)*lambda/2.0)
        sin_sq_theta_over_lambda_sq(n_k) = (norm(k)/2.0)**2
        theta_i = int(n_2theta*(2.0*theta - range_2theta(1))/(range_2theta(2)-range_2theta(1)))+1
        theta_i_of_k(n_k) = theta_i
        ! print *, "k ", n_k, k, theta_i
    endif
  end do
  end do
  end do

  allocate(exp_2_pi_i_k_dot_r(n_k, at%n))
  allocate(k_dot_r(n_k, at%n))
  call matrix_product_sub(k_dot_r, ks, at%pos)
  two_pi_i = 2.0*PI*cmplx(0.0_dp, 1.0_dp)
  !$omp parallel do
  do i_at=1, at%n
      exp_2_pi_i_k_dot_r(:,i_at) = exp(two_pi_i*k_dot_r(:,i_at))
  end do
  !$omp end parallel do
  deallocate(k_dot_r)

  allocate(Z_of_type(at%n))
  allocate(ASF(n_k))
  n_types = 0
  do i_Z=1, maxval(at%Z)
    if (any(at%Z == i_Z)) then
        n_types = n_types + 1
        Z_of_type(n_types) = i_Z
    endif
  end do
  do i_type=1, n_types
      ASF = atomic_structure_factor(Z_of_type(i_type), sin_sq_theta_over_lambda_sq(:))
      do i_at = 1, at%n
        if (at%Z(i_at) == Z_of_type(i_type)) exp_2_pi_i_k_dot_r(:,i_at) = exp_2_pi_i_k_dot_r(:,i_at) * ASF
      end do
  end do
  deallocate(ASF, Z_of_type)

  allocate(F_k(n_k))
  ! F_k = sum(exp_2_pi_i_k_dot_r, dim=2)
  !$omp parallel do
  do i_k=1, n_k
    F_k(i_k) = sum(exp_2_pi_i_k_dot_r(i_k,:))
  end do
  !$omp end parallel do
  !$omp parallel do
  do i_2theta = 1, n_2theta
    xrd(i_2theta) = sum(abs(F_k)**2, mask=(theta_i_of_k == i_2theta)) / at%n
  end do
  !$omp end parallel do

  deallocate(ks, exp_2_pi_i_k_dot_r, theta_i_of_k, F_k, sin_sq_theta_over_lambda_sq)

end subroutine xrd_calc

function atomic_structure_factor(Z, sin_sq_theta_over_lambda_sq) result (ASF)
    integer :: Z
    real(dp) :: sin_sq_theta_over_lambda_sq(:)
    real(dp) :: ASF(size(sin_sq_theta_over_lambda_sq))

    integer :: i

    ! print *, "Z ", Z, ASF_c(Z), ASF_a(:,Z), ASF_b(:,Z)
    ASF = ASF_c(Z)
    do i=1, 4
        ASF(:) = ASF(:) + ASF_a(i,Z)*exp(-ASF_b(i,Z)*sin_sq_theta_over_lambda_sq(:))
    end do
end function

subroutine rdfd_calc(rdfd, at, zone_center, zone_atom_center, bin_width, n_bins, zone_width, n_zones, gaussian_smoothing, gaussian_sigma, &
                     center_mask_str, neighbour_mask_str, bin_pos, zone_pos)
  real(dp), intent(inout) :: rdfd(:,:)
  type(Atoms), intent(inout) :: at
  real(dp), intent(in) :: zone_center(3), bin_width, zone_width
  integer :: zone_atom_center
  integer, intent(in) :: n_bins, n_zones
  logical, intent(in) :: gaussian_smoothing
  real(dp), intent(in) :: gaussian_sigma
  character(len=*), intent(in) :: center_mask_str, neighbour_mask_str
  real(dp), intent(inout), optional :: bin_pos(:), zone_pos(:)

  logical, allocatable :: center_mask_a(:), neighbour_mask_a(:)
  integer :: i_at, j_at, i_bin, i_zone
  integer, allocatable :: n_in_zone(:)
  real(dp) :: r, bin_inner_rad, bin_outer_rad, my_zone_center(3)

  allocate(center_mask_a(at%N))
  allocate(neighbour_mask_a(at%N))
  call is_in_mask(center_mask_a, at, center_mask_str)
  call is_in_mask(neighbour_mask_a, at, neighbour_mask_str)

  allocate(n_in_zone(n_zones))
  n_in_zone = 0

  if (present(zone_pos)) then
    if (zone_width > 0.0_dp) then
      do i_zone=1, n_zones
        zone_pos(i_zone) = (real(i_zone,dp)-0.5_dp)*zone_width
      end do
    else
      zone_pos(1) = -1.0_dp
    endif
  endif
  if (present(bin_pos)) then
    do i_bin=1, n_bins
      if (gaussian_smoothing) then
        bin_pos(i_bin) = (i_bin-1)*bin_width
      else
        bin_pos(i_bin) = (real(i_bin,dp)-0.5_dp)*bin_width
      endif
    end do
  endif

!  if (gaussian_smoothing) then
!    call set_cutoff(at, n_bins*bin_width+5.0_dp*gaussian_sigma)
!    call calc_connect(at)
!  endif

  rdfd = 0.0_dp
  do i_at=1, at%N ! loop over center atoms
    if (.not. center_mask_a(i_at)) cycle

    !calc which zone the atom is in
    if (zone_width > 0.0_dp) then
      if (zone_atom_center > 0) then
	 if (zone_atom_center <= at%N) then
	    my_zone_center = at%pos(:,zone_atom_center)
	 else
	    call system_abort("rdfd_calc got zone_atom_center="//zone_atom_center//" out of range (1.."//at%N//")")
	 endif
      else
	 my_zone_center = zone_center
      endif
      r = distance_min_image(at, my_zone_center, at%pos(:,i_at))
      i_zone = int(r/zone_width)+1
      if (i_zone > n_zones) cycle
    else
      i_zone = 1
    endif

    !count the number of atoms in that zone
    n_in_zone(i_zone) = n_in_zone(i_zone) + 1

    !calc rdfd in each bin for this zone
    if (gaussian_smoothing) then
      call density_sample_radial_mesh_Gaussians(rdfd(:,i_zone), at, center_i=i_at, rad_bin_width=bin_width, n_rad_bins=n_bins, &
        gaussian_sigma=gaussian_sigma, mask_str=neighbour_mask_str, accumulate = .true.)
    else
      !loop over atoms and advance the bins
      do j_at=1, at%N
        if (j_at == i_at) cycle
        if (.not. neighbour_mask_a(j_at)) cycle
        r = distance_min_image(at, i_at, j_at)
        i_bin = int(r/bin_width)+1
        if (i_bin <= n_bins) rdfd(i_bin,i_zone) = rdfd(i_bin,i_zone) + 1.0_dp
      end do ! j_at
    endif ! gaussian_smoothing
  end do ! i_at

  if (.not. gaussian_smoothing) then
    !calculate local density by dividing bins with their volumes
    do i_bin=1, n_bins
      bin_inner_rad = real(i_bin-1,dp)*bin_width
      bin_outer_rad = real(i_bin,dp)*bin_width
      rdfd(i_bin,:) = rdfd(i_bin,:)/(4.0_dp/3.0_dp*PI*bin_outer_rad**3 - 4.0_dp/3.0_dp*PI*bin_inner_rad**3)
    end do
  end if

  ! normalise zones by the number of atoms in that zone
  do i_zone=1, n_zones
    if (n_in_zone(i_zone) > 0) rdfd(:,i_zone) = rdfd(:,i_zone)/real(n_in_zone(i_zone),dp)
  end do
  ! normalise with the global density
  if (count(neighbour_mask_a) > 0) then
    rdfd = rdfd / (count(neighbour_mask_a)/cell_volume(at))
  endif

end subroutine rdfd_calc

subroutine propdf_radial_calc(histograms, at, bin_width, n_bins, &
  zone_center, zone_width, n_zones, gaussian_sigma, mask_str, property, bin_pos, zone_pos)
  real(dp), intent(inout) :: histograms(:,:)
  type(Atoms), intent(in) :: at
  real(dp), intent(in) :: bin_width, zone_width
  integer, intent(in) :: n_bins, n_zones
  real(dp), intent(in) :: gaussian_sigma
  real(dp), intent(in) :: zone_center(3)
  character(len=*), intent(in) :: mask_str, property
  real(dp), intent(inout), optional :: bin_pos(:), zone_pos(:)

  logical, allocatable :: mask_a(:)
  integer :: i_at, i_bin_ctr, i_bin, i_zone
  integer, allocatable :: n_in_zone(:)
  real(dp) :: n, quant, r, bin_ctr
  logical :: has_mass
  logical :: doing_KE
  real(dp), pointer :: prop_a(:), prop_3a(:,:)
  logical :: prop_is_scalar

  allocate(mask_a(at%N))
  call is_in_mask(mask_a, at, mask_str)

  doing_KE = (trim(property) == 'KE') 

  if (len_trim(property) == 0) &
    call system_abort("propdf_radial_calc has no property specified.  Are you sure you passed a propdf_radial_property value?")

  if (doing_KE) then ! KE is a special case
    has_mass = has_property(at, 'mass')
    if (.not. has_property(at, 'velo')) &
      call system_abort("propdf_radial_calc has no 'velo' property in Atoms structure")
  else ! assign pointers for property
    if (.not. has_property(at, property)) &
       call system_abort("propdf_radial_calc has no '"//trim(property)//"' property in Atoms structure")

    prop_is_scalar = .true.
    if (.not. assign_pointer(at, trim(property), prop_a)) then
       prop_is_scalar = .false.
       if (.not. assign_pointer(at, trim(property), prop_3a)) then
	  call system_abort("propdf_radial_calc failed to assign 1- or 3-array real(dp) pointer for property '"//trim(property))
       endif
    endif
  endif

  allocate(n_in_zone(n_zones))
  n_in_zone = 0

  if (present(zone_pos)) then
    if (zone_width > 0.0_dp) then
      do i_zone=1, n_zones
        zone_pos(i_zone) = (real(i_zone,dp)-0.5_dp)*zone_width
      end do
    else
      zone_pos(1) = -1.0_dp
    endif
  endif
  if (present(bin_pos)) then
    do i_bin=1, n_bins
      bin_pos(i_bin) = (i_bin-1)*bin_width
    end do
  endif

  histograms = 0.0_dp
  do i_at=1, at%N ! loop over atoms
    if (.not. mask_a(i_at)) cycle

    !calc which zone the atom is in
    if (zone_width > 0.0_dp) then
      r = distance_min_image(at, zone_center, at%pos(:,i_at))
      i_zone = int(r/zone_width)+1
      if (i_zone > n_zones) cycle
    else
      i_zone = 1
    endif

    !count the number of atoms in that zone
    n_in_zone(i_zone) = n_in_zone(i_zone) + 1

    if (doing_KE) then
      !calc KE for this atom
      if (has_mass) then
	quant = kinetic_energy(at%mass(i_at), at%velo(:,i_at))
      else
	quant = kinetic_energy(ElementMass(at%Z(i_at)), at%velo(:,i_at))
      endif
    else
      if (prop_is_scalar) then
	quant = prop_a(i_at)
      else
	quant = norm(prop_3a(:,i_at))
      endif
    endif
    i_bin_ctr = quant/bin_width
    do i_bin=max(1,i_bin_ctr-floor(6*gaussian_sigma/bin_width)), min(n_bins,i_bin_ctr+floor(6*gaussian_sigma/bin_width))
      bin_ctr = (i_bin-1)*bin_width
      histograms(i_bin,i_zone) = histograms(i_bin,i_zone) + exp(-0.5*(quant-bin_ctr)**2/gaussian_sigma**2)
    end do
  end do ! i_at

  ! normalise zones by the number of atoms in that zone
  n = gaussian_sigma/sqrt(2.0_dp*PI)
  do i_zone=1, n_zones
    if (n_in_zone(i_zone) > 0) histograms(:,i_zone) = n*histograms(:,i_zone)/real(n_in_zone(i_zone),dp)
  end do

end subroutine propdf_radial_calc

subroutine adfd_calc(adfd, at, zone_center, n_angle_bins, dist_bin_width, n_dist_bins, zone_width, n_zones, &
		     center_mask_str, neighbour_1_mask_str, neighbour_1_max_dist, neighbour_2_mask_str, & 
		     dist_bin_rc2, angle_bin_pos, dist_bin_pos, zone_pos)
  real(dp), intent(inout) :: adfd(:,:,:)
  type(Atoms), intent(inout) :: at
  real(dp), intent(in) :: zone_center(3), dist_bin_width, zone_width
  integer, intent(in) :: n_angle_bins, n_dist_bins, n_zones
  character(len=*), intent(in) :: center_mask_str, neighbour_1_mask_str, neighbour_2_mask_str
  real(dp), intent(in) :: neighbour_1_max_dist
  logical, intent(in) :: dist_bin_rc2
  real(dp), intent(inout), optional :: angle_bin_pos(:), dist_bin_pos(:), zone_pos(:)

  logical, allocatable :: center_mask_a(:), neighbour_1_mask_a(:), neighbour_2_mask_a(:)
  integer :: i_at, j_at, k_at, i_dist_bin, i_angle_bin, i_zone, ji, ki
  integer, allocatable :: n_in_zone(:)
  real(dp) :: dist_bin_inner_rad, dist_bin_outer_rad, angle_bin_width, angle_bin_min, angle_bin_max
  real(dp) :: r, r_ij, r_ik, r_jk, jik_angle
  integer :: sj(3), sk(3)
  real(dp) :: c_ij(3), c_ik(3)

  allocate(center_mask_a(at%N), neighbour_1_mask_a(at%N), neighbour_2_mask_a(at%N))
  call is_in_mask(center_mask_a, at, center_mask_str)
  call is_in_mask(neighbour_1_mask_a, at, neighbour_1_mask_str)
  call is_in_mask(neighbour_2_mask_a, at, neighbour_2_mask_str)

  allocate(n_in_zone(n_zones))
  n_in_zone = 0

  if (present(zone_pos)) then
    if (zone_width > 0.0_dp) then
      do i_zone= 1, n_zones
	zone_pos(i_zone) = (real(i_zone,dp)-0.5_dp)*zone_width
      end do
    else
      zone_pos(1) = -1.0_dp
    end if
  end if
  if (present(dist_bin_pos)) then
    do i_dist_bin=1, n_dist_bins
      dist_bin_pos(i_dist_bin) = (real(i_dist_bin,dp)-0.5_dp)*dist_bin_width
    end do
  end if
  angle_bin_width = PI/real(n_angle_bins,dp)
  if (present(angle_bin_pos)) then
    do i_angle_bin=1, n_angle_bins
      angle_bin_pos(i_angle_bin) = (real(i_angle_bin,dp)-0.5_dp)*angle_bin_width
    end do
  end if

  adfd = 0.0_dp
  call set_cutoff(at, max(neighbour_1_max_dist,2*n_dist_bins*dist_bin_width))
  call calc_connect(at)
  do i_at=1, at%N ! center atom
    if (.not. center_mask_a(i_at)) cycle
    if (zone_width > 0.0_dp) then
      r = distance_min_image(at, zone_center, at%pos(:,i_at))
      i_zone = int(r/zone_width)+1
      if (i_zone > n_zones) cycle
    else
      i_zone = 1
    endif
    n_in_zone(i_zone) = n_in_zone(i_zone) + 1
    do ji=1, n_neighbours(at, i_at)
      j_at = neighbour(at, i_at, ji, shift=sj, cosines=c_ij, distance=r_ij)
      if (r_ij == 0.0_dp) cycle
      if (.not. neighbour_1_mask_a(j_at)) cycle
      if (neighbour_1_max_dist > 0.0_dp) then
	if (r_ij > neighbour_1_max_dist) cycle
      else
	if (r_ij > bond_length(at%Z(i_at),at%Z(j_at))*at%nneightol) cycle
      endif
      do ki=1, n_neighbours(at, i_at)
	k_at = neighbour(at, i_at, ki, shift=sk, cosines=c_ik, distance=r_ik)
	if (r_ik == 0.0_dp .or. (k_at == j_at .and. all (sj == sk))) cycle
	if (.not. neighbour_2_mask_a(k_at)) cycle
	r_jk = distance_min_image(at, j_at, k_at)
	if (dist_bin_rc2) then
	   i_dist_bin = int(r_ik/dist_bin_width)+1
	else
	   i_dist_bin = int(r_jk/dist_bin_width)+1
	endif
	if (i_dist_bin > n_dist_bins) cycle
! call print("doing triplet ijk " // i_at // " " // j_at // " "// k_at, PRINT_ALWAYS)
! call print("  Zijk " // at%Z(i_at) // " " // at%Z(j_at) // " " // at%Z(k_at), PRINT_ALWAYS)
! call print("  pi " // at%pos(:,i_at), PRINT_ALWAYS)
! call print("  pj " // at%pos(:,j_at), PRINT_ALWAYS)
! call print("  pk " // at%pos(:,k_at), PRINT_ALWAYS)
! call print("  r_ij " // diff_min_image(at,i_at,j_at) // "     " // r_ij, PRINT_ALWAYS)
! call print("  r_jk " // diff_min_image(at,j_at,k_at) // "     " // distance_min_image(at, j_at, k_at), PRINT_ALWAYS)
! call print("  r_ik " // diff_min_image(at,i_at,k_at) // "     " // distance_min_image(at, i_at, k_at), PRINT_ALWAYS)
        ! jik_angle = angle(diff_min_image(at,i_at,j_at), diff_min_image(at,i_at,k_at))
        jik_angle = angle(c_ij, c_ik)
! call print("  r_ij " // r_ij // " r_jk " // r_jk // " jik_angle " // (jik_angle*180.0/PI), PRINT_ALWAYS)
	i_angle_bin = int(jik_angle/angle_bin_width)+1
	if (i_angle_bin > n_angle_bins) i_angle_bin = n_angle_bins
	adfd(i_angle_bin,i_dist_bin,i_zone) = adfd(i_angle_bin,i_dist_bin,i_zone) + 1.0_dp ! /(2.0_dp*PI*sin(jik_angle))
      end do ! k_at
    end do ! j_at
  end do ! i_at

  do i_angle_bin=1, n_angle_bins
    angle_bin_min = real(i_angle_bin-1,dp)*angle_bin_width
    angle_bin_max = real(i_angle_bin,dp)*angle_bin_width
    do i_dist_bin=1, n_dist_bins
      dist_bin_inner_rad = real(i_dist_bin-1,dp)*dist_bin_width
      dist_bin_outer_rad = real(i_dist_bin,dp)*dist_bin_width
      adfd(i_angle_bin,i_dist_bin,:) = adfd(i_angle_bin,i_dist_bin,:) / (2.0_dp*PI * ((dist_bin_outer_rad**3 - dist_bin_inner_rad**3)/3.0_dp) * (cos(angle_bin_min)-cos(angle_bin_max)))
      ! adfd(i_angle_bin,i_dist_bin,:) = adfd(i_angle_bin,i_dist_bin,:) / (4.0_dp/3.0_dp*PI * (dist_bin_outer_rad**3 - dist_bin_inner_rad**3))
    end do
  end do

  ! normalise zones by the number of atoms in that zone
  do i_zone=1, n_zones
    if (n_in_zone(i_zone) > 0) adfd(:,:,i_zone) = adfd(:,:,i_zone)/real(n_in_zone(i_zone),dp)
  end do
  ! normalise with the global density
  if (count(neighbour_2_mask_a) > 0) then
    adfd = adfd / (count(neighbour_2_mask_a)/cell_volume(at))
  endif

end subroutine adfd_calc

subroutine density_sample_rectilinear_mesh_Gaussians(histogram, at, min_p, sample_dist, n_bins, gaussian_sigma, mask_str, grid_pos, accumulate)
  real(dp), intent(inout) :: histogram(:,:,:)
  type(Atoms), intent(in) :: at
  real(dp), intent(in) :: min_p(3), sample_dist(3)
  integer, intent(in) :: n_bins(3)
  real(dp), intent(in) :: gaussian_sigma
  character(len=*), optional, intent(in) :: mask_str
  real(dp), intent(out), optional :: grid_pos(:,:,:,:)
  logical, optional, intent(in) :: accumulate

  logical :: my_accumulate
  integer :: at_i, i1, i2, i3
  real(dp) :: p(3), w, dist
  logical, allocatable :: mask_a(:)

  my_accumulate = optional_default(.false., accumulate)
  if (.not. my_accumulate) histogram = 0.0_dp

  allocate(mask_a(at%N))
  call is_in_mask(mask_a, at, mask_str)

! ratio of 20/4=5 is bad
! ratio of 20/3=6.66 is bad
! ratio of 20/2.5=8 is borderline (2e-4)
! ratio of 20/2.22=9 is fine  (error 2e-5)
! ratio of 20/2=10 is fine
  if ( (norm(at%lattice(:,1)) < 9.0_dp*gaussian_sigma) .or. &
       (norm(at%lattice(:,2)) < 9.0_dp*gaussian_sigma) .or. &
       (norm(at%lattice(:,3)) < 9.0_dp*gaussian_sigma) ) &
    call print("WARNING: at%lattice may be too small for sigma, errors (noticeably too low a density) may result", PRINT_ALWAYS)

  w = 1.0_dp / (gaussian_sigma*sqrt(2.0_dp*PI))**3

  if (present(grid_pos)) then
    do i1=1, n_bins(1)
      p(1) = min_p(1) + (i1-1)*sample_dist(1)
      do i2=1, n_bins(2)
        p(2) = min_p(2) + (i2-1)*sample_dist(2)
        do i3=1, n_bins(3)
          p(3) = min_p(3) + (i3-1)*sample_dist(3)
          grid_pos(:,i1,i2,i3) = p
        end do
      end do
    end do
  endif

  do at_i=1, at%N
    if (.not. mask_a(at_i)) cycle
    do i1=1, n_bins(1)
      p(1) = min_p(1) + (i1-1)*sample_dist(1)
      do i2=1, n_bins(2)
        p(2) = min_p(2) + (i2-1)*sample_dist(2)
        do i3=1, n_bins(3)
          p(3) = min_p(3) + (i3-1)*sample_dist(3)
          dist = distance_min_image(at,p,at%pos(:,at_i))
          histogram(i1,i2,i3) = histogram(i1,i2,i3) + exp(-0.5_dp*(dist/(gaussian_sigma))**2)*w
        end do ! i3
      end do ! i2
    end do ! i1
  end do ! at_i

  deallocate(mask_a)
end subroutine density_sample_rectilinear_mesh_Gaussians

subroutine density_bin_rectilinear_mesh(histogram, at, min_p, bin_width, n_bins, mask_str, grid_pos, accumulate)
  real(dp), intent(inout) :: histogram(:,:,:)
  type(Atoms), intent(in) :: at
  real(dp), intent(in) :: min_p(3), bin_width(3)
  integer, intent(in) :: n_bins(3)
  character(len=*), optional, intent(in) :: mask_str
  real(dp), intent(out), optional :: grid_pos(:,:,:,:)
  logical, optional, intent(in) :: accumulate

  logical :: my_accumulate
  integer :: i, i1, i2, i3, bin(3)
  logical, allocatable :: mask_a(:)
  real(dp) :: p(3)

  my_accumulate = optional_default(.false., accumulate)
  if (.not. my_accumulate) histogram = 0.0_dp

  if (present(grid_pos)) then
    do i1=1, n_bins(1)
      p(1) = min_p(1) + (real(i1,dp)-0.5_dp)*bin_width(1)
      do i2=1, n_bins(2)
        p(2) = min_p(2) + (real(i2,dp)-0.5_dp)*bin_width(2)
        do i3=1, n_bins(3)
          p(3) = min_p(3) + (real(i3,dp)-0.5_dp)*bin_width(3)
          grid_pos(:,i1,i2,i3) = p
        end do
      end do
    end do
  endif

  allocate(mask_a(at%N))
  call is_in_mask(mask_a, at, mask_str)

  do i=1, at%N
    if (.not. mask_a(i)) cycle
    bin = floor((at%pos(:,i)-min_p)/bin_width)+1
    if (all(bin >= 1) .and. all (bin <= n_bins)) histogram(bin(1),bin(2),bin(3)) = histogram(bin(1),bin(2),bin(3)) + 1.0_dp
  end do

  deallocate(mask_a)

end subroutine density_bin_rectilinear_mesh

!Calculates distances, angles, dihedrals of the given atoms
subroutine geometry_calc(histogram, at, geometry_params, central_atom, geometry_pos, geometry_label)

  real(dp), intent(inout) :: histogram(:)
  type(Atoms), intent(inout) :: at
  type(Table), intent(in) :: geometry_params
  integer, intent(in) :: central_atom
  real(dp), intent(out), optional :: geometry_pos(:)
  character(STRING_LENGTH), intent(out), optional :: geometry_label(:)

  integer :: i, j, geom_type
  integer :: atom1, atom2, atom3, atom4
  real(dp) :: shift(3), bond12(3),bond23(3),bond34(3)
  real(dp) :: a, b, c, dist, min_dist
  logical :: found

  !center around central_atom if requested
  if (central_atom.gt.at%N) call system_abort('central atom is greater than atom number '//at%N)
  if (central_atom.gt.0) then !center around central_atom
     shift = at%pos(1:3,central_atom)
     do j=1,at%N
        at%pos(1:3,j) = at%pos(1:3,j) - shift(1:3)
     enddo
     call map_into_cell(at) !only in this case, otherwise it has been mapped
  endif

  !loop over the parameters to calculate
  do i=1, geometry_params%N
     geom_type=geometry_params%int(1,i)
     atom1 = geometry_params%int(2,i)
     atom2 = geometry_params%int(3,i)
     atom3 = geometry_params%int(4,i)
     atom4 = geometry_params%int(5,i)
     select case (geom_type)
       case(1) !y coord atom1
         if (atom1<1.or.atom1>at%N) call system_abort('atom1 must be >0 and < '//at%N)
         histogram(i) = at%pos(2,atom1)
       case(2) !distance atom1-atom2
         if (atom1<1.or.atom1>at%N) call system_abort('atom1 must be >0 and < '//at%N)
         if (atom2<1.or.atom2>at%N) call system_abort('atom2 must be >0 and < '//at%N)
         !histogram(i) = norm(at%pos(1:3,atom1)-at%pos(1:3,atom2))
         histogram(i) = distance_min_image(at,atom1,atom2)
       case(3) !angle atom1-atom2-atom3
         if (atom1<1.or.atom1>at%N) call system_abort('atom1 must be >0 and < '//at%N)
         if (atom2<1.or.atom2>at%N) call system_abort('atom2 must be >0 and < '//at%N)
         if (atom3<1.or.atom2>at%N) call system_abort('atom3 must be >0 and < '//at%N)
         !histogram(i) = angle(at%pos(1:3,atom1)-at%pos(1:3,atom2), &
         !                     at%pos(1:3,atom3)-at%pos(1:3,atom2))
         histogram(i) = angle(diff_min_image(at,atom2,atom1), &
                              diff_min_image(at,atom2,atom3))
       case(4) !dihedral atom1-(bond12)->atom2-(bond23)->atom3-(bond34)->atom4
         if (atom1<1.or.atom1>at%N) call system_abort('atom1 must be >0 and < '//at%N)
         if (atom2<1.or.atom2>at%N) call system_abort('atom2 must be >0 and < '//at%N)
         if (atom3<1.or.atom2>at%N) call system_abort('atom3 must be >0 and < '//at%N)
         if (atom4<1.or.atom2>at%N) call system_abort('atom4 must be >0 and < '//at%N)
         !bond12(1:3) = at%pos(1:3,atom2)-at%pos(1:3,atom1)
         bond12(1:3) = diff_min_image(at,atom1,atom2)
         !bond23(1:3) = at%pos(1:3,atom3)-at%pos(1:3,atom2)
         bond23(1:3) = diff_min_image(at,atom2,atom3)
         !bond34(1:3) = at%pos(1:3,atom4)-at%pos(1:3,atom3)
         bond34(1:3) = diff_min_image(at,atom3,atom4)
         histogram(i) = atan2(norm(bond23(1:3)) * bond12(1:3).dot.(bond23(1:3).cross.bond34(1:3)), &
                              (bond12(1:3).cross.bond23(1:3)) .dot. (bond23(1:3).cross.bond34(1:3)))
       case(5) !distance of atom3-Z from atom1-atom2 bond
         min_dist=HUGE(1._dp)
         found=.false.
         do j=1, at%N
            if (j==atom1 .or. j==atom2) cycle
            if (at%Z(j)/=atom3) cycle
            if (j<=6) cycle !CH3Cl2- specific!!!
            !only close ones
            a = distance_min_image(at,atom1,atom2)
            b = distance_min_image(at,atom1,j)
            c = distance_min_image(at,atom2,j)
!            if (b>4.0_dp .or. c>4.0_dp) cycle
            if (b>a .or. c>a) cycle
            dist = 0.25_dp*sqrt((a+b+c)*(b+c-a)*(c+a-b)*(a+b-c)) / (0.5_dp*a)
            if (dist<min_dist) then
               found=.true.
               min_dist=dist
            endif
            if (found) then
               histogram(i) = min_dist
            else
               histogram(i) = 0._dp
            endif
         enddo
       case default
         call system_abort("geometry_calc: unknown geometry type "//geom_type)
     end select
     if (present(geometry_pos)) geometry_pos(i) = real(i,dp)
     if (present(geometry_label)) geometry_label(i) = geom_type//'=='//atom1//'--'//atom2//'--'//atom3//'--'//atom4
  enddo

end subroutine geometry_calc


!silica-water inerface: put half-half of the silica slab to the 2 edges of the cell along the axis  normal to the surface
subroutine shift_silica_to_edges(at, axis, silica_center_i, mask_str)
  type(Atoms), intent(inout) :: at
  integer, intent(in) :: axis
  integer, intent(in) :: silica_center_i
  character(len=*), optional, intent(in) :: mask_str

  logical, allocatable :: mask_silica(:)
  integer :: i
  integer :: counter
  integer, allocatable :: Si_atoms(:)
  real(dp) :: com(3), shift(3)
  real(dp), pointer :: mass_p(:)

  if (.not.silica_center_i>0) return

  !find the com of the silica (using the Si atoms) and move that to the edges in this direction
  shift(1:3) = at%pos(1:3,1)
  do i=1,at%N
     at%pos(1:3,i) = at%pos(1:3,i) - shift(1:3)
  enddo
  call map_into_cell(at)
  allocate(mask_silica(at%N))
  call is_in_mask(mask_silica, at, "Si")
  allocate(Si_atoms(count(mask_silica)))
  counter = 0
  do i=1,at%N
     if (mask_silica(i)) then
        counter = counter + 1
        Si_atoms(counter) = i
     endif
  enddo
  !allocate(at%mass(at%N))
  !at%mass = ElementMass(at%Z)
  call add_property(at,"mass",0._dp)
  if (.not.(assign_pointer(at, "mass", mass_p))) call system_abort('??')
  mass_p = ElementMass(at%Z)
  com = centre_of_mass(at,index_list=Si_atoms(1:size(Si_atoms)),origin=1)
  !shift axis to the edge (-0.5*edge_length)
  at%pos(axis,1:at%N) = at%pos(axis,1:at%N) - 0.5_dp * at%lattice(axis,axis) - com(axis)
  call map_into_cell(at) !everyone becomes -0.5b<y<0.5b
  !NO !shift everyone to positive coordinate along axis
  !at%pos(axis,1:at%N) = at%pos(axis,1:at%N) + 0.5_dp * at%lattice(axis,axis)
  deallocate(Si_atoms)
  deallocate(mask_silica)

end subroutine shift_silica_to_edges

subroutine density_axial_calc(histogram, at, axis, silica_center_i,n_bins, gaussian_smoothing, gaussian_sigma, mask_str, axial_pos, accumulate)
  real(dp), intent(inout) :: histogram(:)
  type(Atoms), intent(inout) :: at
  integer, intent(in) :: axis
  integer, intent(in) :: silica_center_i
  integer, intent(in) :: n_bins
  logical, intent(in) :: gaussian_smoothing
  real(dp), intent(in) :: gaussian_sigma
  character(len=*), optional, intent(in) :: mask_str
  real(dp), intent(out), optional :: axial_pos(:)
  logical, optional, intent(in) :: accumulate

  logical :: my_accumulate
  real(dp) :: ax_sample_r, dist, r, exp_arg
  logical, allocatable :: mask_a(:)
  integer at_i, ax_sample_i
  real(dp) :: bin_width

  my_accumulate = optional_default(.false., accumulate)
  if (.not. my_accumulate) histogram = 0.0_dp

  !atommask
  allocate(mask_a(at%N))
  call is_in_mask(mask_a, at, mask_str)

  !bin labels
  bin_width = at%lattice(axis,axis) / n_bins
  if (present(axial_pos)) then
    do ax_sample_i=1, n_bins
      ax_sample_r = (real(ax_sample_i,dp)-0.5_dp)*bin_width !the middle of the bin
      axial_pos(ax_sample_i) = ax_sample_r
    end do
  endif

  if (silica_center_i>0) then ! silica
     call shift_silica_to_edges(at, axis, silica_center_i, mask_str)
  endif

  !now bins are from -axis/2 to axis/2
  !shift everyone to positive coordinate along axis
  at%pos(axis,1:at%N) = at%pos(axis,1:at%N) + 0.5_dp * at%lattice(axis,axis)


  !simply check distances and bin them


  ! ratio of 20/4=5 is bad
  ! ratio of 20/3=6.66 is bad
  ! ratio of 20/2.5=8 is borderline (2e-4)
  ! ratio of 20/2.22=9 is fine  (error 2e-5)
  ! ratio of 20/2=10 is fine
  if ( gaussian_smoothing .and. (at%lattice(axis,axis) < 9.0_dp*gaussian_sigma) ) &
    call print("WARNING: at%lattice may be too small for sigma, errors (noticeably too low a density) may result", PRINT_ALWAYS)

  if (gaussian_smoothing) then

    do at_i=1, at%N
      if (silica_center_i>0 .and. at_i<=silica_center_i) cycle !ignore silica atoms
      if (.not. mask_a(at_i)) cycle
      r = at%pos(axis,at_i)
      do ax_sample_i=1, n_bins

        ax_sample_r = (real(ax_sample_i,dp)-0.5_dp)*bin_width
        dist = abs(r - ax_sample_r)
!Include all the atoms, slow but minimises error
!	  if (dist > 4.0_dp*gaussian_sigma) cycle
          exp_arg = -0.5_dp*(dist/(gaussian_sigma))**2
          if (exp_arg > -20.0_dp) then ! good to about 1e-8
            histogram(ax_sample_i) = histogram(ax_sample_i) + exp(exp_arg)/(gaussian_sigma*sqrt(2.0_dp*PI)) !Gaussian in 1 dimension
          endif

      end do ! ax_sample_i
    end do ! at_i

  else !no gaussian_smoothing

    do at_i=1, at%N
      if (silica_center_i>0 .and. at_i<=silica_center_i) cycle !ignore silica atoms
      if (.not. mask_a(at_i)) cycle
      r = at%pos(axis,at_i)

        histogram(int(r/bin_width)+1) = histogram(int(r/bin_width)+1) + 1

    end do ! at_i

  endif

  deallocate(mask_a)
end subroutine density_axial_calc

!
! Calculates the number of H-bonds along an axis
!  -- 1st center around atom 1
!  -- then calculate COM of Si atoms
!  -- then shift the centre of mass along the y axis to y=0
!  -- calculate H bonds for
!        water - water
!        water - silica
!        silica - water
!        silica - silica
!     interactions
!  -- the definition of a H-bond (O1-H1 - - O2):
!        d(O1,O2) < 3.5 A
!        d(O2,H1) < 2.45 A
!        angle(H1,O1,O2) < 30 degrees
!     source: P. Jedlovszky, J.P. Brodholdt, F. Bruni, M.A. Ricci and R. Vallauri, J. Chem. Phys. 108, 8525 (1998)
!
subroutine num_hbond_calc(histogram, at, axis, silica_center_i,n_bins, gaussian_smoothing, gaussian_sigma, mask_str, num_hbond_pos, num_hbond_type_code, num_hbond_type_label,accumulate)
  real(dp),                          intent(inout) :: histogram(:,:)
  type(Atoms),                       intent(inout) :: at
  integer,                           intent(in)    :: axis
  integer,                           intent(in)    :: silica_center_i
  integer,                           intent(in)    :: n_bins
  logical,                           intent(in)    :: gaussian_smoothing
  real(dp),                          intent(in)    :: gaussian_sigma
  character(len=*),        optional, intent(in)    :: mask_str
  real(dp),                optional, intent(out)   :: num_hbond_pos(:)
  integer,                 optional, intent(out)   :: num_hbond_type_code(:)
  character(STRING_LENGTH), optional, intent(out)   :: num_hbond_type_label(:)
  logical,                 optional, intent(in)    :: accumulate

  logical :: my_accumulate
  real(dp) :: num_hbond_sample_r, dist, r, exp_arg
  logical, allocatable :: mask_a(:)
  integer :: num_hbond_sample_i
  real(dp) :: bin_width
  real(dp), parameter                   :: dist_O2_H1 = 2.45_dp
  real(dp), parameter                   :: dist_O1_O2 = 3.5_dp
  real(dp), parameter                   :: angle_H1_O1_O2 = 30._dp
real(dp) :: min_distance, distance, HOO_angle
integer :: H1, O1, i, j, k, O2, num_atoms, hbond_type

  my_accumulate = optional_default(.false., accumulate)
  if (.not. my_accumulate) histogram = 0.0_dp

  !atommask
  allocate(mask_a(at%N))
  call is_in_mask(mask_a, at, mask_str)

  !bin labels
  bin_width = at%lattice(axis,axis) / n_bins
  if (present(num_hbond_pos)) then
    do num_hbond_sample_i=1, n_bins
      num_hbond_sample_r = (real(num_hbond_sample_i,dp)-0.5_dp)*bin_width !the middle of the bin
      num_hbond_pos(num_hbond_sample_i) = num_hbond_sample_r
    end do
  endif
  if (present(num_hbond_type_label)) then
    num_hbond_type_label(1)="water-water"
    num_hbond_type_label(2)="water-silica"
    num_hbond_type_label(3)="silica-water"
    num_hbond_type_label(4)="silica-silica"
  endif
  if (present(num_hbond_type_code)) then
    num_hbond_type_code(1)=11
    num_hbond_type_code(2)=10
    num_hbond_type_code(3)=01
    num_hbond_type_code(4)=00
  endif

  if (silica_center_i>0) then ! silica
     call shift_silica_to_edges(at, axis, silica_center_i, mask_str)
  endif

  !!calc_connect now, before shifting positions to positive, because it would remap the positions!!
  !call calc_connect including the H-bonds
  call set_cutoff(at,dist_O2_H1)
  call calc_connect(at)

  !now bins are from -axis/2 to axis/2
  !shift everyone to positive coordinate along axis
  at%pos(axis,1:at%N) = at%pos(axis,1:at%N) + 0.5_dp * at%lattice(axis,axis)

  !simply check hbonds and bin them


  ! ratio of 20/4=5 is bad
  ! ratio of 20/3=6.66 is bad
  ! ratio of 20/2.5=8 is borderline (2e-4)
  ! ratio of 20/2.22=9 is fine  (error 2e-5)
  ! ratio of 20/2=10 is fine
  if ( gaussian_smoothing .and. (at%lattice(axis,axis) < 9.0_dp*gaussian_sigma) ) &
    call print("WARNING: at%lattice may be too small for sigma, errors (noticeably too low a density) may result", PRINT_ALWAYS)

!  call set_cutoff(at,dist_O2_H1)
!  call calc_connect(at)

  num_atoms = 0
  do H1=1, at%N
     if(at%Z(H1)/=1) cycle !find H: H1
     !Count the atoms
     call print('Found H'//H1//ElementName(at%Z(H1)),PRINT_ANAL)
     num_atoms = num_atoms + 1

     !find closest O: O1
     min_distance = huge(1._dp)
     O1 = 0
     k = 0
     do i = 1, n_neighbours(at,H1)
        j = neighbour(at,H1,i,distance)
        if (distance<min_distance) then
           min_distance = distance
           k = i  !the closest neighbour is the k-th one
           O1 = j
        endif
     enddo
     if (O1==0) call system_abort('H has no neighbours.')
     !if (at%Z(O1).ne.8) call system_abort('H'//H1//' has not O closest neighbour '//ElementName(at%Z(O1))//O1//'.')
     if (.not.mask_a(O1)) call system_abort('H'//H1//' has not O closest neighbour '//ElementName(at%Z(O1))//O1//'.')

     !loop over all other Os: O2
     do i = 1, n_neighbours(at,H1)
        if (i.eq.k) cycle
        O2 = neighbour(at,H1,i)
        !if (at%Z(O2).ne.8) cycle !only keep O
        if (.not. mask_a(O2)) cycle
        !check O1-O2 distance for definition
        if (distance_min_image(at,O1,O2).gt.dist_O1_O2) cycle
        !check H1-O1-O2 angle for definition
        HOO_angle = angle(diff_min_image(at,O1,H1), &
                          diff_min_image(at,O1,O2)) *180._dp/PI
        call print('HOO_ANGLE '//ElementName(at%Z(H1))//H1//' '//ElementName(at%Z(O1))//O1//' '//ElementName(at%Z(O2))//O2//' '//HOO_angle,PRINT_ANAL)
        if (HOO_angle.gt.angle_H1_O1_O2) cycle

        !We've found a H-bond.

        !Find out the type (what-to-what)
        if (O1>silica_center_i .and. O2>silica_center_i) then  ! water - water
           call print('Found water-water H-bond.',PRINT_ANAL)
           hbond_type = 1
        elseif (O1>silica_center_i .and. O2<=silica_center_i) then  ! water - silica
           call print('Found water-silica H-bond.',PRINT_ANAL)
           hbond_type = 2
        elseif (O1<=silica_center_i .and. O2>silica_center_i) then  ! silica - water
           call print('Found silica-water H-bond.',PRINT_ANAL)
           hbond_type = 3
        elseif (O1<=silica_center_i .and. O2<=silica_center_i) then  ! silica - silica
           call print('Found silica-silica H-bond.',PRINT_ANAL)
           hbond_type = 4
        endif

        !Build histogram
        r = at%pos(axis,H1) !the position of H1

        if (gaussian_smoothing) then !smear the position along axis
           do num_hbond_sample_i=1, n_bins
             num_hbond_sample_r = (real(num_hbond_sample_i,dp)-0.5_dp)*bin_width
             dist = abs(r - num_hbond_sample_r)
!!!!!!Include all the atoms, slow but minimises error
!!!!!!	  if (dist > 4.0_dp*gaussian_sigma) cycle
               exp_arg = -0.5_dp*(dist/(gaussian_sigma))**2
               if (exp_arg > -20.0_dp) then ! good to about 1e-8
                 histogram(num_hbond_sample_i,hbond_type) = histogram(num_hbond_sample_i,hbond_type) + exp(exp_arg)/(gaussian_sigma*sqrt(2.0_dp*PI)) !Gaussian in 1 dimension
               endif
           end do ! num_hbond_sample_i
        else !no gaussian_smoothing
           histogram(int(r/bin_width)+1,hbond_type) = histogram(int(r/bin_width)+1,hbond_type) + 1
        endif
     end do ! i, atom_neighbours

  end do ! H1

  deallocate(mask_a)

end subroutine num_hbond_calc

!
! Calculates the orientation of water molecules along the y axis
!  -- 1st center around atom 1
!  -- then calculate COM of Si atoms
!  -- then shift the centre of mass along the y axis to y=0
!  -- calculate the orientation of the {dipole moment} / {angle half line} of the water: if skip_atoms is set to the last atom of the silica
!
subroutine water_orientation_calc(histogram, at, axis, silica_center_i,n_pos_bins, n_angle_bins, gaussian_smoothing, pos_gaussian_sigma, pos_bin, angle_bin, angle_bin_w, use_dipole_rather_than_angle_bisector, accumulate)
!subroutine water_orientation_calc(histogram, at, axis, silica_center_i,n_pos_bins, n_angle_bins, gaussian_smoothing, pos_gaussian_sigma, angle_gaussian_sigma, pos_bin, angle_bin, angle_bin_w, use_dipole_rather_than_angle_bisector, accumulate)
  real(dp),                          intent(inout) :: histogram(:,:)
  type(Atoms),                       intent(inout) :: at
  integer,                           intent(in)    :: axis
  integer,                           intent(in)    :: silica_center_i
  integer,                           intent(in)    :: n_pos_bins, n_angle_bins
  logical,                           intent(in)    :: gaussian_smoothing
  real(dp),                          intent(in)    :: pos_gaussian_sigma !, angle_gaussian_sigma
  real(dp),                optional, intent(out)   :: pos_bin(:), angle_bin(:)
  real(dp),                optional, intent(inout) :: angle_bin_w(:)
  logical,                 optional, intent(in)    :: use_dipole_rather_than_angle_bisector
  logical,                 optional, intent(in)    :: accumulate

  logical :: my_accumulate
  real(dp) :: sample_r, sample_angle, r
  ! logical, allocatable :: mask_a(:)
  integer :: sample_i
  real(dp) :: pos_bin_width, angle_bin_width
  real(dp) :: sum_w
integer :: n, num_atoms
integer :: O, H1, H2
  real(dp) :: surface_normal(3)
  logical :: use_dipole
real(dp) :: vector_OH1(3), vector_OH2(3)
real(dp) :: bisector_vector(3)
real(dp) :: dipole(3)
real(dp) :: orientation_angle
    real(dp), parameter                   :: charge_O = -0.834_dp
    real(dp), parameter                   :: charge_H = 0.417_dp
real(dp) :: sum_counts
real(dp) :: dist, exp_arg

  my_accumulate = optional_default(.false., accumulate)
  if (.not. my_accumulate) histogram = 0.0_dp

  use_dipole = optional_default(.true.,use_dipole_rather_than_angle_bisector)
  if (use_dipole) then
     call print("Using dipole to calculate angle with the surface normal.",PRINT_VERBOSE)
  else
     call print("Using HOH angle bisector to calculate angle with the surface normal.",PRINT_VERBOSE)
  endif

  !pos bin labels along axis
  pos_bin_width = at%lattice(axis,axis) / real(n_pos_bins,dp)
  if (present(pos_bin)) then
    do sample_i=1, n_pos_bins
      sample_r = (real(sample_i,dp)-0.5_dp)*pos_bin_width !the middle of the bin
      pos_bin(sample_i) = sample_r
    end do
  endif

 !angle bin labels
  angle_bin_width = PI/n_angle_bins
  if (present(pos_bin)) then
    sum_w = 0._dp
    do sample_i=1, n_angle_bins
      sample_angle = (real(sample_i,dp)-0.5_dp)*angle_bin_width !the middle of the bin
      angle_bin(sample_i) = sample_angle
      !the normalised solid angle is (1/4pi) * 2pi * sin((fi_north)-sin(fi_south)) where fi e [-pi/2,pi/2]
      angle_bin_w(sample_i) = 0.5_dp * ( sin(0.5_dp*PI - real(sample_i-1,dp)*angle_bin_width) - &
                                         sin(0.5_dp*PI - real(sample_i, dp)*angle_bin_width) )
      sum_w = sum_w + angle_bin_w(sample_i)
    end do
    angle_bin_w(1:n_angle_bins) = angle_bin_w(1:n_angle_bins) / sum_w
  endif

  !shift silica slab to the edges, water in the middle    || . . . ||
  if (silica_center_i>0) then ! silica
     call shift_silica_to_edges(at, axis, silica_center_i)
  endif

  !!calc_connect now, before shifting positions to positive, because it would remap the positions!!
  !call calc_connect including the H-bonds
  call set_cutoff(at,0._dp)
  call calc_connect(at)

  !now bins are from -axis/2 to axis/2
  !shift everyone to positive coordinate along axis
  at%pos(axis,1:at%N) = at%pos(axis,1:at%N) + 0.5_dp * at%lattice(axis,axis)

  !simply check water orientations and bin them


!  if (gaussian_smoothing) call system_abort('not implemented.')

  surface_normal(1:3) = 0._dp
  surface_normal(axis) = 1._dp

  num_atoms = 0
  do O=silica_center_i+1, at%N !only check water molecules
     if(at%Z(O)==1) cycle !find O
     !Count the atoms
     call print('Found O'//O//ElementName(at%Z(O)),PRINT_ANAL)
     num_atoms = num_atoms + 1

     !find H neighbours
     n = n_neighbours(at,O)
     if (n.ne.2) then ! O with =2 nearest neighbours
        call print("WARNING! water(?) oxygen with "//n//"/=2 neighbours will be skipped!")
        cycle
     endif
     H1 = neighbour(at,O,1)
     H2 = neighbour(at,O,2)
     if ((at%Z(H1).ne.1).or.(at%Z(H2).ne.1)) then !2 H neighbours
        call print("WARNING! water(?) oxygen with non H neighbour will be skipped!")
        cycle
     endif
 
     !We've found a water molecule.

     !Build histogram
     r = at%pos(axis,H1) !the position of O
     !HOH_angle = angle(diff_min_image(at,O,H1), &
     !                  diff_min_image(at,O,H2)) !the H-O-H angle

     if (.not. use_dipole) then
     !VERSION 1.
     !the direction of the HH->O vector, the bisector of the HOH angle
     !vector that is compared to the surface normal:
     !  point from the bisector of the 2 Hs (scaled to have the same bond length)
     !        to the O
         vector_OH1(1:3) = diff_min_image(at, O, H1)
         if (norm(vector_OH1) > 1.2_dp) &
            call system_abort('too long OH bond? '//O//' '//H1//' '//norm(vector_OH1))
         vector_OH2(1:3) = diff_min_image(at, O, H2)
         if (norm(vector_OH2) > 1.2_dp) &
            call system_abort('too long OH bond? '//O//' '//H2//' '//norm(vector_OH1))
         bisector_vector(1:3) = vector_OH1(1:3) / norm(vector_OH1) * norm(vector_OH2)

         ! a.dot.b = |a|*|b|*cos(angle)
         orientation_angle = dot_product((bisector_vector(1:3)),surface_normal(1:3)) / &
                             sqrt(dot_product(bisector_vector(1:3),bisector_vector(1:3))) / &
                             sqrt(dot_product(surface_normal(1:3),surface_normal(1:3)))
     else ! use_dipole

     !VERSION 2.
     !the dipole of the water molecule = sum(q_i*r_i)
     !Calculate the dipole and its angle compared to the surface normal

         !
         dipole(1:3) = ( diff_min_image(at,O,H1)*charge_H + &
                         diff_min_image(at,O,H2)*charge_H )
         if (norm(diff_min_image(at,O,H1)).gt.1.2_dp) call system_abort('too long O-H1 bond (atoms '//O//'-'//H1//'): '//norm(diff_min_image(at,O,H1)))
         if (norm(diff_min_image(at,O,H2)).gt.1.2_dp) call system_abort('too long O-H2 bond (atoms '//O//'-'//H2//'): '//norm(diff_min_image(at,O,H2)))
!call print ('dipole '//dipole(1:3))
    
         ! a.dot.b = |a|*|b|*cos(angle)
         orientation_angle = dot_product((dipole(1:3)),surface_normal(1:3)) / &
                             sqrt(dot_product(dipole(1:3),dipole(1:3))) / &
                             sqrt(dot_product(surface_normal(1:3),surface_normal(1:3)))
     endif

     if (orientation_angle.gt.1._dp) then
        call print('WARNING | correcting cos(angle) to 1.0 = '//orientation_angle)
        orientation_angle = 1._dp
     else if (orientation_angle.lt.-1._dp) then
        call print('WARNING | correcting cos(angle) to -1.0 = '//orientation_angle)
        orientation_angle = -1._dp
     endif
     orientation_angle = acos(orientation_angle)
     if (orientation_angle.lt.0._dp) then
        call print('WARNING | correcting angle to 0.0: '//orientation_angle)
        orientation_angle = 0._dp
     endif
     if (orientation_angle.gt.PI) then
        call print('WARNING | correcting angle to pi : '//orientation_angle)
        orientation_angle = PI
     endif
!call print ('angle '//(orientation_angle*180._dp/pi))

     call print('Storing angle for water '//O//'--'//H1//'--'//H2//' with reference = '//round(orientation_angle,5)//'degrees',PRINT_ANAL)
     call print('   with distance -1/2 b -- '//O//' = '//round(r,5)//'A',PRINT_ANAL)

     if (gaussian_smoothing) then !smear the position along axis
        !call system_abort('not implemented.')
        do sample_i=1, n_pos_bins
          sample_r = (real(sample_i,dp)-0.5_dp)*pos_bin_width
          dist = abs(r - sample_r)
          !Include all the atoms, slow but minimises error
          !	  if (dist > 4.0_dp*gaussian_sigma) cycle
            exp_arg = -0.5_dp*(dist/(pos_gaussian_sigma))**2
            if (exp_arg > -20.0_dp) then ! good to about 1e-8
              histogram(int(orientation_angle/angle_bin_width)+1,sample_i) = histogram(int(orientation_angle/angle_bin_width)+1,sample_i) + exp(exp_arg)/(pos_gaussian_sigma*sqrt(2.0_dp*PI)) !Gaussian in 1 dimension
            endif
        end do ! sample_i

     else !no gaussian_smoothing
        histogram(int(orientation_angle/angle_bin_width)+1,int(r/pos_bin_width)+1) = histogram(int(orientation_angle/angle_bin_width)+1,int(r/pos_bin_width)+1) + 1._dp
     endif

  end do ! O

  !normalise for the number of molecules in each pos_bin
  do sample_i=1,n_pos_bins
     sum_counts = sum(histogram(1:n_angle_bins,sample_i))
     if (sum_counts /= 0) histogram(1:n_angle_bins,sample_i) = histogram(1:n_angle_bins,sample_i) / sum_counts
  enddo

  !normalise for different solid angles of each angle_bin
  do sample_i=1, n_angle_bins
     histogram(sample_i,1:n_pos_bins) = histogram(sample_i,1:n_pos_bins) / angle_bin_w(sample_i)
  enddo 

end subroutine water_orientation_calc

end module structure_analysis_routines_module
