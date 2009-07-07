!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X     libAtoms: atomistic simulation library
!X     
!X     Copyright 2006-2007.
!X
!X     Authors: Gabor Csanyi, Steven Winfield, James Kermode
!X     Contributors: Noam Bernstein, Alessio Comisso
!X
!X     The source code is released under the GNU General Public License,
!X     version 2, http://www.gnu.org/copyleft/gpl.html
!X
!X     If you would like to license the source code under different terms,
!X     please contact Gabor Csanyi, gabor@csanyi.net
!X
!X     When using this software, please cite the following reference:
!X
!X     http://www.libatoms.org
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!X
!X  Atoms module
!X
!% An atoms object contains atomic numbers, all dynamical variables
!% and connectivity information for all the atoms in the simulation cell. 
!% It is initialised like this:
!%> 	  call initialise(MyAtoms,N,lattice)
!% where 'N' is the number of atoms to allocate space for and 'lattice' is a $3\times3$
!% matrix of lattice vectors given as column vectors.
!% 
!% Atoms also contains a Connection object, which stores distance information about
!% the atom neghbours after 'calc_connect' has been called. Rather than using a minimum
!% image convention, all neighbours are stored up to a radius of 'cutoff', including images
!% 
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! $Id: Atoms.f95,v 1.171 2008-07-11 09:58:42 gc121 Exp $
! $Log: not supported by cvs2svn $
! Revision 1.170  2008/06/24 13:49:49  nb326
! atoms_print_cfg_filename() opens file with action=OUTPUT
!
! Revision 1.169  2008/06/23 21:49:58  nb326
! atoms_print_xyz_filename() passes action=OUTPUT for inoutput_initialise
!
! Revision 1.168  2008/06/16 13:59:57  jrk33
! Bug fixed in set_atoms: species was not being set
!
! Revision 1.167  2008/06/16 13:04:32  jrk33
! Overhauled add_atoms. Shortened argument names to match fields in atoms, and now uses blank_rows argument to Table append()
!
! Revision 1.166  2008/05/27 08:12:59  jrk33
! Added string and logical property support to add_atoms - not yet tested.
!
! Revision 1.165  2008/05/24 18:15:44  jrk33
! Incresed some string lengths in print_xyz. Removed cluster printing in create_cluster, since it gets trigged inside metapotential_minim and generates too much output
!
! Revision 1.164  2008/05/19 20:40:35  nb326
! Add optional_argument max_dist to atoms_neighbour
!
! Revision 1.163  2008/05/16 13:29:59  jrk33
! Trivial change in create_cluster
!
! Revision 1.162  2008/05/16 13:27:17  jrk33
! Fixed species bug in create_cluster
!
! Revision 1.161  2008/05/15 01:42:23  nb326
! Fix cell_image_N[abc] for situation where cutoff is large relative to periodicity
!
! Revision 1.160  2008/05/13 16:12:27  jrk33
! Updated docs
!
! Revision 1.159  2008/05/13 15:51:55  jrk33
! Added support for string and logical property types, and extended print_xyz, read_xyz and print_cfg appropriately
!
! Revision 1.158  2008/05/05 13:39:58  jrk33
! Changed table_allocate calls to reflect changes to Table.f95
!
! Revision 1.157  2008/05/02 18:22:06  jrk33
! Incresed default real precision in print_xyz to 5 digits; added complement and difference subroutines
!
! Revision 1.156  2008/04/23 19:35:33  nb326
! Gracefully handle situation where there are no neighbors and MPI is on
!
! Revision 1.155  2008/04/18 15:02:51  jrk33
! Check for an integer valued lattice
!
! Revision 1.154  2008/03/26 14:21:31  jrk33
! Reverted changes to termindex
!
! Revision 1.153  2008/03/19 16:33:05  jrk33
! Fixed bug in create_cluster: even_hydrogens option had stopped working with changes of revision 1.132
!
! Revision 1.152  2008/03/12 09:52:15  jrk33
! Added list_to_property and property_to_list routines
!
! Revision 1.151  2008/03/08 00:02:31  ab686
! simplified implied array construction in remove_atom_single
!
! Revision 1.150  2008/03/07 23:47:02  ab686
! added set_atoms_singlez routine, fixed reshape(optional) bug in add_atoms, fixed bug in add_atoms by adding missing finalise(connect) at the end
!
! Revision 1.149  2008/03/06 17:22:47  jrk33
! changed real_format for parameter writing to f18.6 in print_xyz
!
! Revision 1.148  2008/02/26 17:29:02  jrk33
! Added has_property() convenience function
!
! Revision 1.147  2008/02/26 11:21:52  jrk33
! Remove misc_data from Atoms
!
! Revision 1.146  2008/02/21 11:27:42  jrk33
! Fixed typo in docs
!
! Revision 1.145  2008/02/14 22:57:12  nb326
! this%params=params in read_xyz _after_ atoms_initialise
!
! Revision 1.144  2008/02/12 18:54:40  jrk33
! Updated documentation on atoms%params
!
! Revision 1.143  2008/02/12 18:20:12  jrk33
! Added Atoms%params dictionary, initialised from line 2 of an xyz file
!
! Revision 1.142  2008/02/11 18:44:23  nb326
! Fix spacing
!
! Revision 1.141  2008/02/06 17:22:05  nb326
! In atoms_select, prevent memory leak by finalising to%properties before copying from%properties onto it
!
! Revision 1.140  2008/02/04 14:56:33  jrk33
! atoms_read_xyz: Default to ElementMass(Z) if no mass property present in file
!
! Revision 1.139  2008/02/04 14:07:34  jrk33
! Decoupled atomic number and mass to allow mixed isotope simulations. New mass property must be set
!
! Revision 1.138  2008/02/04 13:30:56  saw44
! Fixed comment typo
!
! Revision 1.137  2007/12/13 17:08:28  nb326
! bfs_grow_list only finalises temp if n >= 1
!
! Revision 1.136  2007/12/05 16:35:49  nb326
! split cell_volume into atoms_cell_volume and lattice_cell_volume
!
! Revision 1.135  2007/11/28 15:34:29  nb326
! More digits in human_readable atomic positions and in connection print
!
! Revision 1.134  2007/11/21 14:55:00  jrk33
! Added per atom average kinetic energy
!
! Revision 1.133  2007/11/14 19:00:23  nb326
! Add optional argument own_neighbour to calc_connect
!
! Revision 1.132  2007/11/12 23:59:41  gc121
! thorough cleanup and debug of create_cluster. introduced do_periodic masking, works now for very narrow slabs. wonder how it ever worked before...
!
! Revision 1.131  2007/11/12 15:00:57  saw44
! Changed my_real_format in print_xyz to have set length of 255 rather than variable length which depended on an optional argument, causing g95 to complain
!
! Revision 1.130  2007/11/08 19:17:11  nb326
! Silly bug in atoms_set_cutoff_factor
!
! Revision 1.129  2007/11/08 19:12:16  nb326
! make set_cutoff_minimum a module  procedure
!
! Revision 1.128  2007/11/08 19:11:27  nb326
! set_cutoff with 0 reverts to default (non uniform, default factor), and add set_cutoff_minimum which only increases cutoff
!
! Revision 1.127  2007/11/07 20:36:55  nb326
! Call wipe_cells in connection_wipe, and call connection_initialise in calc_connect if cutoff .feq. 0
!
! Revision 1.126  2007/11/07 19:06:16  saw44
! Added initialised=true to cutoff=0.0 section of calc_connect
!
! Revision 1.125  2007/11/07 17:20:04  nb326
! Optional real_format passed to atoms_print_xyz
!
! Revision 1.124  2007/11/07 10:11:19  gc121
! map_into_cell now checks to see if connections were initialised before it tries to update the connection table. it updates the travel in any case
!
! Revision 1.123  2007/11/05 21:52:42  nb326
! Add parallel read capability to atoms_read_xyz_filename and atoms_read_xyz_inoutput
!
! Revision 1.122  2007/10/22 18:10:00  gc121
! in create_cluster we zero out certain shifts in directions in which the cluster is periodic. this now fixes a bug, but the code is still a mess with respect to periodicity.
!
! Revision 1.121  2007/10/22 16:20:27  jrk33
! ends of functions and subroutines
!
! Revision 1.120  2007/10/22 16:17:09  jrk33
! ends of functions and subroutines
!
! Revision 1.119  2007/10/17 16:19:45  nb326
! Remove debugging statements
!
! Revision 1.118  2007/10/16 22:26:31  nb326
! atoms_read_xyz can now read from appropriately formatted extendable_str
!
! Revision 1.117  2007/10/11 18:01:40  nb326
! Use mpi_n_procs() and mpi_id() instead of calling get_mpi_size_rank
!
! Revision 1.116  2007/10/11 16:30:40  nb326
! Replace decode_mpi_error with abort_on_mpi_error
!
! Revision 1.115  2007/10/09 16:24:32  nb326
! More fields in read_xyz, and longer properties string in print_xyz
!
! Revision 1.114  2007/10/09 13:34:04  jrk33
! Fixed bug in create_cluster for periodic clusters with cluster larger than atoms structure
!
! Revision 1.113  2007/10/04 15:41:47  nb326
! Use increase_stack, not c_increase_stack
!
! Revision 1.112  2007/09/18 18:27:47  saw44
! Removed printing from parse_atom_mask
!
! Revision 1.111  2007/09/17 15:36:37  nb326
! Try to use c_increase_stack() to increase stack limit if necessary
!
! Revision 1.110  2007/09/14 13:09:48  jrk33
! Applied checkcode, added missing names to ends of functions and subroutines
!
! Revision 1.109  2007/09/06 14:26:13  gc121
! missing finalise statement
!
! Revision 1.108  2007/09/04 14:05:44  nb326
! Finalise props in atoms_print_cfg to prevent mem leaks
!
! Revision 1.107  2007/08/30 14:31:42  gc121
! changed a comment
!
! Revision 1.106  2007/08/30 14:27:35  nb326
! Move diamond, fcc, supercell, make_structure to Structures.f95
!
! Revision 1.105  2007/08/30 14:04:36  nb326
! Fix spacing in human readable output of atoms
!
! Revision 1.104  2007/08/24 10:10:26  nb326
! Finalise local props in atoms_read_xyz to prevent memory leaks
!
! Revision 1.103  2007/08/16 12:48:22  nb326
! Use double precision in matrix3x3_inverse
!
! Revision 1.102  2007/08/16 12:45:05  saw44
! Fixed docs for find_motif
!
! Revision 1.101  2007/08/14 15:05:11  nb326
! Marginally ore debugging output in atoms_initialise
!
! Revision 1.100  2007/08/10 15:35:26  jrk33
! Fixed real properties in atoms_intialise_pointers. Changed connection_print "out" argument to "file". Updated connection_file_write to use new Table array binary_writes for speed
!
! Revision 1.99  2007/08/10 15:27:19  gc121
! now we return immediate form atoms_repoint if this%N==0, but for now I didnt touch the zero-sized pointer assignments, until we conclude the argument about this routine
!
! Revision 1.98  2007/08/10 13:09:16  jrk33
! Removed hardcoded sizes for default_int and default_real in atoms_initialise
!
! Revision 1.97  2007/08/10 10:36:34  jrk33
! Changed atoms_initialise_pointeers to point into zero length arrays correctly if this%N == 0
!
! Revision 1.96  2007/08/09 22:05:08  nb326
! Add space in calc_connect output
!
! Revision 1.95  2007/08/09 17:29:16  gc121
! Add fcc() which creates fcc structure in simple-cubic supercell
!
! Revision 1.94  2007/08/08 17:29:47  gc121
! fixed quotes in a string
!
! Revision 1.93  2007/08/03 00:26:45  ab686
! broken some lines for the broken compiler by Sun
!
! Revision 1.92  2007/08/02 23:39:38  gc121
! removed intent attribute from pointer assignments, several compilers couldnt cope with it
!
! Revision 1.91  2007/07/30 14:19:00  nb326
! Move code that actually adds bond to connect to add_bond(), and create set_bonds()
!
! Revision 1.90  2007/07/27 10:54:19  nb326
! A bit more debugging output on NERD in calc_connect
!
! Revision 1.89  2007/07/24 16:13:40  saw44
! Added find_motif and helper routines. These perform pattern matching based on a paper by J. R. Ullmann and find occurences of a sub structure within an atoms structure
!
! Revision 1.88  2007/07/20 16:45:56  gc121
! swapped operands of a + for clarity in map_into_cell, added _dp
!
! Revision 1.87  2007/07/18 15:02:05  nb326
! Use new verbosity system (intro.tex v 1.13)
!
! Revision 1.86  2007/07/18 13:06:35  nb326
! Use new verbosity system (intro.tex v 1.13)
!
! Revision 1.85  2007/07/17 10:35:35  gc121
! put some ANAL prints inside ifs
!
! Revision 1.84  2007/07/17 08:45:28  gc121
! added printing of atom name in human_readable printing
!
! Revision 1.83  2007/07/16 15:11:18  nb326
! Test for zero cutoff in calc_connect
!
! Revision 1.82  2007/07/14 08:34:32  gc121
! fixed a bug in atoms_initialise: this%N was set early, and in case of no data present, add_property crashed because it tried to reference data%int(1:this%N), even though data was not yet initialised
!
! Revision 1.81  2007/07/13 14:45:08  jrk33
! Fixes for g95: two counts of repeated + split over line break and (/1:3/) -> (/1,2,3/)
!
! Revision 1.80  2007/07/13 14:22:56  jrk33
! Changed print(line,outfile) to print(line,file=outfile) throughout
!
! Revision 1.79  2007/07/12 09:57:47  nb326
! Add remove_property()
!
! Revision 1.78  2007/07/11 15:46:43  nb326
! Don't print 'Lattice vectors:...' or atomic positions if current verbosity is too low
!
! Revision 1.77  2007/07/11 12:40:51  gc121
! added set_cutoff and set_cutoff_factor interfaces. changed default cluster vacuum to 10A
!
! Revision 1.76  2007/07/10 15:41:30  nb326
! assign_pointer is now a logical function, returning true if successful and false if property doesn't exist.  Still aborts for mismatch in pointer type
!
! Revision 1.75  2007/07/10 15:17:18  gc121
! renamed find(mask) to find_indices(mask) and added a find() interface
!
! Revision 1.74  2007/07/10 15:11:26  gc121
! deleted some commented lines of code, and made create_cluster accept an atomlist table with just one int column, in which case the shifts are assumed to be zero
!
! Revision 1.73  2007/07/09 17:05:16  jrk33
! Updated documentation to remove reference to find_neighbour
!
! Revision 1.72  2007/07/09 16:34:26  saw44
! Fixed memory leaks: added finalise(props) to atoms_print_xyz, added finalise(big_data) to supercell
!
! Revision 1.71  2007/07/09 14:45:22  gc121
! added functionality to create_cluster needed by the local energy code
!
! Revision 1.70  2007/07/07 12:36:50  jrk33
! Fixed formatting error in atoms_print for long property names
!
! Revision 1.69  2007/07/03 12:54:27  jrk33
! Added bfs_grow_single, moved bfs_grow to bfs_grow_list and added bfs_grow interface
!
! Revision 1.68  2007/06/22 14:40:07  nb326
! Fix typo
!
! Revision 1.67  2007/06/22 14:27:35  nb326
! Fix long line in comment, and add option to skip in-out-in pattern removal in create_cluster
!
! Revision 1.66  2007/06/12 16:25:10  saw44
! tidied code slightly in distance
!
! Revision 1.65  2007/06/08 13:34:36  jrk33
! Corrected bug in create_cluster: same_lattice -> do_same_lattice to prevent seg fault if arguement not present
!
! Revision 1.64  2007/06/08 13:10:20  jrk33
! Modified atoms_select to accept either a mask or a list of atoms to keep, and removed subset since functionality is now duplicated by select
!
! Revision 1.63  2007/06/07 16:58:48  nb326
! Add atoms_select() using a mask, and fix create_cluster() for periodic cluster for non-rectangular supercells
!
! Revision 1.62  2007/05/16 13:19:56  jrk33
! Downgraded the verbosity of some messages. Fixed bug in create_cluster, map into cell was called before new lattice defined
!
! Revision 1.61  2007/05/16 11:31:24  jrk33
! Simplied create_cluster: it now complains if there are repeated atoms and you want a periodic cluster, rather than removing the images itself. Added call to map_into_cell to ensure final cluster positions are within unit cell.
!
! Revision 1.60  2007/05/15 13:16:45  jrk33
! Changed bfs_grow min_images_only option to first add all the images then remove duplicates to avoid missing atoms
!
! Revision 1.59  2007/05/15 10:20:22  jrk33
! Added min_images_only option to bfs_step and bfs_grow to only include minimum image atoms in lists, and discard_non_min_images subroutine.
!
! Revision 1.58  2007/05/14 15:32:33  saw44
! Added parse_atom_mask routine. Useful for selecting subsets of atoms from command line arguments
!
! Revision 1.57  2007/05/10 16:47:33  jrk33
! bfs_step bug fix: should check if atoms are in output AND input list when removing duplicate images
!
! Revision 1.56  2007/05/10 15:53:35  jrk33
! Added function flatten_images
!
! Revision 1.55  2007/05/10 14:08:52  jrk33
! Removed ptr argument of atoms_add_property subroutines and replaced with atoms_assign_pointer interface to 
! avoid pointers being made invalid by subsequent add_property calls. bfs_step changed to match width of output 
! table to that of input, i.e. only include shifts in output if they were in input. Same for make_convex. 
! create_cluster modified to deal with periodicity correctly
!
! Revision 1.54  2007/05/08 13:46:41  jrk33
! Added optional even_hydrogens argument to create_cluster to remove a hydrogen from cluster if necessary to ensure n_up electrons = n_down electrons
!
! Revision 1.53  2007/05/04 13:24:12  jrk33
! Fixed bugs in bfs_step for case of inputr list provided without shifts
!
! Revision 1.52  2007/05/04 12:33:03  jrk33
! Updated directionality and make_convex to work with multiple image convention
!
! Revision 1.51  2007/05/03 16:30:20  jrk33
! Fixed create_cluster for multiple image convention. New periodic option can be used to specifiy if you want a periodic cluster or one with multiple images in
!
! Revision 1.50  2007/05/03 12:39:14  jrk33
! Fixed a couple of bugs in create_cluster
!
! Revision 1.49  2007/05/03 10:39:26  jrk33
! Various i < j converted to i <= j to allow for self-images
!
! Revision 1.48  2007/05/03 08:20:38  gc121
! corrected bug: forgot to include self images
!
! Revision 1.47  2007/05/02 23:44:44  gc121
! updated connection_print to print index entry of neighbour2 list; put in a distance() function that returns the real distance given two atoms and a shift
!
! Revision 1.46  2007/05/02 22:46:13  gc121
! fixed some bugs in the new calc_connect
!
! Revision 1.45  2007/05/02 17:35:28  jrk33
! Some bugs in calc_connect fixed. bfs_step and create_cluster made multiple image aware, but the latter not working properly yet
!
! Revision 1.44  2007/04/27 15:08:31  jrk33
! Readded an overzealously removed variable
!
! Revision 1.43  2007/04/27 15:02:46  jrk33
! Removed unused variables
!
! Revision 1.42  2007/04/27 14:55:31  jrk33
! Removed symmetry support
!
! Revision 1.41  2007/04/27 11:58:56  gc121
! Massive change to the way we deal with neighbour connections. We no longer follow minimum image conventions, but store all images of neighbouring atoms within a (perhaps huge) cutoff. Much code is expected to be now broken.
!
! Revision 1.40  2007/04/20 17:22:18  jrk33
! Added () within // expression
!
! Revision 1.39  2007/04/18 12:40:49  jrk33
! Updated print routines to new parameter convention
!
! Revision 1.38  2007/04/18 10:06:52  gc121
! added too_few_cells_warning_issued flag to Connection, so that we print the warning only once
!
! Revision 1.37  2007/04/18 09:02:23  jrk33
! Changed function and subroutine names to fit nnew naming convention
!
! Revision 1.36  2007/04/18 08:50:37  jrk33
! Convert function and subroutine names to lower case
!
! Revision 1.35  2007/04/18 01:30:27  gc121
! updated to reflect changes in printing and other naming conventions
!
! Revision 1.34  2007/04/17 09:57:19  gc121
! put copyright statement in each file
!
! Revision 1.33  2007/04/11 15:44:17  saw44
! Updated/Added comments for the documentation generator
!
! Revision 1.32  2007/04/03 16:25:47  jrk33
! Bug fix to Atoms_read_xyz - missing trim() calls
!
! Revision 1.31  2007/04/03 14:02:01  jrk33
! Updated doc comments
!
! Revision 1.30  2007/03/30 16:46:11  jrk33
! Modified print argument order to conform with changes to System
!
! Revision 1.29  2007/03/28 17:42:58  saw44
! Fixed comment about the dictionary in Atoms
!
! Revision 1.28  2007/03/28 16:00:11  jrk33
! Corrected documentation comment. Added optional lookup argument to Add_Property to return column reference of the new property. Changed properties argument of print and read subroutines to be a colon separated list rather than an array
!
! Revision 1.27  2007/03/27 15:12:12  jrk33
! Reformatted long lines in log entries
!
! Revision 1.26  2007/03/27 15:10:11  jrk33
! Moved all atomic properties to data table. The properties dictionary maps field names to 
! columns in this table. move_mask, damp_mask and thermostat_mask are now of type integer to 
! enable them to fit in the table. Add_Property interface is used to add columns to this data 
! table - this is called by Initialise to add pos, travel and Z entries and again by DS_Initialise 
! to add velo, acc, etc. Pointers into data are maintaned for commonly used properties. Addition and 
! removal of atoms has been updated. Print, Print_xyz and print_cfg have been heavily overhauled to
! facilitate specifying which columns to print. Read_xyz has also been updated to read properties. 
! ReadB and WriteB are updated and will not be compatible with files written previosuly. 
! Atoms_Cutoff query function added.                                                                 
!
! Revision 1.25  2007/03/21 16:29:13  jrk33
! Added optional print_avgpos and print_realpos arguments to print_xyz and print_ext_xyz
!
! Revision 1.24  2007/03/21 15:23:55  jrk33
! Moved CLUSTER_VACUUM to be optional argument of Create_Cluster. Moved NNEIGHTOL to be member 
! of Atoms type. Added Atoms_Cutoff function. Added Zero interface to Atoms_Zero.
!
! Revision 1.23  2007/03/13 15:04:33  jrk33
! cosines -> Direction_Cosines
!
! Revision 1.22  2007/03/12 16:52:18  jrk33
! Removed Set_Atoms interface. Edited documentation. Converted all intent(IN/OUT/INTOU) to lowercase
!
! Revision 1.21  2007/03/06 13:40:11  nb326
! Fix typos in documentation
!
! Revision 1.20  2007/03/01 13:51:45  jrk33
! Documentation comments reformatted and edited throughout. Anything starting "!(no space)%" 
! is picked up by the documentation generation script
!
! Revision 1.19  2007/02/21 11:16:02  saw44
! Replaced logical == .false. with .not.logical (keeps g95 happy)
!
! Revision 1.18  2007/02/16 10:31:24  saw44
! Removed use of <> in format strings - not supported by all compilers
!
! Revision 1.17  2007/02/13 15:11:08  saw44
! Replaced multiple finalise routines with one routine which can finalise up to 10 atoms structures. 
! Modified Cells_Initialise to optionally take a number of atoms, which it uses to set the initial 
! length and increment of the cell lists. Removed some code from CalcConnectFast and put into 
! Wipe_Cells and Partition_Atoms, since these are helpful even when CalcConnectFast is no being
!  used (e.g. rdf calculations :-))
!
! Revision 1.16  2007/01/24 11:20:29  saw44
! Removed unused variables
!
! Revision 1.15  2007/01/17 14:23:44  jrk33
! Added optional parallel argument to CalcDists. Defaults to serial operation. Set parallel 
! to true when called from DynamicalSystem.advanceVerlet
!
! Revision 1.14  2007/01/09 16:55:08  saw44
! Added Is_Nearest_Neighbour logical function and fixed error messages in Directionality
!
! Revision 1.13  2007/01/09 16:12:41  jrk33
! Fixed bug for case of zero aux props in Atoms_print_ext_cfg
!
! Revision 1.12  2007/01/03 16:43:21  jrk33
! Removed unused variable local_g from distance8_vec_vec. Check if there is already a Lattice 
! when printing comment to avoid printing twice
!
! Revision 1.11  2007/01/03 10:38:11  nb326
! MPI is a C++ namespace, replace ifdef MPI with ifdef _MPI
!
! Revision 1.10  2006/12/20 18:45:49  jrk33
! Write slightly smaller CFG files by ommitted element mass and name if the same as previous atom
!
! Revision 1.9  2006/12/20 18:34:31  jrk33
! Added Atoms_Print_ext_cfg subroutine to print in AtomEye extended CFG format
!
! Revision 1.8  2006/12/14 16:51:45  saw44
! Added Subset function, which returns the subset of the input atoms given in a table
!
! Revision 1.7  2006/12/12 18:59:31  jrk33
! Fixed bug in Atoms_Read_Ext_xyz - if status was not present always aborted
!
! Revision 1.6  2006/12/12 00:14:44  gc121
! fixed errors, now compiles
!
! Revision 1.5  2006/12/11 23:27:51  gc121
! Moved (:,N) variables from DynamicalSystem to Atoms
!
! Revision 1.4  2006/12/06 18:07:48  nb326
! Correctly handle case when status isn't passed to atoms_read_xyz
!
! Revision 1.3  2006/12/04 17:17:43  gc121
! got changes from LOTF95 that were checked in late
!
! Revision 1.2  2006/12/04 16:21:15  nb326
! diff_at_at sets (optional) shift vector
!
! Revision 1.1.1.1  2006/12/04 11:11:30  gc121
! Imported sources
!
! Revision 1.110  2006/11/28 11:37:40  saw44
! Removed Make_Convex interface: it contained a subroutine and a function, which breaks the 
! Fortran standard and makes g95 complain (although ifort compiled it fine)
!
! Revision 1.109  2006/11/27 14:27:11  saw44
! Added Set_Lattice interface
!
! Revision 1.108  2006/11/27 09:28:41  jrk33
! Create_Cluster now returns atomlist of correct intsize (1)
!
! Revision 1.107  2006/11/24 17:31:41  jrk33
! Check if nneighb_only is present in BFS_grow
!
! Revision 1.106  2006/11/17 13:03:18  saw44
! Added _ to AddAtoms and RemoveAtoms. Added more erro checking to Make_Lattice
!
! Revision 1.105  2006/11/13 23:32:19  gc121
! extended convexity routines
!
! Revision 1.104  2006/11/06 17:37:44  jrk33
! Added Make_Convex function (untested)
!
! Revision 1.103  2006/10/30 20:47:52  gc121
! removed comment no longer relevant
!
! Revision 1.102  2006/10/27 16:38:39  saw44
! Added Set_Lattice routine which changes the lattice, keeping the inverse up-to-date, and 
! optionally recalculates the connectivity or maps atoms back into the cell. The atomic 
! positions are left unchanged. The inverse lattice is now part of the atoms object and 
! doesnt need passing around by TestConnect etc. Corrected typos
!
! Revision 1.101  2006/10/17 09:54:59  jrk33
! NNEIGHTOL reverted to 1.2
!
! Revision 1.100  2006/10/16 16:01:10  gc121
! overloaded print_ext_xyz, added more debugging printint to it
!
! Revision 1.99  2006/10/12 08:31:08  jrk33
! Removed int_part and real_part function calls from Atoms_Read_xyz since they cause seg faults with large files
!
! Revision 1.98  2006/09/04 16:27:50  jrk33
! Parrallelised CalcDists since this was found to be bottleneck when called from advanceVerlet
!
! Revision 1.97  2006/08/09 13:39:26  jrk33
! Moved cluster symmetrising inside NASTINESS II loop so that we iterate until cluster stops
! growing due to either IN-OUT-IN neighbours or missing symmetric image atoms
!
! Revision 1.96  2006/08/02 10:43:13  saw44
! Added Make_Structure function: Create an atoms object from lattice, lattice type (i.e. P,I,F,A,B,C) and motif
!
! Revision 1.95  2006/07/21 13:44:53  saw44
! Changed diff routine to interface: can now compute diff vector from pos to pos, atom to atom, 
! or mixture. Changed Lattice_Inverse to Matrix3x3_Inverse, since this is being used on matrices 
! other than lattice vectors now.
!
! Revision 1.94  2006/07/18 14:02:34  saw44
! Added SVD method to the Directionality calculator
!
! Revision 1.93  2006/07/17 14:23:31  saw44
! Replaced Coplanarity stuff with Directionality stuff and CosAngle_To_Line
!
! Revision 1.92  2006/07/06 21:50:59  jrk33
! Added reading/writing of symmetry data to ReadB and WriteB. We only save tthe real parts of atoms%pos 
! and atoms%Z to reduce file size. Changed Create_Cluster to work better with symmetric clusters - we 
! now make sure clusters themselves are symmetric after dealing with nastiness things by adding missing 
! atoms to make this so
!
! Revision 1.91  2006/07/04 16:36:28  jrk33
! Reformatted log entries
!
! Revision 1.90  2006/07/04 16:27:24  jrk33
! Fixed three problems with remove_images option to Create_Cluster. Check all images 
! of atoms in cluster are also in cluster before terminating to get symmetric 
! termination. Termination atoms were not correctly ordered since Table_Delete moves 
! last item over deleted item: fixed by adding to new table cluster_temp2. Lattice in 
! periodic direction corrected to be original lattice/N_repeat, not maxlen/N_repeat.
!
! Revision 1.89  2006/06/29 17:53:49  gc121
! catch zero volume condition
!
! Revision 1.88  2006/06/23 14:36:22  saw44
! Fixed Lattice_Inverse (returned inverse was transpose of actual inverse), added 
! error checking to print_(ext_)xyz
!
! Revision 1.87  2006/06/21 16:55:36  gc121
! overloaded read_xyz to take a filename, also added checking if conversion of atomic symbol was OK
!
! Revision 1.86  2006/06/20 17:23:18  gc121
! added new copyright notice to include James, Gian, Mike and Alessandro
!
! Revision 1.85  2006/06/14 10:26:43  saw44
! b and c are now optional in Make_Lattice, if not present they are assumed to be the same as a
!
! Revision 1.84  2006/06/08 15:44:00  jrk33
! Connection_Initialise picks sensible initial lengths and increments for the neighbour 
! tables based on density and cutoff of Atoms object
!
! Revision 1.83  2006/06/08 13:55:36  saw44
! Added interface to allow Find_Neighbour to be called with an atoms object. Added 
! Get_lattice_Params subroutine, which is the opposite of Make_Lattice. 
! Added Centre_Of_Mass function.
!
! Revision 1.82  2006/06/08 13:28:55  jrk33
! Bug fix to advanceVerlet for symmetric systems: atoms that had been remapped by 
! May_Into_Cell were not having their displacements correctly taken into account. 
! Now we use realpos instead of pos
!
! Revision 1.81  2006/06/08 11:05:22  jrk33
! Implemented more memory efficient neighbour storing scheme. Full neighbour details are 
! now only stored once, in the neighbour list of the atom with the lower index, and the 
! atom with the large index gets a reference to these details. Changes should be internal 
! to Atoms.f95 providing that you dont make direct access to the connect%neighbour table, 
! which has now been replaced by two tables: connect%neighbour1 and connect%neighbour2
!
! Revision 1.80  2006/06/07 16:22:50  jrk33
! Decreased default table lengths for connect%neighbour and connect%cell
!
! Revision 1.79  2006/06/06 13:05:55  jrk33
! Only print warning about stack size once, the first time Atoms_Initialise is called with a size > 4Mb
!
! Revision 1.78  2006/05/30 13:53:06  jrk33
! if (present(verbosity) .and. verbosity > NORMAL) then ... -> if (present(verbosity)) 
! then if (verbosity > NORMAL) then ...
!
! Revision 1.77  2006/05/30 11:09:44  jrk33
! Removed declarations for unused variables
!
! Revision 1.76  2006/05/25 13:37:51  jrk33
! Removed FIXME in Create_Cluster for non-orthorhombic cells (it's been fixed a while)
!
! Revision 1.75  2006/05/25 11:06:46  jrk33
! Added lattice_inverse: simple function to work out inverse of 3x3 to cut down overhead 
! of calling lapack
!
! Revision 1.74  2006/05/23 15:56:44  jrk33
! Changed size(extras%[int,real],1) to extras%intsize and extras%realsize in Atoms_Print_ext_xyz. 
! On Tru64, size() on unallocated arrays does not appear to return 0
!
! Revision 1.73  2006/05/23 10:29:54  jrk33
! Removed debugging write line from atoms_n_neighbours
!
! Revision 1.72  2006/05/18 10:46:46  jrk33
! Check if image_disp and image are allocated before copying them in Atoms_Assignment
!
! Revision 1.71  2006/05/16 16:46:01  saw44
! Added Cell_Volume function. Added status returning in Read_xyz. Put comment printing
! back into Print_xyz.
!
! Revision 1.70  2006/05/15 17:07:34  jrk33
! Fixed bug in Atoms_print_xyz and Atoms_print_ext_xyz that caused problems when
!  printing unsymmetrised Atoms objects
!
! Revision 1.69  2006/05/12 13:57:42  jrk33
! Periodic clusters now work with non orthorhombic cells (thanks Steve!)
!
! Revision 1.68  2006/05/11 19:59:17  saw44
! removed origin from supercell and diamond. still needs fixing though...
!
! Revision 1.67  2006/05/11 15:44:10  saw44
! Added checks to the N_Neighbours and Nth_Neighbour routines to make sure connectivity 
! data exists. Added the Fit_Box_In_Cell subroutine (see comments) and Make_Lattice 
! (to create a lattice from a,b,c,alpha,beta and gamma). Also added preliminary coplanarity 
! testing for atoms, but this code will change. Do not use for now. I also commented the 
! variables in WriteB and ReadB with numbers so that it is easier to match up the Read/Write order.
!
! Revision 1.66  2006/05/11 10:32:16  jrk33
! Now we check whether clusters should be made periodic due to small cell sizes (only
!  works for orthorhombic cells at the moment)
!
! Revision 1.65  2006/05/10 10:45:48  jrk33
! Code to deal with simple cubic symmetries added. Supercell_Symmetry creates a symmetric
! Atoms object. Create_Cluster, Print_xyz and Print_ext_xyz have been modified to only 
! consider real atoms if there is a symmetry in the system. Functions Is_Real_Atom and
! Is_Image_Atom identify atoms as real or images.
!
! Revision 1.64  2006/05/04 13:08:22  jrk33
! Added warning about stack size if size of atoms positions bigger than 4000k
!
! Revision 1.63  2006/05/02 17:06:44  saw44
! did some optimisation on CalcConnectFast, 8x speed increase. distance8_ interfaces have changed
!
! Revision 1.62  2006/04/26 13:39:00  jrk33
! Modified Create_Cluster to measure dimensions of cluster and set lattice accordingly
!
! Revision 1.61  2006/04/15 16:07:42  jrk33
! Changed length of comment parameter in atoms_read_ext_xyz from 1024 to * to be more general
!
! Revision 1.60  2006/04/13 18:51:45  saw44
! Added check to CalcConnect/Fast to test if the cutoff is too large. Fixed a mangling of
!  the code due to a CVS conflict
!
! Revision 1.59  2006/04/13 11:09:12  jrk33
! Added comment optional parameter to atoms_print_xyz
!
! Revision 1.58  2006/03/31 16:24:02  saw44
! Added extended xyz reader atoms_print_ext_xyz
!
! Revision 1.57  2006/03/24 11:36:24  saw44
! freed tmplist in bfs_grow: possible memory leak
!
! Revision 1.56  2006/03/10 00:34:24  gc121
! fixed bugs in Remove_Atomst_multiple: forgot to update this%N, flipped one loop exit
!  condition, made do loop go from 1,N, it was not terminated before.
!
! Revision 1.55  2006/03/09 16:25:21  gc121
! missed updating this%N and bug in updating this%Z in AddAtoms_multiple. Thanks Gian\!
!
! Revision 1.54  2006/03/03 16:16:54  gc121
! removed warning about nonexistent parameters from init(), as we do a save automatically 
! if we are mapping
!
! Revision 1.53  2006/02/27 17:25:18  gc121
! atoms_print now takes account of the different meanings of atoms%cutoff
!
! Revision 1.52  2006/02/27 17:15:27  gc121
! made the default cutoff = NNEIGHTOL in accordance with its new meaning as a multiplier 
! when use_uniform_cutoff = .false. (which is the default)
!
! Revision 1.51  2006/02/27 17:02:16  gc121
! commented out cutoff and use_uniform_cutoff variables from the Connection type, as 
! they were not used. Also, if atoms%use_uniform_cutoff is true, then cutoff now acts 
! as a multipler for bond_length(Zi,Zj) to determine the cutoff
!
! Revision 1.50  2006/02/27 11:06:46  gc121
! only print message in carve_cluster at VERBOSE level
!
! Revision 1.49  2006/02/21 12:17:36  saw44
! Added Atoms_Set_Atoms_char to take an atomic symbol, created a Set_Atoms interface 
! and made BFS_Step complain if no connectivity data is present
!
! Revision 1.48  2006/02/20 16:52:39  saw44
! Fixed stack overflow bug in create_cluster
!
! Revision 1.47  2006/02/17 12:40:49  saw44
! Updated Print interface to include Connection_Print
!
! Revision 1.46  2006/02/15 14:16:54  saw44
! Removed commas before first entry in write statements
!
! Revision 1.45  2006/02/08 17:53:29  saw44
! Added Atoms_Read_xyz
!
! Revision 1.44  2006/02/06 16:48:22  saw44
! General Code clean-up: routine names changed to match the changes in System and linearalgebra
!
! Revision 1.43  2006/02/02 20:10:11  gc121
!  put in if(allocated(..)) statements into atoms_finalise and connection_finalise()
!
! Revision 1.42  2006/02/02 17:24:19  saw44
! Changed CalcConnectFast to use Table_Wipe for speed
!
! Revision 1.41  2006/01/31 14:28:19  saw44
! Updated ReadB and WriteB argument order
!
! Revision 1.40  2006/01/31 13:57:47  gc121
! inlines some matrix.mult.int-vect calls. if-then on nerdy writes
!
! Revision 1.39  2006/01/31 12:35:33  saw44
! Binary file writing arguments switched back temporarily
!
! Revision 1.38  2006/01/30 13:04:15  gc121
! changed verbosity level of enter/exit messages
!
! Revision 1.37  2006/01/30 11:41:13  gc121
! swapped arguments of is_in_array
!
! Revision 1.36  2006/01/30 10:38:24  gc121
! adjusted verbosity level of some calls
!
! Revision 1.35  2006/01/29 13:12:23  saw44
! Create_Cluster now tests for nearest neighbours
!
! Revision 1.34  2006/01/28 20:47:50  gc121
! changed default to do_terminate=.true. in create_cluster(). there is still a bug, 
! it segfaults when cutoff is much larger than the nearest neighbour distance
!
! Revision 1.33  2006/01/28 16:29:51  saw44
! Fixed BFS_Step: max_length=0 was causing call to Table_Free
!
! Revision 1.32  2006/01/28 15:53:20  saw44
! Fixed CalcConnectFast - was not using NNEIGHTOL
!
! Revision 1.31  2006/01/27 18:28:19  gc121
! fixed BFS_step, added diamond() and supercell()
!
! Revision 1.30  2006/01/27 17:33:47  saw44
! Added output lines to Create_Cluster for debugging
!
! Revision 1.29  2006/01/27 15:19:24  saw44
! Create_Cluster finished - untested
!
! Revision 1.28  2006/01/26 18:18:39  gc121
! added debugging prints
!
! Revision 1.27  2006/01/26 17:55:46  saw44
! Added bond length functions
!
! Revision 1.25  2006/01/26 11:34:31  saw44
! Fixed conflict
!
! Revision 1.24  2006/01/26 01:52:21  gc121
! changed continue statements to cycle. corrected bug: missing NNEIGHTOL factor in 
! connection calc. minor typos
!
! Revision 1.23  2006/01/25 17:35:27  gc121
! changed vector_abs and vector_abs2 to norm() and norm2()
!
! Revision 1.22  2006/01/25 16:07:20  gc121
! enforced intent(in)::this for distance() calcs. corrected some typos.
!
! Revision 1.21  2006/01/25 15:10:28  gc121
! added some query functions, and made comments about routines look nice
!
! Revision 1.20  2006/01/25 14:55:30  saw44
! Periodic commit
!
! Revision 1.19  2006/01/25 11:18:33  saw44
! Moved array searching/sorting from Atoms to linearalgebra
!
! Revision 1.18  2006/01/24 11:49:08  saw44
! Atoms now compiles with the recent changes
!
! Revision 1.17  2006/01/24 11:38:01  saw44
! Periodic commit
!
! Revision 1.16  2006/01/18 11:02:19  saw44
! twobody -> neighbour
!
! Revision 1.15  2006/01/17 13:43:21  gc121
! Removed default_cutoff: user should think about it. Renamed twobody
! to neighbour. Added logical: use_uniform_cutoff. Removed
! default_number_of_neighbours. Removed atoms%lattice_set. Removed
! interface distance8 from atoms. Renamed dist8 to diance8.
! Generalized Atoms_Zero to optionally take any integer array.
! Removed AtomsSetLattice.  
!
! Revision 1.14  2005/12/21 16:07:15  saw44
! Added and tested Atoms_Print_ext_xyz, fixed CalcConnectFast
!
! Revision 1.13  2005/12/16 15:46:55  saw44
! Fixed teleporting atoms bug
!
! Revision 1.12  2005/12/14 17:15:45  saw44
! Added binary file read/write with Binoutput type (until Inoutput is changed).
!
! Revision 1.11  2005/12/05 16:59:55  saw44
! Minor changes to CalcConnectFast, Use*Mask
!
! Revision 1.10  2005/11/30 14:24:43  saw44
! Debugged and Tested. Compiles with ifort, some problems with ifort9 occasionally
!
! Revision 1.9  2005/11/25 11:22:06  saw44
! Debugged. Compiles under ifort. As yet untested
!
! Revision 1.8  2005/11/23 12:23:22  saw44
! Added Connectivity procedures
!
! Revision 1.7  2005/11/17 17:50:44  saw44
! Distance8, Distance27 and Map_Into_Cell added
!
! Revision 1.6  2005/11/17 11:31:31  saw44
! Added Connection type and begun connectivity procedures
!
! Revision 1.5  2005/11/15 14:56:48  saw44
! real(8) -> real(dp) and errors -> System_abort
!
! Revision 1.3  2005/11/11 13:43:09  saw44
! Added DynamicalSystem.f95. Added AddAtoms and RemoveAtom to Atoms.f95
!
! Revision 1.2  2005/11/11 11:18:59  gc121
! forgot to comment cvs magic variables
!
! Revision 1.1.1.1  2005/11/11 10:22:24  gc121
! starting


module  atoms_module

  use system_module
  use extendable_str_module
  use table_module
  use dictionary_module
  use linearalgebra_module
  use periodictable_module
  use minimization_module

  implicit none

  real(dp), parameter :: DEFAULT_NNEIGHTOL = 1.2_dp    !% Default value for 'atoms%nneightol'

  integer,  parameter :: NOT_NEIGHBOUR = 0     !% Returned by 'Find_Neighbour' if unsuccessful

  logical :: printed_stack_warning = .false.   !% OMIT

  type Table_pointer
    type(Table), pointer :: t => null()
  end type Table_pointer

  type Connection

     !% We do not use a minimum image convention, rather, collect all the images of a neigbouring atom
     !% that fall withing the neighbour cutoff. The different images are made distinct in the connection list
     !% by having different 'shift' vectors associated with them. 
     !%
     !% To save storage, the 'neighbour1' table contains all information about the connection
     !% but is only filled in for $i <= j$. 'neighbour2' is just a list of those of $i$'s neighbours 
     !% with $i > j$ together with an index into the 'neighbour1' table of atom $j$.
     !%
     !% In normal use (i.e. outside this module) you don\'t need direct access to the tables 
     !% so you should use the interface functions 'atoms_n_neighbours' and 'atoms_neighbour' which 
     !% hide the distiction between the cases $i <= j$ and $i > j$.
     !%
     !% \begin{tabular}{|c|c|c|c|c|}
     !% \hline
     !% \multicolumn{4}{|c|}{'neighbour1(i)%int'} & 'neighbour1(i)%real' \\
     !% \hline
     !% 1 & 2 & 3 & 4 & 1 \\
     !% \hline
     !% $j$ & 'shift_a' & 'shift_b' & 'shift_c' & $r_{ij}$ \\
     !% \hline
     !% \end{tabular}
     !%
     !% \begin{tabular}{|c|c|}
     !% \hline
     !% \multicolumn{2}{|c|}{'neighbour2(i)%int}' \\
     !% \hline
     !% 1 & 2 \\
     !% \hline
     !% $j$ & $n$ \\
     !% \hline
     !% \end{tabular}
     !%
     !% N.B. If $i$ and $j$ are neighbours with shift 'shift',
     !% then 'norm(atoms%pos(j) - atoms%pos(i) + shift)'
     !% is a minimum. 
     !% Mnemonic: 'shift' is added to $j$ to get closer to $i$.

     ! neighbour1 (i): i <= j for all j in table
     !
     ! int: Nint=4                            real: Nreal=1
     ! -----------------------------------    --------
     ! | 1 |    2    |    3    |    4    |    |   1  |
     ! -----------------------------------    --------
     ! | j | shift_a | shift_b | shift_c |    | r_ij |
     ! -----------------------------------    --------
     !
     !
     ! neighbour2 (i): i > j for all j in table
     !
     ! int: Nint=1   real: Nreal=0
     ! ---------              --
     ! | 1 | 2 |              ||
     ! ---------              --
     ! | j | n |              ||
     ! ---------              --
     !
     !
     ! cell(i,j,k)
     !
     ! int: Nint=1   real: Nreal=0
     ! --------             --
     ! |   1  |             ||
     ! --------             --
     ! | atom |             ||
     ! --------             --

     logical                                    :: initialised = .false.
     logical                                    :: cells_initialised = .false.
     logical                                    :: too_few_cells_warning_issued = .false.

     integer                                    :: cellsNa !% no. of cells in the lattice directions
     integer                                    :: cellsNb, cellsNc

     type(table_pointer), allocatable, dimension(:)     :: neighbour1 !% Neighbour information for pairs $i <= j$. 
                                                              !% Contains full details of $j$, $r_{ij}$ and shift.
     type(table_pointer), allocatable, dimension(:)     :: neighbour2 !% Neighbour information for pairs $i > j$.
                                                              !% Simply contains $j$ and a reference to $j$'s
                                                              !% 'neighbour1' table.

     type(table), allocatable, dimension(:,:,:) :: cell    !% For the linear scaling connection calculator

  end type Connection


  type Atoms

     logical                               :: initialised = .false.
     integer                               :: N = 0 !% The number of atoms held

     logical                               :: use_uniform_cutoff = .false. !% Rather than covalent radii --- 
                                                                           !% default is variable cutoff.
     real(dp)                              :: cutoff = DEFAULT_NNEIGHTOL, cutoff_break = DEFAULT_NNEIGHTOL
     !% if 'use_uniform_cutoff' is true, cutoff is the cutoff distance in \AA{}.
     !% Otherwise, cutoff is a multiplier for 'bond_length(Zi,Zj)'.

     real(dp)                              :: nneightol = DEFAULT_NNEIGHTOL 
                                              !% Count as nearest neighbour if sum of covalent radii
                                              !% times 'this%nneightol' greater than distance between atoms.
                                              !% Used in cluster carving.

     real(dp),              dimension(3,3) :: lattice    !% Lattice vectors, as columns:
     !%\begin{displaymath}
     !%\left(
     !%\begin{array}{ccc}
     !% | & | & | \\
     !% \mathbf{a} & \mathbf{b} & \mathbf{c} \\
     !% | & | & | \\
     !%\end{array}
     !%\right)
     !% = \left(
     !%\begin{array}{ccc}
     !% R_{11} & R_{12} & R_{13} \\
     !% R_{21} & R_{22} & R_{23} \\
     !% R_{31} & R_{32} & R_{33} \\
     !%\end{array}
     !%\right)
     !%\end{displaymath}
     !% i.e. $\mathbf{a} = $ 'lattice(:,1)', $\mathbf{b} = $ 'lattice(:,2)' and
     !% $\mathbf{c} = $ 'lattice(:,3)'.

     !  ( | | | | | | )          ( (1,1) (1,2) (1,3) ) 
     ! (  | | | | | |  )        (                     )
     ! (  |a| |b| |c|  )    =   (  (2,1) (2,2) (2,3)  )
     ! (  | | | | | |  )        (                     )
     !  ( | | | | | | )          ( (3,1) (3,2) (3,3) ) 


     real(dp),              dimension(3,3) :: g          !% Inverse lattice (stored for speed)

     type(Table) :: data !% This table contains all atomic data. 
                         !% Rows correspond to atoms, columns to the various properties.

     type(Dictionary) :: properties !% Mapping of names to column indices for the 'data' table, e.g.
                                    !%>result = get_value(properties,"pos",pos_indices)
                                    !% gives 'pos_indices = (/PROPERTY_REAL,1,3/)', indicating that the 
                                    !% position data consists of real number, beginning at column 1 
                                    !% and ending at column 3.

     type(Dictionary) :: params     !% List of key/value pairs read from comment line of XYZ file


     integer,  pointer, dimension(:)   :: Z      => null()  !% Atomic numbers, dimension is actually $(N)$
     character(TABLE_STRING_LENGTH), pointer, dimension(:) :: species => null() !% Names of elements

     integer, pointer, dimension(:)    :: move_mask => null()  !% Atoms with 'move_mask' set to false are fixed
     integer, pointer, dimension(:)    :: damp_mask => null()  !% Damping is only applied to those atoms with
                                                              !% 'damp_mask' set to 1. 
                                                              !% By default this is set to 1 for all atoms.

     integer, pointer, dimension(:)    :: thermostat_region => null() !% Which thermostat is applied to each atoms. 
                                                              !% By default this is set to 1 for all atoms.

     integer,  pointer, dimension(:,:) :: travel => null()  !% Travel across periodic conditions. Actually $(3,N)$ array.
                                                            !% See 'map_into_cell' below.
     real(dp), pointer, dimension(:,:) :: pos    => null()  !% $(3,N)$ array of atomic positions, in \AA. 
                                                            !% Position of atom $i$ is 'pos(:,i)'
     real(dp), pointer, dimension(:)   :: mass   => null()  !% Atomic masses, dimension is actually $(N)$

     real(dp), pointer, dimension(:,:) :: velo   => null()  !% $(3,N)$ array  of atomic velocities, in \AA/fs.
     real(dp), pointer, dimension(:,:) :: acc    => null()  !% $(3,N)$ array  of accelerations in \AA/fs$^2$
     real(dp), pointer, dimension(:,:) :: avgpos => null()  !% $(3,N)$ array  of time-averaged atomic positions.
     real(dp), pointer, dimension(:,:) :: oldpos => null()  !% $(3,N)$ array  of positions of atoms at previous time step.
     real(dp), pointer, dimension(:) :: avg_ke => null()    !% Time-averaged atomic kinetic energy
     type(Connection)                      :: connect       !% Connectivity object (see above)
     type(Connection)                      :: hysteretic_connect       !% Hysteretic connectivity object (see above)

  end type Atoms

  type atoms_ll
    type (atoms_ll_entry), pointer :: first => null()
    type (atoms_ll_entry), pointer :: last => null()
  end type atoms_ll

  type atoms_ll_entry
    type(atoms) :: at
    type(atoms_ll_entry), pointer :: next => null()
    type(atoms_ll_entry), pointer :: prev => null()
  end type atoms_ll_entry


  interface initialise
     module procedure atoms_initialise, connection_initialise
  end interface initialise

  !% Free up the memory associated with one or more objects.
  interface finalise
     module procedure atoms_finalise,atoms_finalise_multi, connection_finalise
  end interface finalise

  interface wipe
     module procedure connection_wipe
  end interface wipe

  interface zero
     module procedure atoms_zero
  end interface

  !% Overloaded assigment operators for Atoms and Connection objects.
  interface assignment(=)
     module procedure atoms_assignment, connection_assignment
  end interface assignment(=)

  interface set_atoms
     module procedure atoms_set_atoms, atoms_set_atoms_singlez
  end interface

  !% increase cutoff
  interface set_cutoff_minimum
     module procedure atoms_set_cutoff_minimum
  end interface

  !% set (a uniform) cutoff
  interface set_cutoff
     module procedure atoms_set_cutoff
  end interface

  !% set cutoff factor
  interface set_cutoff_factor
     module procedure atoms_set_cutoff_factor
  end interface

  !% Add one or more atoms to an Atoms object.
  interface add_atoms
     module procedure add_atom_single, add_atom_multiple
  end interface add_atoms

  !% Remove one or more atoms from an Atoms object.
  interface remove_atoms
     module procedure remove_atom_single, remove_atom_multiple
  end interface remove_atoms

  !% Add an extra property to this atoms object, as extra columns of 
  !% integers or real numbers in the 'data' table. For example, this interface
  !% is used by the DynamicalSystems module to create the 'velo', 'acc', etc. properties
  !% Optionally, the type and indices of the new property are returned
  !% in the array 'lookup'.
  interface add_property
     module procedure atoms_add_property_int, atoms_add_property_int_a
     module procedure atoms_add_property_real, atoms_add_property_real_a
     module procedure atoms_add_property_str, atoms_add_property_str_a
     module procedure atoms_add_property_logical, atoms_add_property_logical_a 
  end interface

  !% Convenience function to test if a property is present. No checking
  !% of property type is done.
  interface has_property
     module procedure atoms_has_property
  end interface

  !% remove a property from this atoms object
  interface remove_property
    module procedure atoms_remove_property
  end interface

  !% This is a convenience interface to assign pointers to custom properties of this
  !% Atoms object. The pointer is simply directed to the relevant columns of the
  !% 'this%data' table. Adding new properties will invalidate these pointers
  !% because the data table gets moved in memory whenever a new column is added, so
  !% this interface should be recalled whenever properties are added.
  !% Returns true if successful, false if property doesn't exist, and aborts
  !% if property type and type of pointer don't match.
  interface assign_pointer
     module procedure atoms_assign_pointer_int1D, atoms_assign_pointer_int2D
     module procedure atoms_assign_pointer_real1D, atoms_assign_pointer_real2D
     module procedure atoms_assign_pointer_str1D, atoms_assign_pointer_str2D
     module procedure atoms_assign_pointer_logical1D, atoms_assign_pointer_logical2D
  end interface

  !% This interface calculates the distance between the nearest periodic images of two points (or atoms).
  interface distance_min_image
     module procedure distance8_atom_atom, distance8_atom_vec, distance8_vec_atom, distance8_vec_vec
  end interface

  !% This interface calculates the difference vector between the nearest periodic images of two points (or atoms).
  interface diff_min_image
     module procedure diff_atom_atom, diff_atom_vec, diff_vec_atom, diff_vec_vec
  end interface

  !% Write a binary representation of an Atoms or Connection object to a file
  interface write_binary
     module procedure Atoms_File_Write, Connection_File_Write, Table_Pointer_write_binary_a
  end interface

  !% Read an Atoms or Connection object from a binary file
  interface read_binary
     module procedure Atoms_File_Read, Connection_File_Read, Table_Pointer_read_binary_a
  end interface

  !% Print a verbose textual description of an Atoms or Connection object to the default logger or to
  !% a specificied Inoutput object.
  interface print
     module procedure atoms_print, connection_print
  end interface print

  !% Print this Atoms object to an Inouput object or to a file in XYZ format. By default, only
  !% the atomic species and positions are written, but the 'properties' argument
  !% allows you to specify additional properties to be printed,
  !% e.g. 'properties=(/"pos","velo","acc"/)' would print the positions, velocities
  !% and accelerations (9 columns of real numbers). The 'all_properties' option
  !% prints all properties associated with this Atoms object. 
  !% The comment line (line 2) of the XYZ file contains the contents
  !% of the 'params' dictionary which at a minimum includes the lattice
  !% and list of properties printed, but can be used to print other simulation
  !% parameters.
  !%
  !% The lattice is printed in the form:
  !%>    Lattice="R11 R21 R31 R12 R22 R32 R13 R23 R33"
  !% and the list of properties in the file is added in the form:
  !%>    Properties="species:S:1:pos:R:3:velo:R:3:acc:R:3"
  !% with a triplet of colon-separated fields for each property, giving
  !% its name, type ('R' for real,'I' for integer, 'S' for string and 'L' for logical) 
  !% and number of columns.
  !% Other parameters are output in 'key=value' format ---
  !% all types of dictionary entries with the exception, for now, of
  !% complex numbers are supported.
  interface print_xyz
     module procedure atoms_print_xyz, atoms_print_xyz_filename, atoms_ll_print_xyz
  end interface print_xyz

  !% Print in AtomEye extended CFG format. Arguments are as for 'print_xyz' above.
  interface print_cfg
     module procedure atoms_print_cfg, atoms_print_cfg_filename
  end interface

  !% Read atoms from a standard XYZ file into this Atoms object. Properties for
  !% all the fields found in the file are created, and parameters in 'key=value'
  !% form are read from the comment line (line 2) into the 'params'
  !% dictionary within the Atoms object. 
  !%
  !% An overloaded routine skips over a single frame in the input stream.
  interface read_xyz
     module procedure atoms_read_xyz_inoutput, atoms_read_xyz_extendable_str
     module procedure atoms_read_xyz_filename, atoms_read_xyz_skip
  end interface read_xyz

  interface set_lattice
     module procedure atoms_set_lattice
  end interface set_lattice

  !% Select a subset of the atoms in an atoms object, either using a logical 
  !% mask array or a Table, in which case the first 'int' column is taken to 
  !% be a list of the atoms that should be retained.
  interface select
     module procedure atoms_select
  end interface select

  !% calculate volume of unit cell
  interface cell_volume
    module procedure atoms_cell_volume
    module procedure lattice_cell_volume
  end interface

  interface map_into_cell
    module procedure atoms_map_into_cell
    module procedure vec_map_into_cell
  end interface

  !% add new entry to atoms linked list structure
  interface new_entry
     module procedure atoms_ll_new_entry
  end interface new_entry

  !% remove last entry from atoms linked list structure
  interface remove_last_entry
     module procedure atoms_ll_remove_last_entry
  end interface remove_last_entry


contains

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  ! Initialisation and Finalisation
  !
  ! N : number of atoms
  ! lattice: 3x3 matrix of lattice vectors as column vectors
  ! data: optional initial values for data table
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  !% Initialise an Atoms object to store 'N' atoms and specify the initial lattice.
  subroutine atoms_initialise(this,N,lattice,data,properties,params)

    type(Atoms),                intent(inout) :: this
    integer,                    intent(in)    :: N
    real(dp), dimension(3,3),   intent(in)    :: lattice
    type(Table), optional,      intent(in)    :: data
    type(Dictionary), optional, intent(in)    :: properties, params

    integer :: i, stack_size, stack_size_err
    integer, dimension(:), allocatable :: default_int
    real(dp), dimension(:), allocatable :: default_real
    character(TABLE_STRING_LENGTH), dimension(:), allocatable :: default_str


    call atoms_finalise(this)

    if ((present(data) .and. .not. present(properties)) .or. &
         (present(properties) .and. .not. present(data))) &
         call system_abort('Atoms_Initialise: both data and properties must be present')

    if (present(data)) then
       if (data%N /= N) call system_abort('Atoms_Initialise: data%N ('//data%N//') /= N('//N//')')
       if (data%N == 0) call system_abort('Atoms_Initialise: cannot initialise from data with zero length')
       this%data = data
       this%properties = properties
    else
       ! By default, we just add properties for Z, species name and position
       call initialise(this%properties)

       call add_property(this, 'Z', 0)
       call add_property(this, 'pos', 0.0_dp, n_cols=3)
       call add_property(this, 'species', repeat(' ', TABLE_STRING_LENGTH))

       allocate(default_int(this%data%intsize))
       allocate(default_real(this%data%realsize))
       allocate(default_str(this%data%strsize))

       default_int = 0
       default_real = 0.0_dp
       default_str = repeat(' ',TABLE_STRING_LENGTH)

       ! Add N rows to data table
       do i=1,N
          call append(this%data, default_int, default_real, default_str)
       end do

       deallocate(default_int, default_real)
    end if

    if (present(params)) then
       this%params = params
    else
       call initialise(this%params)
    end if

    this%N = N

    call atoms_repoint(this)

    stack_size = 3*N*8/1024
    if (stack_size > 4000 .and. .not. printed_stack_warning) then
       call print('Atoms_Initialise: Stack size must be at least '//stack_size//&
            ' kb, or we get seg faults when functions return arrays.')
       printed_stack_warning = .true.
       call print('Atoms_Initialise: trying to increase stack limit')
       stack_size = int(stack_size/1024)*1024
       stack_size_err = increase_stack(stack_size)
       if (stack_size_err /= 0) then
	call print("Atoms_Initialise: error calling increase_stack err = "// stack_size_err)
       endif
    end if

    this%lattice = lattice
    call matrix3x3_inverse(this%lattice,this%g)

    this%initialised = .true.
  end subroutine atoms_initialise

  !% OMIT
  ! Initialise pointers for convenient access to special columns of this%data
  subroutine atoms_repoint(this)
    type(Atoms), target, intent(inout) :: this
    integer :: i, lookup(3)
    character(len=key_len) :: key

    if(this%N == 0) return

    nullify(this%Z, this%travel, this%pos, this%move_mask, this%damp_mask, &
         this%thermostat_region, this%pos, this%velo, this%acc, this%avgpos, &
         this%oldpos)

    ! Loop over all properties looking for those with special names
    ! which we have pointers for
    do i = 1,this%properties%N

       key = this%properties%keys(i)
       if (.not. get_value(this%properties, key, lookup)) &
            call system_abort('Atoms_repoint: key '//trim(key)//' not found.')

       ! If this%N is zero then point at zero length arrays
       if (this%N == 0) then

          select case(trim(key))

          ! Integer properties
          case('Z')
             this%Z               => this%data%int(lookup(2),:)
          case('travel')
             this%travel          => this%data%int(lookup(2):lookup(3),:)
          case('move_mask')
             this%move_mask       => this%data%int(lookup(2),:)
          case('damp_mask')
             this%damp_mask       => this%data%int(lookup(2),:)
          case('thermostat_region')
             this%thermostat_region => this%data%int(lookup(2),:)

          ! Real properties
          case('mass')
             this%mass            => this%data%real(lookup(2),:)
          case('pos')
             this%pos             => this%data%real(lookup(2):lookup(3),:)
          case('velo')
             this%velo            => this%data%real(lookup(2):lookup(3),:)
          case('acc')
             this%acc             => this%data%real(lookup(2):lookup(3),:)
          case('avgpos')
             this%avgpos          => this%data%real(lookup(2):lookup(3),:)
          case('oldpos')
             this%oldpos          => this%data%real(lookup(2):lookup(3),:)
          case('avg_ke')
             this%avg_ke          => this%data%real(lookup(2),:)

          ! String properties
          case('species')
             this%species         => this%data%str(lookup(2),:)

          end select


       else

          select case(trim(key))

          ! Integer properties
          case('Z')
             this%Z               => this%data%int(lookup(2),1:this%N)
          case('travel')
             this%travel          => this%data%int(lookup(2):lookup(3),1:this%N)
          case('move_mask')
             this%move_mask       => this%data%int(lookup(2),1:this%N)
          case('damp_mask')
             this%damp_mask       => this%data%int(lookup(2),1:this%N)
          case('thermostat_region')
             this%thermostat_region => this%data%int(lookup(2),1:this%N)

          ! Real properties
          case('mass')
             this%mass            => this%data%real(lookup(2),1:this%N)   
          case('pos')
             this%pos             => this%data%real(lookup(2):lookup(3),1:this%N)
          case('velo')
             this%velo            => this%data%real(lookup(2):lookup(3),1:this%N)
          case('acc')
             this%acc             => this%data%real(lookup(2):lookup(3),1:this%N)
          case('avgpos')
             this%avgpos          => this%data%real(lookup(2):lookup(3),1:this%N)
          case('oldpos')
             this%oldpos          => this%data%real(lookup(2):lookup(3),1:this%N)
          case('avg_ke')
             this%avg_ke          => this%data%real(lookup(2),1:this%N)

          ! String properties
          case('species')
             this%species         => this%data%str(lookup(2),1:this%N)


          end select

       end if

    end do

  end subroutine atoms_repoint


  subroutine atoms_finalise(this)
    type(Atoms), intent(inout) :: this

    call finalise(this%data)
    call finalise(this%properties)
    call finalise(this%params)

    ! Nullify pointers
    nullify(this%Z, this%travel, this%mass)
    nullify(this%move_mask, this%thermostat_region, this%damp_mask)
    nullify(this%pos, this%velo, this%acc, this%avgpos, this%oldpos, this%avg_ke)

    call connection_finalise(this%connect)
    call connection_finalise(this%hysteretic_connect)

    this%N = 0

    this%initialised = .false.

  end subroutine atoms_finalise

  !Quick multiple finalisations
  subroutine atoms_finalise_multi(at1,at2,at3,at4,at5,at6,at7,at8,at9,at10)
    type(Atoms),           intent(inout) :: at1,at2
    type(Atoms), optional, intent(inout) :: at3,at4,at5,at6,at7,at8,at9,at10
    call atoms_finalise(at1)
    call atoms_finalise(at2)
    if (present(at3)) call atoms_finalise(at3)
    if (present(at4)) call atoms_finalise(at4)
    if (present(at5)) call atoms_finalise(at5)
    if (present(at6)) call atoms_finalise(at6)
    if (present(at7)) call atoms_finalise(at7)
    if (present(at8)) call atoms_finalise(at8)
    if (present(at9)) call atoms_finalise(at9)
    if (present(at10)) call atoms_finalise(at10)
  end subroutine atoms_finalise_multi

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  ! Overloaded Assignment
  !
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine atoms_assignment(to,from)

    type(Atoms), intent(inout) :: to
    type(Atoms), intent(in)    :: from

    if(.not. from%initialised) &
         call system_abort("Atoms_assignment: 'from' object is not initialised")

    call atoms_initialise(to,from%N,from%lattice, from%data, from%properties, from%params)

    to%use_uniform_cutoff = from%use_uniform_cutoff
    to%cutoff      = from%cutoff
    to%cutoff_break      = from%cutoff_break
    to%nneightol   = from%nneightol

    to%connect     = from%connect
    to%hysteretic_connect     = from%hysteretic_connect

  end subroutine atoms_assignment


  !% Make a copy of the atoms object 'from' without including
  !% connectivity information. Useful for saving the state of a
  !% dynamical simulation without incurring too great a memory
  !% cost. 
  subroutine atoms_copy_without_connect(to, from, properties)

    type(Atoms), intent(inout) :: to
    type(Atoms), intent(in)    :: from
    character(len=*), optional, intent(in) :: properties

    integer :: i, n_properties, n_cols, col, lookup(3)
    character(len=1024) :: tmp_properties(from%properties%N)
    integer :: intcols(from%data%intsize), realcols(from%data%realsize), strcols(from%data%strsize), logicalcols(from%data%logicalsize)
    integer :: last_int_col, last_real_col, last_str_col, last_logical_col

    if(.not. from%initialised) &
         call system_abort("atoms_copy_without_connect: 'from' object is not initialised")

    to%N = from%N
    to%lattice = from%lattice
    to%g = from%g

    if (present(properties)) then
      call finalise(to%data)
      call initialise(to%data, max_length=from%data%N)
      last_int_col=0
      last_real_col=0
      last_str_col=0
      last_logical_col=0
      intcols = -1
      realcols = -1
      strcols = -1
      logicalcols = -1
      call parse_string(properties, ':', tmp_properties, n_properties)
      do i=1, n_properties
	if (.not. get_value(from%properties, tmp_properties(i), lookup)) &
	  call system_abort('ERROR: atoms_copy_without_connect: copying key '//trim(tmp_properties(i))//' not found.')
	n_cols = lookup(3)-lookup(2)+1
	if (lookup(1) == PROPERTY_INT) then
	  do col=1, n_cols
	    intcols(last_int_col+col) = lookup(2)+col-1
	  end do
	  call set_value(to%properties, tmp_properties(i), (/ lookup(1), last_int_col+1, last_int_col+n_cols/) )
	  last_int_col = last_int_col + n_cols
	else if (lookup(1) == PROPERTY_REAL) then
	  do col=1, n_cols
	    realcols(last_real_col+col) = lookup(2)+col-1
	  end do
	  call set_value(to%properties, tmp_properties(i), (/ lookup(1), last_real_col+1, last_real_col+n_cols/) )
	  last_real_col = last_real_col + n_cols
	else if (lookup(1) == PROPERTY_STR) then
	  do col=1, n_cols
	    strcols(last_str_col+col) = lookup(2)+col-1
	  end do
	  call set_value(to%properties, tmp_properties(i), (/ lookup(1), last_str_col+1, last_str_col+n_cols/) )
	  last_str_col = last_str_col + n_cols
	else if (lookup(1) == PROPERTY_LOGICAL) then
	  do col=1, n_cols
	    logicalcols(last_logical_col+col) = lookup(2)+col-1
	  end do
	  call set_value(to%properties, tmp_properties(i), (/ lookup(1), last_logical_col+1, last_logical_col+n_cols/) )
	  last_logical_col = last_logical_col + n_cols
	else
	  call system_abort("ERROR: atoms_copy_without_connect got unknown property type " // lookup(1) // " for property " // trim(tmp_properties(i)))
	end if
      end do
      to%data = subtable(from%data, (/ (i, i=1,from%data%N) /), intcols, realcols, strcols, logicalcols)
    else
      to%data = from%data
      to%properties = from%properties
    endif
    to%params = from%params

    to%use_uniform_cutoff = from%use_uniform_cutoff
    to%cutoff      = from%cutoff
    to%cutoff_break      = from%cutoff_break
    to%nneightol   = from%nneightol

    call atoms_repoint(to)
    to%initialised = .true.

  end subroutine atoms_copy_without_connect

  subroutine atoms_select(to, from, mask, list)
    type(Atoms), intent(inout) :: to
    type(Atoms), intent(in) :: from
    logical, intent(in), optional :: mask(:)
    integer, intent(in), optional :: list(:)


    if ((.not. present(list) .and. .not. present(mask)) .or. (present(list) .and. present(mask))) &
         call system_abort('atoms_select: either list or mask must be present (but not both)')


    if (present(mask)) then
       if (size(mask) /= from%N) &
            call system_abort("atoms_select: mismatched sizes of from " // from%N // " and mask " // size(mask))
       call atoms_initialise(to, count(mask), from%lattice)
    else
       call atoms_initialise(to, size(list), from%lattice)
    end if

    call finalise(to%properties)
    to%properties = from%properties
    call finalise(to%params)
    to%params = from%params
    to%use_uniform_cutoff = from%use_uniform_cutoff
    to%cutoff = from%cutoff
    to%cutoff_break = from%cutoff_break
    to%nneightol = from%nneightol

    if(present(mask)) then
       call select(to%data, from%data, row_mask=mask)
    else
       call select(to%data, from%data, row_list=list)
    end if

    call atoms_repoint(to)
  end subroutine atoms_select

  subroutine atoms_remove_property(this, name)
    type(Atoms), intent(inout) :: this
    character(len=*), intent(in) :: name

    integer :: i, lookup(3), rem_type, rem_start, n_removed
    logical :: dummy

    if (get_value(this%properties, name, lookup)) then
       if (lookup(1) == PROPERTY_INT) then
          call remove_columns(this%data, int_col_min = lookup(2), int_col_max = lookup(3))
       else if (lookup(1) == PROPERTY_REAL) then
          call remove_columns(this%data, real_col_min = lookup(2), real_col_max = lookup(3))
       else if (lookup(1) == PROPERTY_STR) then
          call remove_columns(this%data, str_col_min = lookup(2), str_col_max = lookup(3))
       else if (lookup(1) == PROPERTY_LOGICAL) then
          call remove_columns(this%data, logical_col_min = lookup(2), logical_col_max = lookup(3))
       else
          call system_abort("atoms_remove_property confused by property type " // lookup(1))
      endif
    else
      call system_abort("atoms_remove_property tried to remove non-existent property " // trim(name))
    endif

    call remove_value(this%properties, name)
    rem_type = lookup(1)
    rem_start = lookup(2)
    n_removed = lookup(3)-lookup(2)+1
    
    ! correct columns in subsequent properties
    do i=1,this%properties%N
       dummy = get_value(this%properties, this%properties%keys(i), lookup)
       if (lookup(1) /= rem_type) cycle
       if (lookup(2) > rem_start) then
          lookup(2) = lookup(2) - n_removed
          lookup(3) = lookup(3) - n_removed
          call set_value(this%properties, this%properties%keys(i), lookup)
       end if
    end do

    ! remove_columns will have moved this%data in memory and invalidated pointers
    call atoms_repoint(this)

  end subroutine atoms_remove_property

  function atoms_has_property(this, name)
    type(Atoms), intent(in) :: this
    character(len=*), intent(in) :: name
    logical :: atoms_has_property

    integer :: lookup(3)

    atoms_has_property = get_value(this%properties, name, lookup)
    
  end function atoms_has_property

  subroutine atoms_add_property_int(this, name, value, n_cols, lookup)
    type(Atoms), intent(inout), target :: this
    character(len=*), intent(in) :: name
    integer, intent(in) :: value
    integer, intent(in), optional :: n_cols
    integer, dimension(3), optional, intent(out) :: lookup

    integer :: use_n_cols, new_cols(2), use_lookup(3)

    use_n_cols = 1
    if (present(n_cols)) use_n_cols = n_cols

    ! Check if there's already a property with this name
    if (get_value(this%properties, name, use_lookup)) then
       ! Does it match this type and number of columns?
       if (use_lookup(1) == PROPERTY_INT .and. use_lookup(3)-use_lookup(2)+1 == use_n_cols) then
          call print('Add_Property: property '//trim(name)//' already present', VERBOSE)
          if (present(lookup)) lookup = use_lookup
          return
       else
          call system_abort('Add_Property: incompatible property '//name//' already present')
       end if
    end if

    call append_column(this%data, value, n_cols=use_n_cols, cols=new_cols)
    use_lookup = (/PROPERTY_INT, new_cols(1), new_cols(2)/)
    call set_value(this%properties, name, use_lookup)
    if (present(lookup)) lookup = use_lookup

    ! Append column will have moved this%data in memory and invalidated pointers
    call atoms_repoint(this)
    call print('WARNING: atoms_add_property - pointers invalidated', VERBOSE)

  end subroutine atoms_add_property_int


  subroutine atoms_add_property_int_a(this, name, value, n_cols, lookup)
    type(Atoms), intent(inout),target :: this
    character(len=*), intent(in) :: name
    integer, intent(in), dimension(this%N) :: value
    integer, intent(in), optional :: n_cols
    integer, dimension(3), optional, intent(out) :: lookup

    integer :: use_n_cols, new_cols(2), use_lookup(3)

    use_n_cols = 1
    if (present(n_cols)) use_n_cols = n_cols

    ! Check if there's already a property with this name
    if (Get_Value(this%properties, name, use_lookup)) then
       ! Does it match this type and number of columns?
       if (use_lookup(1) == PROPERTY_INT .and. use_lookup(3)-use_lookup(2)+1 == use_n_cols) then
          call print('Add_Property: property '//trim(name)//' already present', VERBOSE)
          if (present(lookup)) lookup = use_lookup
          return
       else
          call system_abort('Add_Property: incompatible property '//name//' already present')
       end if
    end if

    call append_column(this%data, value, n_cols=use_n_cols, cols=new_cols)
    use_lookup = (/PROPERTY_INT, new_cols(1), new_cols(2)/)
    call set_value(this%properties, name, use_lookup)
    if (present(lookup)) lookup = use_lookup

    ! Append column will have moved this%data in memory and invalidated pointers
    call atoms_repoint(this)
    call print('WARNING: atoms_add_property - pointers invalidated', VERBOSE)

  end subroutine atoms_add_property_int_a


  subroutine atoms_add_property_real(this, name, value, n_cols, lookup)
    type(Atoms), intent(inout),target :: this
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: value
    integer, intent(in), optional :: n_cols
    integer, dimension(3), optional, intent(out) :: lookup

    integer :: use_n_cols, new_cols(2), use_lookup(3)

    use_n_cols = 1
    if (present(n_cols)) use_n_cols = n_cols

    ! Check if there's already a property with this name
    if (Get_Value(this%properties, name, use_lookup)) then
       ! Does it match this type and number of columns?
       if (use_lookup(1) == PROPERTY_REAL .and. use_lookup(3)-use_lookup(2)+1 == use_n_cols) then
          call print('Add_Property: property '//trim(name)//' already present', VERBOSE)
          if (present(lookup)) lookup = use_lookup
          return
       else
          call system_abort('Add_Property: incompatible property '//name//' already present')
       end if
    end if

    call append_column(this%data, value, n_cols=use_n_cols, cols=new_cols)
    use_lookup = (/PROPERTY_REAL, new_cols(1), new_cols(2)/)
    call set_value(this%properties, name, use_lookup)
    if (present(lookup)) lookup = use_lookup

    ! Append column will have moved this%data in memory and invalidated pointers
    call atoms_repoint(this)
    call print('WARNING: atoms_add_property - pointers invalidated', VERBOSE)

  end subroutine atoms_add_property_real


  subroutine atoms_add_property_real_a(this, name, value, n_cols, lookup)
    type(Atoms), intent(inout),target :: this
    character(len=*), intent(in) :: name
    real(dp), intent(in), dimension(this%N) :: value
    integer, intent(in), optional :: n_cols
    integer, dimension(3), optional, intent(out) :: lookup

    integer :: use_n_cols, new_cols(2), use_lookup(3)

    use_n_cols = 1
    if (present(n_cols)) use_n_cols = n_cols

    ! Check if there's already a property with this name
    if (Get_Value(this%properties, name, use_lookup)) then
       ! Does it match this type and number of columns?
       if (use_lookup(1) == PROPERTY_REAL .and. use_lookup(3)-use_lookup(2)+1 == use_n_cols) then
          call print('Add_Property: property '//trim(name)//' already present', VERBOSE)
          if (present(lookup)) lookup = use_lookup
          return
       else
          call system_abort('Add_Property: incompatible property '//name//' already present')
       end if
    end if

    call append_column(this%data, value, n_cols=use_n_cols, cols=new_cols)
    use_lookup = (/PROPERTY_REAL, new_cols(1), new_cols(2)/)
    call set_value(this%properties, name, use_lookup)
    if (present(lookup)) lookup = use_lookup

    ! Append column will have moved this%data in memory and invalidated pointers
    call atoms_repoint(this)
    call print('WARNING: atoms_add_property - pointers invalidated', VERBOSE)

  end subroutine atoms_add_property_real_a


  subroutine atoms_add_property_str(this, name, value, n_cols, lookup)
    type(Atoms), intent(inout), target :: this
    character(len=*), intent(in) :: name
    character(TABLE_STRING_LENGTH), intent(in) :: value
    integer, intent(in), optional :: n_cols
    integer, dimension(3), optional, intent(out) :: lookup

    integer :: use_n_cols, new_cols(2), use_lookup(3)

    use_n_cols = 1
    if (present(n_cols)) use_n_cols = n_cols

    ! Check if there's already a property with this name
    if (get_value(this%properties, name, use_lookup)) then
       ! Does it match this type and number of columns?
       if (use_lookup(1) == PROPERTY_STR .and. use_lookup(3)-use_lookup(2)+1 == use_n_cols) then
          call print('Add_Property: property '//trim(name)//' already present', VERBOSE)
          if (present(lookup)) lookup = use_lookup
          return
       else
          call system_abort('Add_Property: incompatible property '//name//' already present')
       end if
    end if

    call append_column(this%data, value, n_cols=use_n_cols, cols=new_cols)
    use_lookup = (/PROPERTY_STR, new_cols(1), new_cols(2)/)
    call set_value(this%properties, name, use_lookup)
    if (present(lookup)) lookup = use_lookup

    ! Append column will have moved this%data in memory and invalidated pointers
    call atoms_repoint(this)
    call print('WARNING: atoms_add_property - pointers invalidated', VERBOSE)

  end subroutine atoms_add_property_str


  subroutine atoms_add_property_str_a(this, name, value, n_cols, lookup)
    type(Atoms), intent(inout),target :: this
    character(len=*), intent(in) :: name
    character(TABLE_STRING_LENGTH), dimension(this%N) :: value
    integer, intent(in), optional :: n_cols
    integer, dimension(3), optional, intent(out) :: lookup

    integer :: use_n_cols, new_cols(2), use_lookup(3)

    use_n_cols = 1
    if (present(n_cols)) use_n_cols = n_cols

    ! Check if there's already a property with this name
    if (Get_Value(this%properties, name, use_lookup)) then
       ! Does it match this type and number of columns?
       if (use_lookup(1) == PROPERTY_STR .and. use_lookup(3)-use_lookup(2)+1 == use_n_cols) then
          call print('Add_Property: property '//trim(name)//' already present', VERBOSE)
          if (present(lookup)) lookup = use_lookup
          return
       else
          call system_abort('Add_Property: incompatible property '//name//' already present')
       end if
    end if

    call append_column(this%data, value, n_cols=use_n_cols, cols=new_cols)
    use_lookup = (/PROPERTY_STR, new_cols(1), new_cols(2)/)
    call set_value(this%properties, name, use_lookup)
    if (present(lookup)) lookup = use_lookup

    ! Append column will have moved this%data in memory and invalidated pointers
    call atoms_repoint(this)
    call print('WARNING: atoms_add_property - pointers invalidated', VERBOSE)

  end subroutine atoms_add_property_str_a


  subroutine atoms_add_property_logical(this, name, value, n_cols, lookup)
    type(Atoms), intent(inout), target :: this
    character(len=*), intent(in) :: name
    logical, intent(in) :: value
    integer, intent(in), optional :: n_cols
    integer, dimension(3), optional, intent(out) :: lookup

    integer :: use_n_cols, new_cols(2), use_lookup(3)

    use_n_cols = 1
    if (present(n_cols)) use_n_cols = n_cols

    ! Check if there's already a property with this name
    if (get_value(this%properties, name, use_lookup)) then
       ! Does it match this type and number of columns?
       if (use_lookup(1) == PROPERTY_LOGICAL .and. use_lookup(3)-use_lookup(2)+1 == use_n_cols) then
          call print('Add_Property: property '//trim(name)//' already present', VERBOSE)
          if (present(lookup)) lookup = use_lookup
          return
       else
          call system_abort('Add_Property: incompatible property '//name//' already present')
       end if
    end if

    call append_column(this%data, value, n_cols=use_n_cols, cols=new_cols)
    use_lookup = (/PROPERTY_LOGICAL, new_cols(1), new_cols(2)/)
    call set_value(this%properties, name, use_lookup)
    if (present(lookup)) lookup = use_lookup

    ! Append column will have moved this%data in memory and invalidated pointers
    call atoms_repoint(this)
    call print('WARNING: atoms_add_property - pointers invalidated', VERBOSE)

  end subroutine atoms_add_property_logical


  subroutine atoms_add_property_logical_a(this, name, value, n_cols, lookup)
    type(Atoms), intent(inout),target :: this
    character(len=*), intent(in) :: name
    logical, dimension(this%N) :: value
    integer, intent(in), optional :: n_cols
    integer, dimension(3), optional, intent(out) :: lookup

    integer :: use_n_cols, new_cols(2), use_lookup(3)

    use_n_cols = 1
    if (present(n_cols)) use_n_cols = n_cols

    ! Check if there's already a property with this name
    if (Get_Value(this%properties, name, use_lookup)) then
       ! Does it match this type and number of columns?
       if (use_lookup(1) == PROPERTY_LOGICAL .and. use_lookup(3)-use_lookup(2)+1 == use_n_cols) then
          call print('Add_Property: property '//trim(name)//' already present', VERBOSE)
          if (present(lookup)) lookup = use_lookup
          return
       else
          call system_abort('Add_Property: incompatible property '//name//' already present')
       end if
    end if

    call append_column(this%data, value, n_cols=use_n_cols, cols=new_cols)
    use_lookup = (/PROPERTY_LOGICAL, new_cols(1), new_cols(2)/)
    call set_value(this%properties, name, use_lookup)
    if (present(lookup)) lookup = use_lookup

    ! Append column will have moved this%data in memory and invalidated pointers
    call atoms_repoint(this)
    call print('WARNING: atoms_add_property - pointers invalidated', VERBOSE)

  end subroutine atoms_add_property_logical_a


  function atoms_assign_pointer_int1D(this, name, ptr)
    type(Atoms), intent(in) :: this
    character(len=*), intent(in) :: name
    integer, pointer :: ptr(:)
    logical :: atoms_assign_pointer_int1D

    atoms_assign_pointer_int1D = atoms_assign_pointer(this, name, int1D_ptr=ptr)
  end function atoms_assign_pointer_int1D

  function atoms_assign_pointer_int2D(this, name, ptr)
    type(Atoms), intent(in) :: this
    character(len=*), intent(in) :: name
    integer, pointer :: ptr(:,:)
    logical :: atoms_assign_pointer_int2D

    atoms_assign_pointer_int2D = atoms_assign_pointer(this, name, int2D_ptr=ptr)
  end function atoms_assign_pointer_int2D

  function atoms_assign_pointer_real1D(this, name, ptr)
    type(Atoms), intent(in) :: this
    character(len=*), intent(in) :: name
    real(dp), pointer :: ptr(:)
    logical :: atoms_assign_pointer_real1D

    atoms_assign_pointer_real1D = atoms_assign_pointer(this, name, real1D_ptr=ptr)
  end function atoms_assign_pointer_real1D

  function atoms_assign_pointer_real2D(this, name, ptr)
    type(Atoms), intent(in) :: this
    character(len=*), intent(in) :: name
    real(dp), pointer :: ptr(:,:)
    logical :: atoms_assign_pointer_real2D

    atoms_assign_pointer_real2D = atoms_assign_pointer(this, name, real2D_ptr=ptr)
  end function atoms_assign_pointer_real2D

  function atoms_assign_pointer_str1D(this, name, ptr)
    type(Atoms), intent(in) :: this
    character(len=*), intent(in) :: name
    character(TABLE_STRING_LENGTH), pointer :: ptr(:)
    logical :: atoms_assign_pointer_str1D

    atoms_assign_pointer_str1D = atoms_assign_pointer(this, name, str1D_ptr=ptr)
  end function atoms_assign_pointer_str1D

  function atoms_assign_pointer_str2D(this, name, ptr)
    type(Atoms), intent(in) :: this
    character(len=*), intent(in) :: name
    character(TABLE_STRING_LENGTH), pointer :: ptr(:,:)
    logical :: atoms_assign_pointer_str2D

    atoms_assign_pointer_str2D = atoms_assign_pointer(this, name, str2D_ptr=ptr)
  end function atoms_assign_pointer_str2D

  function atoms_assign_pointer_logical1D(this, name, ptr)
    type(Atoms), intent(in) :: this
    character(len=*), intent(in) :: name
    logical, pointer :: ptr(:)
    logical :: atoms_assign_pointer_logical1D

    atoms_assign_pointer_logical1D = atoms_assign_pointer(this, name, logical1D_ptr=ptr)
  end function atoms_assign_pointer_logical1D

  function atoms_assign_pointer_logical2D(this, name, ptr)
    type(Atoms), intent(in) :: this
    character(len=*), intent(in) :: name
    logical, pointer :: ptr(:,:)
    logical :: atoms_assign_pointer_logical2D

    atoms_assign_pointer_logical2D = atoms_assign_pointer(this, name, logical2D_ptr=ptr)
  end function atoms_assign_pointer_logical2D



  !% OMIT
  function atoms_assign_pointer(this, name, int1D_ptr, int2D_ptr, real1D_ptr, real2D_ptr, &
       str1D_ptr, str2D_ptr, logical1D_ptr, logical2D_ptr)
    type(Atoms), intent(in), target :: this
    character(len=*), intent(in) :: name
    integer, optional, pointer :: int1D_ptr(:), int2D_ptr(:,:)
    real(dp), optional, pointer :: real1D_ptr(:), real2D_ptr(:,:)
    character(TABLE_STRING_LENGTH), optional, pointer :: str1D_ptr(:), str2D_ptr(:,:)
    logical, optional, pointer :: logical1D_ptr(:), logical2D_ptr(:,:)
    logical :: atoms_assign_pointer

    integer :: lookup(3)

    if (.not. get_value(this%properties, name, lookup)) then
      atoms_assign_pointer = .false.
      return
    endif

    select case (lookup(1))

       case(PROPERTY_INT)
          if (lookup(2) == lookup(3)) then
             ! it's a rank 1 integer property, have we got a pointer for that
             if (present(int1D_ptr)) then
                int1D_ptr => this%data%int(lookup(2),1:this%N)
             else
                call system_abort('atoms_assign_pointer: property '//name//' needs a rank 1 integer pointer')
             end if
          else
             ! it's a rank 2 int property
             if (present(int2D_ptr)) then
                int2D_ptr => this%data%int(lookup(2):lookup(3),1:this%N)
             else
                call system_abort('atoms_assign_pointer: property '//name//' needs a rank 2 integer pointer')
             end if
          end if

       case(PROPERTY_REAL)
          if (lookup(2) == lookup(3)) then
             ! it's a rank 1 real property, have we got a pointer for that
             if (present(real1D_ptr)) then
                real1D_ptr => this%data%real(lookup(2),1:this%N)
             else
                call system_abort('atoms_assign_pointer: property '//name//' needs a rank 1 real(dp) pointer')
             end if
          else
             ! it's a rank 2 real property
             if (present(real2D_ptr)) then
                real2D_ptr => this%data%real(lookup(2):lookup(3),1:this%N)
             else
		atoms_assign_pointer = .false.
                call system_abort('atoms_assign_pointer: property '//name//' needs a rank 2 real(dp) pointer')
             end if
          end if

       case(PROPERTY_STR)
          if (lookup(2) == lookup(3)) then
             ! it's a rank 1 str property, have we got a pointer for that
             if (present(str1D_ptr)) then
                str1D_ptr => this%data%str(lookup(2),1:this%N)
             else
                call system_abort('atoms_assign_pointer: property '//name//' needs a rank 1 str pointer')
             end if
          else
             ! it's a rank 2 str property
             if (present(str2D_ptr)) then
                str2D_ptr => this%data%str(lookup(2):lookup(3),1:this%N)
             else
		atoms_assign_pointer = .false.
                call system_abort('atoms_assign_pointer: property '//name//' needs a rank 2 str) pointer')
             end if
          end if

       case(PROPERTY_LOGICAL)
          if (lookup(2) == lookup(3)) then
             ! it's a rank 1 logical property, have we got a pointer for that
             if (present(logical1D_ptr)) then
                logical1D_ptr => this%data%logical(lookup(2),1:this%N)
             else
                call system_abort('atoms_assign_pointer: property '//name//' needs a rank 1 logical pointer')
             end if
          else
             ! it's a rank 2 logical property
             if (present(logical2D_ptr)) then
                logical2D_ptr => this%data%logical(lookup(2):lookup(3),1:this%N)
             else
		atoms_assign_pointer = .false.
                call system_abort('atoms_assign_pointer: property '//name//' needs a rank 2 logical pointer')
             end if
          end if

    end select

    atoms_assign_pointer = .true.

  end function atoms_assign_pointer


  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  !% Set the cutoff (uniform or factor) to at least the requested value
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine atoms_set_cutoff_minimum(this, cutoff, cutoff_break)
    type(Atoms),      intent(inout) :: this
    real(dp),         intent(in)    :: cutoff
    real(dp), optional, intent(in)    :: cutoff_break

    if (present(cutoff_break)) then
      this%cutoff_break = max(this%cutoff_break, cutoff_break)
      this%cutoff = max(this%cutoff, cutoff)
    else
      this%cutoff_break = max(this%cutoff_break, cutoff)
      this%cutoff = max(this%cutoff, cutoff)
    endif

  end subroutine atoms_set_cutoff_minimum

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  !% Specify a uniform neighbour cutoff throughout the system.
  !% If zero, revert to default (uniform_cutoff=false, factor=DEFAULT_NNEIGHTOL)
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine atoms_set_cutoff(this, cutoff, cutoff_break)
    type(Atoms),      intent(inout) :: this
    real(dp),         intent(in)    :: cutoff
    real(dp), optional, intent(in)    :: cutoff_break

    if (present(cutoff_break)) then
      if (cutoff .feq. 0.0_dp) then
	this%cutoff = DEFAULT_NNEIGHTOL
      else
	this%cutoff = cutoff
      endif
      if (cutoff_break .feq. 0.0_dp) then
	this%cutoff_break = DEFAULT_NNEIGHTOL
      else
	this%cutoff_break = cutoff_break
      endif
    else
      if (cutoff .feq. 0.0_dp) then
	this%cutoff = DEFAULT_NNEIGHTOL
	this%cutoff_break = DEFAULT_NNEIGHTOL
      else
	this%cutoff = cutoff
	this%cutoff_break = cutoff
      endif
    endif

    if (cutoff .feq. 0.0_dp) then
      this%use_uniform_cutoff = .false.
    else
      this%use_uniform_cutoff = .true.
    endif

  end subroutine atoms_set_cutoff

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  !% Specify the neighbour cutoff to be a mulitple of the bond length
  !% of the two atoms' types.
  !% If zero, revert to default (DEFAULT_NNEIGHTOL)
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine atoms_set_cutoff_factor(this, factor, factor_break)
    type(Atoms),      intent(inout) :: this
    real(dp),         intent(in)    :: factor
    real(dp), optional, intent(in)    :: factor_break

    if (present(factor_break)) then
      if (factor .feq. 0.0_dp) then
	this%cutoff = DEFAULT_NNEIGHTOL
      else
	this%cutoff = factor
      endif
      if (factor_break .feq. 0.0_dp) then
	this%cutoff_break = DEFAULT_NNEIGHTOL
      else
	this%cutoff_break = factor_break
      endif
    else
      if (factor .feq. 0.0_dp) then
	this%cutoff = DEFAULT_NNEIGHTOL
	this%cutoff_break = DEFAULT_NNEIGHTOL
      else
	this%cutoff = factor
	this%cutoff_break = factor
      endif
    endif
    this%use_uniform_cutoff = .false.

  end subroutine atoms_set_cutoff_factor

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  !% Return the actual cutoff in \AA{} used by this atoms object
  !% used to form 'Z1---Z2' bonds. If 'this%use_uniform_cutoff' is
  !% true, then this is simply 'this%cutoff', otherwise the
  !% cutoff is used multiplied by the 'Z1---Z2' bond-length.
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  function atoms_cutoff(this, Z1, Z2)
    type(Atoms), intent(in) :: this
    integer, intent(in) :: Z1, Z2
    real(dp) :: atoms_cutoff

    if (this%use_uniform_cutoff) then
       atoms_cutoff = this%cutoff
    else
       atoms_cutoff = this%cutoff*bond_length(Z1, Z2)
    end if

  end function atoms_cutoff

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  !% Return the actual cutoff in \AA{} used by this atoms object
  !% used to break 'Z1---Z2' bonds. If 'this%use_uniform_cutoff' is
  !% true, then this is simply 'this%cutoff', otherwise the
  !% cutoff is used multiplied by the 'Z1---Z2' bond-length.
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  function atoms_cutoff_break(this, Z1, Z2)
    type(Atoms), intent(in) :: this
    integer, intent(in) :: Z1, Z2
    real(dp) :: atoms_cutoff_break

    if (this%use_uniform_cutoff) then
       atoms_cutoff_break = this%cutoff_break
    else
       atoms_cutoff_break = this%cutoff_break*bond_length(Z1, Z2)
    end if

  end function atoms_cutoff_break


  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  !% Zero data in an Atoms structure ---
  !% this doesn\'t finalise it or change it\'s size. We zero 'this%pos'
  !% and 'this%Z'. 
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine atoms_zero(this,indices)

    type(Atoms), intent(inout)    :: this
    integer, optional, intent(in) :: indices !% Optionally only zero the specified indices.

    if(present(indices)) then
       this%pos(:,indices) = 0.0_dp
       this%Z(indices) = 0
    else
       this%pos = 0.0_dp
       this%Z = 0
    end if

  end subroutine atoms_zero


  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  !% Change the lattice vectors, keeping the inverse lattice vectors
  !% up to date. Optionally map the existing atoms into the new cell
  !% and recalculate connectivity.
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine atoms_set_lattice(this,new_lattice,remap,reconnect)

    type(Atoms),              intent(inout) :: this
    real(dp), dimension(3,3), intent(in)    :: new_lattice
    logical, optional,        intent(in)    :: remap, reconnect

    this%lattice = new_lattice
    call matrix3x3_inverse(this%lattice,this%g)

    if (present(reconnect)) then
       if (reconnect) then
          call calc_connect(this)
          return
       end if           
    end if

    if (present(remap)) then
       if (remap) call map_into_cell(this)
    end if

  end subroutine atoms_set_lattice

  subroutine atoms_set_atoms_singlez(this, Z)
    type(Atoms), intent(inout) :: this
    integer, intent(in)        :: Z
    integer, allocatable, dimension(:) :: Zarray

    allocate(Zarray(this%N))
    Zarray = Z
    call atoms_set_atoms(this, Zarray)
    deallocate(Zarray)
  end subroutine atoms_set_atoms_singlez

  !% Set atomic numbers and optionally masses (if mass property is present)
  !% If 'mass' is not specified then 'ElementMass(Z)' is used.
  subroutine atoms_set_atoms(this, Z, mass)
    type(Atoms), intent(inout) :: this
    integer, dimension(:), intent(in) :: Z
    real(dp), optional, dimension(:), intent(in) :: mass
    
    this%Z = Z
    if (has_property(this, 'mass')) then
       this%mass = ElementMass(Z)
       ! Optionally override with user specified masses
       if (present(mass)) this%mass = mass
    end if
    if (has_property(this, 'species')) then
       this%species = ElementName(Z)
    end if
    
  end subroutine atoms_set_atoms

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  ! 
  ! Simple query functions
  !
  !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  !% Return the number of neighbour that atom $i$ has.
  !% If the optional arguments max_dist or max_factor are present 
  !% then only neighbours closer than this cutoff are included.
  function atoms_n_neighbours(this, i, max_dist, max_factor, alt_connect) result(n)
    type(Atoms), intent(in), target :: this
    integer, intent(in) :: i
    real(dp), optional, intent(in) :: max_dist, max_factor
    type(Connection), optional, intent(in), target :: alt_connect
    integer :: n

    integer :: j, m
    real(dp) :: r_ij
    type(Connection), pointer :: use_connect

    if (present(alt_connect)) then
      use_connect => alt_connect
    else
      use_connect => this%connect
    endif

    if (.not. use_connect%initialised) &
       call system_abort('Atoms_N_Neighbours: Atoms structure has no connectivity data. Call calc_connect first.')

    if (.not. associated(use_connect%neighbour1(i)%t)) then
      n = 0
      return
    endif

    if (.not. present(max_dist) .and. .not. present(max_factor)) then
       ! All neighbours
       n = use_connect%neighbour1(i)%t%N + use_connect%neighbour2(i)%t%N
    else if (present(max_dist)) then
       ! Only count neighbours within max_dist distance of i
       n = 0
       do m=1,use_connect%neighbour1(i)%t%N+use_connect%neighbour2(i)%t%N
          j = atoms_neighbour(this, i, m, distance=r_ij)
          if (r_ij < max_dist) n = n + 1
       end do
    else if (present(max_factor)) then
       ! Only count neighbours within max_factor of i
       n = 0
       do m=1,use_connect%neighbour1(i)%t%N+use_connect%neighbour2(i)%t%N
          j = atoms_neighbour(this, i, m, distance=r_ij)
          if (r_ij < bond_length(this%Z(i),this%Z(j))*max_factor) n = n + 1
       end do
    else
       call system_abort('Atoms_N_Neighbours: optional arguments max_dist and max_factor must not both be present')
    end if

  end function atoms_n_neighbours

  !% Return the index of the $n^{\mbox{\small{th}}}$ neighbour of atom $i$. Together with the
  !% previous function, this facilites a loop over the neighbours of atom $i$. Optionally, we
  !% return other geometric information, such as distance, direction cosines and difference vector,
  !% and also an direct index into the neighbour tables. If $i <= j$, this is an index into 'neighbour1(i)',
  !% if $i > j$, it is an index into 'neighbour1(j)'
  !%
  !%>   do n = 1,atoms_n_neighbours(at, i)
  !%>      j = atoms_neighbour(at, i, n, distance, diff, cosines, shift, index)
  !%>
  !%>      ...
  !%>   end do
  !%
  !% if distance > max_dist, return 0, and do not waste time calculating other quantities
  function atoms_neighbour(this, i, n, distance, diff, cosines, shift, index, max_dist, jn, alt_connect) result(j)
    type(Atoms), intent(in), target :: this
    integer ::i, j, n
    real(dp), optional, intent(out) :: distance
    real(dp), dimension(3), optional, intent(out) :: diff
    real(dp), optional, intent(out) :: cosines(3)
    integer,  optional, intent(out) :: shift(3)
    integer,  optional, intent(out) :: index
    real(dp), optional, intent(in)  :: max_dist
    integer,  optional, intent(out) :: jn
    type(Connection), optional, intent(in), target :: alt_connect

    real(dp)::mydiff(3)
    integer ::myshift(3)
    integer ::i_n1n, j_n1n, i_njn
    type(Connection), pointer :: use_connect

    if (present(alt_connect)) then
      use_connect => alt_connect
    else
      use_connect => this%connect
    endif

    if (.not. associated(use_connect%neighbour1(i)%t)) then
      call system_abort("called atoms_neighbour on atom " // i // " which has no allocated neighbour1 table")
      return
    endif

    ! First we give the neighbour2 entries (i > j) then the neighbour1 (i <= j)
    ! This order chosen to give neighbours in approx numerical order but doesn't matter
    if (use_connect%initialised) then
       i_n1n = n-use_connect%neighbour2(i)%t%N
       if (n <= use_connect%neighbour2(i)%t%N) then
          j = use_connect%neighbour2(i)%t%int(1,n)
          j_n1n = use_connect%neighbour2(i)%t%int(2,n)
          if(present(index)) index = j_n1n
       else if (i_n1n <= use_connect%neighbour1(i)%t%N) then
          j = use_connect%neighbour1(i)%t%int(1,i_n1n)
          if(present(index)) index = i_n1n
       else
          call system_abort('atoms_neighbour: '//n//' out of range for atom '//i//&
               ' Should be in range 1 < n <= '//atoms_n_neighbours(this, i))
       end if
    else
       call system_abort('atoms_neighbour: Atoms structure has no connectivity data. Call calc_connect first.')
    end if

    if(present(jn)) then
       if(i < j) then
          do i_njn = 1, use_connect%neighbour2(j)%t%N
             if( (use_connect%neighbour2(j)%t%int(1,i_njn)==i) .and. &
             & (use_connect%neighbour2(j)%t%int(2,i_njn)==i_n1n) ) jn = i_njn
          enddo
       elseif(i > j) then
          jn = j_n1n + use_connect%neighbour2(j)%t%N
       else
          do i_njn = 1, use_connect%neighbour1(j)%t%N
             if( (use_connect%neighbour1(j)%t%int(1,i_njn) == i) .and. &
             & all(use_connect%neighbour1(j)%t%int(2:4,i_njn) == -use_connect%neighbour1(i)%t%int(2:4,i_n1n))) &
             & jn = i_njn + use_connect%neighbour2(j)%t%N
          enddo
       endif
    endif
          
    ! found neighbour, now check for optional requests
    if(present(distance)) then
       if(i <= j) then 
          distance = use_connect%neighbour1(i)%t%real(1,i_n1n)
       else
          distance = use_connect%neighbour1(j)%t%real(1,j_n1n)
       end if
       if (present(max_dist)) then
	 if (distance > max_dist) then
	   j = 0
	   return
	 endif
       endif
    else
       if (present(max_dist)) then
	 if (i <= j) then
	   if (use_connect%neighbour1(i)%t%real(1,i_n1n) > max_dist) then
	     j = 0
	     return
	   endif
	 else
	   if (use_connect%neighbour1(j)%t%real(1,j_n1n) > max_dist) then
	     j = 0
	     return
	   endif
	 endif
       endif
    end if

    if(present(diff) .or. present(cosines) .or. present(shift)) then
       if (i <= j) then
          myshift = use_connect%neighbour1(i)%t%int(2:4,i_n1n)
       else
          myshift = -use_connect%neighbour1(j)%t%int(2:4,j_n1n)
       end if

       if(present(shift)) shift = myshift

       if(present(diff) .or. present(cosines)) then
          mydiff = this%pos(:,j) - this%pos(:,i) + (this%lattice .mult. myshift)
          if(present(diff)) diff = mydiff
          if(present(cosines)) then
	    if (norm(mydiff) > 0.0_dp) then
	      cosines = mydiff / norm(mydiff)
	    else
	      cosines = 0.0_dp
	    endif
	  endif
       end if
    end if

  end function atoms_neighbour


  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  ! Adding and Removing atoms: adding/removing single or multiple atoms
  ! For the atoms variable in Dynamical System, this should be called
  ! shift there to avoid inconsistencies with DS's data
  ! 
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine add_atom_single(this, pos, Z, mass, travel)

    type(Atoms),       intent(inout)            :: this
    real(dp),          intent(in), dimension(3) :: pos
    integer,           intent(in)               :: Z
    real(dp), optional,  intent(in)             :: mass
    integer, optional, intent(in), dimension(3) :: travel

    if(present(travel)) then
       if(present(mass)) then
          call add_atom_multiple(this, pos=reshape(pos, (/3,1/)), Z=(/Z/), mass=(/mass/), travel=reshape(travel, (/3,1/)))
       else
          call add_atom_multiple(this, pos=reshape(pos, (/3,1/)), Z=(/Z/), travel=reshape(travel, (/3,1/)))
       end if
    else
       if(present(mass)) then
          call add_atom_multiple(this, pos=reshape(pos, (/3,1/)), Z=(/Z/), mass=(/mass/))
       else
          call add_atom_multiple(this, pos=reshape(pos, (/3,1/)), Z=(/Z/))
       end if
    end if
  end subroutine add_atom_single

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine add_atom_multiple(this, pos, Z, mass,  velo, acc, travel, data)

    type(Atoms),       intent(inout)              :: this
    real(dp), optional, intent(in), dimension(:,:) :: pos
    integer,  optional, intent(in), dimension(:)   :: Z
    real(dp), optional, intent(in), dimension(:)  :: mass
    integer, optional, intent(in), dimension(:,:) :: travel
    real(dp),optional, intent(in), dimension(:,:) :: velo, acc

    type(Table), optional, intent(in)             :: data

    integer                                       :: oldN,i, lookup(3)
    logical :: dummy

    if (present(data)) then

       this%N = this%N + data%N

       if (this%data%intsize /= data%intsize) &
            call system_abort('Add_Atoms: this%data%intsize /= data%intsize')
       if (this%data%realsize /= data%realsize) &
            call system_abort('Add_Atoms: this%data%realsize /= data%realsize')

       call append(this%data, data)
       call atoms_repoint(this)
    
    else
       if (.not. present(Z)) call system_abort('Atoms_Add: Z must be present if data is not')
       oldN = this%N
       this%N = this%N + size(Z)

       !Check the sizes of the input arrays for consistency
       call check_size('Pos',pos,(/3,size(Z)/),'Add_Atom')
       if (present(travel)) call check_size('Travel',travel,(/3,size(Z)/),'Add_Atom')
       if (present(velo)) call check_size('Velo', velo, (/3,size(Z)/), 'Add_Atom')
       if (present(acc)) call check_size('Acc', acc, (/3,size(Z)/), 'Add_Atom')

       call append(this%data, blank_rows=size(Z))
       call atoms_repoint(this)

       ! First check the integer properties...
       if (.not. get_value(this%properties, 'Z', lookup)) &
            call system_abort('Atoms_Add: this atoms has no Z property')
       this%data%int(lookup(2),oldN+1:this%N) = Z

       ! set species from Z
       if (.not. get_value(this%properties, 'species', lookup)) &
            call system_abort('Atoms_Add: this atoms has no species property')
       do i=1,size(Z)
          this%data%str(lookup(2),oldN+i) = ElementName(Z(i))
       end do

       if (present(travel)) then
          if (.not. get_value(this%properties, 'travel', lookup)) &
               call system_abort('Atoms_Add: this atoms has no travel property')
          this%data%int(lookup(2):lookup(3),oldN+1:this%N) = travel
       else
          if (get_value(this%properties, 'travel', lookup)) &
          & this%data%int(lookup(2):lookup(3),oldN+1:this%N) = 0
       end if

       ! Set masks to 1 if properties for them exist
       if (get_value(this%properties, 'move_mask', lookup)) &
            this%data%int(lookup(2):lookup(3),oldN+1:this%N) = 1
       if (get_value(this%properties, 'damp_mask', lookup)) &
            this%data%int(lookup(2):lookup(3),oldN+1:this%N) = 1
       if (get_value(this%properties, 'thermostat_region', lookup)) &
            this%data%int(lookup(2):lookup(3),oldN+1:this%N) = 1

       ! ... and now the real properties
       if (get_value(this%properties, 'mass', lookup)) then
          if (present(mass)) then
             this%data%real(lookup(2),oldN+1:this%N) = mass
          else
             this%data%real(lookup(2),oldN+1:this%N) = ElementMass(Z)
          end if
       else if (present(mass)) then
          ! mass specified but property doesn't yet exist, so create it...
          call add_property(this, 'mass', ElementMass(this%Z))
          dummy = get_value(this%properties, 'mass', lookup)
          ! ... and then override for new atoms
          this%data%real(lookup(2),oldN+1:this%N) = mass
       end if

       if (.not. present(pos)) &
            call system_abort('Atoms_Add: pos must be present if data is not')
       if (.not. get_value(this%properties, 'pos', lookup)) &
            call system_abort('Atoms_Add: this atoms has no pos property')
       this%data%real(lookup(2):lookup(3),oldN+1:this%N) = pos

       if (present(velo) .and. get_value(this%properties, 'velo', lookup)) &
            this%data%real(lookup(2):lookup(3),oldN+1:this%N) = velo

       if (present(acc) .and. get_value(this%properties, 'acc', lookup)) &
            this%data%real(lookup(2):lookup(3),oldN+1:this%N) = acc

    end if

    call finalise(this%connect)

  end subroutine add_atom_multiple


  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine remove_atom_single(this, i)

    type(Atoms), intent(inout) :: this
    integer,     intent(in)    :: i

    call remove_atom_multiple(this,(/i/))

  end subroutine remove_atom_single

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


  subroutine remove_atom_multiple(this, atom_indices)

    type(Atoms), intent(inout)                 :: this
    integer,     intent(in), dimension(:)      :: atom_indices

    !Delete the connection data because the atomic indices become mangled
    call connection_finalise(this%connect)

    ! Remove rows from data table
    call delete_multiple(this%data, atom_indices)

    ! update N
    this%N = this%N - size(atom_indices)

    call atoms_repoint(this)

  end subroutine remove_atom_multiple



  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  ! 
  !% Map atomic fractional positions back into the unit cell
  !% $-0.5 \le t_x,t_y,t_z < 0.5$
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine atoms_map_into_cell(this)
    type(Atoms), intent(inout) :: this

    integer      :: i,j, n, m
    integer      :: shift(3)
    logical      :: mapped


    ! Loop over all atoms
    ! Convert cartesians to lattice co-ords
    ! If outside unit cell:
    !    add lattice vectors to get back into cell
    !    update travel
    !    update shifts

    if (.not. has_property(this, 'travel')) &
         call add_property(this, 'travel', 0, n_cols=3)

    do i=1,this%N
      call map_into_cell(this%pos(:,i), this%lattice, this%g, shift, mapped)
      if (mapped) then
	 this%travel(:,i) = this%travel(:,i) - shift
	 if (this%connect%initialised) then
	   do n=1,atoms_n_neighbours(this, i)  ! Loop over all atom i's neighbours

	      j = atoms_neighbour(this, i, n, index=m) ! get neighbour

	      ! now update the data for atom i and the current neighbour
	      if (i <= j) then
		 this%connect%neighbour1(i)%t%int(2:4,m) = this%connect%neighbour1(i)%t%int(2:4,m) + shift
	      else
		 this%connect%neighbour1(j)%t%int(2:4,m) = this%connect%neighbour1(j)%t%int(2:4,m) - shift
	      end if
	   end do
	end if ! this%connect%initialised
      end if ! mapped
    end do ! i=1..N

  end subroutine atoms_map_into_cell

  subroutine vec_map_into_cell(pos, lattice, g, shift, mapped)
    real(dp), intent(inout) :: pos(3)
    real(dp), intent(in) :: lattice(3,3), g(3,3)
    integer, intent(out), optional :: shift(3)
    logical, intent(out), optional :: mapped

    integer n, k
    real(dp) :: lattice_coord(3)
    logical :: my_mapped

   lattice_coord = g .mult. pos(:)
   my_mapped = .false.
   if (present(shift)) shift = 0

   do n=1,3
      if ((lattice_coord(n) < -0.5_dp) .or. (lattice_coord(n) >= 0.5_dp)) then
	 k = floor(lattice_coord(n)+0.5_dp)
	 lattice_coord(n) = lattice_coord(n) - k
	 if (present(shift)) shift(n) = -k
	 my_mapped = .true.
      end if
   end do

   ! if this atom has been mapped then recalculate its position and shifts for its neighbours
   if (my_mapped) then
      pos(:) = lattice .mult. lattice_coord
   end if
   if (present(mapped)) mapped = my_mapped
  end subroutine vec_map_into_cell


  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  !% The subroutine 'calc_dists' updates the stored distance tables using 
  !% the stored connectivity and shifts. This should be called every time
  !% any atoms are moved (e.g. it is called by 'advance_verlet').
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  subroutine calc_dists(this, parallel)

    type(Atoms), intent(inout) :: this
    logical, optional, intent(in) :: parallel
    integer                    :: i, j, n, index
    integer, dimension(3)      :: shift
    real(dp), dimension(3)     :: j_pos
    logical :: do_parallel
#ifdef _MPI
    integer::error, Nelements, mpi_pos, mpi_old_pos
    include "mpif.h"
    real(dp), allocatable :: mpi_send(:), mpi_recv(:)
#endif
    ! Flag to specify whether or not to parallelise calculation.
    ! Only actually run in parallel if parallel==.true. AND
    ! _MPI is #defined. Default to serial mode.
    do_parallel = .false.
    if (present(parallel)) do_parallel = parallel

#ifdef _MPI
    if (do_parallel) then
       ! Nelements = sum(this%connect%neighbour1(i)%t%N)
       Nelements = 0
       do i=1,this%N
          Nelements = Nelements + this%connect%neighbour1(i)%t%N
       end do

       allocate(mpi_send(Nelements))
       allocate(mpi_recv(Nelements))
       if (Nelements > 0) then
	 mpi_send = 0.0_dp
	 mpi_recv = 0.0_dp
	end if
       mpi_pos = 1
    end if
#endif

    if (.not.this%connect%initialised) &
         call system_abort('CalcDists: Connect is not yet initialised')

    do i = 1, this%N

#ifdef _MPI
       if (do_parallel) then
          mpi_old_pos = mpi_pos
          mpi_pos = mpi_pos + this%connect%neighbour1(i)%t%N

          ! cycle loop if processor rank does not match
          if(mod(i, mpi_n_procs()) .ne. mpi_id()) cycle
       end if
#endif

       do n = 1, atoms_n_neighbours(this, i) 

          j = atoms_neighbour(this, i, n, shift=shift, index=index)

          j_pos = this%pos(:,j) + ( this%lattice .mult. shift )

          if (i <= j) then
             this%connect%neighbour1(i)%t%real(1,index) = norm(j_pos - this%pos(:,i))
          else
             this%connect%neighbour1(j)%t%real(1,index) = norm(j_pos - this%pos(:,i))
          end if

       end do

#ifdef _MPI
       if (do_parallel) then
	  if (mpi_old_pos <= mpi_pos-1) then
	    mpi_send(mpi_old_pos:mpi_pos-1) = &
		 this%connect%neighbour1(i)%t%real(1,1:this%connect%neighbour1(i)%t%N)
	  end if
       end if
#endif      

    end do

#ifdef _MPI
    if (do_parallel) then
       ! collect mpi results
       if (Nelements > 0) then
	 call mpi_allreduce(mpi_send, mpi_recv, &
	      size(mpi_send), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, error)
	 call abort_on_mpi_error(error, "Calc_Dists: MPI_ALL_REDUCE()")
       end if

       mpi_pos = 1
       do i=1, this%N
	  if (this%connect%neighbour1(i)%t%N > 0) then
	    this%connect%neighbour1(i)%t%real(1,1:this%connect%neighbour1(i)%t%N) = &
		 mpi_recv(mpi_pos:mpi_pos+this%connect%neighbour1(i)%t%N-1)
	    mpi_pos = mpi_pos + this%connect%neighbour1(i)%t%N
	  endif
       end do

       if (Nelements > 0) then
	 deallocate(mpi_send, mpi_recv)
       end if
    end if
#endif

  end subroutine calc_dists


  !% Difference vector between atoms $i$ and $j$ if they are separated by a shift of 'shift'
  function diff(this, i, j, shift)
    type(Atoms), intent(in)    :: this
    integer,     intent(in)    :: i,j
    integer,  dimension(3)     :: shift
    real(dp), dimension(3) :: diff

    diff = this%pos(:,j) - this%pos(:,i) + (this%lattice .mult. shift)

  end function diff


  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  !  diff_min_image interface
  !  
  !  return relative vector from one position to another, adhering to PBC
  !  and minimum image conventions
  ! 
  !  Flavours are: atom-atom, vector-atom, atom-vector, vector-vector
  !
  !  All are accessible using the 'diff_min_image' interface
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


  function diff_atom_atom(this, i, j, shift)

    type(Atoms), intent(in)    :: this
    integer,     intent(in)    :: i,j
    integer,  dimension(3), optional     :: shift
    real(dp), dimension(3)     :: diff_atom_atom
    real(dp)                   :: dummy

    integer, dimension(3) :: myshift

    dummy = distance_min_image(this,i,j, shift=myshift)

    diff_atom_atom = this%pos(:,j) - this%pos(:,i) + (this%lattice .mult. myshift)

    if (present(shift)) shift = myshift

  end function diff_atom_atom


  function diff_vec_atom(this, v, j)

    type(Atoms), intent(in)    :: this
    real(dp), dimension(3)     :: v
    integer,     intent(in)    :: j
    real(dp), dimension(3)     :: diff_vec_atom
    integer,  dimension(3)     :: shift
    real(dp)                   :: dummy

    dummy = distance_min_image(this,v,j, shift=shift)

    diff_vec_atom = this%pos(:,j) - v + (this%lattice .mult. shift)

  end function diff_vec_atom

  function diff_atom_vec(this, i, w)

    type(Atoms), intent(in)    :: this
    integer,     intent(in)    :: i
    real(dp), dimension(3)     :: w
    real(dp), dimension(3)     :: diff_atom_vec
    integer,  dimension(3)     :: shift
    real(dp)                   :: dummy

    dummy = distance_min_image(this,i,w, shift=shift)

    diff_atom_vec = w - this%pos(:,i) + (this%lattice .mult. shift)

  end function diff_atom_vec

  function diff_vec_vec(this, v, w)

    type(Atoms), intent(in)    :: this
    real(dp), dimension(3)     :: v, w
    real(dp), dimension(3)     :: diff_vec_vec
    integer,  dimension(3)     :: shift
    real(dp)                   :: dummy

    dummy = distance_min_image(this,v,w, shift=shift)

    diff_vec_vec = w - v + (this%lattice .mult. shift)

  end function diff_vec_vec

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  !% Return the real position of the atom, taking into account the
  !% stored travel across the periodic boundary conditions.
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  function realpos(this,i)

    type(Atoms), intent(in) :: this
    integer,     intent(in) :: i
    real(dp), dimension(3)  :: realpos

    realpos = (this%lattice .mult. this%travel(:,i)) + this%pos(:,i)

  end function realpos


  !% Return distance between atoms $i$ and $j$ if they are separated by a shift
  !% of 'shift'.
  function distance(this, i, j, shift)
    type(Atoms), intent(in)::this
    integer,     intent(in)::i, j, shift(3)
    real(dp)::distance

    distance = norm(this%pos(:,j)+(this%lattice .mult. shift)-this%pos(:,i))
  end function distance

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  ! distance_min_image interface
  !
  ! Actual distance computing routines. 
  !
  ! The real work is done in the function that
  ! computes the distance of two general vector positions.  when
  ! atomic indices are specified, they are first converted to vector
  ! positions.
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  function distance8_atom_atom(this,i,j,shift)

    type(Atoms),                        intent(in)   :: this
    integer,                            intent(in)   :: i,j
    integer,  optional, dimension(3),   intent(out)  :: shift
    real(dp)                                         :: distance8_atom_atom

    distance8_atom_atom = distance8_vec_vec(this,this%pos(:,i),this%pos(:,j),shift)

  end function distance8_atom_atom

  function distance8_atom_vec(this,i,v,shift)

    type(Atoms),                        intent(in)  :: this
    integer,                            intent(in)  :: i
    real(dp),           dimension(3),   intent(in)  :: v
    integer,  optional, dimension(3),   intent(out) :: shift
    real(dp)                                        :: distance8_atom_vec

    distance8_atom_vec = distance8_vec_vec(this,this%pos(:,i),v,shift)

  end function distance8_atom_vec

  function distance8_vec_atom(this,v,j,shift)

    type(Atoms),                        intent(in)  :: this
    real(dp),           dimension(3),   intent(in)  :: v
    integer,                            intent(in)  :: j
    integer,  optional, dimension(3),   intent(out) :: shift
    real(dp)                                        :: distance8_vec_atom

    distance8_vec_atom = distance8_vec_vec(this,v,this%pos(:,j),shift)

  end function distance8_vec_atom

  ! This is the general function

  function distance8_vec_vec(this,v,w,shift)

    type(Atoms),                        intent(in)  :: this
    real(dp),           dimension(3),   intent(in)  :: v,w
    integer,  optional, dimension(3),   intent(out) :: shift
    real(dp)                                        :: distance8_vec_vec, dist2, tmp
    real(dp),           dimension(3)                :: dvw, lattice_coord
    integer,            dimension(3)                :: init_val
    integer                                         :: i,j,k

    !get the difference vector and convert to lattice co-ordinates
    !use the precomputed matrix inverse if possible
    dvw = w - v
    lattice_coord = this%g .mult. dvw

    init_val = (/0,0,0/)

    !work out which block of 8 cells we are testing
    where (lattice_coord > 0.0_dp) init_val = -1

    dist2 = huge(1.0_dp) ! effectively +ve infinity

    !now loop over the cells and test
    do k=init_val(3), init_val(3)+1
       do j=init_val(2), init_val(2)+1
          do i=init_val(1), init_val(1)+1

             !construct the shifted vector
             tmp = norm2(dvw + this%lattice(:,1)*i +this%lattice(:,2)*j + this%lattice(:,3)*k)
             !test if it is the smallest so far and store the shift if necessary
             if (tmp < dist2) then
                dist2 = tmp
                if (present(shift)) shift = (/i,j,k/)
             end if

          end do
       end do
    end do

    distance8_vec_vec = sqrt(dist2)

  end function distance8_vec_vec

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  !% Cosine of the angle j--i--k
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  function cosine(this,i,j,k)

    type(Atoms), intent(in) :: this
    integer,     intent(in)    :: i,j,k
    real(dp)                   :: cosine
    real(dp), dimension(3)     :: ij, ik

    if ((i == j) .or. (i == k)) call system_abort('Cosine: i == j or i == k')

    ij = diff_min_image(this,i,j)
    ik = diff_min_image(this,i,k)
    cosine = (ij .dot. ik) / (norm(ij)*norm(ik))

  end function cosine

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  !% Cosine of the angle n--i--m where {$n,m$} are the {$n$th, $m$th} neighbours of i 
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  function cosine_neighbour(this,i,n,m)

    type(Atoms), intent(in) :: this
    integer,     intent(in)    :: i,n,m
    real(dp)                   :: cosine_neighbour
    real(dp), dimension(3)     :: in, im
    integer::j

    if(n == m) then
       cosine_neighbour = 1.0_dp
       return
    end if
    j = atoms_neighbour(this, i, n, diff=in)
    j = atoms_neighbour(this, i, m, diff=im)
    cosine_neighbour = (in .dot. im) / (norm(in)*norm(im))

  end function cosine_neighbour

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  ! Direction cosines
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  !% Given two atoms $i$ and $j$ and a shift returns the direction 
  !% cosines of the differnece vector from $i$ to $j$.
  function direction_cosines(this,i,j,shift)

    type(Atoms), intent(in) :: this
    integer,     intent(in)    :: i,j, shift(3)
    real(dp), dimension(3)     :: direction_cosines, diffv

    diffv = diff(this,i,j,shift)
    direction_cosines = diffv / norm(diffv)

  end function direction_cosines

  !% Direction cosines of the difference vector from $i$ to $j$
  function direction_cosines_min_image(this,i,j)

    type(Atoms), intent(in) :: this
    integer,     intent(in)    :: i,j
    real(dp), dimension(3)     :: direction_cosines_min_image, diffv

    if (i == j) call system_abort('Cosines: i == j')

    diffv = diff_min_image(this,i,j)
    direction_cosines_min_image = diffv / norm(diffv)

  end function direction_cosines_min_image

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  ! Connectivity procedures: Initialise and Finalise
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


  !% Initialise a Connection object for a given number of atoms 'N'.
  !% If the optional Atoms argument is present then we calculate
  !% the atomic density to initialise the default lengths of the neighbour
  !% list for efficient memory usage.
   subroutine connection_initialise(this,N, pos, lattice, g, cutoff, use_uniform_cutoff, origin, extent)
    type(Connection),   intent(inout) :: this
    integer,            intent(in)    :: N    ! No. of atoms
    real(dp), optional, intent(in) :: pos(:,:), lattice(3,3), g(3,3), cutoff
    logical, optional, intent(in) :: use_uniform_cutoff
    real(dp), optional, intent(in) :: origin(3), extent(3,3)

    ! If already initialised, destroy the existing data an start again
    if (this%initialised) call connection_finalise(this)

    call connection_fill(this, N, pos, lattice, g, cutoff, use_uniform_cutoff, origin, extent)
  end subroutine connection_initialise

  subroutine connection_fill(this, N, pos, lattice, g, cutoff, use_uniform_cutoff, origin, extent)
    type(Connection),   intent(inout) :: this
    integer,            intent(in)    :: N    ! No. of atoms
    real(dp), optional, intent(in) :: pos(:,:), lattice(3,3), g(3,3), cutoff
    logical, optional, intent(in) :: use_uniform_cutoff
    real(dp), optional, intent(in) :: origin(3), extent(3,3)

    integer                           :: i, n0
    real(dp)                          :: mycutoff, extent_inv(3,3), subregion_center(3)
    logical :: do_subregion

    ! If we're given at atoms structure then use density N/V to work
    ! out roughtly how many neighbours to expect within sphere of radius
    ! cutoff. Assume half of these belong in neighbour1 and half in 
    ! neighbour2. Increment is set to n0/2.
    if (present(lattice) .and. present(cutoff) .and. present(use_uniform_cutoff)) then
       if(use_uniform_cutoff) then
          mycutoff = cutoff
       else
          mycutoff = 3.0_dp*cutoff
       endif
       n0 = int(0.5_dp*4.0_dp/3.0_dp*PI*mycutoff**3*N/Cell_Volume(lattice))
    else
       n0 = 5
    end if

    if (present(origin) .and. present(extent)) then
      if (.not.present(lattice) .or. .not.present(g)) &
	call system_abort ("connection_fill got origin and extent, so trying to do subregion, but lattice or g are missing")
      do_subregion = .true.
      call matrix3x3_inverse(extent,extent_inv)
      subregion_center = origin + 0.5_dp*sum(extent,2)
    else
      do_subregion = .false.
    endif

    if (.not. allocated(this%neighbour1)) allocate(this%neighbour1(N))
    if (.not. allocated(this%neighbour2)) allocate(this%neighbour2(N))
    do i=1,N
       if (do_subregion) then
	 if (.not. is_in_subregion(pos(:,i), subregion_center, lattice, g, extent_inv)) then
	   if (associated(this%neighbour1(i)%t)) then
	      call connection_remove_atom(this, i)
	      call finalise(this%neighbour1(i)%t)
	      call finalise(this%neighbour2(i)%t)
	      deallocate(this%neighbour1(i)%t)
	      deallocate(this%neighbour2(i)%t)
	   endif
	   cycle
	 endif 
       endif
       if (.not. associated(this%neighbour1(i)%t)) then
	 allocate(this%neighbour1(i)%t)
	 call allocate(this%neighbour1(i)%t,4,1, 0, 0, max(n0, 1))
	 this%neighbour1(i)%t%increment = max(n0/2, 1)

	 allocate(this%neighbour2(i)%t)
	 call allocate(this%neighbour2(i)%t,2,0, 0, 0, max(n0, 1))
	 this%neighbour2(i)%t%increment = max(n0/2, 1)
       endif
    end do

    this%initialised = .true.

  end subroutine connection_fill

  function is_in_subregion(p, center, lattice, lattice_inv, extent_inv)
    real(dp), intent(in) :: p(3), center(3)
    real(dp), intent(in) :: lattice(3,3), lattice_inv(3,3), extent_inv(3,3)
    logical :: is_in_subregion

    real(dp) :: relative_p(3), extent_lattice_relative_p(3)

    relative_p = p - center
    call map_into_cell(relative_p, lattice, lattice_inv)

    extent_lattice_relative_p = extent_inv .mult. relative_p

    if (any(extent_lattice_relative_p < -0.5_dp) .or. any(extent_lattice_relative_p > 0.5_dp)) then
      is_in_subregion = .false.
    else
      is_in_subregion = .true.
    endif

  end function is_in_subregion


  !% OMIT
  subroutine connection_cells_initialise(this,cellsNa,cellsNb,cellsNc,Natoms)

    type(Connection),  intent(inout) :: this
    integer,           intent(in)    :: cellsNa,cellsNb,cellsNc ! No. cells in a,b,c directions
    integer, optional, intent(in)    :: Natoms                  ! Number of atoms
    integer                          :: i,j,k,av_atoms,stdev_atoms,Ncells

    if (this%cells_initialised) call connection_cells_finalise(this)

    !Set length and increment based on binomial statistics (if possible)
    Ncells = cellsNa * cellsNb * cellsNc
    if (present(Natoms)) then
       av_atoms = Natoms / Ncells
       stdev_atoms = int(sqrt(real((Natoms * (Ncells - 1) / (Ncells * Ncells)),dp)))
    else
       av_atoms = 100   !defaults if number of atoms is not given
       stdev_atoms = 10
    end if

    allocate( this%cell(cellsNa,cellsNb,cellsNc) )

    do k=1,cellsNc
       do j=1,cellsNb
          do i=1,cellsNa
             call allocate(this%cell(i,j,k),1,0,0,0,max(1,(av_atoms+2*stdev_atoms))) !Good for 97.5% of cells
             call set_increment(this%cell(i,j,k),max(1,stdev_atoms))
          end do
       end do
    end do

    this%cellsNa = cellsNa
    this%cellsNb = cellsNb
    this%cellsNc = cellsNc

    this%cells_initialised = .true.

  end subroutine connection_cells_initialise


  !% Finalise this connection object
  subroutine connection_finalise(this)

    type(Connection), intent(inout) :: this
    integer                         :: i

    !do nothing if not initialised / already finalised
    if (.not.this%initialised) return

    if (allocated(this%neighbour1)) then
      do i=1,size(this%neighbour1)
	 if (associated(this%neighbour1(i)%t)) then
	   call finalise(this%neighbour1(i)%t)
	   deallocate(this%neighbour1(i)%t)
	 endif
      end do
    endif

    if (allocated(this%neighbour2)) then
      do i=1,size(this%neighbour2)
	 if (associated(this%neighbour2(i)%t)) then
	   call finalise(this%neighbour2(i)%t)
	   deallocate(this%neighbour2(i)%t)
	 endif
      end do
    endif


    if(allocated(this%neighbour1)) deallocate(this%neighbour1)
    if(allocated(this%neighbour2)) deallocate(this%neighbour2)

    call connection_cells_finalise(this)

    this%initialised = .false.

  end subroutine connection_finalise

  !% Wipe the contents of the connection tables, but keep the allocation
  subroutine connection_wipe(this)

    type(Connection), intent(inout) :: this
    integer                         :: i

    !do nothing if not initialised / already finalised
    if (.not.this%initialised) return

    do i=1,size(this%neighbour1)
       call wipe(this%neighbour1(i)%t)
    end do

    do i=1,size(this%neighbour2)
       call wipe(this%neighbour2(i)%t)
    end do

    call wipe_cells(this)

  end subroutine connection_wipe

  !% OMIT
  subroutine connection_cells_finalise(this)

    type(Connection), intent(inout) :: this
    integer                         :: i,j,k

    if (.not.this%cells_initialised) return

    do k=1,this%cellsNc
       do j=1,this%cellsNb
          do i=1,this%cellsNa
             call finalise(this%cell(i,j,k))
          end do
       end do
    end do

    deallocate(this%cell)

    this%cells_initialised = .false.

  end subroutine connection_cells_finalise

  subroutine connection_assignment(to,from)

    type(Connection), intent(inout) :: to
    type(Connection), intent(in)    :: from
    integer                         :: i,j,k

    call connection_finalise(to)

    if (.not.from%initialised) return


    call connection_initialise(to,size(from%neighbour1))


    ! Use Append to append the source tables to our (empty) destination tables
    do i=1,size(from%neighbour1)
       call append(to%neighbour1(i)%t,from%neighbour1(i)%t)
    end do

    do i=1,size(from%neighbour2)
       call append(to%neighbour2(i)%t,from%neighbour2(i)%t)
    end do


    !If cell data is present then copy that too
    if (from%cells_initialised) then
       call connection_cells_initialise(to,from%cellsNa,from%cellsNb,from%cellsNc)
       do k = 1, from%cellsNc
          do j = 1, from%cellsNb
             do i = 1, from%cellsNa
                call append(to%cell(i,j,k),from%cell(i,j,k))
             end do
          end do
       end do
    end if

  end subroutine connection_assignment

  !% OMIT
  subroutine wipe_cells(this)

    type(Connection) :: this
    integer          :: i,j,k

    if (this%cells_initialised) then

       do k = 1, this%cellsNc
          do j = 1, this%cellsNb
             do i = 1, this%cellsNa
                call wipe(this%cell(i,j,k))
             end do
          end do
       end do

    end if

  end subroutine wipe_cells

  !% Test if an atom's $n$th neighbour is one if its nearest neighbours
  function is_nearest_neighbour(this,i,n, alt_connect)

    type(Atoms), intent(in), target :: this
    integer,     intent(in) :: i,n
    type(Connection), intent(in), optional, target :: alt_connect
    logical                 :: is_nearest_neighbour

    real(dp)                :: d
    integer :: j

    is_nearest_neighbour = .false.

    j = atoms_neighbour(this, i, n, distance=d, alt_connect=alt_connect)
    if (d < (bond_length(this%Z(i),this%Z(j))*this%nneightol)) &
         is_nearest_neighbour = .true.

  end function is_nearest_neighbour

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !
  ! Connectivity procedures
  !
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  !% Test if atom $i$ is a neighbour of atom $j$ and update 'this%connect' as necessary.
  !% Called by 'calc_connect'. The 'shift' vector is added to the position of the $j$ atom
  !% to get the correct image position.
  subroutine test_form_bond(this,cutoff, use_uniform_cutoff, Z, pos, lattice, i,j, shift, check_for_dup)

    type(Connection), intent(inout) :: this
    real(dp), intent(in) :: cutoff
    logical, intent(in) :: use_uniform_cutoff
    integer, intent(in) :: Z(:)
    real(dp), intent(in) :: pos(:,:), lattice(3,3)
    integer,     intent(in)    :: i,j
    integer,     intent(in)    :: shift(3)
    logical, intent(in), optional :: check_for_dup

    logical                    :: do_check_for_dup
    integer                    :: index
    real(dp)                   :: d
    real(dp)                   :: use_cutoff

    if (i > j) return

    if (.not. associated(this%neighbour1(i)%t) .or. .not. associated(this%neighbour1(j)%t)) return

    do_check_for_dup = optional_default(.false., check_for_dup)

    if(current_verbosity() >= ANAL) then
       call print('Entering test_form_bond, i = '//i//' j = '//j, ANAL)
       call print('use_uniform_cutoff = '//use_uniform_cutoff, ANAL)
    end if

    !Determine what cutoff distance to use
    if (use_uniform_cutoff) then
       use_cutoff = cutoff
    else
       use_cutoff = bond_length(Z(i),Z(j)) * cutoff
    end if

    d = norm(pos(:,j)+(lattice .mult. shift) - pos(:,i))

    if (do_check_for_dup) then
      index = find(this%neighbour1(i)%t, (/ j, shift /)) 
      if (index /= 0) then ! bond is already in table
	if (current_verbosity() >= ANAL) call print('test_form_bond had check_for_dup=T, found bond already in table', ANAL)
	this%neighbour1(i)%t%real(1,index) = d
	return
      endif
    endif

    if(current_verbosity() >= ANAL)  call print('d = '//d, ANAL)
    if (d < use_cutoff) then
       call add_bond(this, pos, lattice, i, j, shift, d)
    end if

    if(current_verbosity() >= ANAL) call print('Leaving test_form_bond', ANAL)

  end subroutine test_form_bond

  !% Test if atom $i$ is no longer neighbour of atom $j$ with shift $s$, and update 'this%connect' as necessary.
  !% Called by 'calc_connect'. The 'shift' vector is added to the position of the $j$ atom
  !% to get the correct image position.
  function test_break_bond(this,cutoff_break, use_uniform_cutoff, Z, pos, lattice, i,j, shift)
    type(Connection), intent(inout) :: this
    real(dp), intent(in) :: cutoff_break
    logical, intent(in) :: use_uniform_cutoff
    integer, intent(in) :: Z(:)
    real(dp), intent(in) :: pos(:,:), lattice(3,3)
    integer,     intent(in)    :: i,j
    integer, intent(in)        :: shift(3)
    logical test_break_bond

    real(dp)                   :: d
    real(dp)                   :: cutoff

    test_break_bond = .false.
    if (i > j) return

    if (.not. associated(this%neighbour1(i)%t) .or. .not. associated(this%neighbour1(j)%t)) return

    if(current_verbosity() >= ANAL) then
       call print('Entering test_break_bond, i = '//i//' j = '//j, ANAL)
       call print('use_uniform_cutoff = '//use_uniform_cutoff // " cutoff_break = "// cutoff_break, ANAL)
    end if

    !Determine what cutoff distance to use
    if (use_uniform_cutoff) then
       cutoff = cutoff_break
    else
       cutoff = bond_length(Z(i),Z(j)) * cutoff_break
    end if

    d = norm(pos(:,j)+(lattice .mult. shift) - pos(:,i))
    if(current_verbosity() >= ANAL)  call print('d = '//d//' cutoff = '//cutoff//' i = '//i//' j = '//j, ANAL)
    if (d > cutoff) then
	if(current_verbosity() >= ANAL) call print('removing bond from tables', ANAL)
       call remove_bond(this, i, j, shift)
       test_break_bond = .true.
    end if

    if(current_verbosity() >= ANAL) call print('Leaving test_break_bond', ANAL)

  end function test_break_bond

  subroutine set_bonds(this, pairs, shifts)
    type(Atoms), intent(inout) :: this
    integer, intent(in) :: pairs(:,:)
    integer, intent(in) :: shifts(:,:)

    integer i

    if (.not.this%connect%initialised) then
       call connection_initialise(this%connect, this%N)
    else
       ! otherwise just wipe the connection table
       call wipe(this%connect)
    end if

    if (size(pairs,1) /= 2) call system_abort("set_bond pairs not a 2xN array")
    if (size(shifts,1) /= 3) call system_abort("set_bond shifts not a 3xN array")
    if (size(pairs,2) /= size(shifts,2)) call system_abort("set_bonds called with mismatching pairs and shifts sizes")

    do i=1, size(pairs,2)
      call add_bond(this%connect, this%pos, this%lattice, pairs(1,i), pairs(2,i), shifts(:,i))
    end do
  end subroutine set_bonds

  subroutine add_bond(this, pos, lattice, i, j, shift, d)
    type(Connection), intent(inout) :: this
    real(dp), intent(in) :: pos(:,:), lattice(3,3)
    integer,     intent(in)    :: i,j
    integer,     intent(in)    :: shift(3)
    real(dp), intent(in), optional :: d

    real(dp) :: dd
    integer :: ii, jj, index

    if (.not.this%initialised) then
       call system_abort("add_bond called on uninitialized connection")
    endif

    if (.not. associated(this%neighbour1(i)%t) .or. .not. associated(this%neighbour1(j)%t)) then
      call system_abort("tried to add_bond for atoms i " // i // " j " // j // " which have associated(neighbour1()%t " // &
	associated(this%neighbour1(i)%t) // " " // associated(this%neighbour1(j)%t) // " one of which is false")
    endif

    if (i > j) then
      ii = j
      jj = i
    else
      ii = i
      jj = j
    endif

    if (present(d)) then
      dd = d
    else
      dd = norm(pos(:,j)+(lattice .mult. shift) - pos(:,i))
    endif

    ! Add full details to neighbour1 for smaller of i and j
    call append(this%neighbour1(ii)%t, (/jj, sign(1,jj-ii)*shift /), (/ dd /))
    if(ii .ne. jj) then		
       index = this%neighbour1(min(ii,jj))%t%N
       ! Put a reference to this in neighbour2 for larger of i and j
       call append(this%neighbour2(jj)%t, (/ ii, index/))
    end if

  end subroutine add_bond

  subroutine remove_bond(this, i, j, shift)
    type(Connection), intent(inout) :: this
    integer,     intent(in)    :: i,j
    integer,     intent(in), optional    :: shift(3)

    integer :: ii, jj, iii, jjj, jjjj, r_index, n_removed

    if (.not. associated(this%neighbour1(i)%t) .or. .not. associated(this%neighbour1(j)%t)) then
      call system_abort("tried to remove_bond for atoms i " // i // " j " // j // " which have associated(neighbour1()%t " // &
	associated(this%neighbour1(i)%t) // " " // associated(this%neighbour1(j)%t) // " one of which is false")
    endif

    if (i > j) then
      ii = j
      jj = i
    else
      ii = i
      jj = j
    endif
    ! now ii <= jj

    r_index = 1
    n_removed = 0
    do while (r_index /= 0)
      ! remove entry from neighbour1(ii)
      if (present(shift)) then
	r_index = find(this%neighbour1(ii)%t, (/ jj, shift /) )
      else
	r_index = find(this%neighbour1(ii)%t, (/ jj, 0, 0, 0 /), (/ .true., .false., .false., .false./) )
      endif
      if (r_index == 0) then
	if (n_removed == 0) then
	  if (present(shift)) then
	    call print("WARNING: remove bond called for i " // i // " j " // j // " shift " // shift // &
		       " couldn't find a bond to remove", ERROR)
	  else
	    call print("WARNING: remove bond called for i " // i // " j " // j // &
		       " couldn't find a bond to remove", ERROR)
	  endif
	endif
      else ! r_index /= 0
	n_removed = n_removed + 1
	call delete(this%neighbour1(ii)%t, r_index, keep_order = .true.)
	! remove entry from neighbour2(jj)
	if (ii /= jj) then
	  call delete(this%neighbour2(jj)%t, (/ ii, r_index /), keep_order = .true.)
	endif
	! renumber other neighbour2 entries
	do iii=r_index, this%neighbour1(ii)%t%N
	  ! jjj is another neighbour of ii
	  jjj = this%neighbour1(ii)%t%int(1,iii)
	  if (jjj > ii) then
	    ! jjjj is r_index in jjj's neighbour2 of pointer back to ii's neighbour1
	    jjjj = find(this%neighbour2(jjj)%t, (/ ii, iii+1 /) )
	    if (jjjj /= 0) then
	      ! decrement reference in jjj's neighbour2 table
	      this%neighbour2(jjj)%t%int(:,jjjj) = (/ ii, iii /)
	    else
	      call system_abort("Couldn't find neighbor to fix neighbour2 of")
	    endif
	  endif
	end do
      end if ! r_index == 0
    end do ! while r_index /= 0

  end subroutine remove_bond


  !% Fast $O(N)$ connectivity calculation routine. It divides the unit cell into similarly shaped subcells,
  !% of sufficient size that sphere of radius 'cutoff' is contained in a subcell, at least in the directions 
  !% in which the unit cell is big enough. In narrow directions, the unit cell is replicated on the fly to collect
  !% neighbours that are images of the stored atoms. For very small unit cells, there is only one subcell, so the routine
  !% is equivalent to the standard $O(N^2)$ method.
  subroutine calc_connect_hysteretic(this, alt_connect, origin, extent, own_neighbour)
    type(Atoms), intent(inout), target           :: this
    type(Connection), intent(inout), target, optional :: alt_connect
    real(dp), optional :: origin(3), extent(3,3)
    logical, optional, intent(in) :: own_neighbour

    integer                              :: cellsNa,cellsNb,cellsNc,i,j,k,i2,j2,k2,i3,j3,k3,i4,j4,k4,n1,n2,atom1,atom2
    integer                              :: cell_image_Na, cell_image_Nb, cell_image_Nc
    real(dp)                             :: cutoff
    integer :: ji, s_ij(3)
    logical my_own_neighbour
    type(Connection), pointer :: use_connect

    if (present(alt_connect)) then
      use_connect => alt_connect
    else
      use_connect => this%connect
    endif

    my_own_neighbour = optional_default(.false., own_neighbour)

    if (this%cutoff < 0.0_dp .or. this%cutoff_break < 0.0_dp) then
       call system_abort('calc_connect: Negative cutoff radius ' // this%cutoff // ' ' // this%cutoff_break )
    end if

    if (this%cutoff > this%cutoff_break) then
       call system_abort('calc_connect: Negative hysteresis cutoff radius formation ' // this%cutoff // ' > breaking ' // this%cutoff_break )
    end if

    if ((this%cutoff .feq. 0.0_dp) .or. (this%cutoff_break .feq. 0.0_dp)) then
      call wipe(use_connect)
      return
    endif

    !Calculate the cutoff value we should use in dividing up the simulation cell
    if (this%use_uniform_cutoff) then
       !The cutoff has been specified by the user, so use that value
       cutoff = this%cutoff
    else
       !Otherwise, find the maximum covalent radius of all the atoms in the structure, double it
       !and multiple by cutoff. This makes sure we can deal with the worst case scenario of big atom 
       !bonded to big atom
       cutoff = 0.0_dp
       do i = 1, this%N
          if (ElementCovRad(this%Z(i)) > cutoff) cutoff = ElementCovRad(this%Z(i))
       end do
       cutoff = (2.0_dp * cutoff) * this%cutoff
    end if

    call print("calc_connect: cutoff calc_connect " // cutoff, NERD)

    if (present(origin) .and. present(extent)) then
      cellsNa = 1
      cellsNb = 1
      cellsNc = 1
    else
      call divide_cell(this%lattice, cutoff, cellsNa, cellsNb, cellsNc)
    endif

    call print("calc_connect: cells_N[abc] " // cellsNa // " " // cellsNb // " " // cellsNc, NERD)

    ! If the lattice has changed, then the cells need de/reallocating
    if ((cellsNa /= use_connect%cellsNa) .or. &
         (cellsNb /= use_connect%cellsNb) .or. &
         (cellsNc /= use_connect%cellsNc)) call connection_cells_finalise(use_connect)

    ! Allocate space for the connection object if needed
    if (present(origin) .and. present(extent)) then
      if (.not.use_connect%initialised) then
	 call connection_initialise(use_connect, this%N, this%pos, this%lattice, this%g, this%cutoff, this%use_uniform_cutoff, origin, extent)
      else
	 call connection_fill(use_connect, this%N, this%pos, this%lattice, this%g, this%cutoff, this%use_uniform_cutoff, origin, extent)
      end if
    else
      if (.not.use_connect%initialised) then
	 call connection_initialise(use_connect, this%N, this%pos, this%lattice, this%g, this%cutoff, this%use_uniform_cutoff)
      else
	 call connection_fill(use_connect, this%N, this%pos, this%lattice, this%g, this%cutoff, this%use_uniform_cutoff)
      end if
    endif

    if (.not.use_connect%cells_initialised) then
      call connection_cells_initialise(use_connect, cellsNa, cellsNb, cellsNc,this%N)
    endif

    ! Partition the atoms into cells
    call partition_atoms(use_connect, this)


    ! figure out how many unit cell images we will need to loop over in each direction
    call fit_box_in_cell(cutoff, cutoff, cutoff, this%lattice, cell_image_Na, cell_image_Nb, cell_image_Nc)
    ! cell_image_N{a,b,c} apply to a box of side 2*cutoff. Since we loop from -cell_image_N to
    ! +cell_image_N we can reduce them as follows

    cell_image_Na = max(1,(cell_image_Na+1)/2)
    cell_image_Nb = max(1,(cell_image_Nb+1)/2)
    cell_image_Nc = max(1,(cell_image_Nc+1)/2)

    call print('calc_connect: image cells '//cell_image_Na//'x'//cell_image_Nb//'x'//cell_image_Nc, NERD)

    ! look for bonds that have been broken, and remove them
    do i=1, this%N
      ji = 1
      do
	if (ji > atoms_n_neighbours(this, i, alt_connect=use_connect)) exit
	j = atoms_neighbour(this, i, ji, shift = s_ij, alt_connect=use_connect)
	if (.not. test_break_bond(use_connect, this%cutoff_break, this%use_uniform_cutoff, &
	  this%Z, this%pos, this%lattice, i, j, s_ij)) then
	  ji = ji + 1 ! we didn't break this bond, so go to next one
	              ! if we did break a bond, ji now points to a different bond, so don't increment it
	endif
      end do
!      do ji=1, atoms_n_neighbours(this, i, alt_connect=use_connect)
!	j = atoms_neighbour(this, i, ji, shift = s_ij, alt_connect=use_connect)
!	call test_break_bond(use_connect, this%cutoff_break, this%use_uniform_cutoff, &
!	  this%Z, this%pos, this%lattice, i, j, s_ij)
!      end do
    end do

    ! Here is the main loop:
    ! Go through each cell and update the connectivity between atoms in this cell and neighbouring cells
    ! N.B. test_form_bond updates both atoms i and j, so only update if i <= j to avoid doubling processing

    ! defaults for cellsNx = 1
    k3 = 1; k4 = 1; j3 = 1; j4 = 1; i3 = 1; i4 = 1
    ! Loop over all cells
    do k = 1, cellsNc
       do j = 1, cellsNb
          do i = 1, cellsNa

             !Loop over atoms in cell(i,j,k)
             do n1 = 1, use_connect%cell(i,j,k)%N

                atom1 = use_connect%cell(i,j,k)%int(1,n1)

                ! Loop over neighbouring cells, applying PBC
                do k2 = -cell_image_Nc, +cell_image_Nc

                   ! the stored cell we are in 
                   if(cellsNc > 1) k3 = mod(k+k2-1+cellsNc,cellsNc)+1 

                   ! the shift we need to get to the cell image
                   k4 = (k+k2-k3)/cellsNc

                   do j2 = -cell_image_Nb, +cell_image_Nb
                      ! the stored cell we are in                 
                      if(cellsNb > 1) j3 = mod(j+j2-1+cellsNb,cellsNb)+1 

                      ! the shift we need to get to the cell image
                      j4 = (j+j2-j3)/cellsNb

                      do i2 = -cell_image_Na, +cell_image_Na
                         ! the stored cell we are in                 
                         if(cellsNa > 1) i3 = mod(i+i2-1+cellsNa,cellsNa)+1 

                         ! the shift we need to get to the cell image
                         i4 = (i+i2-i3)/cellsNa

                         ! The cell we are currently testing atom1 against is cell(i3,j3,k3)
                         ! with shift (i4,j4,k4)
                         ! loop over it's atoms and test connectivity if atom1 < atom2

                         do n2 = 1, use_connect%cell(i3,j3,k3)%N

                            atom2 = use_connect%cell(i3,j3,k3)%int(1,n2)
                            ! omit atom2 < atom1
                            if (atom1 > atom2) cycle
                            ! omit self in the same cell without shift
                            if (.not. my_own_neighbour .and. (atom1 == atom2 .and. & 
                                 (i4==0 .and. j4==0 .and. k4==0) .and. &
                                 (i==i3 .and. j==j3 .and. k==k3))) cycle

                            call test_form_bond(use_connect, this%cutoff, this%use_uniform_cutoff, &
			      this%Z, this%pos, this%lattice, atom1,atom2, (/i4,j4,k4/), .true.)

                         end do ! n2

                      end do ! i2
		   end do ! j2
                end do ! k2

             end do ! n1

          end do ! i
       end do ! j
    end do ! k

  end subroutine calc_connect_hysteretic

  subroutine connection_remove_atom(this, i)
    type(Connection), intent(inout) :: this
    integer, intent(in) :: i

    integer :: ji, j, jj, s_ij(3), n_entries

    n_entries=this%neighbour1(i)%t%N
    do ji=n_entries, 1, -1
      j = this%neighbour1(i)%t%int(1,ji)
      s_ij = this%neighbour1(i)%t%int(2:4,ji)
      call remove_bond(this, i, j, s_ij)
    end do

    n_entries=this%neighbour2(i)%t%N
    do ji=n_entries, 1, -1
      j = this%neighbour2(i)%t%int(1,ji)
      jj = this%neighbour2(i)%t%int(2,ji)
      s_ij = this%neighbour1(j)%t%int(2:4,jj)
      call remove_bond(this, j, i, s_ij)
    end do

  end subroutine connection_remove_atom

  !% Fast $O(N)$ connectivity calculation routine. It divides the unit cell into similarly shaped subcells,
  !% of sufficient size that sphere of radius 'cutoff' is contained in a subcell, at least in the directions 
  !% in which the unit cell is big enough. In narrow directions, the unit cell is replicated on the fly to collect
  !% neighbours that are images of the stored atoms. For very small unit cells, there is only one subcell, so the routine
  !% is equivalent to the standard $O(N^2)$ method.
  subroutine calc_connect(this, own_neighbour)
    type(Atoms), intent(inout)           :: this
    logical, optional, intent(in) :: own_neighbour

    integer                              :: cellsNa,cellsNb,cellsNc,i,j,k,i2,j2,k2,i3,j3,k3,i4,j4,k4,n1,n2,atom1,atom2
    integer                              :: cell_image_Na, cell_image_Nb, cell_image_Nc
    real(dp)                             :: cutoff
    logical my_own_neighbour

    my_own_neighbour = optional_default(.false., own_neighbour)

    my_own_neighbour = optional_default(.false., own_neighbour)

    if (this%cutoff < 0.0_dp .or. this%cutoff_break < 0.0_dp) then
       call system_abort('calc_connect: Negative cutoff radius ' // this%cutoff // ' ' // this%cutoff_break )
    end if

    if ((this%cutoff .feq. 0.0_dp) .or. (this%cutoff_break .feq. 0.0_dp)) then
      call wipe(this%connect)
      return
    endif

    !Calculate the cutoff value we should use in dividing up the simulation cell
    if (this%use_uniform_cutoff) then
       !The cutoff has been specified by the user, so use that value
       cutoff = this%cutoff
    else
       !Otherwise, find the maximum covalent radius of all the atoms in the structure, double it
       !and multiple by cutoff. This makes sure we can deal with the worst case scenario of big atom 
       !bonded to big atom
       cutoff = 0.0_dp
       do i = 1, this%N
          if (ElementCovRad(this%Z(i)) > cutoff) cutoff = ElementCovRad(this%Z(i))
       end do
       cutoff = (2.0_dp * cutoff) * this%cutoff
    end if

    call print("calc_connect: cutoff calc_connect " // cutoff, NERD)

    call divide_cell(this%lattice, cutoff, cellsNa, cellsNb, cellsNc)

    call print("calc_connect: cells_N[abc] " // cellsNa // " " // cellsNb // " " // cellsNc, NERD)

    ! If the lattice has changed, then the cells need de/reallocating
    if ((cellsNa /= this%connect%cellsNa) .or. &
         (cellsNb /= this%connect%cellsNb) .or. &
         (cellsNc /= this%connect%cellsNc)) call connection_finalise(this%connect)

    ! Allocate space for the connection object if needed
    if (.not.this%connect%initialised) then
       call connection_initialise(this%connect, this%N, this%pos, this%lattice, this%g, this%cutoff, this%use_uniform_cutoff)
    else
       ! otherwise just wipe the connection table
       call wipe(this%connect)
    end if

    if (.not.this%connect%cells_initialised) &
         call connection_cells_initialise(this%connect, cellsNa, cellsNb, cellsNc,this%N)

    ! Partition the atoms into cells
    call partition_atoms(this%connect, this)

    ! figure out how many unit cell images we will need to loop over in each direction
    call fit_box_in_cell(cutoff, cutoff, cutoff, this%lattice, cell_image_Na, cell_image_Nb, cell_image_Nc)
    ! cell_image_N{a,b,c} apply to a box of side 2*cutoff. Since we loop from -cell_image_N to
    ! +cell_image_N we can reduce them as follows

    cell_image_Na = max(1,(cell_image_Na+1)/2)
    cell_image_Nb = max(1,(cell_image_Nb+1)/2)
    cell_image_Nc = max(1,(cell_image_Nc+1)/2)

    call print('calc_connect: image cells '//cell_image_Na//'x'//cell_image_Nb//'x'//cell_image_Nc, NERD)

    ! Here is the main loop:
    ! Go through each cell and update the connectivity between atoms in this cell and neighbouring cells
    ! N.B. test_form_bond updates both atoms i and j, so only update if i <= j to avoid doubling processing

    ! defaults for cellsNx = 1
    k3 = 1; k4 = 1; j3 = 1; j4 = 1; i3 = 1; i4 = 1
    ! Loop over all cells
    do k = 1, cellsNc
       do j = 1, cellsNb
          do i = 1, cellsNa

             !Loop over atoms in cell(i,j,k)
             do n1 = 1, this%connect%cell(i,j,k)%N

                atom1 = this%connect%cell(i,j,k)%int(1,n1)

                ! Loop over neighbouring cells, applying PBC
                do k2 = -cell_image_Nc, +cell_image_Nc

                   ! the stored cell we are in 
                   if(cellsNc > 1) k3 = mod(k+k2-1+cellsNc,cellsNc)+1 

                   ! the shift we need to get to the cell image
                   k4 = (k+k2-k3)/cellsNc

                   do j2 = -cell_image_Nb, +cell_image_Nb
                      ! the stored cell we are in                 
                      if(cellsNb > 1) j3 = mod(j+j2-1+cellsNb,cellsNb)+1 

                      ! the shift we need to get to the cell image
                      j4 = (j+j2-j3)/cellsNb

                      do i2 = -cell_image_Na, +cell_image_Na
                         ! the stored cell we are in                 
                         if(cellsNa > 1) i3 = mod(i+i2-1+cellsNa,cellsNa)+1 

                         ! the shift we need to get to the cell image
                         i4 = (i+i2-i3)/cellsNa

                         ! The cell we are currently testing atom1 against is cell(i3,j3,k3)
                         ! with shift (i4,j4,k4)
                         ! loop over it's atoms and test connectivity if atom1 < atom2

                         do n2 = 1, this%connect%cell(i3,j3,k3)%N

                            atom2 = this%connect%cell(i3,j3,k3)%int(1,n2)
                            ! omit atom2 < atom1
                            if (atom1 > atom2) cycle
                            ! omit self in the same cell without shift
                            if (.not. my_own_neighbour .and. (atom1 == atom2 .and. & 
                                 (i4==0 .and. j4==0 .and. k4==0) .and. &
                                 (i==i3 .and. j==j3 .and. k==k3))) cycle

                            call test_form_bond(this%connect,this%cutoff, this%use_uniform_cutoff, &
			      this%Z, this%pos, this%lattice, atom1,atom2, (/i4,j4,k4/))

                         end do

                      end do
                   end do
                end do

             end do

          end do
       end do
    end do

  end subroutine calc_connect

  !
  !% Spatially partition the atoms into cells. The number of cells in each dimension must already be
  !% set (cellsNa,b,c). Pre-wiping of the cells can be skipped (e.g. if they are already empty).
  !
  subroutine partition_atoms(this, at, dont_wipe)

    type(Connection), intent(inout) :: this
    type(Atoms), intent(inout) :: at
    logical, optional, intent(in)    :: dont_wipe

    logical                          :: my_dont_wipe
    integer                          :: i,j,k,n

    ! Check inputs
    if (.not.this%cells_initialised) call system_abort('Partition_Atoms: Cells have not been initialised')
    my_dont_wipe = .false.
    if (present(dont_wipe)) my_dont_wipe = dont_wipe

    ! Wipe the cells
    if (.not.my_dont_wipe) call wipe_cells(this)

    ! Make sure all atomic positions are within the cell
    call map_into_cell(at)

    do n = 1, at%N
       if (.not. associated(this%neighbour1(n)%t)) cycle ! not in active subregion

       call cell_of_pos(this, at%g, at%pos(:,n), i, j, k)
       !Add the atom to this cell
       call append( this%cell(i,j,k), (/n/) )
    end do

  end subroutine partition_atoms

  subroutine cell_of_pos(this, g, pos, i, j, k)
    type(Connection), intent(in) :: this
    real(dp), intent(in) :: g(3,3)
    real(dp), intent(in) :: pos(3)
    integer, intent(out) :: i, j, k

    real(dp) :: t(3)

     t = g .mult. pos
     i = floor(real(this%cellsNa,dp) * (t(1)+0.5_dp)) + 1
     j = floor(real(this%cellsNb,dp) * (t(2)+0.5_dp)) + 1
     k = floor(real(this%cellsNc,dp) * (t(3)+0.5_dp)) + 1

     ! Very small numerical errors in the lattice inverse can lead to bad i,j,k values.
     ! Test for this:
     if (i < 1) then
	i = 1
     else if (i > this%cellsNa) then
	i = this%cellsNa
     end if

     if (j < 1) then
	j = 1
     else if (j > this%cellsNb) then
	j = this%cellsNb
     end if

     if (k < 1) then
	k = 1
     else if (k > this%cellsNc) then
	k = this%cellsNc
     end if

    end subroutine cell_of_pos


   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   !
   ! Geometry procedures
   !
   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

   !
   !% Given a simulation cell defined by lattice vectors, how many
   !% times can the cell be divided along the lattice vectors into
   !% subcells such that a sphere of radius 'cutoff' with centre in
   !% one subcell does not spill out of the surrounding $3 \times 3$ subcell block?
   !
   subroutine divide_cell(lattice,cutoff,Na,Nb,Nc)

      real(dp), dimension(3,3), intent(in)  :: lattice  !% Box defined by lattice vectors
      real(dp),                 intent(in)  :: cutoff   !% Radius of sphere
      integer,                  intent(out) :: Na,Nb,Nc !% Number of supercells required along $x$, $y$ and $z$
      real(dp)                              :: cellVol
      real(dp), dimension(3)                :: a, b, c

      a = lattice(:,1); b = lattice(:,2); c = lattice(:,3)
      cellVol = abs( scalar_triple_product(a,b,c) )

      Na = max(1,int( cellVol / (cutoff * norm(b .cross. c) ) )) ! round down to nearest int >= 1
      Nb = max(1,int( cellVol / (cutoff * norm(c .cross. a) ) ))
      Nc = max(1,int( cellVol / (cutoff * norm(a .cross. b) ) ))

      call print('divide_cell: '//Na//'x'//Nb//'x'//Nc//' cells.', NERD)

   end subroutine divide_cell

   !
   !% Returns the maximum cutoff radius for 'calc_connect', given the lattice if we want to avoid image neghbours
   !
   function max_cutoff(lattice)

     real(dp), dimension(3,3), intent(in) :: lattice
     real(dp)                             :: Max_Cutoff
     real(dp), dimension(3)               :: a,b,c
     real(dp)                             :: cellVol, ra, rb, rc

     a = lattice(:,1); b = lattice(:,2); c = lattice(:,3)
     cellVol = abs( scalar_triple_product(a,b,c) )

     if(cellVol == 0.0_dp) call system_abort("Max_cutoff(): cell volume is exactly 0.0!")

     ra = cellVol / norm(b .cross. c)
     rb = cellVol / norm(c .cross. a)
     rc = cellVol / norm(a .cross. b)

     Max_Cutoff = 0.5_dp * min(ra,rb,rc)

   end function max_cutoff



   !
   !% Returns the (unsigned) volume of the simulation cell of 'this'
   !
   function atoms_cell_volume(this)
     type(Atoms), intent(in) :: this
     real(dp)                :: atoms_Cell_Volume

     atoms_cell_volume = cell_volume(this%lattice)
   end function atoms_cell_volume

   !% Returns the (unsigned) volume of the simulation cell of lattice
   function lattice_cell_volume(lattice)
     real(dp), intent(in) :: lattice(3,3)
     real(dp)                :: lattice_Cell_Volume
 
     real(dp), dimension(3)  :: a, b, c

     a = lattice(:,1)
     b = lattice(:,2)
     c = lattice(:,3)

     lattice_Cell_Volume = abs(Scalar_Triple_Product(a,b,c))

   end function lattice_cell_volume


   !
   ! Fit_Box_In_Cell
   !
   !% Given an orthogonal box, oriented along the cartesian axes with lengths '2*rx', '2*ry' and '2*rz'
   !% and centred on the origin, what parameters must we pass to supercell to make a system big 
   !% enough from our original cell defined by lattice for the box to fit inside?
   !%
   !% The required supercell parameters are returned in 'Na', 'Nb', 'Nc' respectively. If, e.g. 'Na = 1'
   !% then we don\'t need to supercell in the a direction. 'Na > 1' means we must supercell by 'Na' times
   !% in the $x$ direction.
   !
   subroutine fit_box_in_cell(rx,ry,rz,lattice,Na,Nb,Nc)

     real(dp),                 intent(in)  :: rx, ry, rz
     real(dp), dimension(3,3), intent(in)  :: lattice
     integer,                  intent(out) :: Na, Nb, Nc
     !local variables
     real(dp), dimension(3,3)              :: g                 !inverse of lattice
     real(dp)                              :: maxa, maxb, maxc  !the current maximum in the a,b,c directions
     real(dp), dimension(3)                :: lattice_coords    !the coordinates of a corner in terms of a,b,c
     integer                               :: i,j,k
     !This subroutine works by calculating the coordinates of the corners of the box
     !in terms of the lattice vectors, taking the maximum in the lattice directions and
     !rounding up to the nearest integer.

     call matrix3x3_inverse(lattice,g)
     maxa = 0.0_dp; maxb=0.0_dp; maxc=0.0_dp

     !For the time being, do all eight corners. This is most likely overkill, but it's still quick.
     do k = -1, 1, 2
        do j = -1, 1, 2
           do i = -1, 1, 2
              lattice_coords = g .mult. (/i*rx,j*ry,k*rz/)
              if (abs(lattice_coords(1)) > maxa) maxa = abs(lattice_coords(1))
              if (abs(lattice_coords(2)) > maxb) maxb = abs(lattice_coords(2))
              if (abs(lattice_coords(3)) > maxc) maxc = abs(lattice_coords(3))
           end do
        end do
     end do

     Na = ceiling(2.0_dp*maxa)
     Nb = ceiling(2.0_dp*maxb)
     Nc = ceiling(2.0_dp*maxc)

   end subroutine fit_box_in_cell

   !
   ! Make_Lattice
   !
   !% Make a matrix of lattice vectors from the lengths 'a','b','c'
   !% and the angles 'alpha', 'beta' and 'gamma'.
   !% One length must be supplied. Any missing angle is assumed to be 90 degrees
   !% and any missing length is assumed to be 'a'.
   !% The vectors are created in a right-handed order.
   !
   function make_lattice(a,b,c,alpha,beta,gamma) result(lattice)

     real(dp),           intent(in) :: a
     real(dp), optional, intent(in) :: b,c
     real(dp), optional, intent(in) :: alpha,beta,gamma
     real(dp), dimension(3,3)       :: lattice
     real(dp)                       :: my_b, my_c,            &
                                       cos_alpha, cos2_alpha, &
                                       cos_beta,  cos2_beta,  &
                                       cos_gamma, cos2_gamma, &
                                       sin_gamma, sin2_gamma, &
                                       my_alpha, my_beta,my_gamma
     my_b = a
     my_c = a

     if (present(b)) my_b = b
     if (present(c)) my_c = c

     my_alpha = PI/2.0_dp
     my_beta  = PI/2.0_dp
     my_gamma = PI/2.0_dp

     if (present(alpha)) my_alpha = alpha
     if (present(beta))  my_beta  = beta
     if (present(gamma)) my_gamma = gamma

     if ( (my_alpha <= 0.0_dp) .or. (my_alpha <= 0.0_dp) .or. (my_gamma <= 0.0_dp) ) &
          call system_abort('Make_Lattice: Negative angles are not permitted')

     if ( (my_alpha + my_beta) < my_gamma ) call system_abort('Make_Lattice: alpha + beta < gamma')
     if ( (my_beta + my_gamma) < my_alpha ) call system_abort('Make_Lattice: beta + gamma < alpha')
     if ( (my_gamma + my_alpha) < my_beta ) call system_abort('Make_Lattice: gamma + alpha < beta')

     cos_alpha = cos(my_alpha); cos2_alpha = cos_alpha*cos_alpha
     cos_beta  = cos(my_beta);  cos2_beta  = cos_beta *cos_beta
     cos_gamma = cos(my_gamma); cos2_gamma = cos_gamma*cos_gamma
     sin_gamma = sin(my_gamma); sin2_gamma = sin_gamma*sin_gamma

     ! a
     lattice(1,1) = a
     lattice(2,1) = 0.0_dp
     lattice(3,1) = 0.0_dp
     ! b
     lattice(1,2) = my_b * cos_gamma
     lattice(2,2) = my_b * sin_gamma
     lattice(3,2) = 0.0_dp
     ! c
     lattice(1,3) = my_c * cos_beta
     lattice(2,3) = my_c * (cos_alpha - cos_beta*cos_gamma) / sin_gamma
     lattice(3,3) = my_c * sqrt(1.0_dp - (cos2_alpha + cos2_beta - 2.0_dp*cos_alpha*cos_beta*cos_gamma)/ sin2_gamma)

   end function make_lattice

   !
   !% Opposite of Make_Lattice.
   !% Given a lattice, return a,b,c,alpha,beta and gamma (if needed)
   !
   subroutine get_lattice_params(lattice,a,b,c,alpha,beta,gamma)

     real(dp), dimension(3,3), intent(in)  :: lattice
     real(dp), optional,       intent(out) :: a,b,c,alpha,beta,gamma
     !local variables
     real(dp), dimension(3)                :: a_1,a_2,a_3

     a_1 = lattice(:,1)
     a_2 = lattice(:,2)
     a_3 = lattice(:,3)

     if (present(a)) a = norm(a_1)
     if (present(b)) b = norm(a_2)
     if (present(c)) c = norm(a_3)

     if (present(alpha)) alpha = acos( (a_2 .dot. a_3) / (norm(a_2)*norm(a_3)) )
     if (present(beta))  beta  = acos( (a_3 .dot. a_1) / (norm(a_3)*norm(a_1)) )
     if (present(gamma)) gamma = acos( (a_1 .dot. a_2) / (norm(a_1)*norm(a_2)) )

   end subroutine get_lattice_params

   !
   ! Centre_Of_Mass
   !
   !% Calculate the centre of mass of an atoms object, using the closest images to the origin atom,
   !% or first atom if this is not specified.
   !% If an index_list is present, just calculate it for that subset of atoms (then the origin atom is
   !% the first in this list unless it is specified separately).
   !%
   !% Note: Because the origin can be specified separately it need not be one of the atoms in the 
   !% calculation.
   function centre_of_mass(at,index_list,origin) result(CoM)

     type(atoms),                      intent(in) :: at
     integer,                optional, intent(in) :: origin
     integer,  dimension(:), optional, intent(in) :: index_list
     real(dp), dimension(3)                       :: CoM
     !local variables
     integer                                      :: i, my_origin
     real(dp)                                     :: M_Tot

     if (.not. has_property(at, 'mass')) &
          call system_abort('center_of_mass: Atoms has no mass property')

     if (present(origin)) then
        if (origin > at%N .or. origin < 1) call system_abort('Centre_Of_Mass: Invalid origin atom')
        my_origin = origin
     else
        if (present(index_list)) then
           my_origin = index_list(1)
        else
           my_origin = 1
        end if
     end if

     CoM = 0.0_dp
     M_Tot = 0.0_dp

     if (present(index_list)) then

        do i = 1, size(index_list)
           if (index_list(i) > at%N .or. index_list(i) < 1) &
                call system_abort('Centre_Of_Mass: Invalid atom in index_list')
           CoM = CoM + at%mass(index_list(i)) * diff_min_image(at,my_origin,index_list(i))
           M_Tot = M_Tot + at%mass(index_list(i))
        end do

     else

        do i = 1, at%N
           CoM = CoM + at%mass(i) * diff_min_image(at,my_origin,i)
           M_Tot = M_Tot + at%mass(i)
        end do

     end if

     CoM = CoM / M_Tot
     CoM = CoM + at%pos(:,my_origin)

   end function centre_of_mass

   !
   ! Directionality ellipsoid:
   !
   !% Given an origin atom and a list of other atoms, give information as to whether the other atoms
   !% are distributed roughly linearly, planar or spherically around the origin atom.
   !%
   !% The most notable use is to check that the splines in adjustable potential will be able to reproduce
   !% a randomly oriented force difference well.
   !%
   !% The information returned is the set of eigenvectors and associated eigenvalues of the directionality
   !% ellipsoid. One large e-value suggests roughly linear clustering, two similar and one small e-values suggest
   !% a planar distribution, while three similar e-values suggests almost spherical distribution (when copies of
   !% the atoms reflected through the origin atom are also considered).
   !% 
   !% To acheive a more spherical distribution, atoms along the e-vector(s) with the smallest e-value(s) should be
   !% added to the index list (See 'CosAngle_To_Line' below).
   !%
   !% The matrix which is diagonalised is an average of the outer products of the unit vectors from the origin
   !% atom to the other atoms.
   !%
   !% An outer product has 1 eigenvector which is the vector it was constructed from with
   !% eigenvalue 1 and the other eigenvectors have eigenvalue 0.
   !%
   !% The eigenvalues of the averaged matrix sum to 1.
   !% 
   subroutine directionality(this,origin,list,evalues,evectors,method)

     type(Atoms),                      intent(in)  :: this     !% The input atoms structure
     integer,                          intent(in)  :: origin   !% The origin atom
     type(table),                      intent(in)  :: list     !% Indices and shifts of the other atoms relative to origin
     real(dp), dimension(3),           intent(out) :: evalues  !% Eigenvalues of the directionality matrix
     real(dp), dimension(3,3),         intent(out) :: evectors !% Eigenvectors of the directionality matrix
     integer, optional,                intent(in)  :: method   !% 'METHOD = 1' Directionality ellipsoid method \\
                                                               !% 'METHOD = 2' Singular Value Decomposition method (default)

     !local variables
     integer                                       :: i, j, k, l, n, my_method, lwork, info, jshift(3)
     real(dp), dimension(3)                        :: r_ij,rhat_ij
     real(dp), dimension(3,3)                      :: A, B 
     real(dp), allocatable, dimension(:,:)         :: vectors, u
     real(dp), allocatable, dimension(:)           :: work

     my_method = 2

     if (present(method)) then
        if (method < 1 .or. method > 2) then
           write(line,'(a,i0,a)')'Directionality: Method = ',method,' does not exist'
           call system_abort(line)
        end if
        my_method = method
     end if

     if (list%intsize /= 4) &
          call system_abort('Directionality: list must have 4 int columns for indices and shifts')
     if (list%N == 0) call system_abort('Directionality: list table has no entries')

     i = origin

     select case(my_method)

     case(1)        ! **** Directionality Ellipsoid method ****

        B = 0.0_dp

        !Construct the directionality matrix
        do n = 1, list%N
           j = list%int(1,n)
           jshift = list%int(2:4,n)
           if (j > this%N) then
              write(line,'(a,i0,a,i0,a)')'Directionality: Atom ',j,' is out of range (',this%N,')'
              call system_abort(line)
           end if
           r_ij = diff(this,i,j,jshift)
           rhat_ij = r_ij / norm(r_ij)
           forall(k=1:3,l=1:3) A(k,l) = rhat_ij(k)*rhat_ij(l)
           B = B + A
        end do

        B = B / real(list%N,dp)

        !Find eigenvalues/eigenvectors of the directionality matrix
        call diagonalise(B,evalues,evectors)

     case(2)       ! **** Singular Value Decomposition method ****

        lwork = 2 * max(3*min(3,list%N)+max(3,list%N),5*min(3,list%N))

        allocate(work(lwork), vectors(list%N,3), u(list%N,3))

        !Fill 'vectors' with unit vectors
        do n = 1, list%N
           j = list%int(1,n)
           jshift = list%int(2:4,n)
           if (j > this%N) then
              write(line,'(a,i0,a,i0,a)')'Directionality: Atom ',j,' is out of range (',this%N,')'
              call system_abort(line)
           end if
           r_ij = diff(this,i,j,jshift)
           vectors(n,:) = r_ij / norm(r_ij)
        end do   

        call dgesvd('N','A',list%N,3,vectors,list%N,evalues,u,list%N,evectors,3,work,lwork,info)
        !No left singular vectors, all right singular vectors

        if (info/=0) then
           if (info < 0) then
              write(line,'(a,i0,a)')'Directionality: Problem with argument ',-info,' passed to DGESVD'
              call system_abort(line)
           else
              call system_abort('Directionality: DBDSQR (called from DGESVD) did not converge')
           end if
        end if

        deallocate(work, vectors, u)

        evectors = transpose(evectors)

     end select

   end subroutine directionality

   !
   ! CosAngle_To_Line
   !
   !% For use with the 'Directionality' routine above.
   !% Given an atom ('atom') and a direction ('dir') return the absolute value of the cosine of the angle
   !% between the the line running through 'atom' in direction 'dir' and the line
   !% between 'atom' and 'test_atom'
   !%
   !% The idea is that this will be called with the origin atom from 'Directionality' as 'atom',
   !% the eigenvector associated with the smallest eigenvalue as 'dir' and a potentially new 
   !% atom to connect to 'atom' with a spline as 'test_atom'.
   !% 
   !% If the result is close to 1 then accept the 'test_atom', and reject if close to zero

   function cosangle_to_line(this,atom,dir,test_atom)

     type(atoms),            intent(in) :: this
     integer,                intent(in) :: atom
     real(dp), dimension(3), intent(in) :: dir
     integer,                intent(in) :: test_atom
     real(dp)                           :: CosAngle_To_Line
     !local variables
     real(dp), dimension(3)             :: r_ab

     !Sanity checks
     if (atom > this%N) then
        write(line,'(a,i0,a,i0,a)')'CosAngle_To_Line: Atom ',atom,' out of range (',this%N,')'
        call system_abort(line)
     end if

     if (test_atom > this%N) then
        write(line,'(a,i0,a,i0,a)')'CosAngle_To_Line: Test atom ',test_atom,' out of range (',this%N,')'
        call system_abort(line)
     end if

     if (norm2(dir) .feq. 0.0_dp) call system_abort('CosAngle_To_Line: A non-zero direction is required')

     r_ab = diff_min_image(this,atom,test_atom)

     CosAngle_To_Line = abs(dir .dot. r_ab) / (norm(dir)*norm(r_ab))

   end function cosangle_to_line

   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   !
   ! I/O procedures  
   !
   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

   subroutine atoms_print(this,file,properties)
      type(Atoms),    intent(inout)             :: this
      type(Inoutput), intent(inout),optional, target :: file
      character(*), optional, intent(in)     :: properties 
      type(Inoutput), pointer :: my_out

      if(.not.this%initialised) call system_abort('Atoms_Print: Atoms structure not initialised')

      if(current_verbosity() <= SILENT) return ! for efficiency

      if (present(file)) then
         my_out => file
      else
         my_out => mainlog
      end if

      call print('Atoms Structure: ', NORMAL, my_out)
      call print('Number of atoms = '//this%N, NORMAL, my_out)

      if(this%use_uniform_cutoff) then
         call print('Bond-formation cutoff radius = '//this%cutoff//' Angstroms', NORMAL, my_out)
         call print('Bond-breaking cutoff radius = '//this%cutoff_break//' Angstroms', NORMAL, my_out)
      else
         call print('Bond-formation cutoff radius = '//this%cutoff//' *bond_length', NORMAL, my_out)
         call print('Bond-breaking cutoff radius = '//this%cutoff_break//' *bond_length', NORMAL, my_out)
      end if

      call print('Lattice vectors:', NORMAL, my_out)
      call print('a = ('//this%lattice(:,1)//')', NORMAL, my_out)
      call print('b = ('//this%lattice(:,2)//')', NORMAL, my_out)
      call print('c = ('//this%lattice(:,3)//')', NORMAL, my_out)

      if (present(properties)) then
	 call print_xyz(this, my_out, properties=properties,human_readable=.true., real_format='f10.5')
      else
	 call print_xyz(this, my_out, all_properties=.true.,human_readable=.true., real_format='f10.5')
      end if

      if (this%connect%initialised) then
	call verbosity_push_decrement()
	call connection_print(this%connect, my_out)
	call verbosity_pop()
      end if

      call print('',NORMAL, my_out)
   end subroutine atoms_print


   subroutine atoms_print_xyz(this, xyzfile, comment, properties, all_properties, human_readable, real_format, mask)

      type(Atoms),            intent(inout)    :: this     !% Atoms object to print
      type(Inoutput),         intent(inout) :: xyzfile  !% Inoutput object to write to
      character(*), optional, intent(in)    :: comment  !% Comment line (line #2 of xyz file)
      character(*), optional, intent(in) :: properties  !% Colon separated list of properties from 'this%data' 
                                                        !% to be printed. If not specified, we print only the 
                                                        !% atomic positions, unless 'all_properties' is true.
      logical,      optional, intent(in)    :: all_properties !% Should we print all properties (default is false)
      logical,      optional, intent(in)    :: human_readable !% If set to true, pretty-print table of 
                                                              !% atomic properties.
      character(len=*), optional, intent(in) :: real_format
      logical, optional, intent(in) :: mask(:)

      logical             :: got_species
      integer             :: i,lookup(3), n_properties
      logical             :: do_all_properties, do_human_readable, dummy
      character(len=2048) :: my_properties, prop_names
      character(len=1024) :: tmp_properties(this%properties%N)
      type(Dictionary)    :: use_properties
      character(len=1024) :: my_real_format

      my_real_format = optional_default('f12.5',real_format)

      do_all_properties = .false.
      if (present(all_properties)) do_all_properties = all_properties

      do_human_readable = .false.
      if (present(human_readable)) do_human_readable = human_readable

      if (.not.this%initialised) &
           call system_abort('Atoms_Print_xyz: Atoms structure not initialised')

      ! Try to fill in species from Z
      if (all(this%species == repeat(' ',TABLE_STRING_LENGTH))) then
         call print('atoms_print_xyz: filling in species from Z',VERBOSE)
         do i=1,this%N
            this%species(i) = ElementName(this%Z(i))
         end do
      end if

      if (.not. present(properties)) then

         ! Print all properties, starting with species and position
         if (do_all_properties) then
            use_properties = this%properties

            ! Move 'species' to start of list
            call swap(use_properties, 'species', use_properties%keys(1))

            ! Move 'pos' to second in list
            call swap(use_properties, 'pos', use_properties%keys(2))
         else
            ! just print the species and positions
            use_properties = subset(this%properties,(/'species','pos    '/))
         end if

      else
         ! properties is colon separated list of properties

         my_properties = properties

         ! Backwards compatibility: prepend species tag if it's missing
         call parse_string(my_properties, ':', tmp_properties, n_properties)
	 got_species = .false.
	 do i=1, n_properties
	  if (tmp_properties(i) == 'species') then
	    if (i /= 1) then
	      tmp_properties(i) = tmp_properties(1)
	      tmp_properties(1) = 'species'
	    endif
	    got_species = .true.
	  endif
	 end do
         
	 if (got_species) then
	   my_properties = trim(tmp_properties(1))
	   do i=2, n_properties
	     my_properties = trim(my_properties) // ':' // tmp_properties(i)
	   end do
	 else
	   my_properties = 'species:'//trim(my_properties)
	 endif

         call parse_string(my_properties, ':', tmp_properties, n_properties)

         use_properties = subset(this%properties, tmp_properties(1:n_properties))

      end if

      ! Check 1st property looks like species and 
      ! 2nd property has three real columns (i.e. looks like a position)
      if (.not. do_human_readable) then
         dummy = get_value(use_properties, use_properties%keys(1), lookup)
         if (lookup(1) /= PROPERTY_STR .or. (lookup(3)-lookup(2)+1) /= 1) &
           call system_abort('Atoms_print_xyz: first property must be a species string')

         dummy = get_value(use_properties, use_properties%keys(2), lookup)
         if (lookup(1) /= PROPERTY_REAL .or. (lookup(3)-lookup(2)+1) /= 3) &
           call system_abort('Atoms_print_xyz: seccond property must be a real 3-vector, e.g. pos or avgpos')
      end if

      ! Get list of properties to print
      prop_names = dict_prop_names_string(use_properties,.true.)

      if (.not. do_human_readable) then
         ! First write the number of atoms
	 if (present(mask)) then
	    call print(""//count(mask), file=xyzfile)
	 else
	    call print(""//this%N, file=xyzfile)
	 endif

         ! Then the lattice, property names and other parameters
         call set_value(this%params, 'Lattice', reshape(this%lattice, (/9/)))
         call set_value(this%params, 'Properties', prop_names)
         line = write_string(this%params, real_format='f18.6')

         ! Finally any optional comment
         if (present(comment)) then 
            line = trim(line)//' '//trim(comment)
         end if

         call print(line,file=xyzfile)
      else
         call print('Properties printed are:')
         do i=1,use_properties%N
            dummy = get_value(use_properties, use_properties%keys(i), lookup)
            call print(i//' '//trim(use_properties%keys(i))//' '//(lookup(3)-lookup(2)+1)//' columns')
         end do
      end if

      call print(this%data,file=xyzfile,real_format=my_real_format,&
           properties=use_properties, mask=mask)

      call finalise(use_properties)

    end subroutine atoms_print_xyz


   function prop_names_string(this)
     type(Atoms), intent(in) :: this
     character(len=2048) :: prop_names_string

     prop_names_string=dict_prop_names_string(this%properties, .false.)
   end function prop_names_string

   function dict_prop_names_string(this,with_types)
     type(Dictionary), intent(in) :: this
     logical, intent(in), optional :: with_types
     character(len=2048) :: dict_prop_names_string

     character(len=1) :: prop_types(4) = (/'I', 'R', 'S', 'L'/)
     character(len=1024) :: tmp
     integer :: i
     integer :: lookup(3)
     logical :: dummy
     logical :: my_with_types

     my_with_types = optional_default(.true., with_types)

     dict_prop_names_string = ""
     do i=1,this%N
	dummy = get_value(this, this%keys(i), lookup)
	if (my_with_types) then
	  write(tmp,'(i0)') lookup(3)-lookup(2)+1 ! Number of columns for this property
	  dict_prop_names_string=trim(dict_prop_names_string)//trim(this%keys(i))//':'// &
				 prop_types(lookup(1))//':'//trim(tmp)//':'
	else
	  dict_prop_names_string=trim(dict_prop_names_string)//trim(this%keys(i))//':'
	endif
     end do
     ! Remove trailing ':'
     dict_prop_names_string = dict_prop_names_string(1:len_trim(dict_prop_names_string)-1)
   end function dict_prop_names_string

   subroutine atoms_print_xyz_filename(this, xyzfilename, comment, properties, all_properties, human_readable, real_format)

      type(Atoms),            intent(inout)    :: this     
      character(*),           intent(in)    :: xyzfilename
      character(*), optional, intent(in)    :: comment
      character(*), optional, intent(in)    :: properties 
      logical,      optional, intent(in)    :: all_properties
      logical,      optional, intent(in)    :: human_readable
      character(*), optional, intent(in)    :: real_format

      type(Inoutput)  :: file

      call initialise(file, xyzfilename, action=OUTPUT)
      call atoms_print_xyz(this, file, comment, properties, all_properties, human_readable, real_format)
      call finalise(file)

    end subroutine atoms_print_xyz_filename

    subroutine atoms_read_xyz_filename(this, xyzfilename, comment, status, properties, mpi_comm)
      type(Atoms), intent(inout) :: this
      character(*), intent(in) :: xyzfilename
      character(len=*), optional, intent(out) :: comment
      integer, optional, intent(out) :: status
      character(len=*), optional, intent(out) :: properties 
      integer, optional, intent(in) :: mpi_comm

      type(Inoutput)  :: file
      logical :: master_only

      if (present(mpi_comm)) then
	master_only = .true.
      else
	master_only = .false.
      end if
      call initialise(file, xyzfilename, master_only=master_only)
      call read_xyz(this, file, comment, status, properties, mpi_comm)
      call finalise(file)
    end subroutine atoms_read_xyz_filename

    subroutine atoms_read_xyz_inoutput(this, xyzfile, comment, status, properties, mpi_comm)
      type(Atoms),    intent(inout)           :: this
      type(Inoutput), intent(inout)           :: xyzfile
      character(len=*), optional, intent(out) :: comment
      integer, optional, intent(out)          :: status
      character(len=*), optional, intent(out) :: properties 
      integer, optional, intent(in)           :: mpi_comm

      type(extendable_str) :: es

      if (present(mpi_comm)) then
	call initialise(es)
	call read(es, xyzfile%unit, mpi_comm=mpi_comm, keep_lf=.true.)
	call atoms_read_xyz_extendable_str(this, es, comment, status, properties)
	call finalise(es)
      else
	call atoms_read_xyz(this, xyz_inoutput=xyzfile, comment=comment, status=status, properties=properties)
      endif
    end subroutine atoms_read_xyz_inoutput

    subroutine atoms_read_xyz_extendable_str(this, xyzfile_es, comment, status, properties)
      type(Atoms),    intent(inout)           :: this
      type(extendable_str), intent(inout)     :: xyzfile_es
      character(len=*), optional, intent(out) :: comment
      character(len=*), optional, intent(out) :: properties 
      integer, optional, intent(out)          :: status

      call atoms_read_xyz(this, xyz_extendable_str=xyzfile_es, comment=comment, status=status, properties=properties)
    end subroutine atoms_read_xyz_extendable_str

    subroutine parse_line_io_es(io, es, delimiters, fields, num_fields, status)
      type(Inoutput), intent(inout), optional :: io
      type(extendable_str), intent(inout), optional :: es
      character(*),               intent(in)    :: delimiters
      character(*), dimension(:), intent(inout) :: fields
      integer,                    intent(out)   :: num_fields
      integer, optional,          intent(out)   :: status

      if (present(io) .and. present(es)) &
	call system_abort("parse_line_io_es called with both io and es")
      if (.not. present(io) .and. .not. present(es)) &
	call system_abort("parse_line_io_es called with neither io nor es")

      if (present(io)) then
	call parse_line(io, delimiters, fields, num_fields, status)
      end if
      if (present(es)) then
	call parse_line(es, delimiters, fields, num_fields, status)
      end if
    end subroutine parse_line_io_es

    function read_line_io_es(io, es, status)
      type(Inoutput), intent(inout), optional       :: io
      type(extendable_str), intent(inout), optional :: es
      integer, optional,          intent(out)       :: status
      character(len=1024) :: read_line_io_es

      if (present(io) .and. present(es)) &
	call system_abort("read_line_io_es called with both io and es")
      if (.not. present(io) .and. .not. present(es)) &
	call system_abort("read_line_io_es called with neither io nor es")

      if (present(io)) then
	read_line_io_es = read_line(io, status)
      end if
      if (present(es)) then
	read_line_io_es = read_line(es, status)
      end if
    end function read_line_io_es

    subroutine atoms_read_xyz(this, xyz_inoutput, xyz_extendable_str, comment, status, properties)
      type(Atoms),    intent(inout)                 :: this
      type(Inoutput), intent(inout), optional       :: xyz_inoutput
      type(extendable_str), intent(inout), optional :: xyz_extendable_str
      character(len=*), optional, intent(out)       :: comment
      character(len=*), optional, intent(out)       :: properties 
      integer, optional, intent(out)                :: status

      type(Table) :: props
      integer                             :: i, j, k, N, lookup(3), num_fields, &
           req_num_fields, my_status, n_cols, field_count
      character(len=2048), dimension(200) :: fields
      character(len=1024)                 :: use_comment
      character(len=1024)                  :: prop_names, tmp
      character(3)                        :: delimiters
      real(dp), dimension(9)              :: tmplattice
      integer, dimension(9)               :: tmp_int_lattice
      real(dp), dimension(3,3)            :: lattice
      logical                             :: lattice_found, properties_found, Z_found
      real(dp)                            :: max_dist
      character(len=1024) :: filename

      type(Dictionary) :: params

      delimiters = ' ="'
      max_dist = 0.0_dp

      if (present(xyz_inoutput) .and. present(xyz_extendable_str)) &
	call system_abort("atoms_read_xyz called with both inoutput and extendable_str")

      if (present(xyz_inoutput)) then
	filename = "from " // trim(xyz_inoutput%filename)
      else
	filename = "from extendable_str"
      endif

      !Read the number of atoms
      call parse_line_io_es(xyz_inoutput,xyz_extendable_str,delimiters,fields,num_fields,my_status)

      if(present(status)) then
         status = my_status
         if(status /= 0) return
      else
	 if (my_status > 0) then
	    call system_abort('Atoms_read_xyz: Error reading from '//filename)
	 else if (my_status < 0) then
	    call system_abort('Atoms_read_xyz: End of file when reading from '//filename)
	 end if
      end if

      N = string_to_int(fields(1))

      ! Next, read the comment line into this%params dictionary
      !  and try to find the lattice and properties keys.
      use_comment = read_line_io_es(xyz_inoutput, xyz_extendable_str)
      call read_string(params, use_comment)

      ! Is there a real valued lattice?
      lattice_found = get_value(params, 'Lattice', tmplattice)
      lattice = reshape(tmplattice, (/3,3/))
      if (.not. lattice_found) then
         ! Try looking for an integer lattice instead
         lattice_found = get_value(params, 'Lattice', tmp_int_lattice)
         if (lattice_found) then
            lattice = reshape(real(tmp_int_lattice,dp),(/3,3/))
         end if
      end if

      ! If there's no 'Properties=' tag in parameters then just look for species and positions
      properties_found = get_value(params, 'Properties', prop_names)
      if (.not. properties_found) then
         prop_names = 'species:S:1:pos:R:3'
      end if

      ! Backward compatibility - prepend species tag if it's missing
      if (prop_names(1:7) == 'pos:R:3') prop_names = 'species:S:1:'//prop_names

      if (present(comment)) then
         ! Remove Lattice="..." and Properties=... from comment line
         if (lattice_found) then
            i = index(use_comment,'Lattice')
            tmp = use_comment(i:)
            tmp(index(tmp,'"'):index(tmp,'"')) = ' ' ! kill first " after lattice
            use_comment = adjustl(use_comment(:i-1)//use_comment(i+index(tmp, '"'):))
         end if

         if (properties_found) then
            i = index(use_comment, 'Properties')
            tmp = use_comment(i:)
            use_comment = adjustl(use_comment(:i-1)//use_comment(i+index(tmp, ' '):))
         end if

         comment = use_comment
      end if

      call atoms_initialise(this, N, lattice)
      this%params = params
      call finalise(params)

      req_num_fields = 0
      call allocate(props, 3, 0, 0, 0)

      if (present(properties)) properties = ''

      call parse_string(prop_names,':',fields,num_fields)
      Z_found = .false.
      do i=1,num_fields,3
         n_cols = string_to_int(fields(i+2))
         req_num_fields = req_num_fields + n_cols

         if (trim(fields(i)) == 'Z') Z_found = .true.

         select case(trim(fields(i+1)))
            case('I')
               call add_property(this, fields(i), 0, n_cols=n_cols)

            case('R')
               call add_property(this, fields(i), 0.0_dp, n_cols=n_cols)

            case('S')
               call add_property(this, fields(i), repeat(' ',TABLE_STRING_LENGTH), n_cols=n_cols)

            case('L')
               call add_property(this, fields(i), .false., n_cols=n_cols)
               
            case default
               call system_abort('Atoms_read_xyz: unknown property type '//fields(i+1))   

         end select
         if (present(properties)) then 
            properties = trim(properties)//trim(fields(i))//':'
         end if

         if (.not. get_value(this%properties, fields(i), lookup)) &
              call system_abort('Atoms_read_xyz: key '//fields(i)//' not found')
         call append(props, lookup)
      end do

      ! Remove trailing ':'
      if (present(properties)) properties = trim(properties(1:len(properties)-1))

      ! Now read in the atomic data
      do i = 1, N

         call parse_line_io_es(xyz_inoutput,xyz_extendable_str,delimiters,fields,num_fields)
         if (num_fields < req_num_fields) then
            call print_title('ERROR')
            call print('Fields:')
            do j = 1, num_fields
               call print(j//': '//trim(fields(j)))
            end do
            call print('Required number of fields = '//req_num_fields)
            call system_abort('Atoms_read_xyz: not enough fields in xyz '//filename)
         end if
           
         field_count = 1
         do j=1,props%N
            do k=props%int(2,j),props%int(3,j)
               select case(props%int(1,j))
                  case(PROPERTY_INT)
                     read (fields(field_count),*) this%data%int(k,i)

                  case(PROPERTY_REAL)
                     read (fields(field_count),*) this%data%real(k,i)
                     
                  case(PROPERTY_STR)
                     this%data%str(k,i) = fields(field_count)

                  case(PROPERTY_LOGICAL)
                     read (fields(field_count),*) this%data%logical(k,i)

               end select
               field_count = field_count + 1
            end do
         end do
      end do

      ! Set Z from species if we didn't read a Z property from file
      if (.not. Z_found) then
         call print('Atoms_read_xyz: Setting Z from species', VERBOSE)
         do i=1,this%N
            this%species(i) = ElementFormat(this%species(i))
            this%Z(i) = Atomic_Number(this%species(i))
            if (this%Z(i) == 0) call system_abort('Atoms_read_xyz: Could not parse atom name: '//this%species(i))
         end do
      end if

      call finalise(props)

      ! If we didn't get a lattice from the file then make a safe one
      if (.not.lattice_found) then
         lattice = 0.0_dp

         max_dist = maxval(norm(this%pos,1))
         lattice(1,1) = 5.0_dp*max_dist
         lattice(2,2) = 5.0_dp*max_dist
         lattice(3,3) = 5.0_dp*max_dist

         call atoms_set_lattice(this, lattice)
      end if

    end subroutine atoms_read_xyz



   subroutine atoms_read_xyz_skip(xyzfile, status)

     type(Inoutput),    intent(inout) :: xyzfile
     integer, optional, intent(out)   :: status
     integer                          :: i,n,my_status
     character(80)                    :: text

     if (present(status)) status = 0

     !Read number of atoms
     text = Read_Line(xyzfile,my_status)
     if (my_status /= 0) then
        if (present(status)) then
           status = my_status
           return
        else
           call system_abort('Atoms_Read_xyz_skip: Problem skipping frame in file '//trim(xyzfile%filename))
        end if
     end if
     n = String_To_Int(text)

     !Now read past the comment and that number of atoms
     do i = 1, n+1
        text = Read_Line(xyzfile,my_status)
        if (my_status /= 0) then
           if (present(status)) then
              status = my_status
              return
           else
              call system_abort('Atoms_Read_xyz_skip: Problem skipping frame in file '//trim(xyzfile%filename))
           end if
        end if
     end do

   end subroutine atoms_read_xyz_skip


   subroutine atoms_print_cfg(this, cfgfile, comment, properties, all_properties)

      type(Atoms), target,    intent(in)    :: this     !% Atoms object to print
      type(Inoutput),         intent(inout) :: cfgfile  !% Inoutput object to write to
      character(*), optional, intent(in)    :: comment  !% Comment line (line #2 of xyz file)
      character(*), optional, intent(in)    :: properties !% Colon separated list of columns from 'this%data' to printed.
                                                          !% If not specified, we print only the atom
                                                          !% positions, unless 'all_properties' is true.      
      logical,      optional, intent(in)    :: all_properties !% Should we print all properties (default is false)


      integer                                 :: i,j,k,lookup(3), n_properties, n_cols, n_aux
      character(len=1)                        :: prop_types(2)
      logical                                 :: do_all_properties, dummy
      character(len=1024) :: use_properties(this%properties%N), tmp
      type(Table) :: props
      real(dp), dimension(:,:), pointer :: pos
      real(dp), dimension(:), allocatable :: mass

      do_all_properties = .false.
      if (present(all_properties)) do_all_properties = all_properties

      if (.not.this%initialised) &
           call system_abort('Atoms_Print_cfg: Atoms structure not initialised')

      if (any(this%Z == 0)) call system_abort('Atoms_Print_cfg: Atomic number is zero')

      prop_types = (/'I', 'R'/)
      call table_allocate(props, 3, 0, 0, 0)

      if (.not. present(properties)) then

         ! Print all real and int properties, starting with position
         if (do_all_properties) then
            
            n_properties = 1
            do i=1,this%properties%N
               dummy = get_value(this%properties, this%properties%keys(i), lookup)
               if (lookup(1) == PROPERTY_INT .or. lookup(1) == PROPERTY_REAL) then
                  use_properties(n_properties) = this%properties%keys(i)
                  n_properties = n_properties + 1
               end if
            end do
            n_properties = n_properties - 1

            ! Move 'pos' to beginning of list
            i = lookup_entry_i(this%properties, 'pos')

            tmp = use_properties(1)
            use_properties(1) = use_properties(i)
            use_properties(i) = tmp
         else
            ! Only print the atomic positions
            use_properties(1) = 'pos'
            n_properties = 1
         end if

      else
         call parse_string(properties, ':', use_properties, n_properties)
      end if

      ! Get list of properties to print
      n_cols = 0
      do i=1,n_properties

         if (.not. get_value(this%properties, trim(use_properties(i)),lookup)) &
              call system_abort('Atoms_print_cfg: Property '// &
              trim(use_properties(i))//' not found in this%properties')

         call append(props, lookup)

         if (lookup(1) /= PROPERTY_INT .and. lookup(1) /= PROPERTY_REAL) &
              call system_abort('Atoms_print_cfg: Bad property type for '// &
              trim(use_properties(i)))

         n_cols = n_cols + (lookup(3)-lookup(2)+1)
      end do

      ! Check first property has three real columns (i.e. looks like a position)
      if (props%int(1,1) /= PROPERTY_REAL .or. (props%int(3,1)-props%int(2,1)+1) /= 3) &
           call system_abort('Atoms_print_cfg: first property must be a real 3-vector, e.g. pos or avgpos')

      ! Iniitalise data to 'pos'-like property
      pos => this%data%real(props%int(2,1):props%int(3,1),1:this%N)

      ! Header line
      write (line, '(a,i0)') 'Number of particles = ', this%N
      call print(line,file= cfgfile)

      if (present(comment)) then 
         write (line, '(a, a)') '# ', trim(comment)
         call print(line,file= cfgfile)
      end if

      ! Write lattice vectors - note H0(i,:) = at%lattice(:,i) since
      ! AtomEye uses C array ordering and we use Fortran.
      do i=1,3
         do j=1,3
            write(line, '(a,i0,a,i0,a,f16.8)') 'H0(',i,',',j,') = ',this%lattice(j,i)
            call print(line,file= cfgfile)
         end do
      end do

      ! Number of entries
      call print('.NO_VELOCITY.', file=cfgfile)
      write (line, '(a,i0)') 'entry_count = ', n_cols
      call print(line,file= cfgfile)

      ! Names for auxiliary properties
      n_aux = 0
      do i = 2,n_properties

         do j=props%int(2,i),props%int(3,i)
            if (props%int(2,i) == props%int(3,i)) then
               write (line, '(a,i0,a,a)') 'auxiliary[',n_aux,'] = ',trim(use_properties(i))
            else
               write (line, '(a,i0,a,a,a,i0)') 'auxiliary[',n_aux,'] = ',trim(use_properties(i)),'_',j-props%int(2,i)
            end if
            n_aux = n_aux + 1
            call print(line,file= cfgfile)
         end do
      end do

      ! 3 lines for first atom
      allocate(mass(this%N))
      if (has_property(this, 'mass')) then
         mass = this%mass
      else
         mass = ElementMass(this%Z)
      end if
      write (line, '(f12.4)') mass(1)/MASSCONVERT
      call print(line,file= cfgfile)
      write (line, '(a)') ElementName(this%Z(1))
      call print(line,file= cfgfile)

     ! First the pos-like property...
     write (line, '(3f16.8)') (this%g .mult. pos(:,1))
     do j=2,n_properties ! ...then all the other properties
        do k=props%int(2,j),props%int(3,j)
           select case(props%int(1,j))
              case(PROPERTY_INT)
                 write(tmp,'(a,i8)')' ',this%data%int(k,1)

              case(PROPERTY_REAL)
                 write(tmp,'(a,f16.8)') ' ',this%data%real(k,1)
           end select
           line = trim(line)//trim(tmp)
        end do
     end do
     call print(line,file= cfgfile)

     do i=2,this%N
        ! Only need to write mass and element name if it's different
        ! from that of the previous element
        if (this%Z(i) /= this%Z(i-1)) then
           write (line, '(f12.4)') mass(i)/MASSCONVERT
           call print(line,file= cfgfile)
           write (line, '(a)') ElementName(this%Z(i))
           call print(line,file= cfgfile)
        end if

        write (line, '(3f16.8)') (this%g .mult. pos(:,i))
        do j=2,n_properties
           do k=props%int(2,j),props%int(3,j)
              select case(props%int(1,j))
                 case(PROPERTY_INT)
                    write(tmp,'(a,i8)')' ',this%data%int(k,i)

                 case(PROPERTY_REAL)
                    write(tmp,'(a,f16.8)') ' ',this%data%real(k,i)
              end select
              line = trim(line)//trim(tmp)
           end do
        end do
        call print(line,file= cfgfile)
     end do

     call finalise(props)
     deallocate(mass)
   end subroutine atoms_print_cfg


   subroutine atoms_print_cfg_filename(this, cfgfilename, comment, properties, all_properties)

      type(Atoms),            intent(in)    :: this     
      character(*),           intent(in)    :: cfgfilename
      character(*), optional, intent(in)    :: comment
      character(*), optional, intent(in)    :: properties 
      logical,      optional, intent(in)    :: all_properties

      type(Inoutput)  :: file

      call initialise(file, cfgfilename, action=OUTPUT)
      call atoms_print_cfg(this, file, comment, properties, all_properties)
      call finalise(file)

    end subroutine atoms_print_cfg_filename



   subroutine atoms_file_write(this,outfile)

      type(Atoms),     intent(inout) :: this
      type(Inoutput),  intent(inout) :: outfile

      !Header
      call write_binary('Atoms',outfile)

      !Scalar logicals
      call write_binary(this%initialised,outfile)        !1
      if (.not.this%initialised) return
      call write_binary(this%use_uniform_cutoff,outfile) !2
      !Scalar integers
      call write_binary(this%N,outfile)                  !1

      !Scalar reals
      call write_binary(this%cutoff,outfile)             !1
      call write_binary(this%cutoff_break,outfile)             !1

      !Real arrays
      call write_binary(this%lattice,outfile)            !1

      !Derived types
      call write_binary(this%data, outfile)              !1
      call write_binary(this%properties, outfile)        !2
      call write_binary(this%params, outfile)            !3
      call write_binary(this%connect,outfile)            !4

   end subroutine atoms_file_write


   subroutine atoms_file_read(this,infile)

      type(Atoms),     intent(inout) :: this
      type(Inoutput),  intent(inout) :: infile
      integer                        :: N
      real(dp)                       :: cutoff, cutoff_break
      character(5)                   :: id
      logical                        :: initialised, use_uniform_cutoff
      real(dp), dimension(3,3)       :: lattice            

      call atoms_finalise(this)

      call read_binary(id,infile)
      if (id /= 'Atoms') then
         write(line,'(a,a)') 'Atoms_File_Read: Bad Atoms structure in file ',trim(infile%filename)
         call system_abort(line)
      end if
      !Scalar logicals
      call read_binary(initialised,infile)                !1
      if (.not.initialised) return
      call read_binary(use_uniform_cutoff,infile)         !2
      !Scalar integers
      call read_binary(N,infile)                          !1
      !Scalar reals
      call read_binary(cutoff,infile)                     !1
      call read_binary(cutoff_break,infile)                     !1
      !Real arrays
      call read_binary(lattice,infile)                    !1

      call atoms_initialise(this, N, lattice)

      this%use_uniform_cutoff = use_uniform_cutoff
      this%cutoff = cutoff
      this%cutoff_break = cutoff_break

      !Derived types
      call read_binary(this%data,   infile)               !1
      call read_binary(this%properties, infile)           !2
      call read_binary(this%params, infile)               !3
      call read_binary(this%connect,infile)               !4

      call atoms_repoint(this)

   end subroutine atoms_file_read

  subroutine connection_print(this,file)

    type(Connection), intent(in)    :: this
    type(Inoutput),   optional, intent(inout) :: file
    integer                         :: i,j,k,n,m

    if (.not.this%initialised) call system_abort('Connection_Print: Connection object not initialised')

    call print('Connectivity data:',file=file)
    call print('-------------------------------------------------',file=file)
    call print('|    I    |    J    |    Shift     | Distance   |',file=file)
    call print('-------------------------------------------------',file=file)

    do i = 1, size(this%neighbour1)

       if (.not. associated(this%neighbour1(i)%t)) cycle

       if ((this%neighbour1(i)%t%N + this%neighbour2(i)%t%N) > 0) then
	 write(line,'(a47)')'| Neighbour1 (i <= j)                           |'
	 call print(line, file=file)
       endif

       do j = 1, this%neighbour1(i)%t%N
          if(j == 1) then
             write(line,'(a2,i7,a3,i7,a3,3i4,a3,f10.5,a2)')'| ',    &
                 & i, ' | ', this%neighbour1(i)%t%int(1,j),' | ',this%neighbour1(i)%t%int(2:4,j),&
                 &' | ',this%neighbour1(i)%t%real(1,j),' |'
          else
             write(line,'(a12,i7,a3,3i4,a3,f10.5,a2)')'|         | ',    &
                  this%neighbour1(i)%t%int(1,j),' | ',this%neighbour1(i)%t%int(2:4,j),' | ',this%neighbour1(i)%t%real(1,j),' |'
          end if
          call print(line,file=file)
       end do

       if (this%neighbour2(i)%t%N > 0) then

          write(line,'(a47)')'-----------------------------------------------'
          call print(line, file=file)
          write(line,'(a47)')'| Neighbour2 (i > j)                            |'
          call print(line, file=file)


          do j = 1, this%neighbour2(i)%t%N

             if(j == 1) then
                write(line,'(a2,i7,a3,i7,a3,i7,a)') '| ', i, ' | ', this%neighbour2(i)%t%int(1,j),' | ',&
                       &this%neighbour2(i)%t%int(2,j),'                 |'
             else
                write(line,'(a12,i7,a3,i7,a)') '|         | ', this%neighbour2(i)%t%int(1,j),' | ',&
                       &this%neighbour2(i)%t%int(2,j),'                 |'
             end if
             call print(line,file=file)
          end do

       end if

       if ((this%neighbour1(i)%t%N + this%neighbour2(i)%t%N) > 0) then
	 write(line,'(a47)')'-----------------------------------------------'
	 call print(line,file=file)
       endif

    end do

    write(line,'(a47)')''
    call print(line,file=file)


    !Print the cell lists if they exist
    if (this%cells_initialised .and. current_verbosity() > NORMAL) then
       call verbosity_push_decrement()

       write(line,'(a11)')'Cell Lists:'
       call print(line,file=file)
       write(line,'(a70)')'----------------------------------------------------------------------'
       call print(line,file=file)

       ! Write atomic indices 10 per line
       do k = 1, this%cellsNc
          do j = 1, this%cellsNb
             do i = 1, this%cellsNa
                write(line,'(a7,i0,a1,i0,a1,i0,a3)')'Cell ( ',i,' ',j,' ',k,' ):'
                call print(line,file=file)
                n = this%cell(i,j,k)%N / 10 ! integer division will round down
                do m = 1, n
                   write(line,'(10i7)') this%cell(i,j,k)%int(1,(10*m-9):(10*m))
                   call print(line,file=file)
                end do
                m = this%cell(i,j,k)%N - 10 * n
                if (m /= 0) then
                   write(line,'(10i7)') this%cell(i,j,k)%int(1,(10*n+1):(10*n+m))
                   call print(line,file=file)
                end if
                write(line,'(a70)')'----------------------------------------------------------------------'
                call print(line,file=file)
             end do
          end do
       end do

       call verbosity_pop()

    end if

  end subroutine connection_print

  subroutine table_pointer_write_binary_a(this, outfile)
    type(Table_Pointer), intent(in) :: this(:)
    type(Inoutput), intent(inout) :: outfile

    integer i

    call write_binary(size(this), outfile)
    do i=1, size(this)
      call write_binary(this(i)%t, outfile)
    end do
  end subroutine table_pointer_write_binary_a


  subroutine connection_file_write(this,outfile)
    type(Connection), intent(in) :: this
    type(Inoutput),   intent(inout) :: outfile
    integer                         :: N

    N = size(this%neighbour1)
    call write_binary('Connection',outfile)
    call write_binary(this%initialised,outfile)
    if (.not.this%initialised) return
    call write_binary(this%cells_initialised,outfile)
    call write_binary(this%too_few_cells_warning_issued, outfile)
    call write_binary(N,outfile)
    call write_binary(this%cellsNa,outfile)
    call write_binary(this%cellsNb,outfile)
    call write_binary(this%cellsNc,outfile)

    call write_binary(this%neighbour1,outfile)
    call write_binary(this%neighbour2,outfile)

    if (this%cells_initialised) &
         call write_binary(this%cell,outfile)

  end subroutine connection_file_write

  subroutine table_pointer_read_binary_a(this, outfile)
    type(Table_Pointer), allocatable, intent(inout) :: this(:)
    type(Inoutput), intent(inout) :: outfile

    integer i, N

    call read_binary(N, outfile)
    allocate(this(N))
    do i=1, N
      call read_binary(this(i)%t, outfile)
    end do
  end subroutine table_pointer_read_binary_a

  subroutine connection_file_read(this,infile)
    type(Connection), intent(inout) :: this
    type(Inoutput), intent(inout)   :: infile
    integer                         :: cellsNa,cellsNb,cellsNc,N
    character(10)                   :: id
    logical                         :: initialised, cells_initialised

    call connection_finalise(this)

    call read_binary(id,infile)
    if (id /= 'Connection') then
       write(line,'(a,a)')'Connection_File_Read: Bad Connection structure in file ',trim(infile%filename)
       call system_abort(line)
    end if

    call read_binary(initialised,infile)
    if (.not.initialised) return
    call read_binary(cells_initialised,infile)
    call read_binary(this%too_few_cells_warning_issued, infile)
    call read_binary(N,infile)
    call read_binary(cellsNa,infile)
    call read_binary(cellsNb,infile)
    call read_binary(cellsNc,infile)

    call connection_initialise(this,N)

    call read_binary(this%neighbour1, infile)
    call read_binary(this%neighbour2, infile)

    if (cells_initialised) then
       call read_binary(this%cell, infile)

       this%cellsNa = cellsNa
       this%cellsNb = cellsNb
       this%cellsNc = cellsNc

       this%cells_initialised = .true.
    end if

  end subroutine connection_file_read


  !Bond_Length
  !% Returns the sum of the covalent radii of two atoms
  function bond_length(z1,z2)
    integer, intent(in) :: z1,z2
    real(dp)            :: bond_length
    bond_length = ElementCovRad(z1) + ElementCovRad(z2)
  end function bond_length

  !Termination_Bond_Rescale
  !% Calculates the rescale ratio of a Z1--H bond 
  !% generate from a Z1--Z2 bond.
  function termination_bond_rescale(z1,z2)
    integer,  intent(in) :: z1, z2
    real(dp)             :: termination_bond_rescale
    termination_bond_rescale = (ElementCovRad(z1) + ElementCovRad(1)) &
         / (ElementCovRad(z1) + ElementCovRad(z2))
  end function termination_bond_rescale

  !% Parses an atom_mask, which is string consisting of the '@' symbol followed by a comma separated
  !% list of indices or ranges into a table containing all the indices it represents.
  !% E.g. '@1,37-39,54,99-102' is expanded to a table with 1, 37, 38, 39, 54, 99, 100,
  !% 101, 102 as its first integer column. There must be no spaces in the mask.
  subroutine parse_atom_mask(mask_in,atom_indices)

    character(*),  intent(in)    :: mask_in
    type(Table),   intent(out)   :: atom_indices
    character(len(mask_in))      :: mask
    character(20), dimension(50) :: fields
    integer                      :: Nfields, i, j, n, a, b, c

    call allocate(atom_indices,1,0,0,0,100)

    mask = adjustl(mask_in)

    if (mask(1:1)=='@') then  !Found atom mask
       !Remove @ symbol
       mask = mask(2:)
       !Parse into comma separated numbers or ranges
       call parse_string(mask,',',fields,Nfields)
       if (Nfields==0) call system_abort('parse_atom_mask: Atom mask contained no indices')
       do i = 1, Nfields
          !Is this field a single atom or a contiguous range?
          n = scan(fields(i),'-')
          if (n /= 0) then
             !it's a range. get start number
             a = string_to_int(fields(i)(:n-1))
             !get end number
             b = string_to_int(fields(i)(n+1:))
             if (a > b) then
                c = b
                b = a
                a = c
             end if
             !add all these atoms to the list
             do j = a, b
                if (.not.is_in_array(int_part(atom_indices,1),j)) call append(atom_indices,j)
             end do
          else
             !it's a single atom.
             a = string_to_int(fields(i))
             if (.not.is_in_array(int_part(atom_indices,1),j)) call append(atom_indices,a)
          end if
       end do
    else
       call system_abort('parse_atom_mask: Invalid atom mask: '//mask_in)
    end if

  end subroutine parse_atom_mask

  !% Find atoms which have integer property 'prop' set to
  !% true (i.e. set to 1) and return them in the Table 'list'.
  subroutine property_to_list(at, prop, list)
    type(Atoms), intent(in) :: at
    character(*), intent(in) :: prop
    type(Table), intent(inout) :: list
    
    integer :: i
    integer, dimension(:), pointer :: p

    if (.not. assign_pointer(at, prop, p)) &
         call system_abort('property_to_list: at does not have property "'//trim(prop)//'".')

    call allocate(list, 1, 0, 0, 0)
    call append(list, pack((/ (i, i=1,size(p)) /), p == 1))

  end subroutine property_to_list

  !% Convert the Table 'list' to a single column integer property in 'at',
  !% with atoms in list marked with a 1 and absent 
  subroutine list_to_property(at, list, prop)
    type(Atoms), intent(inout) :: at
    type(Table), intent(in) :: list
    character(*), intent(in) :: prop

    integer, dimension(:), pointer :: p

    ! Add property if necessary
    call add_property(at, prop, 0)

    if (.not. assign_pointer(at, prop, p)) &
         call system_abort('list_to_property: at does not have property "'//trim(prop)//'".')

    p = 0 ! Set all entries to zero
    p(int_part(list,1)) = 1 ! Set entries in 'list' to 1
    
  end subroutine list_to_property


  !% Return the complement of a list, i.e. all those atoms not included
  !% in list. Result is in outlist on exit.
  !% (not yet tested).
  subroutine complement(at, inlist, outlist)
    type(Atoms), intent(in) :: at
    type(Table), intent(in) :: inlist
    type(Table), intent(out) :: outlist

    integer :: i
    integer, allocatable, dimension(:) :: inarray
    
    call table_allocate(outlist,1,0,0,0)
    
    allocate(inarray(inlist%N))
    inarray = int_part(inlist,1)
    
    do i=1,at%N
       if (.not. is_in_array(inarray,i)) &
            call append(outlist,i)
    end do

    deallocate(inarray)
    
  end subroutine complement


  !% Return the difference between list1 and list2 in outlist.
  !% That is, those elements in list1 but not in list2
  !% (not yet tested).
  subroutine difference(list1, list2, outlist)
    type(Table), intent(in) :: list1, list2
    type(Table), intent(out) :: outlist

    integer :: i
    integer, dimension(:), allocatable :: array1, array2

    if (list1%N <= list2%N) &
         call system_abort('difference: list1%N ('//(list1%N)//') <= list2%N ('//(list2%N)//').')

    call table_allocate(outlist, 1, 0, 0, 0)

    allocate(array1(list1%N),array2(list2%N))
    array1 = int_part(list1,1)
    array2 = int_part(list2,1)

    do i=1,list1%N
       if (.not. is_in_array(array2, array1(i))) &
            call append(outlist,list1%int(1,i))
    end do

    deallocate(array1,array2)

  end subroutine difference

  !% move atoms around following neighbor list bonds so that all are in the same periodic image
  !%    (that of 'seed', if present)
  !% poorly tested, especially for situations where not all atoms are in one connected clump
  !% probably needs a better subroutine name
  subroutine coalesce_in_one_periodic_image(this, seed, stat)
    type(Atoms), intent(inout) :: this
    integer, intent(in), optional :: seed
    integer, intent(out), optional :: stat

    integer :: i, ji, jji, j, shift(3), delta_shift(3), jj, k, ki
    integer :: n_neighbours, max_n_neighbours
    integer, allocatable :: shifts(:,:,:)
    logical, allocatable :: touched(:)
    integer :: seed_val

    seed_val=optional_default(1, seed)

    max_n_neighbours = 0
    do i=1, this%N
      n_neighbours = atoms_n_neighbours(this, i)
      if (n_neighbours > max_n_neighbours) max_n_neighbours = n_neighbours
    end do
    allocate(shifts(3,max_n_neighbours,this%N))
    shifts = 0

    do i=1, this%N
      do ji=1, atoms_n_neighbours(this, i)
	j = atoms_neighbour(this, i, ji, shift=shift)
	shifts(:,ji,i) = shift
      end do
    end do

    allocate(touched(this%N))
    touched = .false.
    touched(seed_val) = .true.

    do while (any(shifts /= 0))
      ! look for atoms i that have been touched
      do i=1, this%N
	if (touched(i)) then
	  ! look at neighbors of i
	  do ji=1, atoms_n_neighbours(this,i)
	    j = atoms_neighbour(this, i, ji)
	    ! if neighbor has 0 shift, it doesn't need to move
	    if (any(shifts(:,ji,i) /= 0)) then
	      ! atom j does need to move
	      ! atoms with non-zero shift should not have been touched before
	      if (touched(j)) then
		if (present(stat)) then
		  call print("ERROR: undo_pbcs tried to move atom " // j // " twice", ERROR)
		  stat = 1
		else
		  call system_abort("undo_pbcs tried to move atom " // j // " twice")
		endif
	      endif

	      ! shift atom to zero out shift
	      this%pos(:,j) = this%pos(:,j) + shifts(1,ji,i)*this%lattice(:,1) + &
					  shifts(2,ji,i)*this%lattice(:,2) + &
					  shifts(3,ji,i)*this%lattice(:,3)
	      ! fix shifts of j's neighbours
	      delta_shift = shifts(:,ji,i)
	      do jji=1, atoms_n_neighbours(this, j)
		! fix shifts from j to its neighbors
		jj = atoms_neighbour(this, j, jji)
		shifts(:,jji,j) = shifts(:,jji,j) + delta_shift(:)
		! fix shifts from j's neighbours to it
		do ki=1, atoms_n_neighbours(this, jj)
		  k = atoms_neighbour(this, jj, ki)
		  if (k == j) then
		    shifts(:,ki,jj) = shifts(:,ki,jj) - delta_shift(:)
		  endif
		end do ! ki
	      end do ! jji
	      ! shifts(:,ji,i) = 0
	      touched(j) = .true.
	    else
	      ! atom j doesn't need to move
	      touched(j) = .true.
	    endif ! any(shift) /= 0
	  end do ! ji
	endif ! touched i
      end do ! i
    end do ! some shift isn't zero

    deallocate(touched)
    deallocate(shifts)

  end subroutine coalesce_in_one_periodic_image

  subroutine atoms_ll_new_entry(this, atoms_p)
    type(atoms_ll), target, intent(inout) :: this
    type(atoms), intent(out), pointer :: atoms_p

    type(atoms_ll_entry), pointer :: entry

    allocate(entry)
    if (associated(this%first) .or. associated(this%last)) then
      if (.not. associated(this%first) .or. .not. associated(this%last)) &
	call system_abort("new_entry associated(first) " // associated(this%first) // &
		  " associated(this%lat) " // associated(this%last))
      this%last%next => entry
      entry%prev => this%last
      this%last => entry
    else
      this%first => entry
      this%last => entry
    endif

    atoms_p => entry%at
  end subroutine atoms_ll_new_entry

  subroutine atoms_ll_remove_last_entry(this)
    type(atoms_ll), target, intent(inout) :: this

    if (associated(this%last)) then
      call finalise(this%last%at)
      nullify(this%last%prev%next)
      this%last => this%last%prev
    endif

  end subroutine atoms_ll_remove_last_entry

  subroutine atoms_ll_print_xyz(this, file)
    type(atoms_ll), intent(inout), target :: this
    type(inoutput), optional, intent(inout) :: file

    type(atoms_ll_entry), pointer :: entry
    integer :: i

    entry => this%first
    i = 1
    do while (associated(entry)) 
      call print("atoms_ll entry # " // i, file=file)
      call print_xyz(entry%at, xyzfile=file)
      entry => entry%next
      i = i + 1
    end do
  end subroutine atoms_ll_print_xyz

  function closest_atom(this, r, cell_image_Na, cell_image_Nb, cell_image_Nc, dist)
    type(Atoms), intent(in) :: this
    real(dp), intent(in) :: r(3)
    integer, intent(in) :: cell_image_Na, cell_image_Nb, cell_image_Nc
    real(dp), intent(out), optional :: dist
    integer :: closest_atom

    integer :: i, j, k
    integer i2, j2, k2, i3, j3, k3, i4, j4, k4, n2, atom_i
    integer :: cellsNa, cellsNb, cellsNc
    real(dp) :: pos(3), cur_dist, min_dist

    if (.not. this%connect%initialised) call system_abort("closest_atom must have initialised connection object")

    call cell_of_pos(this%connect, this%g, r, i, j, k)

    cellsNa = this%connect%cellsNa
    cellsNb = this%connect%cellsNb
    cellsNc = this%connect%cellsNc

    min_dist = 1.0e38_dp
    closest_atom = 0
    ! Loop over neighbouring cells, applying PBC
    do k2 = -cell_image_Nc, +cell_image_Nc

       ! the stored cell we are in 
       if(cellsNc > 1) k3 = mod(k+k2-1+cellsNc,cellsNc)+1 

       ! the shift we need to get to the cell image
       k4 = (k+k2-k3)/cellsNc

       do j2 = -cell_image_Nb, +cell_image_Nb
	  ! the stored cell we are in                 
	  if(cellsNb > 1) j3 = mod(j+j2-1+cellsNb,cellsNb)+1 

	  ! the shift we need to get to the cell image
	  j4 = (j+j2-j3)/cellsNb

	  do i2 = -cell_image_Na, +cell_image_Na
	     ! the stored cell we are in                 
	     if(cellsNa > 1) i3 = mod(i+i2-1+cellsNa,cellsNa)+1 

	     ! the shift we need to get to the cell image
	     i4 = (i+i2-i3)/cellsNa

	     ! The cell we are currently testing atom1 against is cell(i3,j3,k3)
	     ! with shift (i4,j4,k4)
	     ! loop over it's atoms and test connectivity if atom1 < atom2

	     do n2 = 1, this%connect%cell(i3,j3,k3)%N
		atom_i = this%connect%cell(i3,j3,k3)%int(1,n2)
		pos = this%pos(:,atom_i) + ( this%lattice .mult. (/ i4, j4, k4 /) )
		cur_dist = norm(pos-r)
		if (cur_dist < min_dist) then
		  min_dist = cur_dist
		  closest_atom = atom_i
		endif

	     end do

	  end do
       end do
    end do

    if (present(dist)) dist = min_dist

  end function closest_atom

end module atoms_module
