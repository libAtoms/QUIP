# H0 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# H0 X
# H0 X   libAtoms+QUIP: atomistic simulation library
# H0 X
# H0 X   Portions of this code were written by
# H0 X     Albert Bartok-Partay, Silvia Cereda, Gabor Csanyi, James Kermode,
# H0 X     Ivan Solt, Wojciech Szlachta, Csilla Varnai, Steven Winfield.
# H0 X
# H0 X   Copyright 2006-2010.
# H0 X
# H0 X   These portions of the source code are released under the GNU General
# H0 X   Public License, version 2, http://www.gnu.org/copyleft/gpl.html
# H0 X
# H0 X   If you would like to license the source code under different terms,
# H0 X   please contact Gabor Csanyi, gabor@csanyi.net
# H0 X
# H0 X   Portions of this code were written by Noam Bernstein as part of
# H0 X   his employment for the U.S. Government, and are not subject
# H0 X   to copyright in the USA.
# H0 X
# H0 X
# H0 X   When using this software, please cite the following reference:
# H0 X
# H0 X   http://www.libatoms.org
# H0 X
# H0 X  Additional contributions by
# H0 X    Alessio Comisso, Chiara Gattinoni, and Gianpietro Moras
# H0 X
# H0 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

MATH_LINKOPTS=-L/usr/local/Cluster-Apps/intel/mkl/10.0.1.014/lib/em64t -Bstatic -lmkl_lapack -lmkl_em64t -lguide -Bdynamic -lpthread
FOX_LIBDIR=/home/kermode/QUIP/build.linux_x86_64_ifort_icc/../FoX-4.0.3/objs.linux_x86_64_ifort_icc/lib
FOX_INCDIR=/home/kermode/QUIP/build.linux_x86_64_ifort_icc/../FoX-4.0.3/objs.linux_x86_64_ifort_icc/finclude
NETCDF_LIBDIR=/home/kermode/lib
NETCDF_INCDIR=/home/kermode/include
NETCDF4=1
EXTRA_LINKOPTS=
OPTIM=${DEFAULT_OPTIM}


