Installation instructions:

Create your on arch.make or copy one of the pre-defined arch.*.

Type "make".


### Old readme.txt

Src.tar.gz contains the entire source tree.
 
You need to set the following environment variables

alias make gmake
set hosttype = linux  (or aix alpha CYGWIN HPUX irix64 OSF1 sun4 sgi Darwin)
setenv SYS `uname`
setenv LIB_PATH ~/Co/${hosttype}Lib
setenv BIN_PATH ~/Co/${hosttype}Bin
setenv INC_PATH ~/Co/Include
setenv MAN_PATH ~/Co/Man

The Makefile is for gmake only.

The sequence of make directories is:

IO/
Scalar/
Min/
Timer/ 
VecMat/ 
VecMat2/ 
VecMat3/ 
Atoms/
AX/
A/

