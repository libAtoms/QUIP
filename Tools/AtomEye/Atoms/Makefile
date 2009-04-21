#################################################
# Atoms: -llapack -lblas -lm                    #
#        -lVecMat3 -lVecMat -lScalar -lIO       #
#                                               #
# Physical constants, macros, and configuration #
# operators for atomistic simulations.          #
#                                               #
# Dec. 12, 1999  Ju Li <liju99@mit.edu>         #
#################################################

LIB = Atoms
HEADERS = NIST.h ${LIB}.h  xyz_netcdf.h # plural
STATIC_LIB_FILE = lib${LIB}.a
PROGS = mul iniT punchcylinder vcut trans scut voronoirize \
	label_immobile c2p p2c cfg2xyz xtal vasp2cfg vasp2out cfg2vasp \
	chemical_randomize perturb_s perturb_x fch2cfg chg2xsf \
	cfg2dacapo annotate_atomic_strain affine_deform linear_path \
	period_reduce infectH concat div
OBJS = ${LIB}.o Config.o Neighborlist.o Gr.o Motion.o Voronoi.o \
	VASP.o Gaussian.o Dacapo.o LeastSquareStrain.o ChainofStates.o xyz_netcdf.o # plural
SRCS = $(addsuffix .c, $(basename ${OBJS}))
GDBTMP = /tmp/${LIB}.gdb
EXESUFFIX =

include Makefile.inc
include Makefile.common
ifeq (${QUIP_ARCH},)
  include Makefile.arch
else
  include Makefile.${QUIP_ARCH}
endif

USE_SYSTEM_LIB = ${ATOMEYE_SYSTEM_LIB}

MY_BIN_DIR = ../bin
MY_LIB_DIR = ../lib
MY_INC_DIR = ../include
MY_BAK_DIR = ../bak

USE_MY_LIB = \
   -lAX -lAtoms -lVecMat3 -lVecMat -lIO -lScalar -lTimer

CFLAGS += -D_$(SYS) ${CCFLAGS} -I${MY_INC_DIR} -DHAVE_ATOMEYE ${ATOMEYE_LIB}
LDFLAGS = ${SYSTEM_LIB_FLAG} ${USE_SYSTEM_LIB} -L${MY_LIB_DIR} ${USE_MY_LIB}

.PHONY: all  i install  t  c clean  b backup  h  r  g  p  a access

all:	install

%.o:    %.c
	${CC} ${CFLAGS} -c $< -o $@
Atoms.o: Atoms.c Mendeleyev.c
	${CC} ${CFLAGS} -c $< -o $@

${STATIC_LIB_FILE}:	${OBJS}
	-@echo "***************************************************"; \
	rm -f $@; \
	${AR} $@ ${OBJS}; \
	${RANLIB} $@; \
	${NM} $@; \
	echo "***************************************************"

i:	install
install: ${STATIC_LIB_FILE} ${HEADERS}
	-@cp -f $< ${MY_LIB_DIR}; \
	chmod a+r ${MY_LIB_DIR}/$<; \
	cp -f ${HEADERS} ${MY_INC_DIR}; \
	dump=`ls -ld ${MY_LIB_DIR}/$< | sed 's/    */  /g'`; \
	month=`echo $$dump | cut -d' ' -f 6`; \
	prev=`echo "$$dump" | sed "s/$$month.*//"`; \
	ww=`echo "$$prev" | wc -c | sed 's/ //g'`; \
	l=`expr $$ww - 1`; \
	machine_name=`hostname | cut -f 1 -d .`:"               "; \
	machine_name=`echo Current time on "$$machine_name" | cut -c1-$$l`; \
	echo "$$machine_name`date '+%b %d %H:%M' | sed 's/    */  /g'`"; \
	a=`echo $$HOME | sed 's/\//~/g'`; \
	for file in ${HEADERS}; do \
	 chmod a+r ${MY_INC_DIR}/$$file; \
	 a=`echo $$HOME | sed 's/\//~/g'`; \
	 ls -l ${MY_INC_DIR}/$$file | sed 's/    */  /g' | sed 's/\//~/g' |\
	 sed "s/$$a//" | sed -e "s/~/\//g" | sed "s/\//~\//"; \
	done; \
	ls -l ${MY_LIB_DIR}/$< | sed 's/    */  /g' | sed 's/\//~/g' | \
	sed "s/$$a//" | sed -e "s/~/\//g" | sed "s/\//~\//"

r:
	-rm -f *_TEST${EXESUFFIX}; \
	rm -f ${STATIC_LIB_FILE}

h:
	cp -f ${HEADERS} ${MY_INC_DIR}

c:	clean
clean:
	@echo
	-rm -f *~; \
	rm -f \#*\#; \
	rm -f .\#*; \
	rm -f *.o; \
	rm -f ${STATIC_LIB_FILE}; \
	rm -f *_TEST; \
	rm -f ${CHECKCONST_TEMPLATE}.c; \
	rm -f a.out; \
	rm -f $(addsuffix ${EXESUFFIX}, ${PROGS}); \
	rm -f *.m *.m~; \
	rm -f Notes/*~; \
	rm -f Notes/[0-9]; \
	rm -f core
	@echo

b:	backup
backup: clean
	@tmp=/tmp/backupto.tar; \
	echo "taring the directory $(PWD) ..."; \
	dir=`echo $(PWD) | sed -e "s/.*\///g"`; \
	cd ..; \
	echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%; \
	tar cvf $$tmp $$dir/* \
	`ls -d $$dir/.* | sed -e "s/$$dir\/\.\.$$//" | \
	sed -e "s/$$dir\/\.$$//"`; \
	echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%; \
	echo compressing $$tmp with gzip ...; \
	gzip --verbose -c $$tmp > $$tmp.gz; \
	echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%; \
	nbytes=`ls -l $$tmp.gz | sed -e 's/  */ /g' | cut -d' ' -f 5`; \
	echo $$tmp.gz = $$nbytes bytes.; \
	for backup in $(MY_BAK_DIRS); do \
	 if [ -d $$backup ]; then \
	  echo ' '; \
	  echo "copying $$tmp.gz to $$backup/$$dir.tar.gz ..."; \
	  cp -i $$tmp.gz $$backup/$$dir.tar.gz; \
	 fi; \
	done; \
	echo

g:
	@export CC=gcc; export LD=g77; export CCFLAGS=-g; \
	make --environment-overrides --no-print-directory \
	`echo "${MAKECMDGOALS}" | sed -e 's/^g //' | sed -e 's/^g$///' `; \
	exit 99

p:	${PROGS}
${PROGS}:${SRCS} ${HEADER}  ${ADD_DEPENDENCE}
	@for obj in $(OBJS); do \
	 base=`echo $$obj | sed -e 's/\.o$$//g'`; \
	 grep _$@ $$base.c > /dev/null; \
	 if [ $$? = "0" ]; then \
	  ${CC} -D_$@ ${CFLAGS} -c $$base.c -o "$$base"_$@".o"; \
	  if [ $$? = "0" ]; then \
	   objs="$$base"_$@".o $$objs"; else exit 99; fi; \
	 else \
	  make --no-print-directory $$base.o > /dev/null; \
	  objs="$$objs $$base.o"; \
	 fi; \
	done; \
	${LD} ${CFLAGS} $$objs -o $@${EXESUFFIX} ${LDFLAGS}; \
	cp -f $@${EXESUFFIX} ${MY_BIN_DIR}; \
	chmod u+x ${MY_BIN_DIR}/$@${EXESUFFIX}; \
	ls -l ${MY_BIN_DIR}/$@${EXESUFFIX} | sed -e 's/   */  /g'



.PRECIOUS:	%_TEST${EXESUFFIX}
Makefile: ;

a:	access
access:
	-@chmod -R 700 *; \
	chmod 744 *.h; \
	if [ -f ${STATIC_LIB_FILE} ]; then chmod 744 ${STATIC_LIB_FILE}; fi; \
	ls -la
