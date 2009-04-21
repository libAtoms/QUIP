###########################################
# Threaded Atomistic Configuration Viewer #
#                                         #
# Fri Mar 31 2000 Ju Li <liju99@mit.edu>  #
###########################################

PROG = A
LIB = AtomEye
STATIC_LIB_FILE = lib${LIB}.a
HEADER = atomeyelib.h
OBJS  = $(PROG).o primitives.o viewport.o utils.o xtal_shift.o info.o \
	rcut_patch.o scratch.o geo.o LeastSquareStrain.o \
	cui.o cui_AX.o #p3dp.o p3dp_info.o p3dp_AX.o P3DCore.o
SRCS = $(addsuffix .c, $(basename ${OBJS}))
INPUT = 1D66.pdb.bz2
USE_SYSTEM_LIB = -lm -L/usr/X11R6/lib/ -lXext -lX11 -lpthread \
	-lmpi -llam -laio -laio -lutil
USE_MY_LIB = -lAX -lpng -lz -ljpeg -lAtoms -lVecMat3 -lVecMat -lIO -lScalar \
	-lTimer -llapack -lblas -lXpm -lreadline -lhistory -ltermcap
SHELL = sh
GDBTMP = /tmp/${PROG}.gdb
#EXESUFFIX =

include Makefile.inc
include Makefile.common
ifeq (${QUIP_ARCH},)
  include Makefile.arch
else
  include Makefile.${QUIP_ARCH}
endif

USE_SYSTEM_LIB = ${ATOMEYE_LIBS}
LD = ${LINKER}

MY_BIN_DIR = ../bin
MY_LIB_DIR = ../lib
MY_INC_DIR = ../include
MY_BAK_DIR = ../bak

USE_MY_LIB = \
   -lAX -lAtoms -lVecMat3 -lVecMat -lIO -lScalar -lTimer

CFLAGS += -D_$(SYS) ${CCFLAGS} -I${MY_INC_DIR} \
	-I.. -DUSE_CUI ${ATOMEYE_LIB}
LDFLAGS += -L${MY_LIB_DIR} ${USE_MY_LIB} ${SYSTEM_LIB_FLAG} ${USE_SYSTEM_LIB}

.PHONY: all  i install  c clean  b backup  t test  g  r relink  a access

all:	i t
# all:	t #

%.o:	%.c
	${CC} ${CFLAGS} -c $< -o $@
A.o:	A.c Icon/icon.c ${MY_LIB_DIR}/libAX.a ${MY_LIB_DIR}/libAtoms.a ${MY_LIB_DIR}/libVecMat3.a ${MY_LIB_DIR}/libVecMat.a ${MY_LIB_DIR}/libIO.a ${MY_LIB_DIR}/libScalar.a ${MY_LIB_DIR}/libTimer.a
	${CC} ${CFLAGS} -c $< -o $@

$(PROG):$(OBJS) Makefile
	${LD} ${CFLAGS} ${OBJS} -o $@${EXESUFFIX} ${LDFLAGS}
	#@if [ "${CCFLAGS}" != "-g" ] ; then strip $@${EXESUFFIX}; fi; 

i:	install
install:$(PROG)
	-@cp -f $< ${MY_BIN_DIR}; \
	chmod a+rx ${MY_BIN_DIR}/$<; \
	dump=`ls -ld ${MY_BIN_DIR}/$< | sed 's/    */  /g'`; \
	month=`echo $$dump | cut -d' ' -f 6`; \
	prev=`echo "$$dump" | sed "s/$$month.*//"`; \
	ww=`echo "$$prev" | wc -c | sed 's/ //g'`; \
	l=`expr $$ww - 1`; \
	machine_name=`hostname | cut -f 1 -d .`:"               "; \
	machine_name=`echo Current time on "$$machine_name" | cut -c1-$$l`; \
	echo "$$machine_name`date '+%b %d %H:%M' | sed 's/    */  /g'`"; \
	a="@"`echo $$HOME | sed 's/\//~/g'`; \
	ls -l ${MY_BIN_DIR}/$< | sed 's/    */  /g' | \
	sed 's/ \// @~/g'| sed 's/\//~/g' | sed "s/$$a//" | \
	sed -e "s/~/\//g" | sed "s/ \// ~\//" | sed "s/@\//\//"

g:
	@export CC=gcc; export LD=g77; export CCFLAGS=-g; \
	make --environment-overrides --no-print-directory \
	`echo "${MAKECMDGOALS}" | sed -e 's/^g //' | sed -e 's/^g$///' `; \
	exit 99

t:	test
test:	$(PROG)${EXESUFFIX}
	-@echo ' '; \
	echo ========== Running \"$(PROG)${EXESUFFIX} $(INPUT)\" ==========; \
	if [ "${CCFLAGS}" = "-g" ] ; then \
	 echo run $(INPUT) > ${GDBTMP}; gdb ./$(PROG)${EXESUFFIX} -x ${GDBTMP}; \
	else \
	 ./$(PROG)${EXESUFFIX} $(INPUT); \
	fi; \
	echo ======== Done running \"$(PROG)${EXESUFFIX} $(INPUT)\" =======; \
	echo ' '

r:	relink
relink:	
	rm -f $(PROG)${EXESUFFIX}

c:      clean
clean:
	@echo
	-rm -f $(PROG)${EXESUFFIX}; \
	rm -f .\#* \#*\#; \
	rm -f a.out; \
	rm -f *~; \
	rm -f *.o *.png *.jpg *.eps; \
	rm -rf scr_anim Jpg; \
	rm -f *_TEST
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


.PRECIOUS:	%_TEST${EXESUFFIX}
Makefile: ;

a:	access
access:
	-@chmod -R 700 *; \
	chmod 744 *.h; \
	if [ -f ${PROG}${EXESUFFIX} ]; \
	then chmod 755 ${PROG}${EXESUFFIX}; fi; \
	ls -la

${STATIC_LIB_FILE}:	${OBJS}
	-@echo "***************************************************"; \
	rm -f $@; \
	${AR} $@ ${OBJS}; \
	${RANLIB} $@; \
	${NM} $@; \
	echo "***************************************************"
	-@cp -f ${STATIC_LIB_FILE} ${MY_LIB_DIR}; \
	chmod a+r ${MY_LIB_DIR}/${STATIC_LIB_FILE}; \
	cp -f ${HEADER} ${MY_INC_DIR}
