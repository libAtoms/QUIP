###########################################
# Threaded Atomistic Configuration Viewer #
#                                         #
# Fri Mar 31 2000 Ju Li <liju99@mit.edu>  #
###########################################

PROG = A
HEADER = $(PROG).h
OBJS  = $(PROG).o primitives.o viewport.o utils.o xtal_shift.o info.o \
	rcut_patch.o scratch.o geo.o LeastSquareStrain.o
SRCS = $(addsuffix .c, $(basename ${OBJS}))
INPUT = 1D66.pdb.bz2
USE_SYSTEM_LIB = -lm -L/usr/X11R6/lib/ -lXext -lX11 -lpthread
USE_MY_LIB = -lAX -lpng -lz -ljpeg -lAtoms -lVecMat3 -lVecMat -lIO -lScalar \
	-lTimer -llapack -lblas -lXpm -lreadline -lhistory -ltermcap
SHELL = sh
GDBTMP = /tmp/${PROG}.gdb
EXESUFFIX =

ifndef SYS
 SYS = Linux
endif
ifdef INC_PATH
 MY_INC_DIR = $(INC_PATH)
else
 MY_INC_DIR = $(HOME)/Co/Include
endif
ifdef LIB_PATH
 MY_LIB_DIR = $(LIB_PATH)
else
 MY_LIB_DIR = $(HOME)/Co/Lib
endif
ifdef BIN_PATH
 MY_BIN_DIR = $(BIN_PATH)
else
 MY_BIN_DIR = $(HOME)/Co/Bin
endif
ifdef BAK_PATH
 MY_BAK_DIRS = $(BAK_PATH)
else
 MY_BAK_DIRS = $(HOME)/Co/Bak
endif

CC = cc
LD = f77 # if we link to f77 compiled libraries
AR = ar -r
RANLIB = ranlib
NM = nm -p
SYSTEM_LIB_FLAG =
ifeq ($(SYS), IRIX64)
 CCFLAGS = -64 -r10000 -TARG:proc=r10000:platform=IP27 -O2 -diag_suppress 1197\
  -LNO:opt=1:fission=2:fusion=2:cs1=32K:cs2=4M:cs2=1024M:ou_further=3:pwr2=ON\
  -OPT:Olimit=0:div_split=ON:roundoff=3:fast_exp=ON:fast_sqrt=ON -s
 LD = f77 -dont_warn_unused
 LDFLAGS = -L${MY_LIB_DIR} ${USE_MY_LIB} -lX11 -lXpm -lXext -lm -lpthread
endif
ifeq ($(SYS), IRIX)
 AR = ar -rc
 RANLIB = echo There is no \"ranlib\" to be used on
 CCFLAGS = -O3
 LD = f77 -dont_warn_unused
 LDFLAGS = -L${MY_LIB_DIR} ${USE_MY_LIB} -lX11 -lXpm -lXext -lm -lpthread
endif
ifeq ($(SYS), SunOS)
 CCFLAGS = -xO4 -fsimple -dalign -fast
 SHELL = bash
endif
ifeq ($(SYS), alpha)
 CCFLAGS = -fast -O3 -fomit-frame-pointer
 CC = ccc
 LD = fort -nofor_main
 USE_MY_LIB = -lAX -lpng -lz -ljpeg -lAtoms -lVecMat3 -lVecMat \
	-lIO -lScalar -lTimer -lcxml
endif
ifeq ($(SYS), HPUX)
 CC = gcc
 LD = gcc
 CCFLAGS = -O6 -ffast-math -funroll-loops -inline-functions
 USE_SYSTEM_LIB = -lm -lpthread \
	-L/opt/fortran90/lib/ -lcl -l:libX11.1 -l:libXext.1 -ldld
 NM = nm --print-file-name --extern-only --defined-only --no-sort --portability
endif
ifeq ($(SYS), OSF1)
 CC = cc -non_shared -pthread
 LD = cc -non_shared -pthread
 CCFLAGS = 
 USE_SYSTEM_LIB = 
 USE_MY_LIB = -lAX -lAtoms -lVecMat3 \
	-lVecMat -lScalar -lIO -lTimer -llapack -lblas \
	-lpthread -lreadline -ljpeg -lXpm -lpng -lz \
	-lm -ldb -lfor -lUfor -lFutil -lX11 -lXext -ldnet_stub -ltermcap
endif
ifeq ($(SYS), CYGWIN)
 CC = gcc
 LD = g77
 CCFLAGS = -O6 -ffast-math -funroll-loops -finline-functions \
	-I/usr/X11R6/include/ -I/usr/local/include/
 NM = nm --print-file-name --extern-only --defined-only --no-sort --portability
 BACKUPDIRS = /asm/home/Moon/Bak /Home/Sun/Bak
 EXESUFFIX = .exe
#  USE_SYSTEM_LIB = -lm -L/usr/X11R6/lib/ -lX11 -L/usr/local/lib/ -lcygipc #
 USE_SYSTEM_LIB = -lm -L/usr/X11R6/lib/ -lX11 
endif
ifeq ($(SYS), Linux)
#  CC = /usr/local/pgcc/bin/gcc #
#  LD = /usr/local/pgcc/bin/g77 -static #
#  CCFLAGS = -Wall -O6 -mpentiumpro -march=pentiumpro -mmx -ffast-math \ #
# 	-funroll-loops -inline-functions -s #
 CC = icc -w
 LD = icc -static
#  CCFLAGS = -g #
 CCFLAGS = -O2 -axW
# --- James stuff? ---
# USE_MY_LIB = -lAX -lpng -lz -ljpeg -lAtoms -lVecMat3 -lVecMat -lIO -lScalar \
#	-lTimer -L/opt/intel/mkl70/lib/32/ -lmkl_lapack -lmkl_ia32 \
#	-lguide -lXpm -lreadline -lhistory -ltermcap -lsvml -ldl
 USE_MY_LIB = -lAX -lpng -lz -ljpeg -lAtoms -lVecMat3 -lVecMat -lIO -lScalar \
       -lTimer -L/opt/intel_mkl_100/lib/em64t -lmkl_lapack -lmkl_em64t \
       -lguide -lXpm -lreadline -lhistory -ltermcap -lsvml -ldl -lnetcdf -lpthread
 USE_SYSTEM_LIB = -lm \
 	-lutil -lXext -lX11 -lxcb -lxcb-xlib -lXau -lpthread
endif
ifeq ($(SYS), Darwin)
 CC = gcc
 LD = cc 
 CCFLAGS = -O3 -fomit-frame-pointer -fno-schedule-insns \
	-fschedule-insns2 -mcpu=750
 USE_MY_LIB = -faltivec -framework vecLib -lAX \
 $(MY_LIB_DIR)/libpng.a $(MY_LIB_DIR)/libz.a $(MY_LIB_DIR)/libjpeg.a -lAtoms \
 -lVecMat3 -lVecMat -lIO -lScalar -lTimer \
 $(MY_LIB_DIR)/libXpm.a $(MY_LIB_DIR)/libreadline.a \
 $(MY_LIB_DIR)/libhistory.a $(MY_LIB_DIR)/libncurses.a
 USE_SYSTEM_LIB = -lm /usr/X11R6/lib/libXext.a /usr/X11R6/lib/libX11.a \
 -lpthread
# one may use 'otool -L' instead of 'ldd' to check dynamic lib #
endif

CFLAGS += -D_$(SYS) ${CCFLAGS} -I${MY_INC_DIR}
LDFLAGS += -L${MY_LIB_DIR} ${USE_MY_LIB} ${SYSTEM_LIB_FLAG} ${USE_SYSTEM_LIB}

.PHONY: all  i install  c clean  b backup  t test  g  r relink  a access

all:	i t
# all:	t #

%.o:	%.c
	${CC} ${CFLAGS} -c $< -o $@
A.o:	A.c Icon/icon.c
	${CC} ${CFLAGS} -c $< -o $@

$(PROG):$(OBJS) Makefile
	${LD} ${CFLAGS} ${OBJS} -o $@${EXESUFFIX} ${LDFLAGS}
	-@if [ "${CCFLAGS}" != "-g" ] ; then strip $@${EXESUFFIX}; fi; 

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

%_TEST:	${SRCS} ${HEADER}
	${CC} -D_$@ ${CFLAGS} -c ${SRCS}
	${LD} ${CFLAGS} ${OBJS} -o $@${EXESUFFIX} ${LDFLAGS}
	chmod u+x $@${EXESUFFIX}

#	Selective test of certain function, say "make Mprint"
%:	%.c
%:	%.o
%.xpm:
%:	%_TEST Makefile
	-@echo ' '; \
	echo ============ Running \"$@_TEST${EXESUFFIX}\" ============; \
	./$@_TEST${EXESUFFIX}; \
	echo ========== Done running \"$@_TEST${EXESUFFIX}\" =========; \
	echo ' '

.PRECIOUS:	%_TEST${EXESUFFIX}
Makefile: ;

a:	access
access:
	-@chmod -R 700 *; \
	chmod 744 *.h; \
	if [ -f ${PROG}${EXESUFFIX} ]; \
	then chmod 755 ${PROG}${EXESUFFIX}; fi; \
	ls -la
