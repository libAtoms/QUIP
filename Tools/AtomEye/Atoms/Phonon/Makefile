###########################################
# libPhonon: -llapack -lblas -lm          #
#            -lVecMat -lScalar -lIO       #
#                                         #
# Force constant, phonon diagonalization, #
# and thermodynamic functions.            #
#                                         #
# Jul 9 2000 Ju Li <liju99@mit.edu>.      #
###########################################

LIB = Phonon
HEADER = ${LIB}.h # singular
STATIC_LIB_FILE = lib${LIB}.a
PROGS = cv dk dgamma
OBJS = ${LIB}.o Debye.o # plural
SRCS = $(addsuffix .c, $(basename ${OBJS}))
ADD_DEPENDENCE = Makefile
USE_SYSTEM_LIB = -lm
USE_MY_LIB = -lVecMat -lScalar -lIO -llapack -lblas
SHELL = sh
GDBTMP = /tmp/${LIB}.gdb
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
 CCFLAGS = -64 -r10000 -TARG:proc=r10000:platform=IP27 -O3 -diag_suppress 1197\
  -LNO:opt=1:fission=2:fusion=2:cs1=32K:cs2=4M:cs2=1024M:ou_further=3:pwr2=ON\
  -OPT:Olimit=0:div_split=ON:roundoff=3:fast_exp=ON:fast_sqrt=ON -s
endif
ifeq ($(SYS), IRIX)
 AR = ar -rc
 RANLIB = echo There is no \"ranlib\" to be used on
 CCFLAGS = -O3
 LD = f77 -dont_warn_unused
endif
ifeq ($(SYS), SunOS)
 CCFLAGS = -xO4 -fsimple -dalign -fast
 SHELL = bash
endif
ifeq ($(SYS), alpha)
 CCFLAGS = -fast -O3 -fomit-frame-pointer
 CC = ccc
 LD = fort -nofor_main
 USE_MY_LIB = -lVecMat -lScalar -lIO -lcxml
endif
ifeq ($(SYS), HPUX)
 CC = gcc
 LD = gcc
 CCFLAGS = -O6 -ffast-math -funroll-loops -inline-functions
 USE_SYSTEM_LIB = -lm -L/opt/fortran90/lib/ -lcl
 NM = nm --print-file-name --extern-only --defined-only --no-sort --portability
endif
ifeq ($(SYS), CYGWIN)
 CC = gcc
 LD = g77
 CCFLAGS = -O6 -ffast-math -funroll-loops -inline-functions
 NM = nm --print-file-name --extern-only --defined-only --no-sort --portability
 BACKUPDIRS = /asm/home/Moon/Bak /Home/Sun/Bak
 EXESUFFIX = .exe
endif
ifeq ($(SYS), Linux)
 CC = /usr/local/pgcc/bin/gcc
 LD = /usr/local/pgcc/bin/g77 
 CCFLAGS = -O6 -mpentiumpro -march=pentiumpro -mmx -ffast-math \
	-funroll-loops -inline-functions -s
 NM = nm --print-file-name --extern-only --defined-only --no-sort --portability
endif

CFLAGS = -D_$(SYS) ${CCFLAGS} -I${MY_INC_DIR}
LDFLAGS = ${SYSTEM_LIB_FLAG} ${USE_SYSTEM_LIB} -L${MY_LIB_DIR} ${USE_MY_LIB}

.PHONY: all  i install  c clean  b backup  r  g  p  h  a access

all:	install

%.o:	%.c ${HEADER}  ${ADD_DEPENDENCE}
	${CC} ${CFLAGS} -c $< -o $@

${STATIC_LIB_FILE}:	${OBJS}
	-@echo "***************************************************"; \
	rm -f $@; \
	${AR} $@ ${OBJS}; \
	${RANLIB} $@; \
	${NM} $@; \
	echo "***************************************************"

i:	install
install:	${STATIC_LIB_FILE}
	-@cp -f $< ${MY_LIB_DIR}; \
	chmod a+rx ${MY_LIB_DIR}/$<; \
	cp -f ${HEADER} ${MY_INC_DIR}; \
	chmod a+rx ${MY_INC_DIR}/${HEADER}; \
	dump=`ls -ld ${MY_INC_DIR}/${HEADER} | sed 's/    */  /g'`; \
	month=`echo $$dump | cut -d' ' -f 6`; \
	prev=`echo "$$dump" | sed "s/$$month.*//"`; \
	ww=`echo "$$prev" | wc -c | sed 's/ //g'`; \
	l=`expr $$ww - 1`; \
	machine_name=`hostname | cut -f 1 -d .`:"               "; \
	machine_name=`echo Current time on "$$machine_name" | cut -c1-$$l`; \
	echo "$$machine_name`date '+%b %d %H:%M' | sed 's/    */  /g'`"; \
	a="@"`echo $$HOME | sed 's/\//~/g'`; \
	ls -l ${MY_INC_DIR}/${HEADER} | sed 's/    */  /g' | \
	sed 's/ \// @~/g'| sed 's/\//~/g' | sed "s/$$a//" | \
	sed -e "s/~/\//g" | sed "s/ \// ~\//" | sed "s/@\//\//"; \
	ls -l ${MY_LIB_DIR}/$< | sed 's/    */  /g' | \
	sed 's/ \// @~/g'| sed 's/\//~/g' | sed "s/$$a//" | \
	sed -e "s/~/\//g" | sed "s/ \// ~\//" | sed "s/@\//\//"

h:
	cp -f ${HEADER} ${MY_INC_DIR}

r:	
	rm -f *_TEST${EXESUFFIX};

c:	clean
clean:
	@echo
	-rm -f *~; \
	rm -f \#*\#; \
	rm -f *.o; \
	rm -f ${STATIC_LIB_FILE}; \
	rm -f *_TEST${EXESUFFIX}; \
	rm -f *.out; \
	rm -f $(addsuffix ${EXESUFFIX}, ${PROGS}); \
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

%_TEST: $(SRCS) ${HEADER}
	for obj in $(OBJS); do \
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
	chmod u+x $@${EXESUFFIX}

#	Selective test of certain function, say "make spectral_sumcorr"
%:	%.f
%:	%.c
%:	%.o
%:	%_TEST
	-@echo ' '; \
	echo ============ Running \"$@_TEST${EXESUFFIX}\" ============; \
	if [ "${CCFLAGS}" = "-g" ] ; then \
	 echo run > ${GDBTMP}; gdb ./$@_TEST${EXESUFFIX} -x ${GDBTMP}; \
	else \
	 ./$@_TEST${EXESUFFIX}; \
	fi; \
	echo ========== Done running \"$@_TEST${EXESUFFIX}\" =========; \
	echo ' '

.PRECIOUS:	%_TEST${EXESUFFIX}
Makefile: ;

a:	access
access:
	-@chmod -R 700 *; \
	chmod 744 *.h; \
	if [ -f ${STATIC_LIB_FILE} ]; then chmod 744 ${STATIC_LIB_FILE}; fi; \
	ls -la
