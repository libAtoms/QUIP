############################################
# libAX: -lX11 -lXext -lpng -lz -ljpeg -lm #
#        -lreadline -lncurses              #
#        -lScalar -lIO -lTimer             #
#                                          #
# Accelerated low-level graphics library   #
# with X Window Shared Memory Extension.   #
#                                          #
# Jan.13, 2000 Ju Li <liju99@mit.edu>      #
############################################

# singular
LIB = AX
HEADER = ${LIB}.h
STATIC_LIB_FILE = lib${LIB}.a
PROGS = cmap
# plural
OBJS = AX.o pixel.o jpg2eps.o Sorter.o Scan.o 3D.o
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


CFLAGS += -D_$(SYS) ${CCFLAGS} -I${MY_INC_DIR}
LDFLAGS = ${SYSTEM_LIB_FLAG} ${USE_SYSTEM_LIB} -L${MY_LIB_DIR} ${USE_MY_LIB}

.PHONY: all  i install  c clean  b backup  r  g  h  p  a access

all:	install

%.o:	%.c
	${CC} ${CFLAGS} -c $< -o $@
Scan3D.o: Scan.c 3D.c
	${CC} ${CFLAGS} -c Scan3D.c -o $@

${STATIC_LIB_FILE}:	${OBJS}
	-@echo "***************************************************"; \
	rm -f $@; \
	${AR} $@ ${OBJS}; \
	${RANLIB} $@; \
	${NM} $@; \
	echo "***************************************************"

i:	install
install:	${STATIC_LIB_FILE}
	-@cp -f $< ${MY_LIB_DIR}/; \
	chmod a+rx ${MY_LIB_DIR}/$<; \
	cp -f ${HEADER} ${MY_INC_DIR}/; \
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
	cp -f ${HEADER} ${MY_INC_DIR}/

r:	
	rm -f *_TEST${EXESUFFIX};

c:	clean
clean:
	@echo
	-rm -f $(PROGS); \
	rm -f *~; \
	rm -f \#*\#; \
	rm -f .\#*; \
	rm -f *.o; \
	rm -f ${STATIC_LIB_FILE}; \
	rm -f *_TEST${EXESUFFIX}; \
	rm -f a.out; \
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
	${MAKE} --environment-overrides --no-print-directory \
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


#	Selective test of certain function, say "make Mprint"
%:	%.f
%:	%.c
%:	%.o
%:	%_TEST ${HEADER}
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
