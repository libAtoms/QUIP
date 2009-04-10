#########################################
# Fastest Double Array Sorting Routines #
#########################################

PROG = match
HEADER = $(PROG).h
OBJ = $(PROG).o nr.o jp.o glibc.o merge.o
SRC = $(addsuffix .c, $(basename ${OBJ}))
USE_SYSTEM_LIB = -lm
USE_MY_LIB = -lTimer -lScalar
SHELL = sh

# powerful but stupid batch execution machine
EXE_MACHINE = mmm
# network source development directory [mount-point at execution machine]
CODE_BASE = /home/Moon/Z/
# reliable local directory to which source will be transfered and compiled
EXE_MACHINE_BASE = /usr/people/liju99/Z/

SYSTEM_LIB_FLAG =
ifeq (${HOST}, $(EXE_MACHINE))
 CC = cc  # SGI O2 options
 CCFLAGS = -64 -r10000 -TARG:proc=r10000:platform=IP27 -O3\
  -LNO:opt=1:fission=2:fusion=2:cs1=32K:cs2=4M:cs2=1024M:ou_further=3:pwr2=ON\
  -OPT:Olimit=0:div_split=ON:roundoff=3:fast_exp=ON:fast_sqrt=ON
 LD = f77
 BACKUPDIRS = /home/Moon/Bak /usr/people/liju99/Bak
else # media-savvy source development workstation
 CC = gcc # gcc options
 CCFLAGS = -O3 -fforce-mem -fforce-addr -finline-functions -ffast-math\
           -fstrength-reduce -fthread-jumps -funroll-loops -fcse-follow-jumps \
           -fcse-skip-blocks -frerun-cse-after-loop -fexpensive-optimizations
#  CCFLAGS = -g #
 LD = g77
 BACKUPDIRS = /home/Moon/Bak /Home/Sun/Bak
endif

CFLAGS = -D_$(SYS) ${CCFLAGS} -I$(INC_PATH)
LDFLAGS = ${SYSTEM_LIB_FLAG} ${USE_SYSTEM_LIB} -L${LIB_PATH} ${USE_MY_LIB}

.PHONY:	all  r run  c clean  b backup  t transfer

all: $(PROG)

$(PROG): $(OBJ) Makefile
	${FC} ${CFLAGS} -o $(PROG) $(OBJ) $(LDFLAGS)

# default make rules
%.o:	%.c ${HEADER} Makefile
	${CC} ${CFLAGS} -c $< -o $@

r:	run
run:	all
	-@echo " "; $(PROG); echo " "

c:	clean
clean:
	@echo
	-rm -f $(PROG); \
	rm -f a.out; \
	rm -f *~; \
	rm -f *.o
	@echo

# http://long-march.mit.edu/liju99/Misc/Shell/but
b:	backup
backup: clean
	@tmp=/tmp/backupto.tar; \
	echo -e "taring the directory $(PWD) ..."; \
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
	for backup in $(BACKUPDIRS); do \
	 if [ -d $$backup ]; then \
	  echo -e "\ncopying $$tmp.gz to $$backup/$$dir.tar.gz ..."; \
	  cp -i $$tmp.gz $$backup/$$dir.tar.gz; \
	 fi; \
	done; \
	echo

t:	transfer
transfer:
ifeq (${HOST}, $(EXE_MACHINE))
	-cd $(CODE_BASE); \
	cp -f Makefile ${SRC} $(CON) $(EXE_MACHINE_BASE); \
	cd $(EXE_MACHINE_BASE); chmod 400 ${SRC};
else
	@-echo -e "\nCurrent machine is ${HOST} ..."; \
	echo "================ trying to connect $(EXE_MACHINE)\
	via RSH ================"; \
	rsh mmm "cd $(CODE_BASE); \
	cp -f Makefile ${SRC} $(CON) $(EXE_MACHINE_BASE); \
	cd $(EXE_MACHINE_BASE); chmod 400 ${SRC}; \
	echo  -----------------------------------------------------------; \
	echo The current time on $(EXE_MACHINE) is'   '`date`; \
	echo And in $(EXE_MACHINE_BASE) ...; \
	ls -l Makefile ${SRC} $(CON)"; \
	echo ===============================================================;\
	echo If the two times do not agree, run the following:;\
	echo % telnet mmm -l $(USER);\
	echo -e "% cd $(CODE_BASE); make t\n"
endif
