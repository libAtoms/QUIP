ifdef BAK_PATH
 MY_BAK_DIRS = $(BAK_PATH)
else
 MY_BAK_DIRS = $(HOME)/Co/Bak
endif

.PHONY: all  i install  t  c clean  b backup  p

all:	install

i:	install
install:
	cd DSpectral; $(MAKE) --no-print-directory install
	cd SSpectral; $(MAKE) --no-print-directory install

c:	clean
clean:
	@echo
	-rm -f *~; \
	rm -f .\#*; \
	rm -f \#*\#; \
	rm -f *.o; \
	rm -f a.out; \
	rm -f Notes/*~; \
	rm -f core
	@echo
	cd DSpectral; $(MAKE) --no-print-directory clean
	cd SSpectral; $(MAKE) --no-print-directory clean

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

p:
	cd DSpectral; $(MAKE) p
	cd SSpectral; $(MAKE) p