PS = main.ps
DVI = $(addsuffix .dvi, $(basename $(PS)))
PDF = $(addsuffix .pdf, $(basename $(PS)))
AUX = $(addsuffix .aux, $(basename $(PS)))
TEX = $(addsuffix .tex, $(basename $(PS)))
BASE = $(basename $(PS))
LATEX = /usr/bin/latex
BIBTEX = bibtex
DVIPS = /usr/bin/dvips
DVIPDFM = dvipdfm
XDVI = /usr/bin/xdvi -s 4 -gamma 1.2 -paper 8.5x11in +statusline
GHOSTVIEW = gv
ACROREAD = acc
BACKUPDIRS = /home/Moon/Bak /Home/Sun/Bak
SHELL = sh
CMDHIST = .cmd
EXTRA = $(wildcard c1/*/main.tex)

CHAP = $(subst /,,$(wildcard */))
PHONY = all  p  g ghostview  x xdvi  c clean  b backup
.PHONY: $(PHONY) $(CHAP) $(CMDHIST)

ifdef chap  # environment variable: setenv chap "c1 c2"
 COMMA :=,
 EMPTY :=
 SPACE := $(EMPTY) $(EMPTY)
 CMD = $(subst $(SPACE),$(COMMA),$(strip $(addsuffix /main,$(chap))))
 SRC = $(addsuffix /main.tex, $(chap))
 FIGURE = $(wildcard $(addsuffix /*.eps, $(chap)))
else        # CMD: "c1/main,c2/main", or RETURN to include all
 CMD =
 SRC = $(wildcard */*.tex) $(addsuffix .tex, $(wildcard ../*/*.m))
 FIGURE = $(wildcard */*.eps)
endif

all:	x

p:	$(PDF)
	@pid=`ps | grep -v grep | grep $(PDF) | cut -d" " -f1,2`; \
	if [ -n "$$pid" ]; then \
	 kill -s SIGHUP $$pid; \
	else \
	 sh -c "$(ACROREAD) $(PDF) &"; \
	fi

g:	ghostview
ghostview: $(PS)
	@pid=`ps | grep -v grep | grep $(PS) | cut -d" " -f1,2`; \
	if [ -n "$$pid" ]; then \
	 kill -s SIGHUP $$pid; \
	else \
	 sh -c "$(GHOSTVIEW) $(PS) &"; \
	fi

x:	xdvi
xdvi:	$(DVI)
	@pid=`ps --cols 300 xw | grep -v grep | grep $$PWD/$(DVI) | sed -e 's/ [\?p].*//g'`; \
	if [ -n "$$pid" ]; then \
	 kill -s USR1 $$pid; \
	else \
	 sh -c "$(XDVI) $$PWD/$(DVI) > /dev/null &"; \
	fi

$(PDF):	$(DVI)
	$(DVIPDFM) -o $(PDF) $(DVI)
$(PS):	$(DVI)
	$(DVIPS) $(DVI) -o $(PS)
%.m.tex: %.m
	lgrind -i -c -d ~/Co/Runtime/LatexStyles/lgrindef -lMATLAB $< > $@
$(DVI): $(CMDHIST) $(AUX);
$(CMDHIST):
	@if [ -e $(CMDHIST) ]; then \
	 if [ "`cat $(CMDHIST)`" != "$(CMD)" ]; then \
	  touch $(TEX); \
	  echo "$(CMD)" > $(CMDHIST); \
	 fi; \
	else \
	 echo "$(CMD)" > $(CMDHIST); \
	fi
$(AUX): $(FIGURE) $(SRC) $(TEX) $(EXTRA)
	@let n=0; \
	while ( ( echo $(CMD) | $(LATEX) $(TEX) || rm -f $(AUX) ) \
	| tee /dev/tty | egrep 'Rerun|undefined' > /dev/null ) do \
	if [ $$n -lt 1 ]; then $(BIBTEX) $(BASE); fi; \
	let n=($$n+1); if [ $$n -gt 4 ]; then break; fi; \
	echo ' '; echo '** RUNNING LATEX AGAIN **'; echo ' '; done

c:	clean
clean:
	@echo
	-rm -f $(DVI); \
	rm -f *~; \
	rm -f \#*\#; \
	rm -f *.aux; \
	rm -f *.log; \
	rm -f *.toc; \
	rm -f *.lof; \
	rm -f *.lot; \
	rm -f *.bbl; \
	rm -f *.blg; \
	rm -f $(CMDHIST)
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
	nbytes=`ls -l $$tmp.gz | sed -e 's/  */ /g' | cut -d" " -f 5`; \
	echo $$tmp.gz = $$nbytes bytes.; \
	for backup in $(BACKUPDIRS); do \
	 if [ -d $$backup ]; then \
	  echo ' '; \
	  echo "copying $$tmp.gz to $$backup/$$dir.tar.gz ..."; \
	  cp -i $$tmp.gz $$backup/$$dir.tar.gz; \
	 fi; \
	done; \
	echo

$(CHAP):
	@export chap="$(filter $(CHAP), $(MAKECMDGOALS))"; \
	make --no-print-directory $(filter $(PHONY), $(MAKECMDGOALS)); \
	exit 99
