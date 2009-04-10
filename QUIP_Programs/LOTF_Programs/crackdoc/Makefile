F90DOC= ../../libAtoms/f90doc.py

.PHONY: clean 

all: crack.ps crack.pdf

crack.tex: ../crack.f95 ../makecrack.f95 ../cracktools.f95 ../crackparams.f95
	${F90DOC} -t '\textsc{crack} Code Documentation' -a 'James Kermode' -n -l ../crack*.f95 ../makecrack.f95 > crack.tex

crack.dvi: crack.tex
	latex crack.tex
	latex crack.tex
	makeindex general
	latex crack.tex

crack.ps: crack.dvi
	dvips crack.dvi -o crack.ps

crack.pdf: crack.dvi
	dvipdfm crack.dvi

clean: 
	-rm crack.tex crack.ps crack.pdf crack.aux crack.dvi crack.log crack.toc general.idx general.ilg general.ind