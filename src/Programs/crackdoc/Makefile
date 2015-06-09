F90DOC= ../../Tools/f90doc.py

.PHONY: clean 
.SUFFIXES: .eps .pdf

all: crack.ps crack.pdf

crack.tex: ../crack.f95 ../../QUIP_Utils/cracktools.f95 ../../QUIP_Utils/crackparams.f95
	${F90DOC} -t '\textsc{crack} Code Documentation' -a 'James Kermode' -n -l ../crack.f95 ../../QUIP_Utils/cracktools.f95 ../../QUIP_Utils/crackparams.f95 > crack.tex

crack.dvi: crack.tex
	latex crack.tex
	latex crack.tex
	makeindex general
	latex crack.tex

crack.ps: crack.dvi
	dvips crack.dvi -o crack.ps

PDF_FIGS = flow-chart.pdf loading.pdf qm-selection-crack.pdf
pdf_figs: ${PDF_FIGS}

.eps.pdf:
	epstopdf $<

crack.pdf: crack.tex pdf_figs
	sed -i -e 's/\[dvips\]{graphicx}/{graphicx}/' crack.tex
	pdflatex crack.tex
	pdflatex crack.tex
	makeindex general
	pdflatex crack.tex
	sed -i -e 's/{graphicx}/[dvips]{graphicx}/' crack.tex

clean: 
	-rm crack.tex crack.ps crack.pdf crack.aux crack.dvi crack.log crack.toc general.idx general.ilg general.ind ${PDF_FIGS}
