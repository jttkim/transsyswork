PSFILES 	= athroot_notes.ps
PDFFILES	= $(PSFILES:%.ps=%.pdf)

all : $(PSFILES) $(PDFFILES)

pics : athroot.r rootmodeltest
	R CMD BATCH athroot.r

%.dvi : %.tex pics
	latex $*
	latex $*

%.ps : %.dvi
	dvips -P amz -P cmz -Z -o $*.ps $*

%.pdf : %.ps
	ps2pdf14 -sPAPERSIZE=a4 $< $@

clean :
	rm -f *.dvi *.bbl *.blg *.log *.aux *~ $(PSFILES) $(PDFFILES)

.PHONY : pics all clean

