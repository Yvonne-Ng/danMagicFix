# Makefile for creating an ATLAS LaTeX document
#------------------------------------------------------------------------------
# By default makes mydocument.pdf using target run_pdflatex
# Replace mydocument with your main filename or add another target set
# Replace BIBTEX = bibtex with BIBTEX = biber if you use biblatex and biber instead of bibtex
# Use "make clean" to cleanup.
# "make cleanall" also deletes the PDF file $(BASENAME).pdf.

LATEX    = latex
PDFLATEX = pdflatex
BIBTEX   = bibtex
# BIBTEX = biber
DVIPS    = dvips
DVIPDF   = dvipdf

BASENAME = BayesianNote

# Default target - make mydocument.pdf with pdflatex
default: run_pdflatex

.PHONY: new draftcover preprintcover auxmat clean help

new:
	sed s/atlas-document/$(BASENAME)/ template/atlas-document.tex >$(BASENAME).tex
	cp template/atlas-document-metadata.tex $(BASENAME)-metadata.tex
	touch $(BASENAME).bib
	touch $(BASENAME)-defs.sty

draftcover:
	cp template/atlas-draft-cover.tex $(BASENAME)-draft-cover.tex
	
preprintcover:
	cp template/atlas-preprint-cover.tex $(BASENAME)-preprint-cover.tex
	
auxmat:
	sed s/atlas-document/$(BASENAME)/ template/atlas-auxmat.tex >$(BASENAME)-auxmat.tex

run_pdflatex: $(BASENAME).pdf
	@echo "Made $<"

run_latex: $(BASENAME).dvi
	$(DVIPDF) -sPAPERSIZE=a4 -dPDFSETTINGS=/prepress ${BASENAME}
	@echo "Made $<"

# Standard Latex targets
%.pdf:	%.tex *.bib
	$(PDFLATEX) $<
	-$(BIBTEX)   $(basename $<)
	$(PDFLATEX) $<
	$(PDFLATEX) $<

%.dvi:	%.tex 
	$(LATEX)    $<
	-$(BIBTEX)   $(basename $<)
	$(LATEX)    $<
	$(LATEX)    $<

%.bbl:	%.tex *.bib
	$(LATEX) $*
	$(BIBTEX) $*
	
help:
	@echo "To create a new document give the commands"
	@echo "make new [BASENAME=mydocument]"
	@echo "make"
	@echo ""
	@echo "If you need a standalone draft cover give the commands:"
	@echo "make draftcover [BASENAME=mydocument]"
	@echo "pdflatex mydocument-draft-cover"
	@echo ""
	@echo "If you need a standalone preprint cover give the commands:"
	@echo "make preprintcover [BASENAME=mydocument]"
	@echo "pdflatex mydocument-preprint-cover"
	@echo ""
	@echo "If you need a document for auxiliary material give the commands:"
	@echo "make auxmat [BASENAME=mydocument]"
	@echo "pdflatex mydocument-auxmat"

%.ps:	%.dvi
	$(DVIPS) $< -o $@

clean:
	-rm *.dvi *.toc *.aux *.log *.out \
		*.bbl *.blg *.brf *.bcf *.run.xml \
		*.cb *.ind *.idx *.ilg *.inx \
		*.synctex.gz *~ ~* spellTmp 
	
cleanall: clean
	-rm $(BASENAME).pdf 
	-rm $(BASENAME)-draft-cover.pdf $(BASENAME)-preprint-cover.pdf
	-rm $(BASENAME)-auxmat.pdf
