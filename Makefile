# latex Makefile
LATEX=pdflatex -interaction nonstopmode
BIBTEX=bibtex
DVIPS=dvips
TARGETS=proj_description.pdf

all: $(TARGETS)

clean: 
	@rm -f *.aux *.bbl *.blg *.dvi *.log

proj_description.pdf: proj_description.tex
	${LATEX} proj_description.tex 
	${BIBTEX} proj_description
	${LATEX} proj_description.tex
	${BIBTEX} proj_description
	${LATEX} proj_description.tex















