.PHONY: all clean

all: manuscript.pdf

manuscript.pdf: manuscript.tex refs.bib
	pdflatex manuscript.tex
	bibtex manuscript
	pdflatex manuscript.tex
	pdflatex manuscript.tex

clean:
	rm -f *.aux *.log *.bbl *.blg *.out *.toc *.pdf