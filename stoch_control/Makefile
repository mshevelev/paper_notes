# Makefile for LaTeX project

MAIN = main
PDF = $(MAIN).pdf
SRC = $(MAIN).tex

all: $(PDF)

$(PDF): $(SRC)
	pdflatex $(MAIN).tex
	pdflatex $(MAIN).tex

clean:
	rm -f *.aux *.log *.out *.pdf *.toc

.PHONY: all clean
