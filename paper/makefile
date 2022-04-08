pdf: paper.pdf

all: pdf

paper.pdf: paper.tex lmtp.bib
	pdflatex paper
	bibtex paper
	pdflatex paper
	pdflatex paper