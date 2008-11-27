pdf :
	pdflatex proposal.tex
	bibtex proposal
	pdflatex proposal.tex
	pdflatex proposal.tex

clean : 
	rm -rfv *.nav *.out *.toc *.aux *.snm *.log *.pdf *.tgz *.bbl *.blg proposal

src :
	which cp
	mkdir proposal ; 
	cp proposal.tex makefile proposal.bib proposal ; 
	tar cvzf proposal.src.tgz proposal/ ; 
	rm -rvf proposal/ ; 

all : src pdf
