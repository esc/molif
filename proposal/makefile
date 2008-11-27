pdf :
	pdflatex proposal.tex
	bibtex proposal
	pdflatex proposal.tex
	pdflatex proposal.tex

clean : 
	rm -v *.nav *.out *.toc *.aux *.snm *.log *.pdf *.tgz *.bbl *.blg

src : 
	mkdir proposal ; 
		cp proposal.tex \
		makefile \
		-t proposal ; 
	tar cvzf proposal.src.tgz proposal/ ; 
	rm -rvf proposal/ ; 

all : src pdf
