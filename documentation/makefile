pdf :
	pdflatex documentation.tex
	bibtex documentation
	pdflatex documentation.tex
	pdflatex documentation.tex

clean : 
	rm -rfv *.nav *.out *.toc *.aux *.snm *.log *.pdf *.tgz *.bbl *.blg documentation

src :
	which cp
	mkdir documentation ; 
	cp documentation.tex makefile documentation.bib documentation ; 
	tar cvzf documentation.src.tgz documentation/ ; 
	rm -rvf documentation/ ; 

all : src pdf
