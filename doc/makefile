all: simulate.pdf

%.pdf: %.tex
	pdflatex -shell-escape $*.tex
	if ( grep -q citation $*.aux ) ; then \
		bibtex $* ; \
		pdflatex -shell-escape $*.tex ; \
	fi
	pdflatex -shell-escape $*.tex

clean:
	rm -f *.aux *.bbl *.blg *.log *~
