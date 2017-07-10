# LaTeX compiler
TEXC := xelatex

default: presentation.md
	echo "Convert Markdown text file to tex file"
	multimarkdown -t beamer presentation.md > presentation.tex
	echo "Compile tex file with xelatex"
	$(TEXC) presentation.tex
	#bibtex presentation 
	#$(TEXC) presentation.tex
	$(TEXC) presentation.tex
	rm -f *.toc *.out *.aux
