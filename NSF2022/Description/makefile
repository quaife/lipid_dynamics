filename=description
compileLatexOptions=--quiet --c-style-errors
compilePdfOptions=

compile= ${compPdftex}
compLatex=latex ${compileLatexOptions}
compPdftex=pdflatex ${compilePdfOptions}

all:
	${compile}  ${filename}
	-bibtex ${filename}
	bibtex ${filename}
	${compile} ${filename}
	${compile} ${filename}
	@echo Done ....!

fast:
	${compile} ${filename}

