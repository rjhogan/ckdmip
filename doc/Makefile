# If you have pdflatex installed, type "make" to create the
# documentation, "make clean" to delete it

BASE=ckdmip_documentation

documentation:
	pdflatex $(BASE).tex
	pdflatex $(BASE).tex
	pdflatex $(BASE).tex

clean:
	rm -f $(BASE).pdf ${BASE}.aux ${BASE}.out ${BASE}.log ${BASE}.toc

.PHONY: documentation clean
