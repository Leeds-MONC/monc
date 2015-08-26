DOCNAME=monc_user_documentation

default: ${DOCNAME}.pdf

${DOCNAME}.pdf:

%.pdf: %.tex
	pdflatex $<
	if fgrep "Rerun to" $*.log; then $(MAKE) --assume-new $^ $@;fi
ifndef NOBIBTEX
	if fgrep "undefined citations"  $*.log;\
	then bibtex $*; \
	$(MAKE) NOBIBTEX=true --assume-new $^ $@;fi
endif

.PHONY:elements

images: 

%.pdf: %.py
	python $<

%.pdf: %.dia	
	dia -e $@ $<

%_tex.pdf: %.fig
	fig2dev -L pdftex $^ $@

%.pdftex_t: %.fig %_tex.pdf
	fig2dev -L pdftex_t -p $*_tex.pdf $^ $@

clean:
	rm -f ${DOCNAME}.pdf *.bbl *.blg *.out *.aux *.log *.toc \
	*.run.xml *-blx.bib
