# Makefile for Sphinx documentation
#
MAIN=SubDyn
INKSCAPE="C:/Program Files/Inkscape/inkscape.exe"

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
PAPER         =

# Internal variables.
PAPEROPT_a4     = -D latex_paper_size=a4
PAPEROPT_letter = -D latex_paper_size=letter
ALLSPHINXOPTS   = -d _build/doctrees $(PAPEROPT_$(PAPER)) $(SPHINXOPTS) .


SVGDIR=figs
SVGOUTDIR=figs
SVGS=$(notdir $(wildcard $(SVGDIR)/*.svg))
SVGS2PDFS=$(patsubst %,$(SVGOUTDIR)/%,$(SVGS:.svg=.pdf))
SVGS2PNGS=$(patsubst %,$(SVGOUTDIR)/%,$(SVGS:.svg=.png))
#INKSCAPE=inkscape


.PHONY: help clean html dirhtml pickle json htmlhelp qthelp latex changes linkcheck doctest

all: pdf

doc2rst:
	pandoc -F pandoc-crossref -F pandoc-citeproc -s -t rst --toc SubDyn_Manual_Rev037.docx -o output.rst   --bibliography=references_SD.bib

# --wrap=preserve
# -F pandoc-crossref
# -F pandoc-eqnos
# --number-section
	

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  pdf       to make pdf from LaTeX files (uses latex)"
	@echo "  html      to make standalone HTML files"
	@echo "  dirhtml   to make HTML files named index.html in directories"
	@echo "  pickle    to make pickle files"
	@echo "  json      to make JSON files"
	@echo "  htmlhelp  to make HTML files and a HTML help project"
	@echo "  qthelp    to make HTML files and a qthelp project"
	@echo "  latex     to make LaTeX files, you can set PAPER=a4 or PAPER=letter"
	@echo "  changes   to make an overview of all changed/added/deprecated items"
	@echo "  linkcheck to check all external links for integrity"
	@echo "  doctest   to run all doctests embedded in the documentation (if enabled)"

clean:
	-rm -rf _build/*

pdf: latex pdf-compile 
# 	diff
# 	cd _build/latex && make all-pdf

pdf-compile: 
	cd _build/latex && pdflatex $(MAIN).tex
# 	cd _build/latex && make all-pdf

pdf-html: pdf html
	cp _build/latex/sampledoc.pdf _build/html

diff:
	cd _build/latex  && latexdiff -p ../../LatexDiffPreamble.tex $(MAIN)-old.tex $(MAIN).tex > $(MAIN)-diff.tex
	cd _build/latex  && cp $(MAIN).aux $(MAIN)-diff.aux
	cd _build/latex  && pdflatex --interaction=nonstopmode $(MAIN)-diff.tex   

# Rule to create pdf or png from svg
$(SVGOUTDIR)/%.pdf:$(SVGDIR)/%.svg 
	$(INKSCAPE) -A -D --file="$<" "$@"
$(SVGOUTDIR)/%.png:$(SVGDIR)/%.svg 
	$(INKSCAPE) -D --export-dpi 300  --file="$<" -e "$@"


html:
	$(SPHINXBUILD) -b html $(ALLSPHINXOPTS) _build/html
	@echo
	@echo "Build finished. The HTML pages are in _build/html."

dirhtml:
	$(SPHINXBUILD) -b dirhtml $(ALLSPHINXOPTS) _build/dirhtml
	@echo
	@echo "Build finished. The HTML pages are in _build/dirhtml."

pickle:
	$(SPHINXBUILD) -b pickle $(ALLSPHINXOPTS) _build/pickle
	@echo
	@echo "Build finished; now you can process the pickle files."

json:
	$(SPHINXBUILD) -b json $(ALLSPHINXOPTS) _build/json
	@echo
	@echo "Build finished; now you can process the JSON files."

htmlhelp:
	$(SPHINXBUILD) -b htmlhelp $(ALLSPHINXOPTS) _build/htmlhelp
	@echo
	@echo "Build finished; now you can run HTML Help Workshop with the" \
	      ".hhp project file in _build/htmlhelp."

qthelp:
	$(SPHINXBUILD) -b qthelp $(ALLSPHINXOPTS) _build/qthelp
	@echo
	@echo "Build finished; now you can run "qcollectiongenerator" with the" \
	      ".qhcp project file in _build/qthelp, like this:"
	@echo "# qcollectiongenerator _build/qthelp/sampledoc.qhcp"
	@echo "To view the help file:"
	@echo "# assistant -collectionFile _build/qthelp/sampledoc.qhc"

latex: $(SVGS2PNGS)
	cp -r figs _build/latex
	$(SPHINXBUILD) -b latex $(ALLSPHINXOPTS) _build/latex
	@echo
	@echo "Build finished; the LaTeX files are in _build/latex."
	@echo "Run \`make all-pdf' or \`make all-ps' in that directory to" \
	      "run these through (pdf)latex."

changes:
	$(SPHINXBUILD) -b changes $(ALLSPHINXOPTS) _build/changes
	@echo
	@echo "The overview file is in _build/changes."

linkcheck:
	$(SPHINXBUILD) -b linkcheck $(ALLSPHINXOPTS) _build/linkcheck
	@echo
	@echo "Link check complete; look for any errors in the above output " \
	      "or in _build/linkcheck/output.txt."

doctest:
	$(SPHINXBUILD) -b doctest $(ALLSPHINXOPTS) _build/doctest
	@echo "Testing of doctests in the sources finished, look at the " \
	      "results in _build/doctest/output.txt."

