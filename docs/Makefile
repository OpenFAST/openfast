
# You can set these variables from the command line.
ALLSPHINXOPTS =  -c conf.py 
SPHINXBUILD   = sphinx-build
PAPER         =
ALLSPHINXOPTS   = -d _build/doctrees  $(SPHINXOPTS) .

BUILDDIR = _build


all: html

$(BUILDDIR):
	mkdir $(BUILDDIR)

html: $(BUILDDIR)
	$(SPHINXBUILD) -b html $(ALLSPHINXOPTS) $(BUILDDIR) 
	@echo
	@echo "Build finished. The HTML pages are in _build/html."
