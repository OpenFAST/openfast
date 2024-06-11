MAIN=BeamDynInput

SVGTEXDIR=svgtex
SVGTEXOUTDIR=figs_svgtex
SVGTEXS=$(notdir $(wildcard $(SVGTEXDIR)/*.svg))
SVGTEXS2PDFS=$(patsubst %,$(SVGTEXOUTDIR)/%,$(SVGTEXS:.svg=.pdf))


INKSCAPE=inkscape
RM=rm -f
MV=mv
SLASH="/"
ifdef OS
	# Windows specific
	#   INKSCAPE="C:/Program Files/Inkscape/inkscape.exe"
   RM=del
   MV=move /Y 
   SLASH="\\"
endif

all: figspdf pdf

pdf: figspdf
	pdflatex -synctex=1 --file-line-error-style --shell-escape --interaction=nonstopmode $(MAIN)
	bibtex $(MAIN).aux

bibtex:
	bibtex $(MAIN).aux

figspdf: $(SVGTEXOUTDIR) $(SVGTEXS2PDFS)

$(SVGTEXOUTDIR)/%.pdf:$(SVGTEXDIR)/%.svg 
	$(INKSCAPE) -z -D -o "$(SVGTEXOUTDIR)/$*.pdf" --export-latex   "$<"
	$(MV) "$(SVGTEXOUTDIR)$(SLASH)$*.pdf_tex" "$(SVGTEXOUTDIR)$(SLASH)$*.tex"
# 	$(INKSCAPE) -z -D --file="$<" --export-pdf="$(SVGTEXOUTDIR)/$*.pdf" --export-latex

$(SVGTEXOUTDIR):
	mkdir $@

clean:
	$(RM) *.log *.aux *.out *.bbl *.blg *.toc
