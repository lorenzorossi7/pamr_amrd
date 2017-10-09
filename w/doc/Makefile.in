.IGNORE:
# Makefile for "PAMR_ref.ps, AMRD_ref.ps and Wiley_chapter.ps"

DOC1=AMRD_ref
DOC2=PAMR_ref
DOC3=Wiley_chapter

# Add figure files here
FIGS1= 
FIGS2=
FIGS3=

all: $(DOC1).ps $(DOC2).ps $(DOC3).ps

$(DOC1).ps: $(DOC1).tex $(FIGS)
	latex $(DOC1)
	latex $(DOC1)
	latex $(DOC1)
	dvips -tletter -o $(DOC1).ps $(DOC1)
	ls -lt $(DOC1).ps

$(DOC2).ps: $(DOC2).tex $(FIGS)
	latex $(DOC2)
	latex $(DOC2)
	latex $(DOC2)
	dvips -tletter -o $(DOC2).ps $(DOC2)
	ls -lt $(DOC2).ps

$(DOC3).ps: $(DOC3).tex $(FIGS)
	latex $(DOC3)
	latex $(DOC3)
	latex $(DOC3)
	dvips -tletter -o $(DOC3).ps $(DOC3)
	ls -lt $(DOC3).ps

clean:
	/bin/rm $(DOC1).ps
	/bin/rm $(DOC1).aux
	/bin/rm $(DOC1).log
	/bin/rm $(DOC1).dvi
	/bin/rm $(DOC1).end
	/bin/rm $(DOC2).ps
	/bin/rm $(DOC2).aux
	/bin/rm $(DOC2).log
	/bin/rm $(DOC2).dvi
	/bin/rm $(DOC2).end
	/bin/rm $(DOC3).ps
	/bin/rm $(DOC3).aux
	/bin/rm $(DOC3).log
	/bin/rm $(DOC3).dvi
	/bin/rm $(DOC3).end
	/bin/rm $(DOC3).inx
	/bin/rm Makefile
