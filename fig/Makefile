
-include target.mk

######################################################################

Sources = Makefile .gitignore README.md LICENSE.md

-include makestuff/static.mk

# include $(ms)/perl.def

##################################################################

## Content

Sources += diagram.tex

diagram.pdf: diagram.tex

# pull: cleanfigs
cleanfigs:
	git checkout master -- diagram.pdf

######################################################################

## Recrunching
## Start building a pipeline for JD; still focused on poster, though

Sources += corrected.R corrected_hist.R corrected_box.R

corrected.Rout: corrected.R
corrected_hist.Rout: corrected.Rout corrected_hist.R
corrected_box.Rout: corrected.Rout corrected_box.R

######################################################################

## Files

Makefile: makestuff

## If this rule doesn't work for you; copy or link makestuff to here
makestuff:
	/bin/ln -s ~/git/makestuff $@

-include $(ms)/git.mk
-include $(ms)/visual.mk

-include $(ms)/newlatex.mk
-include $(ms)/pandoc.mk

-include $(ms)/wrapR.mk
