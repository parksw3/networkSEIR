
-include target.mk

######################################################################

Sources = Makefile .gitignore README.md LICENSE.md

-include makestuff/static.mk

# include $(ms)/perl.def

##################################################################

## Content

Sources += notes.md

Sources += poster.Rnw

poster.pdf: poster.Rnw

poster.pdf: poster.bbl

pull: cleanposter

cleanposter:
	git checkout master -- poster.pdf

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

# -include $(ms)/wrapR.mk
