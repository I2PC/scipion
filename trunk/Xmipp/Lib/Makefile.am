OPTIONAL_LIBS =

if BUILD_QT
OPTIONAL_LIBS += XmippGraphics
endif

SUBDIRS = XmippData Reconstruction Classification XmippInterface \
	OldXmipp $(OPTIONAL_LIBS)

docs:
	for dir in $(SUBDIRS); do cd $$dir; $(MAKE) docs; cd ..; done
