lib_LTLIBRARIES = \
	libXmippInterface.la \
	libXmippInterfaceINRIA.la

libXmippInterface_la_SOURCES = \
	aph3d.cpp \
	aph.cpp \
	aph_origmerg.cpp \
	cc_lattice_io.cpp \
	ccp4.cpp \
	crisp.cpp \
	em.cpp \
	jdl.cpp \
	opendxang.cpp \
	opendx.cpp \
	spider.cpp \
	virus.cpp \
	vrml.cpp

libXmippInterfaceINRIA_la_SOURCES = \
	inria.cpp

EXTRA_DIST = \
	aph3d.h \
	aph.h \
	aph_origmerg.h \
	cc_lattice_io.h \
	ccp4.h \
	crisp.h \
	em.h \
	inria.h \
	jdl.h \
	opendxang.h \
	opendx.h \
	spider.h \
	virus.h \
	vrml.h

INCLUDES = -I$(top_srcdir)/libraries

docs:
	@echo "Creating documentation in XmippInterface ..."
	@doc++ -f -B @XMIPPDIR@/Extra_Docs/banner.html -d Docs documentation
