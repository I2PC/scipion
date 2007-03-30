lib_LTLIBRARIES = libXmippInterface_@host_os@.la \
	libXmippInterfaceINRIA_@host_os@.la

libXmippInterface_@host_os@_la_SOURCES = Src/xmippAPH.cc \
	Src/xmippSpider.cc \
	Src/xmippVirus.cc \
	Src/xmippCrisp.cc \
	Src/xmippVrml.cc \
	Src/xmippOpenDX.cc \
	Src/xmippOpenDXang.cc \
	Src/xmippJDL.cc \
	Src/xmippAPH3D.cc \
	Src/xmippAPHorigmerg.cc \
	Src/xmippCCLattice_IO.cc \
	Src/xmippCCP4.cc \
	Src/xmippEM.cc

libXmippInterfaceINRIA_@host_os@_la_SOURCES = Src/xmippINRIA.cc

EXTRA_DIST = xmippAPH3D.hh        xmippCrisp.hh  xmippOpenDX.hh  xmippVrml.hh \
	xmippAPH.hh          xmippINRIA.hh  xmippSpider.hh \
	xmippAPHorigmerg.hh  xmippJDL.hh    xmippVirus.hh xmippOpenDXang.hh \
	xmippCCLattice_IO.hh xmippCCP4.hh xmippEM.hh

INCLUDES = -I$(top_srcdir)/Lib

docs:
	@echo "Creating documentation in XmippInterface ..."
	@doc++ -f -B @XMIPPDIR@/Extra_Docs/banner.html -d Docs documentation
