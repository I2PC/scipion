lib_LTLIBRARIES = libXmippInterface_@host_os@.la \
	libXmippInterfaceINRIA_@host_os@.la \
	libXmippInterfaceVTK_@host_os@.la

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

libXmippInterfaceVTK_@host_os@_la_SOURCES   = Src/xmippVTK.cc

EXTRA_DIST = xmippAPH3D.hh        xmippCrisp.hh  xmippOpenDX.hh  xmippVrml.hh \
	xmippAPH.hh          xmippINRIA.hh  xmippSpider.hh  xmippVTK.hh \
	xmippAPHorigmerg.hh  xmippJDL.hh    xmippVirus.hh xmippOpenDXang.hh \
	xmippCCLattice_IO.hh xmippCCP4.hh xmippEM.hh


#libXmippInterfaceVTK_la_LIBADD = /usr/local/lib/vtk/libvtkCommon.a \
#	/usr/local/lib/vtk/libvtkHybrid.a    /usr/local/lib/vtk/libvtkpng.a \
#	/usr/local/lib/vtk/libvtkexpat.a      /usr/local/lib/vtk/libvtkImaging.a \
#	/usr/local/lib/vtk/libvtkRendering.a /usr/local/lib/vtk/libvtkFiltering.a \
#	/usr/local/lib/vtk/libvtkIO.a        /usr/local/lib/vtk/libvtktiff.a \
#	/usr/local/lib/vtk/libvtkfreetype.a   /usr/local/lib/vtk/libvtkjpeg.a  \
#	/usr/local/lib/vtk/libvtkzlib.a /usr/local/lib/vtk/libvtkftgl.a  \
#	/usr/local/lib/vtk/libvtkParallel.a /usr/local/lib/vtk/libvtkGraphics.a   /usr/local/lib/vtk/libvtkPatented.a


INCLUDES = -I$(top_srcdir)/Lib

docs:
	@echo "Creating documentation in XmippInterface ..."
	@doc++ -f -B @XMIPPDIR@/Extra_Docs/banner.html -d Docs documentation
