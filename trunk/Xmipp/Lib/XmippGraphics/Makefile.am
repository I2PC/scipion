lib_LTLIBRARIES = libXmippGraphics_@host_os@.la

libXmippGraphics_@host_os@_la_SOURCES = Src/showTools.cc      showTools_moc.cc \
	Src/show2D.cc         show2D_moc.cc \
	Src/showTable.cc      showTable_moc.cc \
	Src/showSel.cc        showSel_moc.cc \
	Src/showVol.cc        showVol_moc.cc \
	Src/showSpectra.cc    showSpectra_moc.cc \
	Src/showSOM.cc        showSOM_moc.cc \
	Src/showSpectraSOM.cc showSpectraSOM_moc.cc \
	Src/showAssignCTF.cc  showAssignCTF_moc.cc

INCLUDES = -I$(top_srcdir)/Lib @QT_CXXFLAGS@
LIBADD = @QT_LIBS@

MOSTLYCLEANFILES = *_moc.cc

showTools_moc.cc : showTools.hh
	@QT_MOC@ showTools.hh -o showTools_moc.cc
show2D_moc.cc : show2D.hh
	@QT_MOC@ show2D.hh -o show2D_moc.cc
showTable_moc.cc : showTable.hh
	@QT_MOC@ showTable.hh -o showTable_moc.cc
showSel_moc.cc : showSel.hh
	@QT_MOC@ showSel.hh -o showSel_moc.cc
showVol_moc.cc : showVol.hh
	@QT_MOC@ showVol.hh -o showVol_moc.cc
showSpectra_moc.cc : showSpectra.hh
	@QT_MOC@ showSpectra.hh -o showSpectra_moc.cc
showSOM_moc.cc : showSOM.hh
	@QT_MOC@ showSOM.hh -o showSOM_moc.cc
showSpectraSOM_moc.cc : showSpectraSOM.hh
	@QT_MOC@ showSpectraSOM.hh -o showSpectraSOM_moc.cc
showAssignCTF_moc.cc : showAssignCTF.hh
	@QT_MOC@ showAssignCTF.hh -o showAssignCTF_moc.cc

EXTRA_DIST = show2D.hh   showSOM.hh      showSpectraSOM.hh  showTools.hh \
	showSel.hh  showSpectra.hh  showTable.hh       showVol.hh \
        showAssignCTF.hh

docs:
	@echo "Creating documentation in XmippInterface ..."
	@doc++ -f -B @XMIPPDIR@/Extra_Docs/banner.html -d Docs documentation
