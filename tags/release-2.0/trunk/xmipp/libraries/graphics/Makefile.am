lib_LTLIBRARIES = libXmippGraphics.la

libXmippGraphics_la_SOURCES = \
	show_2d.cpp \
	show_ctf_estimate.cpp \
	show_plotter.cpp \
	show_selfile.cpp \
	show_som.cpp \
	show_spectra.cpp \
	show_spectra_som.cpp \
	show_table.cpp \
	show_tools.cpp \
	show_vol.cpp

nodist_libXmippGraphics_la_SOURCES = \
	show_2d_moc.cpp \
	show_ctf_estimate_moc.cpp \
	show_plotter_moc.cpp \
	show_selfile_moc.cpp \
	show_som_moc.cpp \
	show_spectra_moc.cpp \
	show_spectra_som_moc.cpp \
	show_table_moc.cpp \
	show_tools_moc.cpp \
	show_vol_moc.cpp

INCLUDES = -I$(top_srcdir)/libraries @QT_CXXFLAGS@
LIBADD = @QT_LIBS@

MOSTLYCLEANFILES = *_moc.cpp

show_2d_moc.cpp : show_2d.h
	@QT_MOC@ show_2d.h -o show_2d_moc.cpp

show_ctf_estimate_moc.cpp : show_ctf_estimate.h
	@QT_MOC@ show_ctf_estimate.h -o show_ctf_estimate_moc.cpp

show_plotter_moc.cpp : show_plotter.h
	@QT_MOC@ show_plotter.h -o show_plotter_moc.cpp

show_selfile_moc.cpp : show_selfile.h
	@QT_MOC@ show_selfile.h -o show_selfile_moc.cpp

show_som_moc.cpp : show_som.h
	@QT_MOC@ show_som.h -o show_som_moc.cpp

show_spectra_moc.cpp : show_spectra.h
	@QT_MOC@ show_spectra.h -o show_spectra_moc.cpp

show_spectra_som_moc.cpp : show_spectra_som.h
	@QT_MOC@ show_spectra_som.h -o show_spectra_som_moc.cpp

show_table_moc.cpp : show_table.h
	@QT_MOC@ show_table.h -o show_table_moc.cpp

show_tools_moc.cpp : show_tools.h
	@QT_MOC@ show_tools.h -o show_tools_moc.cpp

show_vol_moc.cpp : show_vol.h
	@QT_MOC@ show_vol.h -o show_vol_moc.cpp

EXTRA_DIST = \
	show_2d.h \
	show_ctf_estimate.h \
	show_plotter.h \
	show_selfile.h \
	show_som.h \
	show_spectra.h \
	show_spectra_som.h \
	show_table.h \
	show_tools.h \
	show_vol.h

## FIXME Coss...
## images/filesave20.png \
## images/zoomin20.png \
## images/zoomout20.png \
## images/filesave.png \
## images/zoomin.png
## images/zoomout.png
