bin_PROGRAMS = ctfviewer
ctfviewer_SOURCES = ctfviewer.cc ctfviewer_module.cc ctfviewer_module_moc.cc

ctfviewer_LDADD = $(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la \
                  $(top_srcdir)/Lib/XmippGraphics/libXmippGraphics_@host_os@.la \
                  $(top_srcdir)/Lib/Classification/libXmippClassif_@host_os@.la \
	          $(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la

ctfviewer_LDFLAGS = @QT_LIBS@ 

INCLUDES = -I$(top_srcdir)/Lib @QT_CXXFLAGS@

MOSTLYCLEANFILES = *_moc.cc

ctfviewer_module_moc.cc : 
	@QT_MOC@ ctfviewer_module.hh -o ctfviewer_module_moc.cc

EXTRA_DIST = ctfviewer_module.hh
