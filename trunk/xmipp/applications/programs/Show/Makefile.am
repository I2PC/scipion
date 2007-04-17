bin_PROGRAMS = show
show_SOURCES = show.cc
show_LDADD = $(top_srcdir)/Lib/Classification/libXmippClassif_@host_os@.la \
	$(top_srcdir)/Lib/XmippGraphics/libXmippGraphics_@host_os@.la \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la \
	@QT_LIBS@

show_LDFLAGS = 

INCLUDES = -I$(top_srcdir)/Lib @QT_CXXFLAGS@

MOSTLYCLEANFILES = *_moc.cc
