bin_PROGRAMS = show
show_SOURCES = main.cpp
show_LDADD = \
	$(top_srcdir)/libraries/classification/libXmippClassif.la \
	$(top_srcdir)/libraries/graphics/libXmippGraphics.la \
	$(top_srcdir)/libraries/reconstruction/libXmippRecons.la \
	$(top_srcdir)/libraries/data/libXmippData.la \
	@QT_LIBS@

show_LDFLAGS =

INCLUDES = -I$(top_srcdir)/libraries @QT_CXXFLAGS@

MOSTLYCLEANFILES = *_moc.cpp
