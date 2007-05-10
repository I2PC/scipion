bin_PROGRAMS = ctf_view
ctf_view_SOURCES = main.cpp module.cpp
nodist_ctf_view_SOURCES = module_moc.cpp

ctf_view_LDADD = \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/libraries/graphics/libXmippGraphics.la \
	$(top_srcdir)/libraries/classification/libXmippClassif.la \
	$(top_srcdir)/libraries/reconstruction/libXmippRecons.la \
	$(top_srcdir)/external/bilib/libBilib.la \
        $(top_srcdir)/external/cuba/libCuba.la

ctf_view_LDFLAGS = @QT_LIBS@
INCLUDES = -I$(top_srcdir)/libraries @QT_CXXFLAGS@

MOSTLYCLEANFILES = *_moc.cpp

module_moc.cpp :
	@QT_MOC@ module.h -o module_moc.cpp

EXTRA_DIST = module.h
