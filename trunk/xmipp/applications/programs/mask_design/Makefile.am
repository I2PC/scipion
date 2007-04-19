bin_PROGRAMS = mask_design
mask_design_SOURCES = main.cpp maskimg.cpp
nodist_mask_design_SOURCES = maskimg_moc.cpp

mask_design_LDADD = \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/libraries/classification/libXmippClassif.la \
	$(top_srcdir)/external/bilib/libBilib.la \
	@QT_LIBS@

mask_design_LDFLAGS =

INCLUDES = -I$(top_srcdir)/libraries @QT_CXXFLAGS@

maskimg_moc.cpp :
	@QT_MOC@ maskimg.h -o maskimg_moc.cpp

MOSTLYCLEANFILES = *_moc.cpp

EXTRA_DIST = maskimg.h
