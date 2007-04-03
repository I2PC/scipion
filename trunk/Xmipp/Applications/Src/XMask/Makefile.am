bin_PROGRAMS = design_mask
design_mask_SOURCES = xmask.cc \
	maskimg.cc
nodist_design_mask_SOURCES = maskimg_moc.cc

design_mask_LDADD = $(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la \
	$(top_srcdir)/Lib/Classification/libXmippClassif_@host_os@.la \
	@QT_LIBS@

design_mask_LDFLAGS =

INCLUDES = -I$(top_srcdir)/Lib @QT_CXXFLAGS@

maskimg_moc.cc :
	@QT_MOC@ maskimg.hh -o maskimg_moc.cc

MOSTLYCLEANFILES = *_moc.cc

EXTRA_DIST = maskimg.hh
