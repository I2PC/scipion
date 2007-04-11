bin_PROGRAMS = mask_design
mask_design_SOURCES = xmask.cc \
	maskimg.cc
nodist_mask_design_SOURCES = maskimg_moc.cc

mask_design_LDADD = \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la \
	$(top_srcdir)/Lib/Classification/libXmippClassif_@host_os@.la \
	@QT_LIBS@

mask_design_LDFLAGS =

INCLUDES = -I$(top_srcdir)/Lib @QT_CXXFLAGS@

maskimg_moc.cc :
	@QT_MOC@ maskimg.hh -o maskimg_moc.cc

MOSTLYCLEANFILES = *_moc.cc

EXTRA_DIST = maskimg.hh
