bin_PROGRAMS = xmask
xmask_SOURCES = xmask.cc \
	maskimg.cc \
	maskimg_moc.cc

xmask_LDADD = $(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la \
	$(top_srcdir)/Lib/Classification/libXmippClassif_@host_os@.la \
	@QT_LIBS@

INCLUDES = -I$(top_srcdir)/Lib @QT_CXXFLAGS@

maskimg_moc.cc :
	@QT_MOC@ maskimg.hh -o maskimg_moc.cc
