bin_PROGRAMS = idr_art
idr_art_SOURCES = idr_art.cc
idr_art_LDADD = $(top_srcdir)/Lib/XmippInterface/libXmippInterface_@host_os@.la \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la
	

INCLUDES = -I$(top_srcdir)/Lib
