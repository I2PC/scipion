bin_PROGRAMS = idr
idr_SOURCES = idr_art.cc
idr_LDADD = \
	$(top_srcdir)/Lib/XmippInterface/libXmippInterface_@host_os@.la \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la
idr_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/Lib
