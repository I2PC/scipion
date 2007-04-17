bin_PROGRAMS = crystal_unbend
crystal_unbend_SOURCES = unbend.cc
crystal_unbend_LDADD = \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_Interface_@host_os@.la \
	$(top_srcdir)/Lib/XmippInterface/libXmippInterface_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la

crystal_unbend_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/Lib
