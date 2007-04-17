bin_PROGRAMS = crystal_skew
crystal_skew_SOURCES = skew.cc
crystal_skew_LDADD = \
	$(top_srcdir)/Lib/XmippInterface/libXmippInterface_@host_os@.la \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la

crystal_skew_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/Lib
