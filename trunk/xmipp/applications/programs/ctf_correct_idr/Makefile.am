bin_PROGRAMS = ctf_correct_idr
ctf_correct_idr_SOURCES = idr_art.cc
ctf_correct_idr_LDADD = \
	$(top_srcdir)/Lib/XmippInterface/libXmippInterface_@host_os@.la \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la

ctf_correct_idr_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/Lib
