bin_PROGRAMS = evaluate_FSCs
evaluate_FSCs_SOURCES = evaluate_FSCs.cc
evaluate_FSCs_LDADD = $(top_srcdir)/Lib/XmippInterface/libXmippInterface_@host_os@.la \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la

evaluate_FSCs_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/Lib

