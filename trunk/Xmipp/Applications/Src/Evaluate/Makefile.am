bin_PROGRAMS = foms_evaluate
foms_evaluate_SOURCES = evaluate.cc
foms_evaluate_LDADD = \
	$(top_srcdir)/Lib/XmippInterface/libXmippInterface_@host_os@.la \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la

foms_evaluate_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/Lib
