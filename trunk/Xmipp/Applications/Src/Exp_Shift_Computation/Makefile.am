bin_PROGRAMS = exp_shift_computation
exp_shift_computation_SOURCES = exp_shift_computation.cc
exp_shift_computation_LDADD = $(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_Interface_@host_os@.la \
	$(top_srcdir)/Lib/XmippInterface/libXmippInterface_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la 

exp_shift_computation_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/Lib
EXTRA_DIST = exp_shift_computation.hh

