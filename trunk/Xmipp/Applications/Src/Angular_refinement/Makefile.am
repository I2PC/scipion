bin_PROGRAMS = angular_refinement
angular_refinement_SOURCES = angular_refinement.cc
angular_refinement_LDADD = \
	$(top_srcdir)/Lib/XmippInterface/libXmippInterface_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la

angular_refinement_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/Lib

