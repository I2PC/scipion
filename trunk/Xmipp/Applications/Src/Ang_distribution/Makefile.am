bin_PROGRAMS = angular_distribution_show
angular_distribution_show_SOURCES = ang_distribution.cc
angular_distribution_show_LDADD = \
	$(top_srcdir)/Lib/XmippInterface/libXmippInterface_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la

angular_distribution_show_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/Lib
