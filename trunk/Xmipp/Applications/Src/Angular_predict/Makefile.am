bin_PROGRAMS = angular_discrete_assign
angular_discrete_assign_SOURCES = angular_predict.cc
angular_discrete_assign_LDADD = \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/Classification/libXmippClassif_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la
angular_discrete_assign_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/Lib
