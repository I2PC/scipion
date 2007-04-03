bin_PROGRAMS = angular_continuous_assign
angular_continuous_assign_SOURCES = angular_predict_continuous.cc
angular_continuous_assign_LDADD = \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/Classification/libXmippClassif_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la
angular_continuous_assign_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/Lib
