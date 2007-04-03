bin_PROGRAMS = angular_assign_from_tomogram
angular_assign_from_tomogram_SOURCES = angular_predict_tomography.cc
angular_assign_from_tomogram_LDADD = \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/Classification/libXmippClassif_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la
angular_assign_from_tomogram_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/Lib
