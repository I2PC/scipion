bin_PROGRAMS = angular_predict_tomography
angular_predict_tomography_SOURCES = angular_predict_tomography.cc
angular_predict_tomography_LDADD = $(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/Classification/libXmippClassif_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la 

angular_predict_tomography_LDFLAGS = @EXEC_TYPE@
	
INCLUDES = -I$(top_srcdir)/Lib
