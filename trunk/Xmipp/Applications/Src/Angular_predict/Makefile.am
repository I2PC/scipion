bin_PROGRAMS = angular_predict
angular_predict_SOURCES = angular_predict.cc
angular_predict_LDADD = $(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_Classif_@host_os@.la \
	$(top_srcdir)/Lib/Classification/libXmippClassif_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la 
	
INCLUDES = -I$(top_srcdir)/Lib
