bin_PROGRAMS = classify_fcmeans
classify_fcmeans_SOURCES = fcmeans.cc
classify_fcmeans_LDADD = \
	$(top_srcdir)/Lib/Classification/libXmippClassif_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la

classify_fcmeans_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/Lib
