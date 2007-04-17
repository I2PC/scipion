bin_PROGRAMS = classify_fcmeans
classify_fcmeans_SOURCES = main.cpp
classify_fcmeans_LDADD = \
	$(top_srcdir)/libraries/classification/libXmippClassif.la \
	$(top_srcdir)/libraries/data/libXmippData.la

classify_fcmeans_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/libraries
