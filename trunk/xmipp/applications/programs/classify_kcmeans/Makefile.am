bin_PROGRAMS = classify_kcmeans
classify_kcmeans_SOURCES = main.cpp
classify_kcmeans_LDADD = \
	$(top_srcdir)/libraries/classification/libXmippClassif.la \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la

classify_kcmeans_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/libraries
