bin_PROGRAMS = classify_som
classify_som_SOURCES = main.cpp
classify_som_LDADD = \
	$(top_srcdir)/libraries/classification/libXmippClassif.la \
	$(top_srcdir)/libraries/data/libXmippData.la

classify_som_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/libraries
