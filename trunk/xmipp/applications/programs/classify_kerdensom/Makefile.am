bin_PROGRAMS = classify_kerdensom
classify_kerdensom_SOURCES = main.cpp
classify_kerdensom_LDADD = \
	$(top_srcdir)/libraries/classification/libXmippClassif.la \
	$(top_srcdir)/libraries/data/libXmippData.la

classify_kerdensom_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/libraries
