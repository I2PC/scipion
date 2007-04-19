bin_PROGRAMS = classify_fsom
classify_fsom_SOURCES = main.cpp
classify_fsom_LDADD = \
	$(top_srcdir)/libraries/classification/libXmippClassif.la \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la

classify_fsom_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/libraries
