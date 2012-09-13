bin_PROGRAMS = classify_pca
classify_pca_SOURCES = main.cpp
classify_pca_LDADD = \
	$(top_srcdir)/libraries/classification/libXmippClassif.la \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la \
        $(top_srcdir)/external/cuba/libCuba.la

classify_pca_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/libraries

