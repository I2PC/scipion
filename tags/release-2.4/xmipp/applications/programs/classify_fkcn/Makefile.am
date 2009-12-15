bin_PROGRAMS = classify_fkcn
classify_fkcn_SOURCES = main.cpp
classify_fkcn_LDADD = \
	$(top_srcdir)/libraries/classification/libXmippClassif.la \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la \
        $(top_srcdir)/external/cuba/libCuba.la

classify_fkcn_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/libraries
