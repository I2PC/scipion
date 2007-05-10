bin_PROGRAMS = classify_sammon_projection
classify_sammon_projection_SOURCES = main.cpp
classify_sammon_projection_LDADD = \
	$(top_srcdir)/libraries/classification/libXmippClassif.la \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la \
        $(top_srcdir)/external/cuba/libCuba.la

classify_sammon_projection_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/libraries
