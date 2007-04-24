bin_PROGRAMS = statistics
statistics_SOURCES = main.cpp
statistics_LDADD = \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la

statistics_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/libraries
