bin_PROGRAMS = micrograph_window
micrograph_window_SOURCES = main.cpp
micrograph_window_LDADD = \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la

micrograph_window_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/libraries
