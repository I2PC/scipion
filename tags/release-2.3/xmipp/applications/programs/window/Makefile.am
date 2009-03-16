bin_PROGRAMS = window
window_SOURCES = main.cpp
window_LDADD = \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la \
        $(top_srcdir)/external/cuba/libCuba.la

window_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/libraries
