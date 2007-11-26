bin_PROGRAMS = mirror
mirror_SOURCES = main.cpp
mirror_LDADD = \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la \
        $(top_srcdir)/external/cuba/libCuba.la

mirror_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/libraries
