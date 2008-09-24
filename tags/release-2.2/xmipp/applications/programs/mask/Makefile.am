bin_PROGRAMS = mask
mask_SOURCES = main.cpp
mask_LDADD = \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la \
        $(top_srcdir)/external/cuba/libCuba.la

mask_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/libraries

