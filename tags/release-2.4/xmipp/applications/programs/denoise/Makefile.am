bin_PROGRAMS = denoise
denoise_SOURCES = main.cpp
denoise_LDADD = \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la \
        $(top_srcdir)/external/cuba/libCuba.la

denoise_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/libraries
