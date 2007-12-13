bin_PROGRAMS = add_noise
add_noise_SOURCES = main.cpp
add_noise_LDADD = \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la \
        $(top_srcdir)/external/cuba/libCuba.la

add_noise_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/libraries

