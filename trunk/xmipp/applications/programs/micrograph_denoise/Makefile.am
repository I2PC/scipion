bin_PROGRAMS = micrograph_denoise
micrograph_denoise_SOURCES = main.cpp
micrograph_denoise_LDADD = \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la

micrograph_denoise_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/libraries
