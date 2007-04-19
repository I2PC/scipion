bin_PROGRAMS = fourier_transform
fourier_transform_SOURCES = main.cpp
fourier_transform_LDADD = \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la

fourier_transform_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/libraries
