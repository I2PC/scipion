bin_PROGRAMS = make_spectra
make_spectra_SOURCES = main.cpp
make_spectra_LDADD = \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la

make_spectra_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/libraries
