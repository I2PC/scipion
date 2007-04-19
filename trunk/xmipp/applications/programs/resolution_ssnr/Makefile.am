bin_PROGRAMS = resolution_ssnr
resolution_ssnr_SOURCES = main.cpp
resolution_ssnr_LDADD = \
	$(top_srcdir)/libraries/reconstruction/libXmippRecons.la \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la

resolution_ssnr_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/libraries
