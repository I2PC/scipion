bin_PROGRAMS = resolution_ssnr
resolution_ssnr_SOURCES = SSNR.cc
resolution_ssnr_LDADD = \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la
resolution_ssnr_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/Lib
