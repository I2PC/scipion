bin_PROGRAMS = ssnr
ssnr_SOURCES = SSNR.cc
ssnr_LDADD = $(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la 
	
	
ssnr_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/Lib
