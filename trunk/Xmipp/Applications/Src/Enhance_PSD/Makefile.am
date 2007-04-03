bin_PROGRAMS = psd_enhance
psd_enhance_SOURCES = enhance_psd.cc
psd_enhance_LDADD = \
	$(top_srcdir)/Lib/XmippInterface/libXmippInterface_@host_os@.la \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la
psd_enhance_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/Lib
