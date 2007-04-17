bin_PROGRAMS = psd_enhance
psd_enhance_SOURCES = main.cpp
psd_enhance_LDADD = \
	$(top_srcdir)/libraries/interface/libXmippInterface.la \
	$(top_srcdir)/libraries/reconstruction/libXmippRecons.la \
	$(top_srcdir)/libraries/data/libXmippData.la

psd_enhance_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/libraries
