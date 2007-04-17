bin_PROGRAMS = convert_vol2pseudo
convert_vol2pseudo_SOURCES = volsampler.cc
convert_vol2pseudo_LDADD = \
	$(top_srcdir)/Lib/Classification/libXmippClassif_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la

convert_vol2pseudo_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/Lib
