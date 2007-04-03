bin_PROGRAMS = convert_data2vol
convert_data2vol_SOURCES = data2vol.cc
convert_data2vol_LDADD = \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la \
	$(top_srcdir)/Lib/Classification/libXmippClassif_@host_os@.la
convert_data2vol_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/Lib
