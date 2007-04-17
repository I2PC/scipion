bin_PROGRAMS = convert_data2descr
convert_data2descr_SOURCES = toPhantom.cc
convert_data2descr_LDADD = \
	$(top_srcdir)/Lib/Classification/libXmippClassif_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la
convert_data2descr_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/Lib
