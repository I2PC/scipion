bin_PROGRAMS = convert_data2img
convert_data2img_SOURCES = data2img.cc
convert_data2img_LDADD = \
	$(top_srcdir)/Lib/Classification/libXmippClassif_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la
convert_data2img_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/Lib
