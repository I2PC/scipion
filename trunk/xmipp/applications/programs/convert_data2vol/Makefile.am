bin_PROGRAMS = convert_data2vol
convert_data2vol_SOURCES = main.cpp
convert_data2vol_LDADD = \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/libraries/classification/libXmippClassif.la

convert_data2vol_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/libraries
