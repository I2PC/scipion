bin_PROGRAMS = convert_data2img
convert_data2img_SOURCES = main.cpp
convert_data2img_LDADD = \
	$(top_srcdir)/libraries/classification/libXmippClassif.la \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la \
        $(top_srcdir)/external/cuba/libCuba.la

convert_data2img_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/libraries
