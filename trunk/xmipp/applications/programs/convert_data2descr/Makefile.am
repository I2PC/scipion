bin_PROGRAMS = convert_data2descr
convert_data2descr_SOURCES = main.cpp
convert_data2descr_LDADD = \
	$(top_srcdir)/libraries/classification/libXmippClassif.la \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la \
        $(top_srcdir)/external/cuba/libCuba.la

convert_data2descr_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/libraries
