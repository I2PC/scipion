bin_PROGRAMS = convert_vol2pseudo
convert_vol2pseudo_SOURCES = main.cpp
convert_vol2pseudo_LDADD = \
	$(top_srcdir)/libraries/classification/libXmippClassif.la \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la \
        $(top_srcdir)/external/cuba/libCuba.la

convert_vol2pseudo_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/libraries
