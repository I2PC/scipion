bin_PROGRAMS = convert_vol2data
convert_vol2data_SOURCES = main.cpp
convert_vol2data_LDADD = \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la \
        $(top_srcdir)/external/cuba/libCuba.la

convert_vol2data_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/libraries
