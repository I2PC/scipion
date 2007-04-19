bin_PROGRAMS = convert_spi22em
convert_spi22em_SOURCES = main.cpp
convert_spi22em_LDADD =  \
	$(top_srcdir)/libraries/reconstruction/libXmippRecons.la \
	$(top_srcdir)/libraries/reconstruction/libXmippRecons_Interface.la \
	$(top_srcdir)/libraries/interface/libXmippInterface.la \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la

INCLUDES = \
	-I$(top_srcdir)/libraries \
	-I$(top_srcdir)/external
