bin_PROGRAMS = spi22ccp4
spi22ccp4_SOURCES = main.cpp
spi22ccp4_LDADD =  \
	$(top_srcdir)/libraries/reconstruction/libXmippRecons.la \
	$(top_srcdir)/libraries/reconstruction/libXmippRecons_Interface.la \
	$(top_srcdir)/libraries/interface/libXmippInterface.la \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la

INCLUDES = \
	-I$(top_srcdir)/libraries \
	-I$(top_srcdir)/external
