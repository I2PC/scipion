bin_PROGRAMS = convert_spi22em
convert_spi22em_SOURCES = spi22em.cc
conevrt_spi22em_LDADD =  \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_Interface_@host_os@.la \
	$(top_srcdir)/Lib/XmippInterface/libXmippInterface_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la

INCLUDES = -I$(top_srcdir)/Lib
