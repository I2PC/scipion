bin_PROGRAMS =art
art_SOURCES = art.cc
art_LDADD = \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la

art_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/Lib
