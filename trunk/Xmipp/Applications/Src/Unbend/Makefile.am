bin_PROGRAMS = unbend
unbend_SOURCES = unbend.cc
unbend_LDADD = $(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_Interface_@host_os@.la \
	$(top_srcdir)/Lib/XmippInterface/libXmippInterface_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la 

unbend_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/Lib
