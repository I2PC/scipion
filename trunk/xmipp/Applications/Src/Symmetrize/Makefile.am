bin_PROGRAMS =symmetrize
symmetrize_SOURCES = symmetrize.cc
symmetrize_LDADD = $(top_srcdir)/Lib/XmippInterface/libXmippInterface_@host_os@.la \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la 

symmetrize_LDFLAGS = @EXEC_TYPE@
	
INCLUDES = -I$(top_srcdir)/Lib
