bin_PROGRAMS = align2D
align2D_SOURCES = align2d.cc
align2D_LDADD = $(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la 
	

align2D_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/Lib
