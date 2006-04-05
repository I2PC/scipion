bin_PROGRAMS = MLrefine3D
MLrefine3D_SOURCES = MLrefine3D.cc
MLrefine3D_LDADD = $(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la 

MLrefine3D_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/Lib

