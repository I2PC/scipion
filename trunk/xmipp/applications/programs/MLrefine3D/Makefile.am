bin_PROGRAMS = ml_refine3d
ml_refine3d_SOURCES = MLrefine3D.cc
ml_refine3d_LDADD = $(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la
ml_refine3d_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/Lib
