bin_PROGRAMS = ml_align2d
ml_align2d_SOURCES = MLalign2D.cc
ml_align2d_LDADD = \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la
ml_align2d_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/Lib
