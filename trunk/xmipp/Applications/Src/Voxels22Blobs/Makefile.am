bin_PROGRAMS = convert_voxels22blobs
convert_voxels22blobs_SOURCES = voxels22blobs.cc
convert_voxels22blobs_LDADD = \
	$(top_srcdir)/Lib/XmippInterface/libXmippInterface_@host_os@.la \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la

convert_voxels22blobs_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/Lib
