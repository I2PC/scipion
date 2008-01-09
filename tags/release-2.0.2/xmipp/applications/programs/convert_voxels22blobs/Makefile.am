bin_PROGRAMS = convert_voxels22blobs
convert_voxels22blobs_SOURCES = main.cpp
convert_voxels22blobs_LDADD = \
	$(top_srcdir)/libraries/interface/libXmippInterface.la \
	$(top_srcdir)/libraries/reconstruction/libXmippRecons.la \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la \
        $(top_srcdir)/external/cuba/libCuba.la

convert_voxels22blobs_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/libraries
