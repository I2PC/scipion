bin_PROGRAMS = phantom_transform
phantom_transform_SOURCES = main.cpp
phantom_transform_LDADD = \
	$(top_srcdir)/libraries/reconstruction/libXmippRecons.la \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la

phantom_transform_LDFLAGS = @EXEC_TYPE@

INCLUDES = \
	-I$(top_srcdir)/libraries \
	-I$(top_srcdir)/external
