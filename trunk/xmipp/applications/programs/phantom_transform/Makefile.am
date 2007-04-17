bin_PROGRAMS = phantom_transform
phantom_transform_SOURCES = main.cpp
phantom_transform_LDADD = \
	$(top_srcdir)/libraries/reconstruction/libXmippRecons.la \
	$(top_srcdir)/libraries/data/libXmippData.la

phantom_transform_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/libraries
