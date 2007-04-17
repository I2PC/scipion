bin_PROGRAMS = phantom_create
phantom_create_SOURCES = main.cpp
phantom_create_LDADD = \
	$(top_srcdir)/libraries/reconstruction/libXmippRecons.la \
	$(top_srcdir)/libraries/data/libXmippData.la

phantom_create_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/libraries
