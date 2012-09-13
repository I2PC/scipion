bin_PROGRAMS = find_center3d
find_center3d_SOURCES = main.cpp
find_center3d_LDADD = \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la \
        $(top_srcdir)/external/cuba/libCuba.la

find_center3d_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/libraries
