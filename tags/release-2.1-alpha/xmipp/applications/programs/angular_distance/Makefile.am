bin_PROGRAMS = angular_distance
angular_distance_SOURCES = main.cpp
angular_distance_LDADD = \
	$(top_srcdir)/libraries/reconstruction/libXmippRecons.la \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la \
        $(top_srcdir)/external/cuba/libCuba.la

angular_distance_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/libraries

