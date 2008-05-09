bin_PROGRAMS = angular_projection_matching
angular_projection_matching_SOURCES = main.cpp
angular_projection_matching_LDADD = \
	$(top_srcdir)/libraries/reconstruction/libXmippRecons.la \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la \
        $(top_srcdir)/external/cuba/libCuba.la

angular_projection_matching_LDFLAGS = @EXEC_TYPE@

INCLUDES = \
	-I$(top_srcdir)/libraries \
	-I$(top_srcdir)/external
