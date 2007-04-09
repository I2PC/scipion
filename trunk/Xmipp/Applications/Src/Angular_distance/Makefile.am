bin_PROGRAMS = angular_distance
angular_distance_SOURCES = angular_distance.cc
angular_distance_LDADD = \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la

angular_distance_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/Lib

