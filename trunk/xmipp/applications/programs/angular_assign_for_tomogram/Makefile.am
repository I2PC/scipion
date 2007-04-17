bin_PROGRAMS = angular_assign_from_tomogram
angular_assign_from_tomogram_SOURCES = main.cpp
angular_assign_from_tomogram_LDADD = \
	$(top_srcdir)/libraries/reconstruction/libXmippRecons.la \
	$(top_srcdir)/libraries/classification/libXmippClassif.la \
	$(top_srcdir)/libraries/data/libXmippData.la

angular_assign_from_tomogram_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/libraries
