bin_PROGRAMS = angular_discrete_assign
angular_discrete_assign_SOURCES = main.cpp
angular_discrete_assign_LDADD = \
	$(top_srcdir)/libraries/reconstruction/libXmippRecons.la \
	$(top_srcdir)/libraries/classification/libXmippClassif.la \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la

angular_discrete_assign_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/libraries
