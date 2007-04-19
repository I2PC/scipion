bin_PROGRAMS = angular_distribution_show
angular_distribution_show_SOURCES = main.cpp
angular_distribution_show_LDADD = \
	$(top_srcdir)/libraries/interface/libXmippInterface.la \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la

angular_distribution_show_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/libraries
