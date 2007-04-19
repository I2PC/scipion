bin_PROGRAMS = angular_neighbourhood
angular_neighbourhood_SOURCES = main.cpp
angular_neighbourhood_LDADD = \
	$(top_srcdir)/libraries/interface/libXmippInterface.la \
	$(top_srcdir)/libraries/reconstruction/libXmippRecons.la \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la

angular_neighbourhood_LDFLAGS = @EXEC_TYPE@	

INCLUDES = -I$(top_srcdir)/libraries
