bin_PROGRAMS = angular_neighbourhood
angular_neighbourhood_SOURCES = projection_neighbourhood.cc
angular_neighbourhood_LDADD = \
	$(top_srcdir)/Lib/XmippInterface/libXmippInterface_@host_os@.la \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la
angular_neighbourhood_LDFLAGS = @EXEC_TYPE@	
INCLUDES = -I$(top_srcdir)/Lib
