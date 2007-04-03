bin_PROGRAMS = angular_projection_matching
angular_projection_matching_SOURCES = projection_matching.cc
angular_projection_matching_LDADD = \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la
angular_projection_matching_LDFLAGS = @EXEC_TYPE@	
INCLUDES = -I$(top_srcdir)/Lib
