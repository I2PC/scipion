bin_PROGRAMS = projection_matching
projection_matching_SOURCES = projection_matching.cc
projection_matching_LDADD = $(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la 

projection_matching_LDFLAGS = @EXEC_TYPE@	

INCLUDES = -I$(top_srcdir)/Lib
