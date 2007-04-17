bin_PROGRAMS = projection_matching_crystal
projection_matching_crystal_SOURCES = projection_matching_crystal.cc
projection_matching_crystal_LDADD = $(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la 

projection_matching_crystal_LDFLAGS = @EXEC_TYPE@	

INCLUDES = -I$(top_srcdir)/Lib
