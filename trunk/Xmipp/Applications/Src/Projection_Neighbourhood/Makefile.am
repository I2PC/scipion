bin_PROGRAMS = projection_neighbourhood
projection_neighbourhood_SOURCES = projection_neighbourhood.cc
projection_neighbourhood_LDADD = $(top_srcdir)/Lib/XmippInterface/libXmippInterface_@host_os@.la \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la 

projection_neighbourhood_LDFLAGS = @EXEC_TYPE@	

INCLUDES = -I$(top_srcdir)/Lib
