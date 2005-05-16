bin_PROGRAMS = wbp
wbp_SOURCES = wbp.cc
wbp_LDADD = $(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la 

wbp_LDFLAGS = @EXEC_TYPE@	

INCLUDES = -I$(top_srcdir)/Lib
