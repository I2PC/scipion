bin_PROGRAMS = random_mapping
random_mapping_SOURCES = rm.cc
random_mapping_LDADD = $(top_srcdir)/Lib/Classification/libXmippClassif_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la 

random_mapping_LDFLAGS = @EXEC_TYPE@
	
INCLUDES = -I$(top_srcdir)/Lib
