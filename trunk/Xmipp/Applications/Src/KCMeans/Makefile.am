bin_PROGRAMS = kcmeans
kcmeans_SOURCES = kcmeans.cc
kcmeans_LDADD = $(top_srcdir)/Lib/Classification/libXmippClassif_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la 

kcmeans_LDFLAGS = @EXEC_TYPE@
	
INCLUDES = -I$(top_srcdir)/Lib
