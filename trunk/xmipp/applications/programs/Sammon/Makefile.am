bin_PROGRAMS = classify_sammon_projection
classify_sammon_projection_SOURCES = sammon.cc
classify_sammon_projection_LDADD = \
	$(top_srcdir)/Lib/Classification/libXmippClassif_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la
classify_sammon_projection_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/Lib
