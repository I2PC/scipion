bin_PROGRAMS = batchsom
batchsom_SOURCES = batchsom.cc
batchsom_LDADD = $(top_srcdir)/Lib/Classification/libXmippClassif_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la
batchsom_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/Lib

