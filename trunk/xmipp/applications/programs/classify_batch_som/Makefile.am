bin_PROGRAMS = classify_batch_som
classify_batch_som_SOURCES = batchsom.cc
classify_batch_som_LDADD = \
	$(top_srcdir)/Lib/Classification/libXmippClassif_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la

classify_batch_som_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/Lib
