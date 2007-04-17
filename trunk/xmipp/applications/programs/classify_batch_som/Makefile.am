bin_PROGRAMS = classify_batch_som
classify_batch_som_SOURCES = main.cpp
classify_batch_som_LDADD = \
	$(top_srcdir)/libraries/classification/libXmippClassif.la \
	$(top_srcdir)/libraries/data/libXmippData.la

classify_batch_som_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/libraries
