bin_PROGRAMS = classify_som
classify_som_SOURCES = som.cc
classify_som_LDADD = \
	$(top_srcdir)/Lib/Classification/libXmippClassif_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la
classify_som_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/Lib
