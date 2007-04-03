bin_PROGRAMS = classify_kerdensom
classify_kerdensom_SOURCES = kerdensom.cc
classify_kerdensom_LDADD = \
	$(top_srcdir)/Lib/Classification/libXmippClassif_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la
classify_kerdensom_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/Lib
