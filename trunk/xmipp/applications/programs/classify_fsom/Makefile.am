bin_PROGRAMS = classify_fsom
classify_fsom_SOURCES = fsom.cc
classify_fsom_LDADD = \
	$(top_srcdir)/Lib/Classification/libXmippClassif_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la
classify_fsom_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/Lib
