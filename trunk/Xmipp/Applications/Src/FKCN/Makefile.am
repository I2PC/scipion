bin_PROGRAMS = classify_fkcn
classify_fkcn_SOURCES = fkcn.cc
classify_fkcn_LDADD = \
	$(top_srcdir)/Lib/Classification/libXmippClassif_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la
classify_fkcn_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/Lib
