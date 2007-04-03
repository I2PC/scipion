bin_PROGRAMS = classify_kcmeans
classify_kcmeans_SOURCES = kcmeans.cc
classify_kcmeans_LDADD = \
	$(top_srcdir)/Lib/Classification/libXmippClassif_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la
classify_kcmeans_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/Lib
