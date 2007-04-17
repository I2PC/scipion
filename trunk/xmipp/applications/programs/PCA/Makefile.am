bin_PROGRAMS = classify_pca
classify_pca_SOURCES = pca.cc
classify_pca_LDADD = \
	$(top_srcdir)/Lib/Classification/libXmippClassif_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la
classify_pca_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/Lib

