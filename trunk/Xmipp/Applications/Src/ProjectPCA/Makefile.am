bin_PROGRAMS = classify_project_pca
classify_project_pca_SOURCES = projectpca.cc
classify_project_pca_LDADD = \
	$(top_srcdir)/Lib/Classification/libXmippClassif_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la
classify_project_pca_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/Lib
