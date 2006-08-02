bin_PROGRAMS = pdbphantom
pdbphantom_SOURCES = pdbphantom.cc
pdbphantom_LDADD = $(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la

pdbphantom_LDFLAGS = @EXEC_TYPE@
	
INCLUDES = -I$(top_srcdir)/Lib
