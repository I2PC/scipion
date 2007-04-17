bin_PROGRAMS = convert_data2pdb
convert_data2pdb_SOURCES = toPDB.cc
convert_data2pdb_LDADD = \
	$(top_srcdir)/Lib/Classification/libXmippClassif_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la
convert_data2pdb_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/Lib
