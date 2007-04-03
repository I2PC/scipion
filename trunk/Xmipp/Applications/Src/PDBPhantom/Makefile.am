bin_PROGRAMS = convert_pdb2vol
convert_pdb2vol_SOURCES = pdbphantom.cc
convert_pdb2vol_LDADD = \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la
convert_pdb2vol_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/Lib
