bin_PROGRAMS = convert_pdb2vol
convert_pdb2vol_SOURCES = main.cpp
convert_pdb2vol_LDADD = \
	$(top_srcdir)/libraries/reconstruction/libXmippRecons.la \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la \
        $(top_srcdir)/external/cuba/libCuba.la

convert_pdb2vol_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/libraries
