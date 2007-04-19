bin_PROGRAMS = break_symmetry
break_symmetry_SOURCES = main.cpp
break_symmetry_LDADD = \
	$(top_srcdir)/libraries/reconstruction/libXmippRecons.la \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la

break_symmetry_LDFLAGS = @EXEC_TYPE@
INCLUDES = \
	-I$(top_srcdir)/libraries \
	-I$(top_srcdir)/external
