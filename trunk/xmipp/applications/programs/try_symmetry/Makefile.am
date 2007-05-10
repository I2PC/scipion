bin_PROGRAMS = try_symmetry
try_symmetry_SOURCES = main.cpp
try_symmetry_LDADD = \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la \
        $(top_srcdir)/external/cuba/libCuba.la

try_symmetry_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/libraries
