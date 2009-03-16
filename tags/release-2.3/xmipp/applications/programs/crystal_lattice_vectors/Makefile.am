bin_PROGRAMS = crystal_lattice_vectors
crystal_lattice_vectors_SOURCES = main.cpp
crystal_lattice_vectors_LDADD = \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la \
        $(top_srcdir)/external/cuba/libCuba.la

crystal_lattice_vectors_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/libraries
