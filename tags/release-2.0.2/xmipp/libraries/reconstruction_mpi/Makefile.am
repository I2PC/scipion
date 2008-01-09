lib_LTLIBRARIES = libXmippReconsMPI.la

libXmippReconsMPI_la_SOURCES = mpi_run.cpp

CC=@MPICXX@
CXX=@MPICXX@
CXXLINK = $(LIBTOOL) --tag=CXX --mode=link $(CXXLD) $(AM_CXXFLAGS) $(CXXFLAGS) $(AM_LDFLAGS) $(LDFLAGS) -o $@

## mpi_run_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/libraries

EXTRA_DIST = mpi_run.h
