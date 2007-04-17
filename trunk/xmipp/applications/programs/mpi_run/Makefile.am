bin_PROGRAMS = mpi_run
mpi_run_SOURCES = main.cpp
mpi_run_LDADD = \
	$(top_srcdir)/libraries/reconstructionMPI/libXmippReconsMPI.la \
	$(top_srcdir)/libraries/data/libXmippData.la \
	@MPILIBS@

CC=@MPICXX@
CXX=@MPICXX@
CXXLINK = $(LIBTOOL) --tag=CXX --mode=link $(CXXLD) $(AM_CXXFLAGS) $(CXXFLAGS) $(AM_LDFLAGS) $(LDFLAGS) -o $@

mpi_run_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/libraries
