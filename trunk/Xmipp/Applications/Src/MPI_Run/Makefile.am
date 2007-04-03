bin_PROGRAMS = mpi_run
mpi_run_SOURCES = mpi_run.cc
mpi_run_LDADD = $(top_srcdir)/Lib/ReconstructionMPI/libXmippReconsMPI_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la \
	@MPILIBS@

CC=@MPICXX@
CXX=@MPICXX@
CXXLINK = $(LIBTOOL) --tag=CXX --mode=link $(CXXLD) $(AM_CXXFLAGS) $(CXXFLAGS) $(AM_LDFLAGS) $(LDFLAGS) -o $@

mpi_run_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/Lib
