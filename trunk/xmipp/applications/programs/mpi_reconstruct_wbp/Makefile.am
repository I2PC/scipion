bin_PROGRAMS = mpi_reconstruct_wbp
mpi_reconstruct_wbp_SOURCES = main.cpp
mpi_reconstruct_wbp_LDADD = \
	$(top_srcdir)/libraries/reconstruction/libXmippRecons.la \
	$(top_srcdir)/libraries/data/libXmippData.la \
	@MPILIBS@

CC=@MPICXX@
CXX=@MPICXX@
CXXLINK = $(LIBTOOL) --tag=CXX --mode=link $(CXXLD) $(AM_CXXFLAGS) $(CXXFLAGS) $(AM_LDFLAGS) $(LDFLAGS) -o $@

mpi_reconstruct_wbp_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/libraries
