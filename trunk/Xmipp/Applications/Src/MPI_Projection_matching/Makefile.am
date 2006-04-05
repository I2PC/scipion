bin_PROGRAMS = mpi_projection_matching
mpi_projection_matching_SOURCES = mpi_projection_matching.cc
mpi_projection_matching_LDADD = $(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la \
	@MPILIBS@

CC=@MPICXX@
CXX=@MPICXX@
CXXLINK = $(LIBTOOL) --tag=CXX --mode=link $(CXXLD) $(AM_CXXFLAGS) $(CXXFLAGS) $(AM_LDFLAGS) $(LDFLAGS) -o $@

mpi_projection_matching_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/Lib
