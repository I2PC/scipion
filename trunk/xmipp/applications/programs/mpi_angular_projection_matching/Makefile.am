bin_PROGRAMS = mpi_angular_projection_matching
mpi_angular_projection_matching_SOURCES = main.cpp
mpi_angular_projection_matching_LDADD = \
	$(top_srcdir)/libraries/reconstruction/libXmippRecons.la \
	$(top_srcdir)/libraries/data/libXmippData.la \
	@MPILIBS@

CC=@MPICXX@
CXX=@MPICXX@
CXXLINK = $(LIBTOOL) --tag=CXX --mode=link $(CXXLD) $(AM_CXXFLAGS) $(CXXFLAGS) $(AM_LDFLAGS) $(LDFLAGS) -o $@

mpi_angular_projection_matching_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/libraries
