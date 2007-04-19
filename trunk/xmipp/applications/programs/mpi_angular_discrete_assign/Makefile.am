bin_PROGRAMS = mpi_angular_discrete_assign
mpi_angular_discrete_assign_SOURCES = main.cpp
mpi_angular_discrete_assign_LDADD = \
	$(top_srcdir)/libraries/reconstruction/libXmippRecons.la \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la \
	@MPILIBS@

CC=@MPICXX@
CXX=@MPICXX@
CXXLINK = $(LIBTOOL) --tag=CXX --mode=link $(CXXLD) $(AM_CXXFLAGS) $(CXXFLAGS) $(AM_LDFLAGS) $(LDFLAGS) -o $@

mpi_angular_discrete_assign_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/libraries
