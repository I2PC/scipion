bin_PROGRAMS = mpi_MLrefine3D
mpi_MLrefine3D_SOURCES = mpi_MLrefine3D.cc
mpi_MLrefine3D_LDADD = $(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la \
	@MPILIBS@

CC=@MPICXX@
CXX=@MPICXX@
CXXLINK = $(LIBTOOL) --tag=CXX --mode=link $(CXXLD) $(AM_CXXFLAGS) $(CXXFLAGS) $(AM_LDFLAGS) $(LDFLAGS) -o $@

mpi_MLrefine3D_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/Lib
