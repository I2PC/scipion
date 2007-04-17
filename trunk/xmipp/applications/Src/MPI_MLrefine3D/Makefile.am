bin_PROGRAMS = mpi_ml_refine3d
mpi_ml_refine3d_SOURCES = mpi_MLrefine3D.cc
mpi_ml_refine3d_LDADD = \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la \
	@MPILIBS@

CC=@MPICXX@
CXX=@MPICXX@
CXXLINK = $(LIBTOOL) --tag=CXX --mode=link $(CXXLD) $(AM_CXXFLAGS) $(CXXFLAGS) $(AM_LDFLAGS) $(LDFLAGS) -o $@

mpi_ml_refine3d_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/Lib
