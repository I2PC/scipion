bin_PROGRAMS = mpi_ml_align2d
mpi_ml_align2d_SOURCES = mpi_MLalign2D.cc
mpi_ml_align2d_LDADD = \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la \
	@MPILIBS@

CC=@MPICXX@
CXX=@MPICXX@
CXXLINK = $(LIBTOOL) --tag=CXX --mode=link $(CXXLD) $(AM_CXXFLAGS) $(CXXFLAGS) $(AM_LDFLAGS) $(LDFLAGS) -o $@

mpi_ml_align2d_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/Lib
