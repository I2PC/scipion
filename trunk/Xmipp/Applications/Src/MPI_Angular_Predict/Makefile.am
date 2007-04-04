bin_PROGRAMS = mpi_angular_predict
mpi_angular_predict_SOURCES = MPI_angular_predict.cc
mpi_angular_predict_LDADD = \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la \
	@MPILIBS@

CC=@MPICXX@
CXX=@MPICXX@
CXXFLAGS=-w $(AM_CXXFLAGS)
CXXLINK = $(LIBTOOL) --tag=CXX --mode=link $(CXXLD) $(AM_CXXFLAGS) $(CXXFLAGS) $(AM_LDFLAGS) $(LDFLAGS) -o $@

mpi_angular_predict_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/Lib
