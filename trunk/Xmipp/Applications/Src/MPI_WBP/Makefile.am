bin_PROGRAMS = mpi_reconstruct_wbp
mpi_reconstruct_wbp_SOURCES = mpi_wbp.cc
mpi_reconstruct_wbp_LDADD = \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la \
	@MPILIBS@

CC=@MPICXX@
CXX=@MPICXX@
CXXLINK = $(LIBTOOL) --tag=CXX --mode=link $(CXXLD) $(AM_CXXFLAGS) $(CXXFLAGS) $(AM_LDFLAGS) $(LDFLAGS) -o $@

mpi_reconstruct_wbp_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/Lib
