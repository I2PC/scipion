bin_PROGRAMS = MPI_WBP
MPI_WBP_SOURCES = mpi_wbp.cc
MPI_WBP_LDADD = $(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la \
	@MPILIBS@

CC=@MPICXX@
CXX=@MPICXX@
	
MPI_WBP_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/Lib
