bin_PROGRAMS = MPI_Art
MPI_Art_SOURCES = mpi_art.cc
MPI_Art_LDADD = $(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la \
	@MPILIBS@

CC=@MPICXX@
CXX=@MPICXX@
	
#MPI_Art_LDFLAGS = -all-static

INCLUDES = -I$(top_srcdir)/Lib
