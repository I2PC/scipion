bin_PROGRAMS = mpi_art
mpi_art_SOURCES = mpi_art.cc
mpi_art_LDADD = $(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la \
	@MPILIBS@

CC=@MPICXX@
CXX=@MPICXX@
CXXLINK = $(LIBTOOL) --tag=CXX --mode=link $(CXXLD) $(AM_CXXFLAGS) $(CXXFLAGS) $(AM_LDFLAGS) $(LDFLAGS) -o $@
	
mpi_art_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/Lib
