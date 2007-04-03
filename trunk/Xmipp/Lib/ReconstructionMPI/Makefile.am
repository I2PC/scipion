lib_LTLIBRARIES = libXmippReconsMPI_@host_os@.la 

libXmippReconsMPI_@host_os@_la_SOURCES = Programs/Src/Prog_MPI_Run.cc


INCLUDES = -I$(top_srcdir)/Lib

EXTRA_DIST = Programs/Prog_MPI_Run.hh

docs:
	@echo "Creating documentation in Reconstruction ..."
	@doc++ -f -B @XMIPPDIR@/Extra_Docs/banner.html -d Docs documentation
