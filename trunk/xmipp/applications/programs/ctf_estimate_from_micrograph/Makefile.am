bin_PROGRAMS = ctf_estimate_from_micrograph
ctf_estimate_from_micrograph_SOURCES = main.cpp
ctf_estimate_from_micrograph_LDADD = \
	$(top_srcdir)/libraries/reconstruction/libXmippRecons.la \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la \
        $(top_srcdir)/external/cuba/libCuba.la

ctf_estimate_from_micrograph_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/libraries
