bin_PROGRAMS = ctf_estimate_from_psd
ctf_estimate_from_psd_SOURCES = main.cpp
ctf_estimate_from_psd_LDADD = \
	$(top_srcdir)/libraries/reconstruction/libXmippRecons.la \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la \
        $(top_srcdir)/external/cuba/libCuba.la

ctf_estimate_from_psd_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/libraries
