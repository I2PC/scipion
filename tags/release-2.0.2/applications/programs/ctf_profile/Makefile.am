bin_PROGRAMS = ctf_profile
ctf_profile_SOURCES = main.cpp
ctf_profile_LDADD = \
	$(top_srcdir)/libraries/reconstruction/libXmippRecons.la \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la \
        $(top_srcdir)/external/cuba/libCuba.la

ctf_profile_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/libraries
