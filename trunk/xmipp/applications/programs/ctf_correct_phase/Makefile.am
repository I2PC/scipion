bin_PROGRAMS = ctf_correct_phase
ctf_correct_phase_SOURCES = main.cpp
ctf_correct_phase_LDADD = \
	$(top_srcdir)/libraries/reconstruction/libXmippRecons.la \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la

ctf_correct_phase_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/libraries

