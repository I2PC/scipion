bin_PROGRAMS = ctf_correct_phase
ctf_correct_phase_SOURCES = correctphase.cc
ctf_correct_phase_LDADD = \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la

ctf_correct_phase_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/Lib

