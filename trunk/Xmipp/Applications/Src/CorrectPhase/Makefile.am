bin_PROGRAMS = correct_phase
correct_phase_SOURCES = correctphase.cc
correct_phase_LDADD = \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la

correct_phase_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/Lib

