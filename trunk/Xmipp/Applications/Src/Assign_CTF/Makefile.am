bin_PROGRAMS = ctf_estimate_from_micrograph
ctf_estimate_from_micrograph_SOURCES = assign_CTF.cc
ctf_estimate_from_micrograph_LDADD = \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la

ctf_estimate_from_micrograph_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/Lib
