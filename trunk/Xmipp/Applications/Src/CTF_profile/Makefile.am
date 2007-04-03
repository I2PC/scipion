bin_PROGRAMS = ctf_profile
ctf_profile_SOURCES = CTF_profile.cc
ctf_profile_LDADD = \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la
ctf_profile_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/Lib
