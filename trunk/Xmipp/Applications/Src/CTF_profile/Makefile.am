bin_PROGRAMS = CTF_profile
CTF_profile_SOURCES = CTF_profile.cc
CTF_profile_LDADD = $(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la

CTF_profile_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/Lib
