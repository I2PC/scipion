bin_PROGRAMS = makespectra

makespectra_SOURCES = makespectra.cc
makespectra_LDADD = \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la \
	$(top_srcdir)/Lib/OldXmipp/libOldXmipp_@host_os@.la
makespectra_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/Lib -I$(top_srcdir)/Lib/OldXmipp
