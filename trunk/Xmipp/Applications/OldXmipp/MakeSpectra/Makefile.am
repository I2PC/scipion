bin_PROGRAMS = old_make_spectra

old_make_spectra_SOURCES = makespectra.cc
old_make_spectra_LDADD = \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la \
	$(top_srcdir)/Lib/OldXmipp/libOldXmipp_@host_os@.la
old_make_spectra_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/Lib -I$(top_srcdir)/Lib/OldXmipp
