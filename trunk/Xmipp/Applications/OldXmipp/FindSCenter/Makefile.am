bin_PROGRAMS = findcenter

findcenter_SOURCES = findcenter.cc findc.cc
findcenter_LDADD = $(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la \
	$(top_srcdir)/Lib/OldXmipp/libOldXmipp_@host_os@.la

findcenter_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/Lib
