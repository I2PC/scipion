bin_PROGRAMS = find_center2d

find_center2d_SOURCES = findcenter.cc findc.cc
find_center2d_LDADD = $(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la \
	$(top_srcdir)/Lib/OldXmipp/libOldXmipp_@host_os@.la

find_center2d_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/Lib
