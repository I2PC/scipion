bin_PROGRAMS = data2vol
data2vol_SOURCES = data2vol.cc
data2vol_LDADD = $(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la \
	$(top_srcdir)/Lib/Classification/libXmippClassif_@host_os@.la

data2vol_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/Lib
