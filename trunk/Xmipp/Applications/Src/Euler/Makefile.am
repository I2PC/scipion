bin_PROGRAMS = euler
euler_SOURCES = euler.cc
euler_LDADD =$(top_srcdir)/Lib/XmippInterface/libXmippInterface_@host_os@.la \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la

INCLUDES = -I$(top_srcdir)/Lib

