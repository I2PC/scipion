bin_PROGRAMS = assign_multiple_CTFs
assign_multiple_CTFs_SOURCES = assign_multiple_CTFs.cc
assign_multiple_CTFs_LDADD = $(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippInterface_@host_os@.la

assign_multiple_CTFs_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/Lib
