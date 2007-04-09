bin_PROGRAMS = break_sym
break_sym_SOURCES = break_sym.cc
break_sym_LDADD = \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la

break_sym_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/Lib
