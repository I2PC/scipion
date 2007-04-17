bin_PROGRAMS = break_symmetry
break_symmetry_SOURCES = break_sym.cc
break_symmetry_LDADD = \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la

break_symmetry_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/Lib
