bin_PROGRAMS = crystal_shear
crystal_shear_SOURCES = deform.cc
crystal_shear_LDADD = \
	$(top_srcdir)/Lib/XmippInterface/libXmippInterface_@host_os@.la \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la
crystal_shear_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/Lib
