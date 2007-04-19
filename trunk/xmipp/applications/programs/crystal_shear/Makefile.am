bin_PROGRAMS = crystal_shear
crystal_shear_SOURCES = main.cpp
crystal_shear_LDADD = \
	$(top_srcdir)/libraries/interface/libXmippInterface.la \
	$(top_srcdir)/libraries/reconstruction/libXmippRecons.la \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la

crystal_shear_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/libraries
