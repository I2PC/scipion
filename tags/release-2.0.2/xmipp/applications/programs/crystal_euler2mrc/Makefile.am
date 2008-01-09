bin_PROGRAMS = crystal_euler2mrc
crystal_euler2mrc_SOURCES = main.cpp
crystal_euler2mrc_LDADD = \
	$(top_srcdir)/libraries/reconstruction/libXmippRecons.la \
	$(top_srcdir)/libraries/reconstruction/libXmippRecons_Interface.la \
	$(top_srcdir)/libraries/interface/libXmippInterface.la \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la \
        $(top_srcdir)/external/cuba/libCuba.la

crystal_euler2mrc_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/libraries
