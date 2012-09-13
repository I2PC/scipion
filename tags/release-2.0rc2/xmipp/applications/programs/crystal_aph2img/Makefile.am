bin_PROGRAMS = crystal_aph2img
crystal_aph2img_SOURCES = main.cpp
crystal_aph2img_LDADD = \
	$(top_srcdir)/libraries/reconstruction/libXmippRecons.la \
	$(top_srcdir)/libraries/reconstruction/libXmippRecons_Interface.la \
	$(top_srcdir)/libraries/interface/libXmippInterface.la \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la \
        $(top_srcdir)/external/cuba/libCuba.la

crystal_aph2img_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/libraries
