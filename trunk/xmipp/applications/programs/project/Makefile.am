bin_PROGRAMS = project
project_SOURCES = main.cpp
project_LDADD = \
	$(top_srcdir)/libraries/interface/libXmippInterface.la \
	$(top_srcdir)/libraries/reconstruction/libXmippRecons.la \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la \
        $(top_srcdir)/external/cuba/libCuba.la

project_LDFLAGS = @EXEC_TYPE@

INCLUDES = \
	-I$(top_srcdir)/libraries \
	-I$(top_srcdir)/external
