bin_PROGRAMS = phantom_create_random
phantom_create_random_SOURCES = main.cpp
phantom_create_random_LDADD = \
	$(top_srcdir)/libraries/interface/libXmippInterface.la \
	$(top_srcdir)/libraries/reconstruction/libXmippRecons.la \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la \
        $(top_srcdir)/external/cuba/libCuba.la

phantom_create_random_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/libraries
