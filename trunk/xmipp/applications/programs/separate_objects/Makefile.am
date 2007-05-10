bin_PROGRAMS = separate_objects
separate_objects_SOURCES = main.cpp
separate_objects_LDADD = \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la \
        $(top_srcdir)/external/cuba/libCuba.la

separate_objects_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/libraries
