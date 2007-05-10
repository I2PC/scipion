bin_PROGRAMS = create_surface
create_surface_SOURCES = main.cpp
create_surface_LDADD = \
	$(top_srcdir)/libraries/reconstruction/libXmippRecons.la \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la \
        $(top_srcdir)/external/cuba/libCuba.la

create_surface_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/libraries


