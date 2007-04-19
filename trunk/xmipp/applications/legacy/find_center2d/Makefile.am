bin_PROGRAMS = find_center2d

find_center2d_SOURCES = findc.cpp main.cpp
find_center2d_LDADD = \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la \
	$(top_srcdir)/libraries/legacy/libXmippLegacy.la

find_center2d_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/libraries
