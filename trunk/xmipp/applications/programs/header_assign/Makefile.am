bin_PROGRAMS = header_assign
header_assign_SOURCES = main.cpp
header_assign_LDADD = \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la

header_assign_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/libraries
