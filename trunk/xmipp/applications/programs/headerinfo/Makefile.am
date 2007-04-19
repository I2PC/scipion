bin_PROGRAMS = headerinfo
headerinfo_SOURCES = main.cpp
headerinfo_LDADD = \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la

headerinfo_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/libraries
