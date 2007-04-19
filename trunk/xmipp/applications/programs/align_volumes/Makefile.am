bin_PROGRAMS = align_volumes
align_volumes_SOURCES = main.cpp
align_volumes_LDADD = \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la

align_volumes_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/libraries
