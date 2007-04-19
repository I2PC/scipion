bin_PROGRAMS = header_apply
header_apply_SOURCES = main.cpp
header_apply_LDADD = \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la

header_apply_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/libraries
