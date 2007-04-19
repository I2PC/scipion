bin_PROGRAMS = morphology
morphology_SOURCES = main.cpp
morphology_LDADD = \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la

morphology_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/libraries
