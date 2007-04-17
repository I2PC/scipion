bin_PROGRAMS = reconstruct_art
reconstruct_art_SOURCES = main.cpp
reconstruct_art_LDADD = \
	$(top_srcdir)/libraries/reconstruction/libXmippRecons.la \
	$(top_srcdir)/libraries/data/libXmippData.la

reconstruct_art_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/libraries
