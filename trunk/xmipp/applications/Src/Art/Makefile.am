bin_PROGRAMS = reconstruct_art
reconstruct_art_SOURCES = art.cc
reconstruct_art_LDADD = \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la

reconstruct_art_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/Lib
