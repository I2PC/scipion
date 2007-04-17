bin_PROGRAMS = reconstruct_wbp
reconstruct_wbp_SOURCES = main.cpp
reconstruct_wbp_LDADD = \
	$(top_srcdir)/libraries/reconstruction/libXmippRecons.la \
	$(top_srcdir)/libraries/data/libXmippData.la

reconstruct_wbp_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/libraries
