bin_PROGRAMS = reconstruct_wbp
reconstruct_wbp_SOURCES = wbp.cc
reconstruct_wbp_LDADD = $(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la

reconstruct_wbp_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/Lib
