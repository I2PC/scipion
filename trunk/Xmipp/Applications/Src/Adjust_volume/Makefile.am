bin_PROGRAMS = adjust_volume
adjust_volume_SOURCES = adjust_volume.cc
adjust_volume_LDADD = $(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
        $(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la

adjust_volume_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/Lib

