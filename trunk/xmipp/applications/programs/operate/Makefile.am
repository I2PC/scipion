bin_PROGRAMS = operate
operate_SOURCES = main.cpp
operate_LDADD = \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la

operate_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/libraries
