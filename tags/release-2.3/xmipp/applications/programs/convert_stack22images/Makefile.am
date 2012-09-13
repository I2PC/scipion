bin_PROGRAMS = stack22images
stack22images_SOURCES = main.cpp
stack22images_LDADD = \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la \
        $(top_srcdir)/external/cuba/libCuba.la

stack22images_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/libraries
