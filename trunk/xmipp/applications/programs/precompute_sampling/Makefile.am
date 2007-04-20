bin_PROGRAMS = precompute_sampling
precompute_sampling_SOURCES = main.cpp
precompute_sampling_LDADD = \
	$(top_srcdir)/libraries/reconstruction/libXmippRecons.la \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la

precompute_sampling_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/libraries

