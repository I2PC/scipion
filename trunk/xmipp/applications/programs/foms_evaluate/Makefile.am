bin_PROGRAMS = foms_evaluate
foms_evaluate_SOURCES = main.cpp
foms_evaluate_LDADD = \
	$(top_srcdir)/libraries/interface/libXmippInterface.la \
	$(top_srcdir)/libraries/reconstruction/libXmippRecons.la \
	$(top_srcdir)/libraries/data/libXmippData.la

foms_evaluate_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/libraries
