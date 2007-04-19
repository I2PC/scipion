bin_PROGRAMS = convert_img2data
convert_img2data_SOURCES = main.cpp
convert_img2data_LDADD = \
	$(top_srcdir)/libraries/data/libXmippData.la \
	$(top_srcdir)/external/bilib/libBilib.la

convert_img2data_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/libraries
