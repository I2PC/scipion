bin_PROGRAMS = deconvoluteCTF
deconvoluteCTF_SOURCES = deconvoluteCTF.cc
deconvoluteCTF_LDADD = \
	$(top_srcdir)/Lib/Reconstruction/libXmippRecons_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la

deconvoluteCTF_LDFLAGS = @EXEC_TYPE@
INCLUDES = -I$(top_srcdir)/Lib
