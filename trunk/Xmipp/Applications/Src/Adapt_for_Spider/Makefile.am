bin_PROGRAMS = adapt_for_spider
adapt_for_spider_SOURCES = adapt_for_spider.cc
adapt_for_spider_LDADD = $(top_srcdir)/Lib/XmippInterface/libXmippInterface_@host_os@.la \
	$(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la

adapt_for_spider_LDFLAGS = @EXEC_TYPE@

INCLUDES = -I$(top_srcdir)/Lib 
