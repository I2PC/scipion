bin_PROGRAMS = edit
edit_SOURCES = edit.cc textedit.cc textedit_moc.cc

edit_LDADD = $(top_srcdir)/Lib/XmippData/libXmippData_@host_os@.la 
	
edit_LDFLAGS = @QT_LIBS@ 

INCLUDES = -I$(top_srcdir)/Lib @QT_CXXFLAGS@

MOSTLYCLEANFILES = *_moc.cc

textedit_moc.cc : 
	@QT_MOC@ textedit.hh -o textedit_moc.cc

EXTRA_DIST = textedit.hh
