bin_PROGRAMS = edit
edit_SOURCES = main.cpp textedit.cpp
nodist_edit_SOURCES = textedit_moc.cpp

edit_LDADD = \
	$(top_srcdir)/libraries/data/libXmippData.la  \
	$(top_srcdir)/external/bilib/libBilib.la

edit_LDFLAGS = @QT_LIBS@

INCLUDES = -I$(top_srcdir)/libraries @QT_CXXFLAGS@

MOSTLYCLEANFILES = *_moc.cpp

textedit_moc.cpp :
	@QT_MOC@ textedit.h -o textedit_moc.cpp

EXTRA_DIST = textedit.h
