# -*- Mode : makefile -*-
#
# Copyright (C) 2003 Christophe Prud'homme (prudhomm@debian.org)
#
# public domain software
#  

SUFFIXES		= .cpp .hpp .idl .c .h .f .F .o .moc 

INCLUDES		= -I$(top_srcdir)/include

lib_LTLIBRARIES		= libarpack++.la

libarpack___la_SOURCES	= debug.cpp arrseig.cpp arerror.cpp
libarpack___la_LIBADD	= $(all_libs)
libarpack___la_LDFLAGS	= -version-info 2:0:0
