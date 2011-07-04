# -*- makefile -*-
#
# Author: Christophe Prud'homme(prudhomm@debian.org)
#
# Copyright (C) 2004 Christophe Prud'homme
#
# Distributed under the GPL(GNU Public License):
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#

SUFFIXES		= .cpp .hpp .idl .c .h .f .F .o .moc .cc

INCLUDES		= -I$(top_srcdir)/include -I$(srcdir)/.. -I$(top_srcdir)/examples/matrices/ -I$(top_srcdir)/examples/matrices/sym 
LDADD			= $(top_builddir)/src/libarpack++.la

check_PROGRAMS		= dsymreg dsymshf dsymgreg dsymgshf dsymgbkl dsymgcay

dsymreg_SOURCES		= dsymreg.cc
dsymshf_SOURCES		= dsymshf.cc
dsymgreg_SOURCES	= dsymgreg.cc
dsymgshf_SOURCES	= dsymgshf.cc
dsymgbkl_SOURCES	= dsymgbkl.cc
dsymgcay_SOURCES	= dsymgcay.cc

TESTS			= $(check_PROGRAMS)