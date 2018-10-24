# **************************************************************************
# *
# * Authors:    Yunior C. Fonseca Reyna (cfonseca@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os.path
import stat

import Tkinter as tk
import gui
from tree import BoundTree, TreeProvider
from text import TaggedText, openTextFileEditor
from widgets import Button, HotButton


class PluginManagerWindow(gui.Window):
    """ Windows to hold a plugin manager frame inside. """

    def __init__(self, title, master=None, **kwargs):
        if 'minsize' not in kwargs:
            kwargs['minsize'] = (800, 400)
        gui.Window.__init__(self, title, master, **kwargs)

    def setBrowser(self, browser, row=0, column=0):
        self.browser = browser
        browser.grid(row=row, column=column, sticky='news')
        self.itemConfig = browser.tree.itemConfig


class PluginManager(PluginManagerWindow):
    """ Windows to hold a file browser frame inside. """

    def __init__(self, title, master=None, path=None,
                 onSelect=None, shortCuts=None, **kwargs):
        PluginManagerWindow.__init__(self, title, master, **kwargs)
