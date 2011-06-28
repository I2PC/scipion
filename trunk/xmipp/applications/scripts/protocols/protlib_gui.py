#!/usr/bin/env python
'''
#/***************************************************************************
# * Authors:     J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *
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
# *  e-mail address 'xmipp@cnb.csic.es'
# ***************************************************************************
 '''
 
import sys
import os
import string
from Tkinter import *
import tkFont
import tkMessageBox
from protlib_utils import *

class ProtocolStyle():
    ''' Class to define some style settings like font, colors, etc '''
    def __init__(self, configModuleName=None):        
        #Font
        self.FontName = "Helvetica"
        self.FontSize = 10
        self.ButtonFontSize = self.FontSize
        #TextColor
        self.CitationTextColor = "dark olive green"
        self.LabelTextColor = "black"
        self.SectionTextColor = "blue4"
        #Background Color
        self.BgColor = "white"
        self.LabelBgColor = self.BgColor
        self.HighlightBgColor = self.BgColor
        self.ButtonBgColor = "LightBlue"
        self.ButtonActiveBgColor = "LightSkyBlue"
        self.EntryBgColor = "lemon chiffon" 
        self.ExpertLabelBgColor = "light salmon"
        #Color
        self.ListSelectColor = "DeepSkyBlue4"
        self.BooleanSelectColor = "DeepSkyBlue4"
        #Dimensions limits
        self.MaxHeight = 600
        self.MaxWidth = 800
        self.MaxFontSize = 14
        self.MinFontSize = 6
        self.WrapLenght = self.MaxWidth / 2
        
        if configModuleName:
            self.load(configModuleName)
                    
    def load(self, configModuleName):
        mod = loadModule(configModuleName,False)
        if mod:
            modDir = dir(mod)
            selfDir = dir(self)
            for a in modDir:
                if a in selfDir and not a.startswith('_'):
                        self.__dict__[a] = mod.__dict__[a]
                    
    def createFonts(self):
        self.Font = tkFont.Font(family=self.FontName, size=self.FontSize, weight=tkFont.BOLD)

class AutoScrollbar(Scrollbar):
    '''A scrollbar that hides itself if it's not needed.'''
    def set(self, lo, hi):
        if float(lo) <= 0.0 and float(hi) >= 1.0:
            self.tk.call("grid", "remove", self)
        else:
            self.grid()
        Scrollbar.set(self, lo, hi)
       

class BasicGUI(): 
    def __init__(self):
        pass
        
    def createScrollableCanvas(self):
        """Create an scrollable canvas
        It will set .frame and .canvas new properties
        """
        vscrollbar = AutoScrollbar(self.master)
        vscrollbar.grid(row=0, column=1, sticky=N + S)
        hscrollbar = AutoScrollbar(self.master, orient=HORIZONTAL)
        hscrollbar.grid(row=1, column=0, sticky=E + W)
        self.canvas = Canvas(self.master, background=self.style.BgColor,
                        yscrollcommand=vscrollbar.set,
                        xscrollcommand=hscrollbar.set)
        self.canvas.grid(row=0, column=0, sticky=N + S + E + W)
        vscrollbar.config(command=self.canvas.yview)
        hscrollbar.config(command=self.canvas.xview)
        self.master.grid_rowconfigure(0, weight=1)
        self.master.grid_columnconfigure(0, weight=1)
        self.frame = Frame(self.canvas, background=self.style.BgColor)
        self.frame.rowconfigure(0, weight=1)
        self.frame.columnconfigure(0, weight=1)
    
    def createBasicGUI(self):
        """Perform some basic GUI initializations.
        - create the main window
        - create style and fonts
        - create scrollable canvas
        """
        self.master = Tk()
        self.style = ProtocolStyle('config_gui')
        self.style.createFonts()
        
    def resize(self):
        height = self.frame.winfo_reqheight() + 25
        width = self.frame.winfo_reqwidth() + 25
        if height > self.style.MaxHeight:
           height = self.style.MaxHeight
        if width > self.style.MaxWidth:
           width = self.style.MaxWidth
        self.master.geometry("%dx%d%+d%+d" % (width, height, 0, 0))
        
    def launchCanvas(self):
        # Launch the window
        self.canvas.create_window(0, 0, anchor=NW, window=self.frame)
        self.frame.update_idletasks()
        self.canvas.config(scrollregion=self.canvas.bbox("all")) 
    
    def launchGUI(self):
        self.fillGUI()
        self.launchCanvas() 
        self.resize()      
        self.master.mainloop() 
    