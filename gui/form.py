# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************
"""
This modules implements the automatic
creation of protocol form GUI from its
params definition.
"""

import Tkinter as tk
import tkFont

from gui import configureWeigths, Window
from text import TaggedText
from pyworkflow.protocol.params import *


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
                    
    def createFonts(self):
        self.Font = tkFont.Font(family=self.FontName, size=self.FontSize, weight=tkFont.BOLD)


class SectionFrame(tk.Frame):
    def __init__(self, master, frameParam, **args):
        tk.Frame.__init__(self, master, **args)
        self.headerFrame = tk.Frame(self, bd=2, relief=tk.RAISED, bg='#7D0709')
        self.headerFrame.grid(row=0, column=0, sticky='new')
        self.headerLabel = tk.Label(self.headerFrame, text=frameParam.label.get(), fg='white', bg='#7D0709')
        self.headerLabel.grid(row=0, column=0, sticky='nw')
        self.contentFrame = tk.Frame(self, bg='white', bd=0)
        self.contentFrame.grid(row=1, column=0, sticky='news', padx=5, pady=5)
        configureWeigths(self.contentFrame)
        self.columnconfigure(0, weight=1)

    
class FormWindow(Window): 
    def __init__(self, protocol, master=None, **args):
        Window.__init__(self, "Project: test_project", master, weight=False,
                        icon='scipion_bn.xbm', **args)

        self.fontBig = tkFont.Font(size=12, family='verdana', weight='bold')
        self.font = tkFont.Font(size=10, family='verdana')#, weight='bold')
        self.fontBold = tkFont.Font(size=10, family='verdana', weight='bold')
        
        headerFrame = tk.Frame(self.root)
        headerFrame.grid(row=0, column=0, sticky='new')
        headerLabel = tk.Label(headerFrame, text='Protocol: Import micrographs', font=self.fontBig)
        headerLabel.grid(row=0, column=0, padx=5, pady=5)
        
        text = TaggedText(self.root, width=40, height=15)
        text.grid(row=1, column=0, sticky='news')
        text.config(state=tk.DISABLED)
        self.protocol = protocol
        self.createSections(protocol._definition, text)
        
        btnFrame = tk.Frame(self.root)
        btnFrame.columnconfigure(0, weight=1)
        btnFrame.grid(row=2, column=0, sticky='sew')
        # Add create project button
        btn = tk.Button(btnFrame, text='Execute', fg='white', bg='#7D0709', font=self.font, 
                        activeforeground='white', activebackground='#A60C0C')
        btn.grid(row=0, column=0, padx=(5, 100), pady=5, sticky='se')
        
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(1, weight=1)
        
    
    def fillSection(self, sectionParam, sectionFrame):
        parent = sectionFrame.contentFrame
        r = 0
        for paramName, param in sectionParam.getChildParams():
            # Create the label
            if param.isImportant:
                f = self.fontBold
            else:
                f = self.font
            label = tk.Label(parent, text=param.label.get(), bg='white', font=f)
            label.grid(row=r, column=0, sticky='ne', padx=2, pady=2)
            # Create widgets for each type of param
            t = type(param)
            if t is BooleanParam:
                content = tk.Frame(parent, bg='white')
                rb1 = tk.Radiobutton(content, text='Yes', bg='white')
                rb1.grid(row=0, column=0, padx=2)
                rb2 = tk.Radiobutton(content, text='No', bg='white')
                rb2.grid(row=0, column=1, padx=2)
            else:
                v = tk.StringVar()
                content = tk.Entry(parent, width=25, textvariable=v)
                value = getattr(self.protocol, paramName, None)
                if value is not None:
                    v.set(value.get(''))
            
            content.grid(row=r, column=1, padx=2, pady=2, sticky='w')
            r += 1
            
    def createSections(self, formDef, text):
        """Load the list of projects"""
        r = 0
        parent = tk.Frame(text)    
        parent.columnconfigure(0, weight=1)
        
        for s in formDef:
            frame = SectionFrame(parent, s, bg='white')
            frame.grid(row=r, column=0, padx=10, pady=5, sticky='new')
            frame.columnconfigure(0, minsize=300)
            self.fillSection(s, frame)
            r += 1
        
        # with Windows OS
        self.root.bind("<MouseWheel>", lambda e: text.scroll(e))
        # with Linux OS
        self.root.bind("<Button-4>", lambda e: text.scroll(e))
        self.root.bind("<Button-5>", lambda e: text.scroll(e)) 
        text.window_create(tk.INSERT, window=parent)
        
  
if __name__ == '__main__':
    # Just for testing
    from pyworkflow.em import ProtImportMicrographs
    p = ProtImportMicrographs()
    p.sphericalAberration.set(2.3)
    p.samplingRate.set('5.4')
    w = FormWindow(p)
    w.show()
    
   


