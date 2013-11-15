#!bin/xmipp_python
'''
/***************************************************************************
 * Authors:     J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
 *
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/
 '''

import os
import Tkinter as tk
import tkFont as font
import ttk  

from protlib_gui_ext import centerWindows, XmippButton, registerCommonFonts, showInfo, showError, OutputText,\
    getGeometry, Fonts, TaggedText
from protlib_filesystem import getXmippPath, getXmippVersion
        
class AboutFrame(tk.Frame):

    def __init__(self, master):
        tk.Frame.__init__(self, master)
        self.master = master
        self.tabs = {}
        self.createPanels()
        
    def addTab(self, text, file, processLine=None):
        tab = tk.Frame(self.nb)
        self.nb.add(tab, text=text)
        self.tabs[text] = tab
        txt = TaggedText(tab, width=60, height=20, bg='lightgray');
        txt.frame.grid(column=0, row=0, sticky='nsew', padx=10, pady=10)
        f = open(getXmippPath(file))
        for line in f:
            if processLine:
                line = processLine(line)
            if line:
                txt.addLine(line)
        f.close()
        txt.setReadOnly(True)
        return txt
        
    def createPanels(self):
        root = self.master
        root.columnconfigure(0, minsize=100, weight=1)
        root.columnconfigure(1, minsize=360, weight=4)
        root.rowconfigure(0, weight=1)
        registerCommonFonts()
        #left panel
        bgColor = 'white'
        leftFrame = tk.Frame(root, bg=bgColor)
        leftFrame.columnconfigure(0, weight=1)
        leftFrame.rowconfigure(0, minsize=100)
        imgPath = getXmippPath('resources', 'xmipp_logo.gif')
        self.img = tk.PhotoImage(file=imgPath)
        tk.Label(leftFrame, image=self.img, bg=bgColor).grid(column=0, row=0, sticky='we')
        tk.Label(leftFrame, text='Xmipp '+ getXmippVersion(),  font=Fonts['button'], bg=bgColor).grid(column=0, row=1, sticky='we')
#       TODO: insert revision extracting it from git repository
#        tk.Label(leftFrame, text='r12.4.3.11834', bg=bgColor).grid(column=0, row=2, sticky='we')
        leftFrame.grid(column=0, row=0, sticky='nsew', padx=5, pady=5, rowspan=2)
        self.grid(column=1, row=0, sticky='nsew', padx=5, pady=5)
        
        self.nb = ttk.Notebook(self.master)
        self.nb.grid(row=0, column=1, sticky='nsew', padx=5, pady=5)
        
        self.addTab("About", "README")
        self.addTab("Authors", "AUTHORS", getAuthorLine)
        self.addTab("Software", "SOFTWARE")
        self.addTab("License", "COPYING")
        
        self.btn = XmippButton(self.master, text='Close')
        self.btn.grid(row=1, column=1, sticky='se', padx=5, pady=5)        
        self.btn.config(command=self.close)
        self.btn.bind('<Return>', func=self.close)
        self.btn.focus_set()
    
    def close(self, e=None):
        self.master.destroy()

# Helper function to parse text file lines   
def getAuthorLine(inputLine):
    if '(' in inputLine:
        name, mail = inputLine.split('(')
        return  "<%(name)s> (%(mail)s" % locals()
    return None

def createAboutDialog():
    root = tk.Toplevel()
    root.withdraw()
    root.title("Xmipp About")
    root.minsize(width=460, height=350)
    root.columnconfigure(0, weight=1)
    root.rowconfigure(0, weight=1)
    
    AboutFrame(root).grid(row=0, column=0, sticky='nsew')    
    centerWindows(root)
    root.deiconify()
    root.mainloop()
    
