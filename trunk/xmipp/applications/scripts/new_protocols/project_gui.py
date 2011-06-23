#!/usr/bin/env python
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

from Tkinter import *
import tkFont

sections = [\
('Preprocessing', ['Preprocess Micrograph', 'Particles picking', 'Preprocess Particles']),\
('2D', ['Align+Classify', 'Align', 'Classify']),\
('3D', ['Initial Model', 'Model Refinement']),\
('Other', ['Browse'])]

#Font
FontName = "Helvetica"
FontSize = 10

#TextColor
CitationTextColor = "dark olive green"
LabelTextColor = "black"
SectionTextColor = "blue4"

#Background Color
BgColor = "white"
LabelBgColor = BgColor
HighlightBgColor = BgColor
ButtonBgColor = "LightBlue"
ButtonActiveBgColor = "LightSkyBlue"
EntryBgColor = "lemon chiffon" 
ExpertLabelBgColor = "light salmon"

#Color
ListSelectColor = "DeepSkyBlue4"
BooleanSelectColor = "DeepSkyBlue4"

#Dimensions limits
MaxHeight = 800
MaxWidth = 800
MaxFontSize = 14
MinFontSize = 6


        
class XmippProject(Frame):
  
    def __init__(self):
        self.root = Tk()
        #self.root.geometry("250x150+300+300")
        self.frame = Frame(self.root)
        self.createGUI()
        self.root.mainloop()
        
    def dragWindows(self, event):
        btn, menu = self.lastPair
        if btn:
            self.postMenu(btn, menu) #Move the popup menu with the windows
        
    def createMainMenu(self):
        self.menubar = Menu(self.root)
        self.fileMenu = Menu(self.root, tearoff=0)
        self.fileMenu.add_command(label="Exit", command=self.onExit)
        self.menubar.add_cascade(label="File", menu=self.fileMenu)
        self.lastPair = (None, None)
        self.lastIndex = None
        self.root.bind('<Configure>', self.dragWindows)
        
    def addTbLabel(self, text, row):
        '''Add a label to left toolbar'''
        Font = tkFont.Font(family=FontName, size=FontSize+4, weight=tkFont.BOLD)
        label = Label(self.toolbar, text=text, font=Font, fg=SectionTextColor)
        label.grid(row = row, column=0)
        
    def addTbButton(self, text, row):
        '''Add a button to left toolbar'''
        Font = tkFont.Font(family=FontName, size=FontSize, weight=tkFont.BOLD)
        
        menu = Menu(self.frame, bg=ButtonBgColor, activebackground=ButtonBgColor, font=Font, tearoff=0)
        
        
        btn = Button(self.toolbar, bd = 1, text=text, font=Font, relief=RAISED,
                         bg=ButtonBgColor, activebackground=ButtonBgColor)
        btn.grid(row = row, column = 0, sticky=W+E, pady=2, padx=5)
        
        cmd = lambda:self.showPopup(btn, menu)
        menu.add_command(label="option1", command=lambda:self.menuPick(btn, menu, 0))
        menu.add_command(label="option2", command=lambda:self.menuPick(btn, menu, 1))
        btn.config(command=cmd)
        
    def postMenu(self, btn, menu):
        x, y, w = btn.winfo_x(), btn.winfo_y(), btn.winfo_width()
        xroot, yroot = self.root.winfo_x(), self.root.winfo_y()
        menu.post(xroot + x + w + 10, yroot + y)
        btn.config(bg=ButtonActiveBgColor, activebackground=ButtonActiveBgColor)
        
        
    def menuPick(self, btn, menu, index):
        self.postMenu(btn, menu)
        print "index %s, self.lastIndex %s" % (index, self.lastIndex)
        if self.lastIndex:
            menu.entryconfig(self.lastIndex, background=ButtonBgColor, activebackground=ButtonBgColor)
        self.lastIndex = index
        menu.entryconfig(index, background=ButtonActiveBgColor, activebackground=ButtonActiveBgColor)
        
    def showPopup(self, btn, menu):
        lastBtn, lastMenu = self.lastPair
        if lastBtn and lastBtn != btn:
            lastBtn.config(bg=ButtonBgColor)
            lastMenu.unpost()
            
        if lastBtn != btn:
            self.lastPair = (btn, menu)
            self.postMenu(btn, menu)
            
        

    def createGUI(self):
      
        self.root.title("Xmipp Protocols")
        self.createMainMenu()
        
        self.toolbar = Frame(self.root, bd=1, relief=RAISED)

        i = 1
        for k, v in sections:
            self.addTbLabel(k, i)
            i += 1
            for btnName in v:
                self.addTbButton(btnName, i)
                i += 1
            
        self.toolbar.pack(side=LEFT, fill=Y)
        text = Text(self.root)
        text.pack(side=RIGHT, fill=BOTH)
        self.root.config(menu=self.menubar)
        self.frame.pack()
        
       
    def onExit(self):
        self.root.destroy()


if __name__ == '__main__':
    XmippProject()