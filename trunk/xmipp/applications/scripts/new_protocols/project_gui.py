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
('3D', ['Initial Model', 'Model Refinement']) ]

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
        
    def createMainMenu(self):
        self.menubar = Menu(self.root)
        self.fileMenu = Menu(self.root, tearoff=0)
        self.fileMenu.add_command(label="Exit", command=self.onExit)
        self.menubar.add_cascade(label="File", menu=self.fileMenu)
        
    def addTbLabel(self, text, row):
        '''Add a label to left toolbar'''
        Font = tkFont.Font(family=FontName, size=FontSize+4, weight=tkFont.BOLD)
        label = Label(self.toolbar, text=text, font=Font, fg=SectionTextColor)
        label.grid(row = row, column=0)
        
    def addTbButton(self, text, row):
        '''Add a button to left toolbar'''
        Font = tkFont.Font(family=FontName, size=FontSize, weight=tkFont.BOLD)
        btn = Menubutton(self.toolbar, bd = 1, text=text, font=Font, relief=RAISED,
                         bg=ButtonBgColor, activebackground=ButtonActiveBgColor)
        btn.grid(row = row, column = 0, sticky=W+E, pady=2, padx=5)
        menu = Menu(btn, bg=ButtonBgColor, activebackground=ButtonActiveBgColor, font=Font)
        btn['menu'] = menu
        menu.add_command(label="option1")
        menu.add_command(label="option2")
#        self.bGet = Menubutton(self.frame, text=label, relief=RAISED,  
#                                       width=30, bg=ButtonBackgroundColour, 
#                                        activebackground=ButtonActiveBackgroundColour,
#                                        highlightbackground=HighlightBackgroundColour); 
#                self.bMenu = Menu(self.bGet)
#                self.bGet["menu"] = self.bMenu          
#                if not 'childs' in button:
#                    print "ERROR: button '", label , "' doesn't have script neither childs."
#                    sys.exit() 
#                childs = button['childs'].split(", ")
#                for child in childs:
#                    option = setup.LaunchButtons[child];
#                    label = option['title']
#                    self.bMenu.add_radiobutton(label=label, variable=self.which_setup, 
#                                             value = option['script'], command=self.GuiLaunchSetup,
#                                             activebackground=ButtonActiveBackgroundColour,                                             
#                                             selectcolor=ButtonBackgroundColour);

    def createGUI(self):
      
        self.root.title("Popup menu")
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