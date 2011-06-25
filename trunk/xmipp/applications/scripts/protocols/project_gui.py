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

sections = [
('Preprocessing', 
   [['Preprocess Micrograph'], 
    ['Particles picking'], 
    ['Preprocess Particles']]),
('2D', 
   [['Align+Classify', 'ML2D', 'CL2D'], 
    ['Align', 'ML2D', 'CL2D'], 
    ['Classify', 'KerDenSOM', 'Rotational Spectra']]),
('3D', 
   [['Initial Model', 'Common Lines', 'Random Conical Tilt'], 
    ['Model Refinement']]),
('Other', [['Browse']])]

launchDict = {
              'Preprocess Micrograph': 'xmipp_protocol_preprocess_micrographs.py',
              'Particles picking':     'xmipp_protocol_particle_pick.py', 
              'Preprocess Particles':  'xmipp_protocol_preprocess_particles.py', 
              'ML2D':                  'xmipp_protocol_ml2d.py',
              'CL2D':                  'xmipp_protocol_cl2d.py',
              'KerDenSOM':             'xmipp_protocol_kerdensom.py',
              'Rotational Spectra':    'xmipp_protocol_rotspectra.py',
              'Common Lines':          'xmipp_protocol_rct.py',
              'Random conical Tilt':    'xmipp_protocol_commonlines.py'
              }
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
        self.createGUI()
        self.root.mainloop()
        
    def createMainMenu(self):
        self.menubar = Menu(self.root)
        self.fileMenu = Menu(self.root, tearoff=0)
        self.fileMenu.add_command(label="Exit", command=self.onExit)
        self.menubar.add_cascade(label="File", menu=self.fileMenu)
        
   
    def addTbLabel(self, text, row):
        '''Add a label to left toolbar'''
        Font = tkFont.Font(family=FontName, size=FontSize+1, weight=tkFont.BOLD)
        label = Label(self.toolbar, text=text, font=Font, fg=SectionTextColor)
        label.grid(row = row, column=0)
        
    def addTbButton(self, row, text, opts=[]):
        '''Add a button to left toolbar'''
        
        btn = Button(self.toolbar, bd = 1, text=text, font=self.ButtonFont, relief=RAISED,
                         bg=ButtonBgColor, activebackground=ButtonBgColor)
        btn.grid(row = row, column = 0, sticky=W+E, pady=2, padx=5)
        
        btnFrame = Frame(self.frame, width=100, bd=1, relief=GROOVE)
        btnFrame.columnconfigure(0, minsize=130)
        label = Label(btnFrame, text='Protocols:')
        label.grid(row=0, column=0, pady=5)
        if len(opts) > 0:
            i = 0
            for o in opts:
                i += 1
                label = Label(btnFrame, text=o)
                label.grid(row=2*i, column=0, sticky=W)
                btnLaunch = Button(btnFrame, text='Launch')
                btnLaunch.grid(row=2*i+1, column=0, pady=5, sticky=E)
                
        btn.config(command=lambda:self.menuPick(btn, btnFrame))

    def menuPick(self, btn, frame):
        if self.lastSelectedFrame:
            self.lastSelectedFrame.grid_remove()
        self.lastSelectedFrame = frame
        frame.grid(row=0, column=1, sticky=W+E, padx=5)

    def createGUI(self):
        self.root.title("Xmipp Protocols")
        self.createMainMenu()
        self.lastSelectedFrame = None
        #Create a main frame that contains all other widgets
        self.frame = Frame(self.root)
        self.frame.pack(fill=BOTH)
        self.frame.columnconfigure(0, minsize=150, weight=1)
        self.frame.columnconfigure(1, minsize=150, weight=1)
        self.frame.columnconfigure(2, minsize=300, weight=2)
        
        #Configure toolbar frame
        self.toolbar = Frame(self.frame, bd=1, relief=RAISED)
        self.toolbar.grid(row=0, column=0, sticky=N+W+S)#side=LEFT, fill=Y)
        # Create buttons
        self.ButtonFont = tkFont.Font(family=FontName, size=FontSize, weight=tkFont.BOLD)
        i = 1
        for k, v in sections:
            self.addTbLabel(k, i)
            i += 1
            for btn in v:
                self.addTbButton(i, btn[0], btn[1:])
                i += 1
            
        
        text = Text(self.frame)
        text.grid(row=0, column=2)#pack(side=RIGHT, fill=BOTH)
        self.root.config(menu=self.menubar)
        
       
    def onExit(self):
        self.root.destroy()


if __name__ == '__main__':
    XmippProject()