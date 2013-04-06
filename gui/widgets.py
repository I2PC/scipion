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
Some basic GUI widgets are implemented in this module.
The widgets here are suppose to be used to build more complex 
elements.
"""

import Tkinter as tk
import ttk
import gui

class Button(tk.Button):
    def __init__(self, master, text, imagePath=None, tooltip=None, **opts):
        defaults = {'activebackground': gui.cfgButtonActiveBgColor, 'bg': gui.cfgButtonBgColor,
                    'fg': gui.cfgButtonFgColor, 'activeforeground': gui.cfgButtonActiveFgColor,
                    'font': gui.fontButton}
        defaults.update(opts)
        btnImage = gui.getImage(imagePath)
        
        if btnImage:
            #height=28, width=28,
            if not opts.has_key('bg'):
                del defaults['bg']
            tk.Button.__init__(self, master, image=btnImage, bd=0,  **defaults)
            self.image = btnImage
        else:
            tk.Button.__init__(self, master, text=text, **defaults)
            
        if tooltip:
            from tooltip import ToolTip
            ToolTip(self, tooltip, 500)
            
    def setImage(self, imagePath):
        self.image = gui.getImage(imagePath)
        self.config(image=self.image)

class AutoScrollbar(tk.Scrollbar):
    """"A scrollbar that hides itself if it's not needed."""
    
    def set(self, lo, hi):
        if float(lo) <= 0.0 and float(hi) >= 1.0:
            self.grid_remove()
            #self.tk.call("grid", "remove", self)
        else:
            self.grid()
        tk.Scrollbar.set(self, lo, hi)
        

class Scrollable(object):
    """This is a base class for all scrollable widgets.
    If it is enabled, it will wrap the widget with a frame
    and will add vertical and horizontal AutoScrollbar"""
    
    def __init__(self, master, WidgetClass, frame=True, **opts):
        if frame:
            self.frame = tk.Frame(master)
            self.frame.rowconfigure(0, weight=1)
            self.frame.columnconfigure(0, weight=1)
            self.vscroll = AutoScrollbar(self.frame)
            self.vscroll.grid(row=0, column=1, sticky='ns')
            self.hscroll = AutoScrollbar(self.frame, orient=tk.HORIZONTAL)
            self.hscroll.grid(row=1, column=0, sticky='ew')
            WidgetClass.__init__(self, self.frame, 
                                 yscrollcommand=self.vscroll.set,
                                 xscrollcommand=self.hscroll.set, **opts)
            self.vscroll.config(command=self.yview)
            self.hscroll.config(command=self.xview)
            self.grid(row=0, column=0, sticky='news')
            self.grid = self.frame.grid
            self.grid_remove = self.frame.grid_remove
        else:
            WidgetClass.__init__(self, master, **opts)
        
        
class Tree(ttk.Treeview, Scrollable):
    """ This widget acts as a wrapper around the ttk.Treeview"""
    
    def __init__(self, master, frame=True, **opts):
        """Create a new Tree, if frame=True, a container
        frame will be created and an scrollbar will be added"""
        Scrollable.__init__(self, master, ttk.Treeview, frame, **opts)
        
    def getFirst(self):
        ''' Return first selected item or None if selection empty'''
        selection = self.selection()
        if len(selection):
            return selection[0]
        return None
    
    def _moveSelection(self, moveFunc):
        item = self.selection_first()
        if item:
            item = moveFunc(item)
            if item != '':
                self.selection_set(item)
        
    def moveSelectionUp(self, e=None):
        ''' change selection to previous item '''
        self._moveSelection(self.prev)
    
    def moveSelectionDown(self, e=None):
        ''' change selection to to next item '''
        self._moveSelection(self.next)
        
    def moveItemUp(self, e=None):
        '''if selected item is not the first move up one position'''
        item = self.selection_first()
        if item:
            index = self.index(item)
            if index > 0:
                self.move(item, '', index-1)
                
    def moveItemDown(self, e=None):
        '''if selected item is not the first move up one position'''
        item = self.selection_first()
        if item:
            index = self.index(item)
            if self.next(item) != '':
                self.move(item, '', index+1)
                
    def clear(self):
        ''' remove all items '''
        childs = self.get_children('')
        for c in childs:
            self.delete(c)
            
class LabelSlider(ttk.Frame):
    """ Create a personalized frame that contains label, slider and label value
        it also keeps a variable with the value """

    def __init__(self, master, label, from_=0, to=100, value=50, callback=None, step=0.01):
        self.var = tk.DoubleVar()
        self.var.set(float(value))
        ttk.Frame.__init__(self, master)
        ttk.Label(self, text=label).pack(side=tk.LEFT, padx=2, pady=2, anchor='s')
        self.slider = tk.Scale(self, from_=from_, to=to, variable=self.var, 
                                bigincrement=step, resolution=step, orient=tk.HORIZONTAL)
        if callback:
            self.var.trace('w', callback)
        self.slider.pack(side=tk.LEFT, padx=2)
        
    def getValue(self):
        return self.var.get()