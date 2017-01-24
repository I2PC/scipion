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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

from os.path import basename
import Tkinter as tk

import pyworkflow.gui as gui
from pyworkflow.utils.properties import Icon
from pyworkflow.gui.widgets import Button, HotButton

from pyworkflow.em.packages.xmipp3.nma.data import Point, PathData
from pyworkflow.em.packages.xmipp3.nma.gui import PointPath
from pyworkflow.em.packages.xmipp3.nma.plotter import XmippNmaPlotter

 
    
class TrajectoriesWindow(gui.Window):
    """ This class creates a Window that will display some Point's
    contained in a Data object.
    It will allow to draw and adjust trajectories along 2D axes.
    """
    def __init__(self, **kwargs):
        gui.Window.__init__(self,  minsize=(420, 200), **kwargs)
        
        self.dim = kwargs.get('dim')
        self.data = kwargs.get('data')
        self.pathData = PathData(dim=self.dim)
        self.callback = kwargs.get('callback', None)
        self.loadCallback = kwargs.get('loadCallback', None)
        self.numberOfPoints = kwargs.get('numberOfPoints', 10)
        
        self.plotter = None
         
        content = tk.Frame(self.root)
        self._createContent(content)
        content.grid(row=0, column=0, sticky='news')
        content.columnconfigure(0, weight=1)
        #content.rowconfigure(1, weight=1)
        
    def _createContent(self, content):
        self._createFigureBox(content)
        self._createTrajectoriesBox(content)
        
    def _addLabel(self, parent, text, r, c):
        label = tk.Label(parent, text=text, font=self.fontBold)
        label.grid(row=r, column=c, padx=5, pady=5, sticky='ne')
        return label
    
    def _createFigureBox(self, content):
        frame = tk.LabelFrame(content, text='Figure')
        frame.columnconfigure(0, minsize=50)
        frame.columnconfigure(1, weight=1)#, minsize=30)
        # Create the 'Axes' label
        self._addLabel(frame, 'Axes', 0, 0)
        
        # Create a listbox with x1, x2 ...
        listbox = tk.Listbox(frame, height=5, 
                             selectmode=tk.MULTIPLE, bg='white')        
        for x in range(1, self.dim+1):
            listbox.insert(tk.END, 'x%d' % x)        
        listbox.grid(row=0, column=1, padx=5, pady=5, sticky='nw')
        self.listbox = listbox
        
        # Selection controls
        self._addLabel(frame, 'Rejection', 1, 0)  
        # Selection label
        self.selectionVar = tk.StringVar()
        self.clusterLabel = tk.Label(frame, textvariable=self.selectionVar)
        self.clusterLabel.grid(row=1, column=1, sticky='nw', padx=5, pady=(10, 5))
        self._updateSelectionLabel()
        # --- Expression
        expressionFrame = tk.Frame(frame)
        expressionFrame.grid(row=2, column=1, sticky='news')
        tk.Label(expressionFrame, text='Expression').grid(row=0, column=0, sticky='ne')
        self.expressionVar = tk.StringVar()
        expressionEntry = tk.Entry(expressionFrame, textvariable=self.expressionVar, 
                                   width=30, bg='white')
        expressionEntry.grid(row=0, column=1, sticky='nw')
        helpText = 'e.g. x1>0 and x1<100 or x3>20'
        tk.Label(expressionFrame, text=helpText).grid(row=1, column=1, sticky='nw')
        
        # Buttons    
        buttonFrame = tk.Frame(frame)
        buttonFrame.grid(row=5, column=1, sticky='sew', pady=(10, 5))
        buttonFrame.columnconfigure(0, weight=1)    
        resetBtn = Button(buttonFrame, text='Reset', command=self._onResetClick)
        resetBtn.grid(row=0, column=0, sticky='ne', padx=(5, 0))
        updateBtn = Button(buttonFrame, text='Update Plot', imagePath='fa-refresh.png',
                           command=self._onUpdateClick)
        updateBtn.grid(row=0, column=1, sticky='ne', padx=5)
       
        frame.grid(row=0, column=0, sticky='new', padx=5, pady=(10, 5))

    def _createTrajectoriesBox(self, content):
        frame = tk.LabelFrame(content, text='Trajectories')
        frame.columnconfigure(0, minsize=50)
        frame.columnconfigure(1, weight=1)#, minsize=30)

        # Animation name
        self._addLabel(frame, 'Name', 0, 0)
        self.animationVar = tk.StringVar()
        clusterEntry = tk.Entry(frame, textvariable=self.animationVar, 
                                   width=30, bg='white')
        clusterEntry.grid(row=0, column=1, sticky='nw', pady=5)
        
        buttonsFrame = tk.Frame(frame)
        buttonsFrame.grid(row=1, column=1, 
                          sticky='se', padx=5, pady=5)
        buttonsFrame.columnconfigure(0, weight=1)

        self.generateBtn = HotButton(buttonsFrame, text='Generate Animation', state=tk.DISABLED,
                              tooltip='Select trajectory points to generate the animations',
                              imagePath='fa-plus-circle.png', command=self._onCreateClick)
        self.generateBtn.grid(row=0, column=1, padx=5)  
        
        self.loadBtn = Button(buttonsFrame, text='Load', imagePath='fa-folder-open.png',
                              tooltip='Load a generated animation.',command=self._onLoadClick)
        self.loadBtn.grid(row=0, column=2, padx=5)   
                  
        self.closeBtn = Button(buttonsFrame, text='Close', imagePath=Icon.ACTION_CLOSE,
                              tooltip='Close window', command=self.close)
        self.closeBtn.grid(row=0, column=3, padx=(5, 10)) 
               
        frame.grid(row=1, column=0, sticky='new', padx=5, pady=(5, 10))
        
    def _onResetClick(self, e=None):
        """ Clean the expression and the current selection. """
        self.expressionVar.set('')
        self.pathData.clear()
        for point in self.data.iterAll():
            point.setState(Point.NORMAL)
        self._onUpdateClick()
        self.generateBtn.config(state=tk.DISABLED)
        
    def _onCreateClick(self, e=None):
        if self.callback:
            self.callback()
        
    def _onLoadClick(self, e=None):
        if self.loadCallback:
            self.loadCallback()
        
    def setPathData(self, data):
        self.pathData = data
        
    def _evalExpression(self):
        """ Evaluate the input expression and add 
        matching points to the selection.
        """
        value = self.expressionVar.get().strip()
        if value:
            for point in self.data:
                if point.eval(value):
                    point.setState(Point.DISCARDED)
                   
    def setDataIndex(self, indexName, value):
        """ Set which point data index will be used as X, Y or Z. """
        setattr(self.data, indexName, value)
        setattr(self.pathData, indexName, value)
                             
    def _onUpdateClick(self, e=None):
        components = self.listbox.curselection()
        dim = len(components)
           
        if not dim:
            self.showWarning("Please select some Axis before update plots.")
        else: 
            modeList = components
            modeNameList = ['x%d' % (m+1) for m in components]
            missingList = []
                    
            if missingList:
                return [self.errorMessage("Invalid mode(s) *%s*\n." % (', '.join(missingList)), 
                              title="Invalid input")]
            
            if self.plotter is None or self.plotter.isClosed():
                self.plotter = XmippNmaPlotter(data=self.data)
                doShow = True
                #self.plotter.useLastPlot = True
            else:
                self.plotter.clear()
                doShow = False
            
            # Actually plot
            baseList = [basename(n) for n in modeNameList]
            
            self.setDataIndex('XIND', modeList[0])
            self.ps = None
            
            if dim == 1:
                self.plotter.plotArray1D("Histogram for %s" % baseList[0],  
                                    "Deformation value", "Number of images")
            else:
                self.setDataIndex('YIND', modeList[1])
                if dim == 2:
                    self._evalExpression()
                    self._updateSelectionLabel()
                    ax = self.plotter.createSubPlot("Click and drag to add points to the Cluster",
                                                    *baseList)
                    self.ps = PointPath(ax, self.data, self.pathData, 
                                        callback=self._checkNumberOfPoints)
                elif dim == 3:
                    #del self.ps # Remove PointSelector
                    self.setDataIndex('ZIND', modeList[2])
                    self.plotter.plotArray3D("%s %s %s" % tuple(baseList), *baseList)

            if doShow:
                self.plotter.show()
            else:
                self.plotter.draw()

    def _updateSelectionLabel(self):
        self.selectionVar.set('%d / %d points' % (self.data.getDiscardedSize(),
                                                  self.data.getSize()))
        
    def _checkNumberOfPoints(self):
        """ Check that if the number of points was selected
        and add new ones if needed.
        """
        while (self.pathData.getSize() < self.numberOfPoints):
            self.pathData.splitLongestSegment()
        self._onUpdateClick()
        self.generateBtn.config(state=tk.NORMAL)
        
    def getAnimationName(self):
        return self.animationVar.get().strip()
    
    def setAnimationName(self, value):
        self.animationVar.set(value)
        
    def _onClosing(self):
        if self.plotter:
            self.plotter.close()
        gui.Window._onClosing(self)
        
