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
This module implement the wrappers aroung Xmipp CL2D protocol
visualization program.
"""

import os
from os.path import basename
import numpy as np
import Tkinter as tk

import pyworkflow.gui as gui
from pyworkflow.gui.widgets import Button, HotButton
from pyworkflow.utils.path import cleanPath
from pyworkflow.viewer import (ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO)
from pyworkflow.protocol.params import StringParam, BooleanParam
from pyworkflow.em.data import SetOfParticles
from protocol_nma_dimred import XmippProtDimredNMA
from data import Point, Data
import xmipp

from plotter import XmippNmaPlotter, plotArray2D


        
class XmippDimredNMAViewer(ProtocolViewer):
    """ Visualization of results from the NMA protocol
    """
    _label = 'viewer nma dimred'
    _targets = [XmippProtDimredNMA]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    
    def __init__(self, **kwargs):
        ProtocolViewer.__init__(self, **kwargs)
        self.data = self.loadData()
        
    def _defineParams(self, form):
        form.addSection(label='Visualization')
        form.addParam('displayRawDeformation', StringParam, default='1',
                      label='Display raw deformation',
                      help='Type 1 to see the histogram of raw deformation number 1; \n'
                           'type 2 to see the histogram of raw deformation number 2, etc.\n'
                           'Type 1 2 to see the 2D plot of raw deformations number 1 vs 2.\n'
                           'Type 1 2 3 to see the 3D plot of raw deformations 1, 2 and 3; etc.'
                           )
        
        form.addParam('displayClustering', BooleanParam, default=False,
                      label='Open clustering tool?',
                      help='Open a GUI to visualize the images as points'
                           'and select some of them to create new clusters.')
        
        
    def _getVisualizeDict(self):
        return {'displayRawDeformation': self._viewRawDeformation,
                'displayClustering': self._displayClustering,
                } 
                        
    def _viewRawDeformation(self, paramName):
        components = self.displayRawDeformation.get()
        return self._doViewRawDeformation(components)
        
    def _doViewRawDeformation(self, components):
#        components = map(int, self.displayRawDeformation.get().split())
        components = map(int, components.split())
        dim = len(components)
        views = []
        
        if dim > 0:
            modeList = [m - 1 for m in components]
            modeNameList = ['Mode %d' % m for m in components]
            missingList = []
                    
            if missingList:
                return [self.errorMessage("Invalid mode(s) *%s*\n." % (', '.join(missingList)), 
                              title="Invalid input")]
            
            # Actually plot
            plotter = XmippNmaPlotter(data=self.data) 
            baseList = [basename(n) for n in modeNameList]
            
            if dim == 1:
                Point.XIND = modeList[0]
                plotter.plotArray1D("Histogram for %s" % baseList[0], 
                                    "Deformation value", "Number of images")
            else:
                Point.YIND = modeList[1]
                if dim == 2:
                    plotter.plotArray2D("%s vs %s" % tuple(baseList), *baseList)
                elif dim == 3:
                    Point.ZIND = modeList[2]
                    plotter.plotArray3D("%s %s %s" % tuple(baseList), *baseList)
            views.append(plotter)
            
        return views
    
    def _displayClustering(self, paramName):
        return [self.tkWindow(ClusteringWindow, 
                              dim=self.protocol.reducedDim.get(),
                              data=self.data,
                              callback=self._createCluster
                              )]
        
    def _createCluster(self):
        """ Create the cluster with the selected particles
        from the cluster. This method will be called when
        the button 'Create Cluster' is pressed.
        """
        # Write the particles
        prot = self.protocol
        project = prot.getProject()
        inputSet = prot.getInputParticles()
        fnSqlite = prot._getTmpPath('cluster_particles.sqlite')
        cleanPath(fnSqlite)
        partSet = SetOfParticles(filename=fnSqlite)
        partSet.copyInfo(inputSet)
        for point in self.data:
            if point.getState() == Point.SELECTED:
                particle = inputSet[point.getId()]
                partSet.append(particle)
        partSet.write()
        partSet.close()
                
        from protocol_batch_cluster import BatchProtNMACluster
        newProt = project.newProtocol(BatchProtNMACluster)
        newProt.inputNmaDimred.set(prot)
        newProt.sqliteFile.set(fnSqlite)
        
        project.launchProtocol(newProt)
        
    def loadData(self):
        """ Iterate over the images and the output matrix txt file
        and create a Data object with theirs Points.
        """
        matrix = np.loadtxt(self.protocol.getOutputMatrixFile())
        particles = self.protocol.getInputParticles()
        
        data = Data()
        for i, particle in enumerate(particles):
            data.addPoint(Point(pointId=particle.getObjId(),
                                data=matrix[i, :],
                                weight=particle._xmipp_cost.get()))
            
        return data
    
    
class ClusteringWindow(gui.Window):
    def __init__(self, **kwargs):
        gui.Window.__init__(self,  minsize=(420, 200), **kwargs)
        
        self.dim = kwargs.get('dim')
        self.data = kwargs.get('data')
        self.callback = kwargs.get('callback', None)
        self.plotter = None
         
        content = tk.Frame(self.root)
        self._createContent(content)
        content.grid(row=0, column=0, sticky='news')
        content.columnconfigure(0, weight=1)
        #content.rowconfigure(1, weight=1)
        
    def _createContent(self, content):
        self._createFigureBox(content)
        self._createClusteringBox(content)
        
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
        self._addLabel(frame, 'Selection', 1, 0)  
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

    def _createClusteringBox(self, content):
        frame = tk.LabelFrame(content, text='Cluster')
        frame.columnconfigure(0, minsize=50)
        frame.columnconfigure(1, weight=1)#, minsize=30)

        # Cluster line
        self._addLabel(frame, 'Cluster name', 0, 0)
        self.clusterVar = tk.StringVar()
        clusterEntry = tk.Entry(frame, textvariable=self.clusterVar, 
                                   width=30, bg='white')
        clusterEntry.grid(row=0, column=1, sticky='nw', pady=5)
        
        buttonsFrame = tk.Frame(frame, bg='green')
        buttonsFrame.grid(row=1, column=1, 
                          sticky='se', padx=5, pady=5)
        buttonsFrame.columnconfigure(0, weight=1)

        createBtn = HotButton(buttonsFrame, text='Create Cluster', 
                              imagePath='fa-plus-circle.png', command=self._onCreateClick)
        createBtn.grid(row=0, column=1)       
       
        frame.grid(row=1, column=0, sticky='new', padx=5, pady=(5, 10))
        
    def _onResetClick(self, e=None):
        """ Clean the expression and the current selection. """
        self.expressionVar.set('')
        for point in self.data:
            point.setState(Point.NORMAL)
        self._onUpdateClick()
        
    def _onCreateClick(self, e=None):
        if self.callback:
            self.callback()
        
    def _evalExpression(self):
        """ Evaluate the input expression and add 
        matching points to the selection.
        """
        value = self.expressionVar.get().strip()
        if value:
            for point in self.data:
                if point.eval(value):
                    point.setState(Point.SELECTED)
                                        
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
                #self.plotter.useLastPlot = True
            else:
                self.plotter.clear()
            
            # Actually plot
            baseList = [basename(n) for n in modeNameList]
            
            Point.XIND = modeList[0]
            
            if dim == 1:
                self.plotter.plotArray1D("Histogram for %s" % baseList[0],  
                                    "Deformation value", "Number of images")
            else:
                Point.YIND = modeList[1]
                if dim == 2:
                    self._evalExpression()
                    self._updateSelectionLabel()
                    ax = self.plotter.createSubPlot("Click and drag to add points to the Cluster",
                                                    *baseList)
                    self.ps = PointSelector(ax, self.data, callback=self._updateSelectionLabel)
                elif dim == 3:
                    del self.ps # Remove PointSelector
                    Point.ZIND = modeList[2]
                    self.plotter.plotArray3D("%s %s %s" % tuple(baseList), *baseList)

            self.plotter.draw()

    def _updateSelectionLabel(self):
        self.selectionVar.set('%d / %d points' % (self.data.getSelectedSize(),
                                                  self.data.getSize()))
        
    def _onClosing(self):
        if self.plotter:
            self.plotter.close()
        gui.Window._onClosing(self)
        
        
class PointSelector():
    """ Graphical manager base on Matplotlib to handle mouse
    events of click, drag and release and mark some point
    from input Data as 'selected'.
    """
    def __init__(self, ax, data, callback=None):
        self.ax = ax
        self.data = data
        self.createPlots(ax)
        self.press = None
        self.callback = callback
        # connect to all the events we need 
        self.cidpress = self.rectangle_selection.figure.canvas.mpl_connect(
            'button_press_event', self.onPress)
        self.cidrelease = self.rectangle_selection.figure.canvas.mpl_connect(
            'button_release_event', self.onRelease)
        self.cidmotion = self.rectangle_selection.figure.canvas.mpl_connect(
            'motion_notify_event', self.onMotion)
    
    def createPlots(self, ax):
        plotArray2D(ax, self.data)
        self.createSelectionPlot(ax)
    
    def getSelectedData(self):
        xs, ys = [], []
        for point in self.data:
            if point.getState() == 1: # point selected
                xs.append(point.getX())
                ys.append(point.getY())
        return xs, ys
        
    def createSelectionPlot(self, ax):
        xs, ys = self.getSelectedData()
        self.plot_selected, = ax.plot(xs, ys, 'o', ms=8, alpha=0.4,
                                  color='yellow')
         
        self.rectangle_selection, = ax.plot([0], [0], ':')  # empty line      
        
    def onPress(self, event):
        if event.inaxes != self.rectangle_selection.axes: 
            return
        # ignore click event if toolbar is active
        if self.ax.figure.canvas.manager.toolbar._active is not None: 
            return
        
        self.press = True
        self.originX = event.xdata
        self.originY = event.ydata
        
    def onMotion(self, event):
        if self.press is None: return
        self.update(event, addSelected=False)
        
    def onRelease(self, event):
        self.press = None
        self.update(event, addSelected=True)
        
    def inside(self, x, y, xmin, xmax, ymin, ymax):
        return (x >= xmin and x <= xmax and
                y >= ymin and y <= ymax)
        
    def update(self, event, addSelected=False):
        """ Update the plots with selected points.
        Take the selection rectangle and include points from 'event'.
        If addSelected is True, update the data selection.
        """
        # ignore click event if toolbar is active
        if self.ax.figure.canvas.manager.toolbar._active is not None: 
            return
        ox, oy = self.originX, self.originY
        ex, ey = event.xdata, event.ydata
        xs1, ys1 = self.getSelectedData()
        
        x1 = min(ox, ex)
        x2 = max(ox, ex)
        y1 = min(oy, ey)
        y2 = max(oy, ey)
        
        if addSelected:
            xs, ys = [0], [0] # rectangle selection
        else:
            xs = [x1, x2, x2, x1, x1]
            ys = [y1, y1, y2, y2, y1]
            
        for point in self.data:
            x, y = point.getX(), point.getY()
            if self.inside(x, y, x1, x2, y1, y2):
                xs1.append(x)
                ys1.append(y)
                if addSelected:
                    point.setState(Point.SELECTED)
        
        if self.callback: # Notify changes on selection
            self.callback()
        self.plot_selected.set_data(xs1, ys1)        
        self.rectangle_selection.set_data(xs, ys)
        
        self.ax.figure.canvas.draw()
        
        
