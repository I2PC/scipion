# **************************************************************************
# *
# * Authors:     Josue Gomez Blanco (jgomez@cnb.csic.es)
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

import pyworkflow.gui as gui
from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO
from pyworkflow.em import *
from pyworkflow.gui.form import FormWindow
from protocol_resolution3d import *
from plotter import XmippPlotter
from xmipp import *
import numpy as np


FREQ_LABEL = 'frequency (1/A)'


class XmippResolution3DViewer(ProtocolViewer):
    """ Wrapper to visualize different type of data objects
    with the Xmipp program xmipp_showj
    """
    _label = 'viewer resolution3D'
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [XmippProtResolution3D]
    
    def _defineParams(self, form):
        form.addSection(label='Results')
        form.addParam('doShowFsc', BooleanParam, default=True, 
                      label="Display Fourier Shell Correlation?")
        form.addParam('doShowDpr', BooleanParam, default=True, 
                      label="Display Differential Phase Residual?")
        form.addParam('doShowStructureFactor', BooleanParam, default=True, 
                      label="Display Structure factor?")

    def _getVisualizeDict(self):
        return {'doShowFsc': self._viewFsc,
                'doShowDpr': self._viewDpr,
                'doShowStructureFactor': self._viewStructureFactor,
                }
    
    def _viewFsc(self, e=None):
        fscFn = self.protocol._defineFscName()
        md = MetaData(fscFn)
        plotter = self._createPlot("Fourier Shell Correlation", FREQ_LABEL, 'FSC', 
                               md, MDL_RESOLUTION_FREQ, MDL_RESOLUTION_FRC, color='r')
        return [plotter, DataView(fscFn)]
        
    def _viewDpr(self, e=None):
        fscFn = self.protocol._defineFscName()
        md = MetaData(fscFn)
        return [self._createPlot("Differential Phase Residual", FREQ_LABEL, 'DPR', 
                               md, MDL_RESOLUTION_FREQ, MDL_RESOLUTION_DPR),
                DataView(fscFn)]
    
    def _createPlot(self, title, xTitle, yTitle, md, mdLabelX, mdLabelY, color='g'):        
        xplotter = XmippPlotter(1, 1, figsize=(4,4), windowTitle="Plot")
        xplotter.createSubPlot(title, xTitle, yTitle)
        xplotter.plotMdFile(md, mdLabelX, mdLabelY, color)
        return xplotter

    def _adjustPoints(self, data):
        x=data.getXData()
        if x[0]>x[1]:
            aux=x[1]
            x[1]=x[0]
            x[0]=aux
        X=[]
        Y=[]
        for objId in self.md:
            f=self.md.getValue(xmipp.MDL_RESOLUTION_FREQ2,objId)
            if f>=x[0] and f<=x[1]:
                X.append(f)
                Y.append(self.md.getValue(xmipp.MDL_RESOLUTION_LOG_STRUCTURE_FACTOR,objId))
        X=np.array(X)
        Y=np.array(Y)
        A=np.array([np.ones(X.size), X.T])
        beta=np.linalg.lstsq(A.T,Y)[0]
        y = [beta[0]+beta[1]*xi for xi in x]
        Bfactor = -4*beta[1]

        # Update data
        data.getPoint(0).setX(x[0])
        data.getPoint(0).setY(y[0])
        data.getPoint(1).setX(x[1])
        data.getPoint(1).setY(y[1])
        data.bfactor=Bfactor
            
        f = open(self.protocol._getPath('bfactor.txt'), 'w')
        print >> f, x[0], y[0], x[1], y[1], Bfactor
        f.close()
            
    def _loadData(self):
        from pyworkflow.em.packages.xmipp3.nma.data import PathData
        data = PathData(dim=2)
        bfactorFile = self.protocol._getPath('bfactor.txt')
        if os.path.exists(bfactorFile):
            f = open(bfactorFile)
            values = map(float, f.readline().split())
            p1 = data.createEmptyPoint()
            p1.setX(values[0])
            p1.setY(values[1])
            data.addPoint(p1)
            p2 = data.createEmptyPoint()
            p2.setX(values[2])
            p2.setY(values[3])
            data.addPoint(p2)
            data.bfactor=values[4]
        else:
            data.bfactor=0
            
        return data
        
    def _viewStructureFactor(self, e=None):
        strFactFn = self.protocol._defineStructFactorName()
        self.md = MetaData(strFactFn)
        plotter = self._createPlot("Structure Factor", 'frequency^2 (1/A^2)', 'log(Structure Factor)', 
                               self.md, xmipp.MDL_RESOLUTION_FREQ2, xmipp.MDL_RESOLUTION_LOG_STRUCTURE_FACTOR)
        self.path = PointPath(plotter.getLastSubPlot(), self._loadData(), 
                              callback=self._adjustPoints,
                              tolerance=0.1)
        return [plotter]        


from math import sqrt


STATE_NO_POINTS = 0 # on points have been selected, double-click will add first one
STATE_DRAW_POINTS = 1 # still adding points, double-click will set the last one
STATE_ADJUST_POINTS = 2 # no more points will be added, just adjust the current ones


class PointPath():
    """ Graphical manager based on Matplotlib to handle mouse
    events to create a path of points. 
    It also allow to modify the point positions on the path.
    """
    def __init__(self, ax, pathData, callback=None, tolerance=3, maxPoints=2):
        self.ax = ax
        self.callback = callback

        self.dragIndex = None
        self.tolerance = tolerance
        self.maxPoints = maxPoints
        
        self.cidpress = ax.figure.canvas.mpl_connect('button_press_event', self.onClick)
        self.cidrelease = ax.figure.canvas.mpl_connect('button_release_event', self.onRelease)
        self.cidmotion = ax.figure.canvas.mpl_connect('motion_notify_event', self.onMotion)
    
        self.pathData = pathData
        
        if pathData.getSize() == maxPoints: # this means there is a path
            self.setState(STATE_ADJUST_POINTS)
            self.plotPath()
        else:
            self.setState(STATE_DRAW_POINTS)
            self.path_line = None
            self.path_points = None
            
    def setState(self, state, notify=False):
        self.drawing = state
        
        if state == STATE_DRAW_POINTS:
            self.ax.set_title('Click to add two points.')
        elif state == STATE_ADJUST_POINTS:
            self.ax.set_title('Drag points to adjust line, current Bfactor = %0.3f' % self.pathData.bfactor)
        else:
            raise Exception("Invalid PointPath state: %d" % state)
        
        if notify and self.callback:
            self.callback(self.pathData)
        
    def onClick(self, event):
        if event.inaxes!=self.ax: 
            return
        # ignore click event if toolbar is active
        if self.ax.figure.canvas.manager.toolbar._active is not None: 
            return
        ex = event.xdata
        ey = event.ydata
        
        if self.drawing == STATE_DRAW_POINTS:
            point = self.pathData.createEmptyPoint()
            point.setX(ex)
            point.setY(ey)
            self.pathData.addPoint(point)
            
            if self.pathData.getSize() == 1: # first point is added
                self.plotPath()
            else:
                xs, ys = self.getXYData()
                self.path_line.set_data(xs, ys)
                self.path_points.set_data(xs, ys)
            
            if self.pathData.getSize() == self.maxPoints:
                self.setState(STATE_ADJUST_POINTS)
                
            self.ax.figure.canvas.draw()
        
        elif self.drawing == STATE_ADJUST_POINTS: # Points moving state
            self.dragIndex = None
            for i, point in enumerate(self.pathData):
                x = point.getX()
                y = point.getY()
                if sqrt((ex - x)**2 + (ey - y)**2) < self.tolerance:
                    self.dragIndex = i
                    break
          
    def getXYData(self):
        xs = self.pathData.getXData()
        ys = self.pathData.getYData()
        return xs, ys
    
    def plotPath(self):
        xs, ys = self.getXYData()
        self.path_line, = self.ax.plot(xs, ys, alpha=0.75, color='blue')
        self.path_points, = self.ax.plot(xs, ys, 'o', color='red')  # 5 points tolerance, mark line points
        
    def onMotion(self, event):
        if self.dragIndex is None or self.drawing < 2: 
            return
        # ignore click event if toolbar is active
        if self.ax.figure.canvas.manager.toolbar._active is not None: 
            return
        
        ex, ey = event.xdata, event.ydata
        point = self.pathData.getPoint(self.dragIndex)
        point.setX(ex)
        point.setY(ey)
        self.update()
        
    def onRelease(self, event):
        self.dragIndex = None
        if self.drawing == STATE_ADJUST_POINTS:
            self.setState(STATE_ADJUST_POINTS, notify=True)
        self.update()    
        
    def update(self):
        xs, ys = self.getXYData()
        self.path_line.set_data(xs, ys)
        self.path_points.set_data(xs, ys)
        self.ax.figure.canvas.draw()
        
