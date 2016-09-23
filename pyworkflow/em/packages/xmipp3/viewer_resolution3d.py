# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Carlos Oscar Sorzano   (coss@cnb.csic.es)
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

import os
import numpy as np
import Tkinter as tk
from math import sqrt

from pyworkflow.utils.properties import Icon
import pyworkflow.gui as gui
from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO
from pyworkflow.protocol.params import LabelParam
import pyworkflow.em as em
from pyworkflow.gui.widgets import Button, HotButton
import pyworkflow.gui.dialog as dialog
from plotter import XmippPlotter

import xmipp

from convert import getImageLocation
from protocol_resolution3d import XmippProtResolution3D


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
        form.addParam('doShowFsc', LabelParam,
                      label="Display Fourier Shell Correlation?")
        form.addParam('doShowDpr', LabelParam,
                      label="Display Differential Phase Residual?")
        form.addParam('doShowStructureFactor', LabelParam,
                      label="Display B-factor?")

    def _getVisualizeDict(self):
        return {'doShowFsc': self._viewFsc,
                'doShowDpr': self._viewDpr,
                'doShowStructureFactor': self._viewStructureFactor,
                }
    
    def _loadPlots(self, title, plotLabel, resolutionLabel, **kwargs):
        """ Check if the FSC metadata is generated and if so, 
        read the plots and the metadata.
        *args and **kwargs will be passed to self._createPlot function.
        """
        fscFn = self.protocol._defineFscName()
        
        if not os.path.exists(fscFn):
            return [self.errorMessage('The FSC metadata was not produced\n'
                                      'Execute again the protocol with FSC\n'
                                      'and/or DPR options enabled.',
                                      title='Missing result file')]
        md = xmipp.MetaData(fscFn)
        plotter = self._createPlot(title, FREQ_LABEL, plotLabel, md, 
                                   xmipp.MDL_RESOLUTION_FREQ, resolutionLabel,
                                   **kwargs)
        return [plotter, em.DataView(fscFn)]
        
    def _viewFsc(self, e=None):
        return self._loadPlots("Fourier Shell Correlation", 'FSC', 
                               xmipp.MDL_RESOLUTION_FRC, color='r')

    def _viewDpr(self, e=None):
        return self._loadPlots("Differential Phase Residual", 'DPR', 
                       xmipp.MDL_RESOLUTION_DPR)
        
    def _createPlot(self, title, xTitle, yTitle, md, mdLabelX, mdLabelY, color='g', figure=None):        
        xplotter = XmippPlotter(figure=figure)
        xplotter.plot_title_fontsize = 11
        xplotter.createSubPlot(title, xTitle, yTitle, 1, 1)
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
            f.close()
            p1 = data.createEmptyPoint()
            p1.setX(values[0])
            p1.setY(values[1])
            data.addPoint(p1)
            p2 = data.createEmptyPoint()
            p2.setX(values[2])
            p2.setY(values[3])
            data.addPoint(p2)
            data.bfactor = values[4]
        else:
            data.bfactor = 0
            
        return data
        
    def _viewStructureFactor(self, e=None):
        strFactFn = self.protocol._defineStructFactorName()
        self.md = xmipp.MetaData(strFactFn)
        win = self.tkWindow(BfactorWindow, 
                           title='Clustering Tool',
                           callback=self._applyBfactor
                            )
        plotter = self._createPlot("Structure Factor", 'frequency^2 (1/A^2)', 
                                   'log(Structure Factor)',
                                   self.md, xmipp.MDL_RESOLUTION_FREQ2, 
                                   xmipp.MDL_RESOLUTION_LOG_STRUCTURE_FACTOR,
                                   figure=win.figure)
        self.path = PointPath(plotter.getLastSubPlot(), self._loadData(), 
                              callback=self._adjustPoints,
                              tolerance=0.1)
        
        return [win]
    
    def _applyBfactor(self, e=None):
        bFactorFile = self.protocol._getPath('bfactor.txt')
        f = open(bFactorFile)
        values = map(float, f.readline().split())
        f.close()
        self.protocol.setStepsExecutor() # set default execution
        vol = self.protocol.inputVolume.get()
        volPath = getImageLocation(vol)
        maxres = 1. / sqrt(values[2])
        args = '-i %s ' % volPath
        pixelSize = vol.getSamplingRate()
        args += '--sampling %f ' % pixelSize
        args += '--maxres %d ' % maxres
        args += '--adhoc %f ' % -values[4]
        volName = os.path.basename(volPath)
        volOut = self.protocol._getPath(volName) 
        args += '-o %s ' % volOut
        self.protocol.runJob('xmipp_volume_correct_bfactor', args)
        
        #args = '-i %s -o %s' % (volPath, self.protocol._getPath(volName))
        #self.protocol.runJob('xmipp_image_convert', args)
        volSet = self.protocol._createSetOfVolumes()
        volSet.setSamplingRate(pixelSize)
        newVol = vol.clone()
        newVol.setObjId(None)
        newVol.setLocation(volOut)
        volSet.append(newVol)
        volSet.write()
        
        self.objectView(volSet).show()
        

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

    
class BfactorWindow(gui.Window):
    """ This class creates a Window that will display Bfactor plot
    to adjust two points to fit B-factor.
    It will also contain a button to apply the B-factor to 
    the volume and produce a new volumen that can be registered.
    """
    def __init__(self, **kwargs):
        gui.Window.__init__(self, **kwargs)
        
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
        
    def _createFigureBox(self, content):
        from pyworkflow.gui.matplotlib_image import FigureFrame
        figFrame = FigureFrame(content, figsize=(6, 6))
        figFrame.grid(row=0, column=0, padx=5, columnspan=2)
        self.figure = figFrame.figure
        
        applyBtn = HotButton(content, text='Apply B-factor',
                           command=self._onApplyBfactorClick)
        applyBtn.grid(row=1, column=0, sticky='ne', padx=5, pady=5)
        
        closeBtn = Button(content, text='Close', imagePath=Icon.ACTION_CLOSE,
                           command=self.close)
        closeBtn.grid(row=1, column=1, sticky='ne', padx=5, pady=5)
       
    def _onApplyBfactorClick(self, e=None):
        #self._runBeforePreWhitening(self.prot)
        dialog.FlashMessage(self.root, "Applying B-factor...", 
                            func=self.callback)
        
    def _onClosing(self):
        if self.plotter:
            self.plotter.close()
        gui.Window._onClosing(self)
        