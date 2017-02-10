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

from math import sqrt
from pyworkflow.em.packages.xmipp3.nma.plotter import plotArray2D


STATE_NO_POINTS = 0 # on points have been selected, double-click will add first one
STATE_DRAW_POINTS = 1 # still adding points, double-click will set the last one
STATE_ADJUST_POINTS = 2 # no more points will be added, just adjust the current ones


class PointPath():
    """ Graphical manager based on Matplotlib to handle mouse
    events to create a path of points. 
    It also allow to modify the point positions on the path.
    """
    def __init__(self, ax, data, pathData, callback=None, tolerance=3, maxPoints=10):
        self.ax = ax
        self.data = data
        plotArray2D(ax, self.data)
        self.callback = callback

        self.dragIndex = None
        self.tolerance = tolerance
        self.maxPoints = maxPoints
        
        self.cidpress = ax.figure.canvas.mpl_connect('button_press_event', self.onClick)
        self.cidrelease = ax.figure.canvas.mpl_connect('button_release_event', self.onRelease)
        self.cidmotion = ax.figure.canvas.mpl_connect('motion_notify_event', self.onMotion)
        self.cidkey = ax.figure.canvas.mpl_connect('key_press_event', self.onKeyPress)
    
        self.pathData = pathData
        
        if pathData.getSize(): # this means there is a path
            self.setState(STATE_ADJUST_POINTS)
            self.plotPath()
        else:
            self.setState(STATE_NO_POINTS)
            self.path_line = None
            self.path_points = None
            
    def setState(self, state, notify=False):
        self.drawing = state
        if state == STATE_NO_POINTS:
            self.ax.set_title('Double click to select first point of trajectory')
        elif state == STATE_DRAW_POINTS:
            self.ax.set_title('Click to add points and double click to end trajectory.')
        elif state == STATE_ADJUST_POINTS:
            self.ax.set_title('Click and drag points to adjust trajectory.')
        else:
            raise Exception("Invalid PointPath state: %d" % state)
        
        if notify and self.callback:
            self.callback()
        
    def onClick(self, event):
        if event.inaxes!=self.ax: 
            return
        # ignore click event if toolbar is active
        if self.ax.figure.canvas.manager.toolbar._active is not None: 
            return
        
        doubleClick = event.dblclick
        ex = event.xdata
        ey = event.ydata
        
        if self.drawing == STATE_NO_POINTS:
            if doubleClick:
                self.setState(STATE_DRAW_POINTS)
                doubleClick = False 
        
        if self.drawing == STATE_DRAW_POINTS:
            if doubleClick:
                self.setState(STATE_ADJUST_POINTS, notify=True)
            else:
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
                self.ax.figure.canvas.draw()
                
                if self.pathData.getSize() == self.maxPoints:
                    self.setState(STATE_ADJUST_POINTS, notify=True)
                    
        
        if self.drawing == STATE_ADJUST_POINTS and not doubleClick: # Points moving state
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
        self.path_line, = self.ax.plot(xs, ys, alpha=0.75, color='green')
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
        
    def onKeyPress(self, event):
        n = self.pathData.getSize()
        
        if event.key == 'ctrl+z' and n > 1:
            self.pathData.removeLastPoint()
            
        if n == 1: # Removed last point, change state appropriately
            self.setState(STATE_NO_POINTS)
            
        self.update()
        
    def update(self):
        xs, ys = self.getXYData()
        self.path_line.set_data(xs, ys)
        self.path_points.set_data(xs, ys)
        self.ax.figure.canvas.draw()
        
        
        