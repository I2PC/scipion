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

from math import sqrt
from pyworkflow.em.packages.xmipp3.nma.plotter import plotArray2D


        
class PointPath():
    """ Graphical manager based on Matplotlib to handle mouse
    events to create a path of points. 
    It also allow to modify the point positions on the path.
    """
    def __init__(self, ax, data, callback=None):
        self.ax = ax
        self.data = data
        plotArray2D(ax, self.data)
        self.path_line, = ax.plot([0], [0])  # empty line
        self.path_points = None
        self.drawing = 0
        self.dragIndex = None
        self.xs = []
        self.ys = []
        self.cidpress = ax.figure.canvas.mpl_connect('button_press_event', self.onClick)
        self.cidrelease = ax.figure.canvas.mpl_connect('button_release_event', self.onRelease)
        self.cidmotion = ax.figure.canvas.mpl_connect('motion_notify_event', self.onMotion)
        self.cidkey = ax.figure.canvas.mpl_connect('key_press_event', self.onKeyPress)
    
        self.ax.set_title('Double click to select first point of trajectory')
        
    def onClick(self, event):
        if event.inaxes!=self.path_line.axes: return
        doubleClick = event.dblclick
        ex = event.xdata
        ey = event.ydata
        if self.drawing == 0:
            if doubleClick:
                self.drawing += 1 # Pass to drawing points state
                doubleClick = False 
                self.ax.set_title('Click to add points and double click to end trajectory.')
        if self.drawing == 1:
            if doubleClick:
                self.drawing += 1 # Stop drawing points
                self.ax.set_title('Click and drag points to adjust trajectory.')
            else:
                self.xs.append(ex)
                self.ys.append(ey)
                if len(self.xs) == 1: # first point is added
                    self.path_line, = self.ax.plot(self.xs, self.ys, alpha=0.75)
                    self.path_points, = self.ax.plot(self.xs, self.ys, 'o', color='red')  # 5 points tolerance, mark line points
                else:
                    self.path_line.set_data(self.xs, self.ys)
                    self.path_points.set_data(self.xs, self.ys)
                self.ax.figure.canvas.draw()
        
        if self.drawing == 2 and not doubleClick: # Points moving state
            self.dragIndex = None
            for i, (x, y) in enumerate(zip(self.xs, self.ys)):
                if sqrt((ex - x)**2 + (ey - y)**2) < 0.1:
                    self.dragIndex = i
                    break
          
    def onMotion(self, event):
        if self.dragIndex is None or self.drawing < 2: return
        ex, ey = event.xdata, event.ydata
        self.xs[self.dragIndex] = ex
        self.ys[self.dragIndex] = ey
        self.update()
        
    def onRelease(self, event):
        self.dragIndex = None
        
    def onKeyPress(self, event):
        if event.key == 'ctrl+z':
            del self.xs[-1]
            del self.ys[-1]
        self.update()
        
    def update(self):
        self.path_line.set_data(self.xs, self.ys)
        self.path_points.set_data(self.xs, self.ys)
        self.ax.figure.canvas.draw()
        
        
        