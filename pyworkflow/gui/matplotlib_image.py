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
Several Image tools using Matplotlib.
"""

import Tkinter as tk
import matplotlib
matplotlib.use('TkAgg')
import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.cm as cm
from matplotlib.patches import Wedge


class Preview(tk.Frame):
    #def __init__(self, parent, dim, dpi=36, label=None):
    def __init__(self, parent, dim, dpi=36, label=None, col=0, row=0):
        tk.Frame.__init__(self, parent)
        self.dim = dim
        self.bg = np.zeros((dim, dim), float)
        ddim = dim/dpi
        self.figure = Figure(figsize=(ddim, ddim), dpi=dpi, frameon=False)
        self.canvas = FigureCanvasTkAgg(self.figure, master=self)
        self.canvas.get_tk_widget().grid(column=0, row=0)#, sticky=(N, W, E, S))
        if label:
            tk.Label(self, text=label).grid(column=0, row=1)
        self._createAxes()
        
    def setWindowTitle(self,title):
        """ Set window title"""
        self.canvas.set_window_title(title)

    def _createAxes(self):
        """ Should be implemented in subclasses. """
        pass
    
    def _update(self):
        """ Should be implemented in subclasses. """
    
    def clear(self):
        self._update(self.bg)
        
    def updateData(self, Z):
        self.clear()
        self._update(Z)
    
    
class ImagePreview(Preview):
    def __init__(self, parent, dim, dpi=36, label=None, col=0):
        Preview.__init__(self, parent, dim, dpi, label, col)
            
    def _createAxes(self):
        ax = self.figure.add_axes([0,0,1,1], frameon=False)       
        self.figureimg = ax.imshow(self.bg, cmap=cm.gray)#, extent=[-h, h, -h, h])
        ax.set_axis_off()
        self.ax = ax
        
    def _update(self, Z, *args):
        self.figureimg.set_data(Z)
        self.figureimg.autoscale()
        self.figureimg.set(extent=[0, Z.shape[1], 0, Z.shape[0]])
        self.canvas.draw()
        
        
class PsdPreview(Preview):
    def __init__(self, master, dim, lf, hf, dpi=72, label="PSD", Z=None):
        Preview.__init__(self, master, dim, dpi, label)
        self.lf = lf
        self.hf = hf
        if self.ring:
            self.createRing()
        else:
            self.canvas.draw()
                            
    def _createAxes(self):
        #axdef = SubplotZero(self.figure, 111)
        #ax = self.figure.add_subplot(axdef)
        ax = self.figure.add_axes([0.1,0.1,0.8,0.8], frameon=False)
        #ax.xaxis.set_offset_position(0.5)       
        #ax.set_position([0, 0, 1, 1])
        h = 0.5
        ax.set_xlim(-h, h)
        ax.set_ylim(-h, h)
        ax.grid(True)
        self.ring = None
        self.img = ax.imshow(self.bg, cmap=cm.gray, extent=[-h, h, -h, h])
        self.ax = ax
        
    def createRing(self):
        radius = float(self.hf)
        width = radius - float(self.lf)
        self.ring = Wedge((0,0), radius, 0, 360, width=width, alpha=0.15) # Full ring
        self.ax.add_patch(self.ring)
        self.canvas.draw()
        
    def updateFreq(self, lf, hf):
        self.lf = lf
        self.hf = hf
        if self.ring:
            self.ring.remove()
            self.ring = None
        if self.hf:
            self.createRing()
    
    def _update(self, Z):
        if self.ring:
            self.ring.remove()
            self.ring = None
        if self.hf:
            self.createRing()
        self.img.set_data(Z)
        self.img.autoscale()
        self.canvas.draw()
        
        
class MaskPreview(ImagePreview):
    def __init__(self, parent, dim, dpi=36, label=None, col=0, outerRadius=None, innerRadius=0):
        ImagePreview.__init__(self, parent, dim, dpi, label, col)
        if outerRadius is None:
            outerRadius = dim / 2
        self.ring = None
        self.updateMask(outerRadius, innerRadius)
            
    def updateMask(self, outerRadius, innerRadius=0):
        if self.ring is not None:
            self.ring.remove()
        center = self.dim / 2
        width = outerRadius - innerRadius
        self.ring = Wedge((center, center), outerRadius, 0, 360, width=width, alpha=0.15) # Full ring
        self.ax.add_patch(self.ring)
        self.canvas.draw()
        
        
def getPngData(filename):  
    import matplotlib.image as mpimg
    return mpimg.imread(filename)

def createBgImage(dim):
    return np.ones((dim, dim, 3))