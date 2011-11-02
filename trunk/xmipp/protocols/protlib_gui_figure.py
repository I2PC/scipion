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

import Tkinter as tk
import ttk
#import matplotlib
#matplotlib.use('TkAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.cm as cm
from mpl_toolkits.axes_grid.axislines import SubplotZero
from matplotlib.patches import Wedge

# Some utilities to use matplotlib
def createImageFigure(parent, dim, dpi=36):
    from numpy import zeros
    Z = zeros((dim, dim), float)
    ddim = dim/dpi
    figure = Figure(figsize=(ddim, ddim), dpi=dpi, frameon=False)
    # a tk.DrawingArea
    canvas = FigureCanvasTkAgg(figure, master=parent)
    canvas.get_tk_widget().grid(column=0, row=0)#, sticky=(N, W, E, S))
    #figureimg = figure.figimage(Z, cmap=cm.gray)#, origin='lower')
    #axdef = SubplotZero(figure, 111)
    ax = figure.add_axes([0,0,1,1], frameon=False)
#    ax.set_aspect('equal','datalim')
    #h = 0.5
    figureimg = ax.imshow(Z, cmap=cm.gray)#, extent=[-h, h, -h, h])
    #ax.set_aspect(aspect='equal',adjustable='datalim')
    ax.set_axis_off()
    #ax.set_axis_bgcolor('yellow')
    #ax.get_xaxis().set_visible(False)
    #ax.get_yaxis().set_visible(False)
    return [canvas, figure, figureimg]
    
w = None
def showImage(filename, dim=512, dpi=96):
    import xmipp
    from pylab import axes, Slider
    from protlib_xmipp import getImageData
    
    h = 0.5
    lf0 = 0.15
    hf0 = 0.35
    axcolor = 'lightgoldenrodyellow'
    
    img = xmipp.Image()
    img.readPreview(filename, dim)
    xdim, ydim, zdim, n = img.getDimensions()
    Z = getImageData(img)
    xdim += 10
    ydim += 10
    figure = Figure(figsize=(xdim/dpi, ydim/dpi), dpi=dpi, frameon=False)
    # a tk.DrawingArea
    root = tk.Tk()
    root.title(filename)
    canvas = FigureCanvasTkAgg(figure, master=root)
    canvas.get_tk_widget().grid(column=0, row=0)#, sticky=(N, W, E, S))
    #figureimg = figure.figimage(Z, cmap=cm.gray)#, origin='lower')
    ax = SubplotZero(figure, 111)
    ax = figure.add_subplot(ax)
    #axes([0.1, 0.1, 0.9, 0.9])
    #ax.set_aspect(1.0)
    ax.imshow(Z, cmap=cm.gray, extent=[-h, h, -h, h])
    for direction in ["xzero", "yzero"]:
        ax.axis[direction].set_visible(True)
    for direction in ["left", "right", "bottom", "top"]:
        ax.axis[direction].set_visible(False)
    global w
    w = Wedge((0,0), hf0, 0, 360, width=lf0, alpha=0.15) # Full ring
    ax.add_patch(w)
    
    def update(hf, lf):
        global w
        w.remove()
        w2 = Wedge((0,0), hf, 0, 360, width=lf, alpha=0.2)
        ax.add_patch(w2)
        w = w2
        canvas.draw()
        
    root.mainloop()
    
def getPngData(filename):  
    import matplotlib.image as mpimg
    return mpimg.imread(filename)

class PsdFigure():
    def __init__(self, master, dim, lf, hf, dpi=36, Z=None):
        self.dim = dim
        dd = dim/dpi
        figure = Figure(figsize=(dd, dd), dpi=dpi, frameon=False)
        self.canvas = FigureCanvasTkAgg(figure, master=master)
        self.canvas.get_tk_widget().grid(column=0, row=0)#, sticky=(N, W, E, S))
        axdef = SubplotZero(figure, 111)
        ax = figure.add_subplot(axdef)
        h = 0.5
        self.ring = True
        if Z is None:
            from numpy import zeros
            Z = zeros((dim, dim), float)
            self.ring = None
        self.img = ax.imshow(Z, cmap=cm.gray, extent=[-h, h, -h, h])
        
        for direction in ["xzero", "yzero"]:
            ax.axis[direction].set_visible(True)
        
        for direction in ["left", "right", "bottom", "top"]:
            ax.axis[direction].set_visible(False)

        self.figure = figure
        self.ax = ax
        self.lf = lf
        self.hf = hf
        if self.ring:
            self.createRing()
        else:
            self.canvas.draw()
        
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
        self.createRing()
    
    def updateData(self, Z):
        if self.ring:
            self.ring.remove()
            self.ring = None
        if Z is None:
            from numpy import zeros
            Z = zeros((self.dim, self.dim), float)
        else:
            self.createRing()
        self.img.set_array(Z)
        self.img.autoscale()
        self.canvas.draw()
        