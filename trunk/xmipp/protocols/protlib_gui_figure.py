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
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.cm as cm

# Some utilities to use matplotlib
def createImageFigure(parent, dim):
    from numpy import zeros
    Z = zeros((dim, dim), float)
    dpi = 36
    ddim = dim/dpi
    figure = Figure(figsize=(ddim, ddim), dpi=dpi)
    # a tk.DrawingArea
    canvas = FigureCanvasTkAgg(figure, master=parent)
    canvas.get_tk_widget().grid(column=0, row=0)#, sticky=(N, W, E, S))
    figureimg = figure.figimage(Z, cmap=cm.gray)#, origin='lower')
    return (canvas, figureimg)
    
w = None
def showImage(filename):
    import xmipp
    from mpl_toolkits.axes_grid.axislines import SubplotZero
    from pylab import axes, Slider
    from matplotlib.patches import Wedge
    from protlib_xmipp import getImageData
    
    dim = 512
    dpi = 96
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
    axlfreq = figure.add_axes([0.25, 0.06, 0.5, 0.03], axisbg=axcolor)
    axhfreq  = figure.add_axes([0.25, 0.02, 0.5, 0.03], axisbg=axcolor)
    low_freq = Slider(axlfreq, 'Low Freq', 0.01, 0.5, valinit=lf0)
    high_freq = Slider(axhfreq, 'High Freq', 0.01, 0.5, valinit=hf0)
    
    def update(val):
        #ax.remo
        global w
        w.remove()
        hf = high_freq.val
        lf = low_freq.val
        w2 = Wedge((0,0), hf, 0, 360, width=lf, alpha=0.2)
        #w.set(r=high_freq.val)
        #print w.properties()
        ax.add_patch(w2)
        w = w2
        #w.set_clip_path(w2.get_clip_path())
        canvas.draw()
    low_freq.on_changed(update)
    high_freq.on_changed(update)
    #canvas.show()
    root.mainloop()
    
def getPngData(filename):  
    import matplotlib.image as mpimg
    return mpimg.imread(filename)