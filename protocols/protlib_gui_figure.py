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
from protlib_utils import reportError
matplotlib.use('TkAgg')
import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.cm as cm
from mpl_toolkits.axes_grid.axislines import SubplotZero
from matplotlib.patches import Wedge

def createBgImage(dim):
    from numpy import ones
    return ones((dim, dim, 3))

class Preview():
    def __init__(self, parent, dim, dpi=36, label=None, col=0, row=0):
        self.dim = dim
        self.bg = np.zeros((dim, dim), float)
        ddim = dim/dpi
        self.figure = Figure(figsize=(ddim, ddim), dpi=dpi, frameon=False)
        self.frame = tk.Frame(parent)
        self.canvas = FigureCanvasTkAgg(self.figure, master=self.frame)
        self.canvas.get_tk_widget().grid(column=0, row=0)#, sticky=(N, W, E, S))
        if label:
            tk.Label(self.frame, text=label).grid(column=0, row=1)
        self.frame.grid(column=col, row=row)
        self._create_axes()
        
    def setWindowTitle(self,title):
        """ Set window title"""
        self.canvas.set_window_title(title)

    def _create_axes(self):
        pass
    
    def _update(self):
        pass
    
    def clear(self):
        self._update(self.bg)
        
    def updateData(self, Z):
        self.clear()
        self._update(Z)
    
    
class ImagePreview(Preview):
    def __init__(self, parent, dim, dpi=36, label=None, col=0):
        Preview.__init__(self, parent, dim, dpi, label, col)
    
    def _create_axes(self):
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
    def __init__(self, master, dim, lf, hf, dpi=72, Z=None, col=0):
        Preview.__init__(self, master, dim, dpi, label="PSD", col=col)
        self.lf = lf
        self.hf = hf
        if self.ring:
            self.createRing()
        else:
            self.canvas.draw()
                            
    def _create_axes(self):
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
           
               
    
#w = None
class ImageWindow():
    def __init__(self, filename=None, dim=512, dpi=96, image=None, label=None):
        
        dpi = min(dpi, dim)
        
        if filename is None and image is None:
            reportError("You should provide image or filename")
    
        import xmipp
        #from pylab import axes, Slider
        from protlib_xmipp import getImageData
        
        
        h = 0.5
        lf0 = 0.15
        hf0 = 0.35
        axcolor = 'lightgoldenrodyellow'
        
        if image is None:
            image = xmipp.Image()
            image.readPreview(filename, dim)
        if filename is None:
            filename = "No filename"
            
        xdim, ydim, zdim, n = image.getDimensions()
        Z = getImageData(image)
        xdim += 10
        ydim += 10
        figure = Figure(figsize=(xdim/dpi, ydim/dpi), dpi=dpi, frameon=False)
        # a tk.DrawingArea
        self.root = tk.Tk()
        self.root.title(filename)
        self.imagePreview=ImagePreview(self.root,dim, dpi)
        self.imagePreview.updateData(Z)
        
    def show(self):
        self.root.mainloop()
    
    def updateData(self, Z):
        self.imagePreview.updateData(Z)
    
    def updateImage(self, image):
        from protlib_xmipp import getImageData
        Z = getImageData(image)
        self.imagePreview.updateData(Z)
            

def showImage(filename=None, dim=512, dpi=96, image=None, label=None):
    ImageWindow(filename, dim, dpi, image, label).show()
    
    
    
    
    
#    canvas = FigureCanvasTkAgg(figure, master=root)
#    canvas.get_tk_widget().grid(column=0, row=0)#, sticky=(N, W, E, S))
#    ax = figure.add_axes([0.2,0.2,0.6,0.6], frameon=False)
#    ax.set_xlim(-0.5, 0.5)
#    ax.set_ylim(-0.5, 0.5)
#    self.figureimg = ax.imshow(Z, cmap=cm.gray, extent=[-h, h, -h, h])
#    
#    def update(hf, lf):
#        global w
#        w.remove()
#        w2 = Wedge((0,0), hf, 0, 360, width=lf, alpha=0.2)
#        ax.add_patch(w2)
#        w = w2
#        canvas.draw()
        

    
def getPngData(filename):  
    import matplotlib.image as mpimg
    return mpimg.imread(filename)

def showDependencyTree(runsDict):
    ''' This function will create a figure with a dependency 
    tree between the runs of the project '''
    from protlib_gui_figure import XmippPlotter
    import matplotlib.lines as mlines
    XSIZE, YSIZE = 8, 6
    DPI = 100
    XDIM, YDIM = XSIZE*DPI, YSIZE*DPI
    DY = 56
    DX = 50
    FONT = "sans-serif"
    FONTSIZE = 9
    colors = ['#D9F1FA', '#D9F1FA', '#FCCE62', '#D2F5CB', '#F5CCCB', '#F3F5CB', '#416FF0']
    xplotter = XmippPlotter(figsize=(XSIZE, YSIZE), dpi=100, windowTitle='Runs dependencies TREE')
    # Store the right horizontal x position for better packaging of
    # the graph, assuming max of 100 y-levels
    from numpy import zeros
    hLimits = zeros(100)
    #a = xplotter.createSubPlot("Test", "Particle number", "Zscore")
    a = xplotter.createCanvas()
    a.set_xlim(0, XDIM)
    a.set_ylim(top=0, bottom=YDIM)
    
    def showNode(dd, x, y):
        if dd.prot is None:
            nodeText = dd.extRunName
        else:
            nodeText = "%s\n%s" % (dd.protName, dd.runName)
        
        t = a.text(x, y, nodeText, family=FONT, size=FONTSIZE, 
                       bbox = dict(boxstyle="round", fc=colors[dd.state]))
        xplotter.draw()
        box = t.get_bbox_patch()
        dd.width = box.get_width()
        dd.height = box.get_height()
        dd.start = x
        
        return t
        
    def showLevel(dd, level):
        y = level * DY
        
        if len(dd.deps):
            #width = (xmax - xmin) / n
            childs = [runsDict[rn] for rn in dd.deps]
            for c in childs:
                showLevel(c, level + 1)
                
            firstChild = childs[0]
            lastChild = childs[-1]
            
            t = showNode(dd, 0, y)
            dd.start = (lastChild.start + lastChild.width + firstChild.start - dd.width) / 2
            dd.start = max(dd.start, hLimits[level] + DX)
            t.set_position((dd.start, y))            
            hLimits[level] = dd.start + dd.width     
            xx = [dd.start + dd.width/2, 0]
            yy = [y + dd.height - 18, 0]
            for c in childs:
                xx[1] = c.start + c.width/2
                yy[1] = y + DY - c.height
                a.add_line(mlines.Line2D(xx, yy, lw=2., alpha=0.4))
        else:
            t = showNode(dd, hLimits[level] + DX, y)
            hLimits[level] = dd.start + dd.width
    root = runsDict['PROJECT']
    showLevel(root, 0)    
    xplotter.show()



import numpy as np
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt

''' Class to create several plots'''
class XmippPlotter():
    def __init__(self, x=1, y=1, mainTitle="", figsize=None, dpi=100, windowTitle=""):
        
        if figsize is None: # Set some defaults values
            if x == 1 and y == 1:
                figsize = (6, 5)
            elif x == 1 and y == 2:
                figsize = (4, 6)
            elif x == 2 and y == 1:
                figsize = (6, 4)
            else:
                figsize = (8, 6)
        
        # Create grid
        self.grid = gridspec.GridSpec(x, y)#, height_ratios=[7,4])
        self.grid.update(left=0.15, right=0.95, hspace=0.25, wspace=0.4)#, top=0.8, bottom=0.2)  
        self.gridx = x
        self.gridy = y  
        self.figure = plt.figure(figsize=figsize, dpi=dpi)
        #self.mainTitle = mainTitle
        #self.windowTitle = windowTitle
        if (mainTitle):
            self.figure.suptitle(mainTitle)
        if (windowTitle):
            self.figure.canvas.set_window_title(windowTitle) 
        self.plot_count = 0
        self.plot_axis_fontsize = 10
        self.plot_text_fontsize = 8
        self.plot_yformat = '%1.2e'

    def showLegend(self, labels, loc='best'):
        leg = self.last_subplot.legend(tuple(labels), loc=loc)
        for t in leg.get_texts():
            t.set_fontsize(self.plot_axis_fontsize)    # the legend text fontsize
        
    def createSubPlot(self, title, xlabel, ylabel, xpos=None, ypos=None, 
                      yformat=False, projection='rectilinear'):
        '''
        Create a subplot in the figure. 
        You should provide plot title, and x and y axis labels.
        yformat True specified the use of global self.plot_yformat
        Posibles values for projection are: 
            'aitoff', 'hammer', 'lambert', 'mollweide', 'polar', 'rectilinear' 
        
        '''
        if xpos is None:
            self.plot_count += 1
            pos = self.plot_count
        else:
            pos = xpos + (ypos - 1) * self.gridx
        a = self.figure.add_subplot(self.gridx, self.gridy, pos, projection=projection)
        #a.get_label().set_fontsize(12)
        a.set_title(title)
        a.set_xlabel(xlabel)
        a.set_ylabel(ylabel)
            
        if yformat:
            formatter = ticker.FormatStrFormatter(self.plot_yformat)
            a.yaxis.set_major_formatter(formatter)
        a.xaxis.get_label().set_fontsize(self.plot_axis_fontsize)
        a.yaxis.get_label().set_fontsize(self.plot_axis_fontsize)
        labels = a.xaxis.get_ticklabels() + a.yaxis.get_ticklabels()
        for label in labels:
            label.set_fontsize(self.plot_text_fontsize) # Set fontsize
            label.set_text('aa')
                #print label.
                #label.set_visible(False)
        self.last_subplot = a
        self.plot = a.plot
        self.hist = a.hist
        return a
    
    def createCanvas(self):
        a = self.figure.add_subplot(111, axisbg='g')
        a.set_axis_off()
        self.figure.set_facecolor('white')
        return a
    
    def plotAngularDistribution(self, title, md, color='blue'):
        '''Create an special type of subplot, representing the angular
        distribution of weight projections. A metadata should be provided containing
        labels: MDL_ANGLE_ROT, MDL_ANGLE_TILT, MDL_WEIGHT '''
        from math import radians
        from xmipp import MDL_ANGLE_ROT, MDL_ANGLE_TILT, MDL_WEIGHT
        
        max_p = 40
        min_p = 5
        max_w = 2
        min_w = 1
        rot = [radians(md.getValue(MDL_ANGLE_ROT, objId)) for objId in md]
        tilt = [md.getValue(MDL_ANGLE_TILT, objId) for objId in md]
        weight = [md.getValue(MDL_WEIGHT, objId) for objId in md]
        
        if (len(weight) > 0):
            max_w = max(weight)
            min_w = min(weight)
            a = self.createSubPlot(title, 'Min weight=%(min_w)f, Max weight=%(max_w)f' % locals(), '', projection='polar')
        else:
            a = self.createSubPlot(title, 'Empty plot', '', projection='polar')
      
        for i, objId in enumerate(md):
            if (len(weight) > 0):
                pointsize = int((weight[i] - min_w)/(max_w - min_w + 0.001) * (max_p - min_p) + min_p)
            else:
                pointsize = 1
            a.plot(rot[i], tilt[i], markerfacecolor=color, marker='.', markersize=pointsize)
    
    def plotMd(self, md, mdLabelX, mdLabelY, color='g',**args):
        """ plot metadata columns mdLabelX and mdLabelY
            if nbins is in args then and histogram over y data is made
        """
        if mdLabelX:
            xx = []
        else:
            xx = range(1, md.size() + 1)
        yy = []
        for objId in md:
            if mdLabelX:
                xx.append(md.getValue(mdLabelX, objId))
            yy.append(md.getValue(mdLabelY, objId))
        
        nbins = args.pop('nbins', None)
        marker = args.pop('marker', None)
        linestyle = args.pop('linestyle', None)
        if nbins is None:
            if not marker is None:
                args['marker'] = marker     
            if not linestyle is None:
                args['linestyle'] = linestyle
            self.plot(xx, yy, color, **args) #no histogram
        else:
            self.hist(yy,nbins, facecolor=color, **args)
        
    def plotMdFile(self, mdFilename, mdLabelX, mdLabelY, color='g', **args):
        """ plot metadataFile columns mdLabelX and mdLabelY
            if nbins is in args then and histogram over y data is made
        """
        from xmipp import MetaData
        md = MetaData(mdFilename)
        self.plotMd(md, mdLabelX, mdLabelY, color='g',**args)
        
    def show(self):
        plt.tight_layout()
        plt.show()

    def draw(self):
        plt.tight_layout()
        plt.draw()
        
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.path as path
import mpl_toolkits.mplot3d.axes3d as p3

class XmippArrayPlotter1D():
    # Plot column number col from fnArray
    def __init__(self, fnArray, col, mainTitle="", xlabel="", ylabel="", figsize=None, dpi=100):
        import numpy
        a=numpy.loadtxt(fnArray)
        data=a[:,col]
        
        if figsize is None: # Set some defaults values
            figsize = (6, 5)
        fig = plt.figure(figsize=figsize, dpi=dpi)

        ax = fig.add_subplot(111)
        ax.set_title(mainTitle)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        
        # histogram our data with numpy
        n, bins = np.histogram(data, 50)
        
        # get the corners of the rectangles for the histogram
        left = np.array(bins[:-1])
        right = np.array(bins[1:])
        bottom = np.zeros(len(left))
        top = bottom + n
        
        # we need a (numrects x numsides x 2) numpy array for the path helper
        # function to build a compound path
        XY = np.array([[left,left,right,right], [bottom,top,top,bottom]]).T
        
        # get the Path object
        barpath = path.Path.make_compound_path_from_polys(XY)
        
        # make a patch out of it
        patch = patches.PathPatch(barpath, facecolor='blue', edgecolor='gray', alpha=0.8)
        ax.add_patch(patch)
        
        # update the view limits
        ax.set_xlim(left[0], right[-1])
        ax.set_ylim(bottom.min(), top.max())
        
        # plt.show()
        fig.show()

class XmippArrayPlotter2D():
    # 2D plot of columns X and Y from fnArray
    def __init__(self, fnArray, colX, colY, mainTitle="", xlabel="", ylabel="", figsize=None, dpi=100):
        import numpy
        a=numpy.loadtxt(fnArray)
        
        if figsize is None: # Set some defaults values
            figsize = (5, 5)
        fig = plt.figure(figsize=figsize, dpi=dpi)

        ax = fig.add_subplot(111)
        ax.set_title(mainTitle)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        
        ax.plot(a[:,colX], a[:,colY], '.')
        fig.show()

class XmippArrayPlotter3D():
    # 2D plot of columns X and Y from fnArray
    def __init__(self, fnArray, colX, colY, colZ, mainTitle="", xlabel="", ylabel="", zlabel="", figsize=None, dpi=100):
        import numpy
        a=numpy.loadtxt(fnArray)
        
        if figsize is None: # Set some defaults values
            figsize = (5, 5)
        fig = plt.figure(figsize=figsize, dpi=dpi)

        ax = p3.Axes3D(fig)
        ax.set_title(mainTitle)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_zlabel(zlabel)
        
        ax.scatter3D(a[:,colX], a[:,colY], a[:,colZ])
        fig.show()