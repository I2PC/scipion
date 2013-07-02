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
This module implement the classes to create plots on xmipp.
"""
import Tkinter as tk
import ttk
import matplotlib
matplotlib.use('TkAgg')
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
        
