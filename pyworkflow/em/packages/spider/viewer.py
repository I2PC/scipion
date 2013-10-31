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
This module implement the wrappers around xmipp_showj
visualization program.
"""
import Tkinter as tk
from pyworkflow.protocol.params import *
from pyworkflow.viewer import Viewer, ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO
from pyworkflow.utils.graph import Graph
from pyworkflow.gui.graph import LevelTree
from pyworkflow.gui.canvas import Canvas, ImageBox
from pyworkflow.em.packages.xmipp3.viewer import XmippViewer, runShowJ
from pyworkflow.gui.text import showTextfileViewer

from spider import PcaFile
from protocol_filters import SpiderProtFilter
from protocol_custommask import SpiderProtCustomMask
from protocol_align_apsr import SpiderProtAlignAPSR
from protocol_ca_pca import SpiderProtCAPCA
from protocol_ward import SpiderProtClassifyWard


    
class SpiderViewerGeneric(Viewer):
    """ Wrapper to visualize different type of objects
    with the Xmipp program xmipp_showj
    """
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [PcaFile, SpiderProtCustomMask, SpiderProtFilter, SpiderProtAlignAPSR]

    def visualize(self, obj, **args):
        
        if isinstance(obj, PcaFile):
            showTextfileViewer("PCA file", [obj.filename.get()])
        
        elif isinstance(obj, SpiderProtFilter):
            XmippViewer().visualize(obj.outputParticles)
            
        elif isinstance(obj, SpiderProtAlignAPSR):
            xv = XmippViewer()
            xv.visualize(obj.outputAverage)
            xv.visualize(obj.outputParticles)
            
        elif isinstance(obj, SpiderProtCustomMask):
            mask = obj.outputMask
            XmippViewer().visualize(mask)
            # Remove location to visualize the whole stack
            runShowJ(mask.getFileName())
         
         
class SpiderViewerCAPCA(ProtocolViewer):
    """ Visualization of CA PCA Protocol.
    """       
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [SpiderProtCAPCA]
    
    def _defineParams(self, form):
        form.addSection(label='Visualization')
        form.addParam('doShowEigenImages', BooleanParam, label="Show eigenimages?", default=True)
        form.addParam('doShowReconsImages', BooleanParam, 
                      label="Show reconstitued images?", default=True)
        form.addParam('doShowHistogram', BooleanParam, 
                      label="Show eigenvalue histogram?", default=True)
        form.addParam('doShowPcaFile', BooleanParam, #expertLevel=LEVEL_ADVANCED,
                      label="Show IMC file?", default=True)        
                      
        form.addSection(label='Factor maps')
        form.addParam('doShowFactorMaps', BooleanParam, #expertLevel=LEVEL_ADVANCED,
                      label="Show factor maps?", default=True)
        form.addParam('firstFactor', IntParam, #expertLevel=LEVEL_ADVANCED,
                      label="First factor", default=1,
                      help='') 
        form.addParam('secondFactor', IntParam, #expertLevel=LEVEL_ADVANCED,
                      label="Second factor", default=2,
                      help='')
                      


    def _getVisualizeDict(self):
        return {'doShowEigenImages': self._viewParam,
                'doShowReconsImages': self._viewParam,
                'doShowHistogram': self._plotHistogram,
                'doShowFactorMaps': self._plotFactorMaps,
                'doShowPcaFile': self._viewParam,
                }

    def _viewParam(self, param=None):
        if param == 'doShowEigenImages':
            runShowJ(self.protocol._getFileName('eigenimages'))
        elif param == 'doShowReconsImages':
            runShowJ(self.protocol._getFileName('reconstituted'))
        elif param == 'doShowPcaFile':
            showTextfileViewer("PCA files", [self.protocol.imcFile.filename.get()])
            
    def _plotHistogram(self, param=None):
        """ First we parse the cas_EIG file and we read:
        first line: take the number of eigen values.
        then one line per factor and we read the percent and cumulative percent.
        """
        from numpy import arange
        from matplotlib.ticker import FormatStrFormatter
        from pyworkflow.em.packages.xmipp3.plotter import XmippPlotter
        
        fn = self.protocol._getFileName('eigFile')
        f = open(fn)
        values = f.readline().split()
        n = int(values[0]) # Number of factors
        factors = arange(1, n+1)
        percents = []
        cumPercents = []
        
        for i in factors:
            values = f.readline().split()
            percents.append(float(values[1]))
            cumPercents.append(float(values[2]))
            
        f.close()
        
        width = 0.85
        xplotter = XmippPlotter(1,1)
        a = xplotter.createSubPlot('Eigenvalues histogram', 'Eigenvalue number', '%')
        a.set_xticks(factors + 0.45)
        a.xaxis.set_major_formatter(FormatStrFormatter('%1.0f'))
        bars = a.bar(factors, percents, width, color='b')
        
        for i, rect in enumerate(bars):
            h = rect.get_height()
            a.text(rect.get_x()+rect.get_width()/2., h+0.3, '%d' % cumPercents[i],
                ha='center', va='bottom')
        a.set_ylim([0, percents[0] + 5])
        #a.set_xlim([0.8, n + 1])
        
        xplotter.show()
        
        #showTextfileViewer("PCA files", [])
        
    def _plotFactorMaps(self, param=None):
        from pyworkflow.em.packages.xmipp3.plotter import XmippPlotter
        
        # Parse the file
        fn = self.protocol._getFileName('imcFile')
        f = open(fn)
        values = f.readline().split()
        n = int(values[0]) # Number of images
        nf = int(values[1]) # Number of factors
        
        x = self.firstFactor.get()
        y = self.secondFactor.get()
        xFactors = []
        yFactors = []
        i = 0
        while i < n:
            imgFactors = []
            while len(imgFactors) < nf:
                values = f.readline().split()
                imgFactors += [float(v) for v in values]
            xFactors.append(imgFactors[x-1])
            yFactors.append(imgFactors[y-1])
            i += 1
        f.close() 
        
        # Create the plot
        xplotter = XmippPlotter(1,1)
        a = xplotter.createSubPlot("Factor %d vs %d" % (x, y), 
                                   "Factor %d" % x, "Factor %d" % y)
        a.plot(xFactors, yFactors, 'o')
        xplotter.show()
        
    
class SpiderViewerWard(ProtocolViewer):
    """ Visualization of Classify Ward.
    """       
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [SpiderProtClassifyWard]
    
    def _defineParams(self, form):
        form.addSection(label='Visualization')
        form.addParam('doShowDendrogram', BooleanParam, label="Show dendrogram?", default=True,
                      help='')
        form.addParam('minHeight', FloatParam, default=0.5,
                      label='Minimum height',
                      help='The dendrogram will be show until that height')
        form.addParam('doShowClasses', BooleanParam, label="Visualize class averages?", default=True, 
                      help='')
        form.addParam('maxLevel', IntParam, default=5,
                      label='Maximum level',
                      help='Maximum level of classes to show')

    def _getVisualizeDict(self):
        return {'doShowDendrogram': self.visualizeDendrogram,
                'doShowClasses': self.visualizeClasses
                }
        
    def _viewAll(self, *args):
        if self.doShowClasses:
            self.visualizeClasses()
        if self.doShowDendrogram:
            self.visualizeDendrogram()
            
    def visualizeDendrogram(self, e=None):
        import matplotlib.pyplot as plt
        self.plt = plt
        self.step = 0.25
        self.rightMost = 0.0 # Used to arrange leaf nodes at the bottom
        
        node = self.protocol.buildDendrogram()
        self.plotNode(node, self.minHeight.get())    
        plt.xlim([0., self.rightMost + self.step])
        plt.ylim([-0.1, 105])
        plt.show()
    
    def plotNode(self, node, minHeight=-1):
        childs = node.get('childs', [])
        h = node.get('height')
        if h > minHeight and len(childs) > 1:
            x1, y1 = self.plotNode(childs[0], minHeight)
            x2, y2 = self.plotNode(childs[1], minHeight)
            length = node.get('length')
            index = node.get('index')
            xm = (x1 + x2)/2
            x = [x1, x1, x2, x2]
            y = [y1, h, h, y2]
            self.plt.plot(x, y, color='b')
            self.plt.plot(xm, h, 'ro')
            self.plt.annotate("%d(%d)" % (index, length), (xm, h), xytext=(0, -8),
                         textcoords='offset points', va='top', ha='center')
            point = (xm, h)
        else:
            self.rightMost += self.step
            point = (self.rightMost, 0.)
        return point
    
    def _createNode(self, canvas, node, y):
        node.box = ImageBox(canvas, node.path, text=node.getName(), y=y)
        return node.box

    def visualizeClasses(self, e=None):
        classTemplate = "class_%03d"
        averages = '%03d@' + self.protocol._getFileName('averages')
        
        def getInfo2(level, classNo):
            return classTemplate % classNo, averages % classNo
        
        n, p = getInfo2(0, 1)
        g = Graph(rootName=n)
        a = g.getRoot()
        a.path = p
        maxLevel = self.maxLevel.get()
        
        def addChilds(node, nodeNumber, level):
            imgNo = nodeNumber * 2 
            for off in [0, 1]:
                n, p = getInfo2(level+1, imgNo + off)
                b = g.createNode(n)
                b.path = p
                node.addChild(b)
                if level < maxLevel - 2:
                    addChilds(b, imgNo + off, level + 1)
               
        root = tk.Toplevel()
        canvas = Canvas(root, width=600, height=500)
        canvas.grid(row=0, column=0, sticky='nsew')
        root.grid_columnconfigure(0, weight=1)
        root.grid_rowconfigure(0, weight=1)         
            
        addChilds(a, 1, 0)
            
        lt = LevelTree(g)
        lt.DY = 135 # TODO: change in percent of the image size
        lt.setCanvas(canvas)
        lt.paint(self._createNode)
        canvas.updateScrollRegion()
