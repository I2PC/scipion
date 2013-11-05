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

from protocol_ward import SpiderProtClassifyWard


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
        return {'doShowDendrogram': self._plotDendrogram,
                'doShowClasses': self.visualizeClasses
                }
        
    def _viewAll(self, *args):
        if self.doShowClasses:
            self.visualizeClasses()
        if self.doShowDendrogram:
            self._plotDendrogram()
            
    def _plotDendrogram(self, e=None):
        from pyworkflow.em.packages.xmipp3.plotter import XmippPlotter
        xplotter = XmippPlotter()
        self.plt = xplotter.createSubPlot("Dendrogram", "", "")
        self.step = 0.25
        self.rightMost = 0.0 # Used to arrange leaf nodes at the bottom
        
        node = self.protocol.buildDendrogram()
        self.plotNode(node, self.minHeight.get())    
        self.plt.xlim([0., self.rightMost + self.step])
        self.plt.ylim([-0.1, 105])
        
        return self._showOrReturn(xplotter)
    
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
        node.box = SpiderImageBox(canvas, node.path, text=node.getName(), y=y)
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
        
    def getVisualizeDictWeb(self):
        return {'doShowDendrogram': "doVisualizeDendrogram",
                'doShowClasses': "doVisualizeClasses"}
        
    @classmethod
    def getView(cls):
        """ This function will notify the web viewer for this protocol"""
        return "viewerForm"
    
    @classmethod
    def getViewFunction(cls):
        """ This will return the name of the function to view
        in web one (or all) params of the protocol"""
        return "viewerWARD"
        
        
class SpiderImageBox(ImageBox):
    def _onClick(self, e=None):
        self.selected = getattr(self, 'selected', False)
        # On click change the selection state
        self.selected = not self.selected
        if self.selected:
            self.label.config(bd=2, bg='green')
            
        else:
            self.label.config(bd=0, bg='grey')
            