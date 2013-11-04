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
from protocol_ca_pca import SpiderProtCAPCA



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
        
    def getVisualizeDictWeb(self):
        return {'doShowEigenImages': 'doShowEigenImages',
                'doShowReconsImages': 'doShowReconsImages',
                'doShowHistogram': self._plotHistogram,
                'doShowFactorMaps': self._plotFactorMaps,
                'doShowPcaFile': 'doShowPcaFile' }
        
    @classmethod
    def getView(cls):
        """ This function will notify the web viewer for this protocol"""
        return "viewerForm"
    
    @classmethod
    def getViewFunction(cls):
        """ This will return the name of the function to view
        in web one (or all) params of the protocol"""
        return "viewerCAPCA"