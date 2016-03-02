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
This module implements visualization program for Spider CA-PCA protocol.
"""

import Tkinter as tk
from pyworkflow.protocol.params import *
from pyworkflow.viewer import Viewer, ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO,\
    TextView
from pyworkflow.utils.graph import Graph
from pyworkflow.gui.graph import LevelTree
from pyworkflow.gui.canvas import Canvas, ImageBox
from pyworkflow.em.packages.xmipp3.viewer import XmippViewer
from pyworkflow.em.viewer import DataView
from pyworkflow.em.plotter import EmPlotter

from spider import PcaFile
from protocol.protocol_ca_pca import SpiderProtCAPCA



class SpiderViewerCAPCA(ProtocolViewer):
    """ Visualization of CA PCA Protocol. """
           
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [SpiderProtCAPCA]
    _label = 'viewer CAPCA'
    
    def _defineParams(self, form):
        form.addSection(label='Visualization')
        form.addParam('doShowEigenImages', LabelParam, label="Show eigenimages?", default=True,
                      help='Display eigenimages produced by CA/PCA')
        form.addParam('doShowReconsImages', LabelParam,
                      label="Show reconstitued images?", default=True,
                      help='Re-creation of images from eigenvectors can eliminate the noise in the '
                      'reconstituted image, this also results in large data compression. Here we '
                      'reconstitute positive and negative images')
        form.addParam('doShowHistogram', LabelParam,
                      label="Show eigenvalue histogram?", default=True,
                      help='One of the methods to determine what eigenvalues are useful, '
                      'and which are from noise is to view a histogram showing the percentage '
                      'of eigenvalue variance accounted for by each factor.')
        form.addParam('doShowPcaFile', LabelParam, #expertLevel=LEVEL_ADVANCED,
                      label="Show IMC file?", default=True,
                      help='This file contains coordinates of each image in the new vector space.')
        form.addParam('doShowFactorMaps', LabelParam, #expertLevel=LEVEL_ADVANCED,
                      label="Show factor maps?", default=True,
                      help='Once you know which eigenvectors have some meaning and which are from noise, '
                      'you can display 2D factor maps of selected pairs of factors to visualize clustering (if any).')        
        line = form.addLine('Factors')
        line.addParam('firstFactor', IntParam, #expertLevel=LEVEL_ADVANCED,
                      label="First", default=1) 
        line.addParam('secondFactor', IntParam, #expertLevel=LEVEL_ADVANCED,
                      label="Second", default=2)
                      
    def _getVisualizeDict(self):
        return {'doShowEigenImages': self._viewParam,
                'doShowReconsImages': self._viewParam,
                'doShowHistogram': self._plotHistogram,
                'doShowFactorMaps': self._plotFactorMaps,
                'doShowPcaFile': self._viewParam,
                }
        
    def _viewParam(self, param=None):
        if param == 'doShowEigenImages':
            view = DataView(self.protocol._getFileName('eigenimages'))
        elif param == 'doShowReconsImages':
            view = DataView(self.protocol._getFileName('reconstituted'))
        elif param == 'doShowPcaFile':
            view = self.textView([self.protocol.imcFile.filename.get()], "PCA file")
            
        return [view]
            
    def _plotHistogram(self, param=None):
        """ First we parse the cas_EIG file and we read:
        first line: take the number of eigen values.
        then one line per factor and we read the percent and cumulative percent.
        """
        from numpy import arange
        from matplotlib.ticker import FormatStrFormatter
        
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
        xplotter = EmPlotter()
        a = xplotter.createSubPlot('Eigenvalues histogram', 'Eigenvalue number', '%')
        a.set_xticks(factors + 0.45)
        a.xaxis.set_major_formatter(FormatStrFormatter('%1.0f'))
        bars = a.bar(factors, percents, width, color='b')
        
        for i, rect in enumerate(bars):
            h = rect.get_height()
            a.text(rect.get_x()+rect.get_width()/2., h+0.3, '%d' % cumPercents[i],
                ha='center', va='bottom')
        a.set_ylim([0, percents[0] + 5])
        
        return [xplotter]
        
    def _plotFactorMaps(self, param=None):
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
        xplotter = EmPlotter(1,1)
        a = xplotter.createSubPlot("Factor %d vs %d" % (x, y), 
                                   "Factor %d" % x, "Factor %d" % y)
        a.plot(xFactors, yFactors, 'o')
        
        return [xplotter]
        
