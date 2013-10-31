# **************************************************************************
# *
# * Authors:     L. del Cano (ldelcano@cnb.csic.es)
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
This module implement the wrappers aroung Xmipp ML2D protocol
visualization program.
"""
from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO, createPlots
from pyworkflow.em import *
from pyworkflow.gui.form import FormWindow
from protocol_ml2d import XmippProtML2D
from viewer import runShowJ

import numpy as np



class XmippML2DViewer(ProtocolViewer):
    """ Wrapper to visualize different type of data objects
    with the Xmipp program xmipp_showj
    """
    _targets = [XmippProtML2D]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    
    _label = 'Xmipp Viewer ML2D'
    _plotVars = ['doShowLL', 'doShowPmax', 'doShowSignalChange', 'doShowMirror'] 
    
    def _defineParams(self, form):
        form.addSection(label='Visualization')
        form.addParam('doShowClasses', BooleanParam, label="Visualize last iter references", default=True, 
                      help='Visualize last iteration references.')
        form.addParam('doShowPlots', BooleanParam, label="Show all plots per iteration?", default=True,
                      help='Visualize several plots.')
        
        form.addSection(label='Iteration plots')    
        form.addParam('doShowLL', BooleanParam, label="Show Log-Likehood over iterations?", default=False, 
                      help='The Log-Likelihood value should increase.')      
        form.addParam('doShowPmax', BooleanParam, label="Show maximum model probability?", default=False, 
                      help='Show the maximum probability for a model, this should tend to be a deltha function.')      
        form.addParam('doShowSignalChange', BooleanParam, label="Show plot for signal change?", default=False, 
                      help='Should approach to zero when convergence.')      
        form.addParam('doShowMirror', BooleanParam, label="Show mirror fraction for last iteration?", default=False, 
                      help='he the mirror fraction of each reference in last iteration.')      
        
    
    def _getVisualizeDict(self):
        return {'doShowClasses': self._viewIterRefs,
                'doShowPlots': self._viewAllPlots,
                'doShowLL': self._viewPlot,
                'doShowPmax': self._viewPlot,
                'doShowSignalChange': self._viewPlot,
                'doShowMirror': self._viewPlot}
        
    def _viewAll(self, *args):
        if self.doShowClasses:
            self._viewIterRefs()
        if self.doShowPlots:
            self._viewAllPlots()
        else:
            plots = [p for p in  self._plotVars if self.getAttributeValue(p)]
            print plots
            xplotter = createPlots(self.protocol, plots)
            if xplotter:
                xplotter.show()
        
    def _viewAllPlots(self, e=None):
        xplotter = createPlots(self.protocol, self._plotVars)
        if xplotter:
            xplotter.show()
        
    def _viewPlot(self, paramName):
        xplotter = createPlots(self.protocol, [paramName])
        if xplotter:
            xplotter.show()
        
    def _viewIterRefs(self, e=None):
        lastIter = self.protocol._lastIteration()
        runShowJ('classes@' + self.protocol._getIterClasses())
        print "lastIter: ", lastIter
    
    def getVisualizeDictWeb(self):
        return {'doShowClasses': "doShowClasses",
                'doShowPlots': "doAllPlotsML2D",
                'doShowLL': "doShowLL",
                'doShowPmax': "doShowPmax",
                'doShowSignalChange': "doShowSignalChange",
                'doShowMirror': "doShowMirror"}

    @classmethod
    def getView(cls):
        """ This function will notify the web viewer for this protocol"""
        return "viewerForm"
    
    @classmethod
    def getViewFunction(cls):
        """ This will return the name of the function to view
        in web one (or all) params of the protocol"""
        return "viewerML2D"
        
    