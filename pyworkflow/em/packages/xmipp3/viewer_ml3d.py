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
This module implement the wrappers aroung Xmipp ML3D protocol
visualization program.
"""
from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO
from pyworkflow.em import *
from pyworkflow.gui.form import FormWindow
from protocol_ml3d import XmippProtML3D
from viewer import runShowJ

import numpy as np

LAST_ITER = 0
ALL_ITER = 1
SELECTED_ITERS = 2

class XmippML3DViewer(ProtocolViewer):
    """ Wrapper to visualize different type of data objects
    with the Xmipp program xmipp_showj
    """
    _targets = [XmippProtML3D]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    
    _label = 'Xmipp Viewer ML3D'

    
    def _defineParams(self, form):
        form.addSection(label='Preparation')
        form.addParam('doShowGreyScaleRefVol', BooleanParam, label="Visualize the grey-scale corrected reference volume?", default=True)
        form.addParam('doShowFilterRefVol', BooleanParam, label="Visualize the low-pass filtered reference volume?", default=True)
        form.addParam('doShowSeedsVols', BooleanParam, label="Visualize the generated seeds volumes?", default=True)
        
        form.addSection(label='Results per Iteration')
        form.addParam('iterToShow', EnumParam, label="Which iteration do you want to visualize?", default=0, 
                      choices=['last','all','selection'], display=EnumParam.DISPLAY_LIST)
        form.addParam('selectedIters', StringParam, default='',
              label='Selected iterations', condition='iterToShow==2',
              help='Which iteration do you want to visualize.')     
        form.addParam('doShow2DAvgs', BooleanParam, label="Visualize weighted 2D-averages?.", default=True)
        form.addParam('doShow3DRefsVolumes', BooleanParam, label="Visualize the 3D-references volumes?", default=True)
        form.addParam('doShowAngDist', BooleanParam, label="Plot angular distribution of 3D-references?", default=True)
        form.addParam('doShowDataDist', BooleanParam, label="Plot data distribution over 3d-references?", default=True)     
        
        form.addSection(label='Overall Results')
        form.addParam('doShowStatistics', BooleanParam, label="Plot overall convergence statistics?", default=True)
    
    def _getVisualizeDict(self):
        return {'doShowGreyScaleRefVol': self._viewCorrectedVols,
                'doShowFilterRefVol': self._viewFilteredVols,
                'doShowSeedsVols': self._viewGeneratedVols,
#                'doShow2DAvgs': self._view2DAvgs,
#                'doShow3DRefsVolumes': self._view3DRefsVolumes,
#                'doShowAngDist': self._viewAngDist,
#                'doShowDataDist': self._viewDataDist
                }
        
        
    def _viewCorrectedVols(self, e=None):
        runShowJ(self.protocol._getExtraPath("corrected_volumes.stk"))
        
    def _viewFilteredVols(self, e=None):
        runShowJ(self.protocol._getExtraPath("filtered_volumes.stk"))
        
    def _viewGeneratedVols(self, e=None):
        runShowJ(self.protocol._getExtraPath("generated_volumes.stk"))
        
        
    def setVisualizeIterations(self):
        '''Validate and set the set of iterations to visualize.
        If not set is selected, only use the last iteration'''
        self.lastIter = self.protocol._lastIteration(self)
        
        if self.iterToShow.get() == LAST_ITER:
            self.visualizeIters = [self.lastIter]
            
        elif self.iterToShow.get() == ALL_ITER:
            self.visualizeIters = range(1, self.lastIter+1)
        elif self.iterToShow.get() == SELECTED_ITERS:
            selection = self._getListFromRangeString(self.selectedIters.get())
            self.visualizeIters = [it for it in selection if (it > 0 and it <= self.lastIter)]
            
    def _viewPlot(self, paramName):
        xplotter = self.createPlots(self.protocol, [paramName])
        if xplotter:
            xplotter.show()
        

    def getVisualizeDictWeb(self):
        return {'doShowGreyScaleRefVol': "viewCorrectedVols",
                'doShowFilterRefVol': "viewFilteredVols",
                'doShowSeedsVols': "viewGeneratedVols",
#                'doShow2DAvgs': self._view2DAvgs,
#                'doShow3DRefsVolumes': self._view3DRefsVolumes,
#                'doShowAngDist': self._viewAngDist,
#                'doShowDataDist': self._viewDataDist
                }

    @classmethod
    def getView(cls):
        """ This function will notify the web viewer for this protocol"""
        return "viewerForm"
    
    @classmethod
    def getViewFunction(cls):
        """ This will return the name of the function to view
        in web one (or all) params of the protocol"""
        return "viewerML3D"
    