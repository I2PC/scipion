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
Visualization of the ResMap outputs.
"""

import os
import sys

from pyworkflow.protocol.params import BooleanParam
from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO
from pyworkflow.gui.plotter import Plotter
from protocol_resmap import ProtResMap



class ResMapViewer(ProtocolViewer):
    """
    Visualization tools for ResMap results.
    
    ResMap is software tool for computing the local resolution of 3D
    density maps studied in structural biology, primarily by cryo-electron
    microscopy (cryo-EM).
     
    Please find the manual at http://resmap.sourceforge.net 
    """
           
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [ProtResMap]
    _label = 'viewer resmap'
    
    def __init__(self, *args, **kwargs):
        ProtocolViewer.__init__(self, *args, **kwargs)
        sys.path.append(os.environ['RESMAP_HOME'])
    
    def _defineParams(self, form):
        form.addSection(label='Visualization')
        group = form.addGroup('2D Plots')
        group.addParam('doShowVolumeSlices', BooleanParam, default=True, 
                      label="Show volume slices?")
        group.addParam('doShowResMapSlices', BooleanParam, default=True, 
                      label="Show ResMap slices?")               
        group.addParam('doShowResHistogram', BooleanParam, 
                      label="Show resolution histogram?", default=True)
        
        form.addParam('doShowChimera', BooleanParam, 
                      label="Show Chimera animation?", default=True)
        
        
    def _getVisualizeDict(self):
        return {'doShowVolumeSlices': self._showVolumeSlices,
                'doShowResMapSlices': self._showResMapSlices,
                'doShowResHistogram': self._plotHistogram,
                'doShowChimera': self._showChimera,
                }
        
    def _getVolumeMatrix(self, volName):
        from ResMap_fileIO import MRC_Data
        volPath = self.protocol._getPath(volName)
        return MRC_Data(volPath, 'ccp4').matrix
    
    def _showVolumeSlices(self, param=None):
        from ResMap_visualization import plotOriginalVolume
        fig = plotOriginalVolume(self._getVolumeMatrix('volume1.map'))
        return [Plotter(figure=fig)]
        
    def _showResMapSlices(self, param=None):
        from ResMap_visualization import plotResMapVolume
        fig = plotResMapVolume(self._getVolumeMatrix('volume1.map'),
                                  minRes=self.protocol.minRes.get(),
                                  maxRes=self.protocol.maxRes.get())
        return [Plotter(figure=fig)]
             
    def _plotHistogram(self, param=None):
        """ First we parse the cas_EIG file and we read:
        first line: take the number of eigen values.
        then one line per factor and we read the percent and cumulative percent.
        """
        from ResMap_visualization import plotResolutionHistogram
        from cPickle import loads
        histogramData = loads(self.protocol.histogramData.get())
        fig = plotResolutionHistogram(histogramData)
        return [Plotter(figure=fig)]
    
    def _showChimera(self, param=None):
        os.system('chimera "%s" &' % self.protocol._getPath('volume1_resmap_chimera.cmd'))
        
