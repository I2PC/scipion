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

import os
import sys

from pyworkflow.protocol.params import LabelParam
from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO
from pyworkflow.em.viewer import ImageView, ChimeraView
from protocol_resmap import ProtResMap



class ResMapViewer(ProtocolViewer):
    """
    Visualization tools for ResMap results.
    
    ResMap is software tool for computing the local resolution of 3D
    density maps studied in structural biology, primarily by cryo-electron
    microscopy (cryo-EM).
     
    Please find the manual at http://resmap.sourceforge.net 
    """
           
    _environments = [DESKTOP_TKINTER]
    _targets = [ProtResMap]
    _label = 'viewer resmap'
    
    def __init__(self, *args, **kwargs):
        ProtocolViewer.__init__(self, *args, **kwargs)
        sys.path.append(os.environ['RESMAP_HOME'])
    
    def _defineParams(self, form):
        form.addSection(label='Visualization')
        group = form.addGroup('2D Plots')
        group.addParam('doShowVolumeSlices', LabelParam,
                      label="Show volume slices?")
        group.addParam('doShowResMapSlices', LabelParam,
                      label="Show ResMap slices?")               
        group.addParam('doShowResHistogram', LabelParam,
                      label="Show resolution histogram?")
        
        form.addParam('doShowChimera', LabelParam,
                      label="Show Chimera animation?", default=True)
        
        
    def _getVisualizeDict(self):
        return {'doShowVolumeSlices': self._showVolumeSlices,
                'doShowResMapSlices': self._showResMapSlices,
                'doShowResHistogram': self._plotHistogram,
                'doShowChimera': self._showChimera,
                }
        
    def _showVolumeSlices(self, param=None):
        return [self.protocol._plotVolumeSlices()]
        
    def _showResMapSlices(self, param=None):
        return [self.protocol._plotResMapSlices()]
             
    def _plotHistogram(self, param=None):
        return [self.protocol._plotHistogram()]
    
    def _showChimera(self, param=None):
        #os.system('chimera "%s" &' % self.protocol._getPath('volume1_resmap_chimera.cmd'))
        cmdFile = self.protocol._getPath('volume1_resmap_chimera.cmd')
        view = ChimeraView(cmdFile)
        return [view]

        
        
class ResMapViewerWeb(ResMapViewer):
    """
    Same viewer for ResMap web, but using saved images of the plots.
    """
           
    _environments = [WEB_DJANGO]
    
    def _showVolumeSlices(self, param=None):
        return [ImageView(self.protocol._getExtraPath('volume1.map.png'))]
        
    def _showResMapSlices(self, param=None):
        return [ImageView(self.protocol._getExtraPath('volume1_resmap.map.png'))]
    
    def _plotHistogram(self, param=None):
        return [ImageView(self.protocol._getExtraPath('histogram.png'))]
    
        