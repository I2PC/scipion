# **************************************************************************
# *
# * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
# *              J.M. de la Rosa Trevin  (jmdelarosa@cnb.csic.es)
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
# *  e-mail address 'jgomez@cnb.csic.es'
# *
# **************************************************************************

import os
import sys

import pyworkflow.protocol.params as params
import pyworkflow.viewer as pwviewer
import pyworkflow.em.showj as showj
from pyworkflow.em.viewer import DataView

from protocol_autofocus_classify import ProtAutofocusClassify
from protocol_frm_alignment import ProtFrmAlign


VOLUME_SLICES = 0
VOLUME_CHIMERA = 1


class PyTomViewer(pwviewer.ProtocolViewer):
    """
    Base class to visualize results from PyTom protocols
    (right now for ProtFrmAlign and ProtAutofocusClassify
    """
    _environments = [pwviewer.DESKTOP_TKINTER, pwviewer.WEB_DJANGO]
    
    def __init__(self, *args, **kwargs):
        pwviewer.ProtocolViewer.__init__(self, *args, **kwargs)
    
    def _defineParams(self, form):
        form.addSection(label='Visualization')
        group = form.addGroup('General')
        group.addParam('displayVol', params.EnumParam, 
                       choices=['slices', 'chimera'],
                       default=VOLUME_CHIMERA, 
                       display=params.EnumParam.DISPLAY_LIST,
                       label='Display volume with',
                       help='*slices*: display volumes as 2D slices along z axis.\n'
                            '*chimera*: display volumes as surface with Chimera.')

        group = form.addGroup('Iterations')
        self._defineIterationParams(form, group)
        
    def _defineIterationParams(self, form, group):
        """ Overwrite in base classes to define the 
        visuliazation parameters depending on selected iteration.
        """
        pass
        
    def _viewVol(self, volPath):
        if self.displayVol == VOLUME_CHIMERA:
            view = pwviewer.CommandView('%s/bin/chimera %s' % (os.environ['CHIMERA_HOME'], volPath))
        else:
            view = self.getObjectView(volPath)
        return view
    
    def getIter(self):
        return int(self.displayIter.get())
    

class PyTomFrmAlignmentViewer(PyTomViewer):
    """
    Visualization tools for PyTom autofocus classification.
    """
    _targets = [ProtFrmAlign]
    _label = 'viewer autofocus'
    
    def __init__(self, *args, **kwargs):
        PyTomViewer.__init__(self, *args, **kwargs)
    
    def _defineIterationParams(self, form, group):
        group.addParam('displayIter', params.StringParam, default='0',
                       label='Display iteration')
        group.addParam('showAverageMap', params.LabelParam,
                      label='Show average map', important=True)
        group.addParam('showAlignedMaps', params.LabelParam,
                      label='Show aligned particles (volumes)', important=True)        
        
    def _getVisualizeDict(self):
        return {'showAverageMap': self._showAverageMap,
                'showAlignedMaps': self._showAlignedMaps
                }
        
    def _showAverageMap(self, param=None):
        #view = CommandView('xmipp_chimera_client --input "%s" --mode projector 256 &' % volPath)
        volPath = self.protocol.getAverageMap(self.getIter())
        return [self._viewVol(volPath)]
        
    def _showAlignedMaps(self, param=None):
        #view = CommandView('xmipp_chimera_client --input "%s" --mode projector 256 &' % volPath)
        volsSqlite = self.protocol.getVolumesSqlite(self.getIter())
        viewParams = {showj.MODE: showj.MODE_MD}
        
        return [self.getObjectView(volsSqlite, 
                                   other=self.protocol.inputVolumes.get().getObjId(),
                                   viewParams=viewParams)]   
    

class PyTomAutofocusViewer(PyTomViewer):
    """
    Visualization tools for PyTom autofocus classification.
    """
    _environments = [pwviewer.DESKTOP_TKINTER, pwviewer.WEB_DJANGO]
    _targets = [ProtAutofocusClassify]
    _label = 'viewer autofocus'
    
    def __init__(self, *args, **kwargs):
        PyTomViewer.__init__(self, *args, **kwargs)
    
    def _defineParams(self, form):
        form.addSection(label='Visualization')
        group = form.addGroup('General')
        group.addParam('displayVol', params.EnumParam, 
                       choices=['slices', 'chimera'],
                       default=VOLUME_CHIMERA, 
                       display=params.EnumParam.DISPLAY_LIST,
                       label='Display volume with',
                       help='*slices*: display volumes as 2D slices along z axis.\n'
                            '*chimera*: display volumes as surface with Chimera.')

        group = form.addGroup('Iterations')
        group.addParam('displayIter', params.StringParam, default='0',
                       label='Display iteration')
        group.addParam('showDiffMap', params.LabelParam,
                      label='Show difference maps between classes')
        line = group.addLine('Classes')
        line.addParam('classA', params.IntParam, label='A', default=0)
        line.addParam('classB', params.IntParam, label='B', default=1)
        
        form.addParam('showClassMap', params.IntParam, default=0,
                      label='Show class map?')
        
        
    def _getVisualizeDict(self):
        return {'showDiffMap': self._showDiffMaps,
                'showClassMap': self._showClassMap
                }
        
    def _viewVol(self, volPath):
        if self.displayVol == VOLUME_CHIMERA:
            view = pwviewer.CommandView('%s/bin/chimera %s' % (os.environ['CHIMERA_HOME'], volPath))
        else:
            view = DataView(volPath)
        return view
        
        
    def _showDiffMaps(self, param=None):
        #view = CommandView('xmipp_chimera_client --input "%s" --mode projector 256 &' % volPath)
        volPath = self.protocol.getDiffMap(int(self.displayIter.get()),
                                           self.classA.get(),
                                           self.classB.get())
        return [self._viewVol(volPath)]
    
    def _showClassMap(self, param=None):
        volPath = self.protocol.getClassMap(int(self.displayIter.get()),
                                            int(self.showClassMap.get()))
        return [self._viewVol(volPath)]
        
