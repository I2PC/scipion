# **************************************************************************
# *
# * Authors:  Carlos Oscar Sanchez Sorzano (coss@cnb.csic.es)
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

from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO
from pyworkflow.em.viewer import DataView, ChimeraView
from pyworkflow.em.packages.xmipp3.viewer import XmippViewer
from pyworkflow.protocol.params import PointerParam, LabelParam

from protocol_normalize_strain import XmippProtNormalizeStrain

class XmippNormalizeStrainViewer(ProtocolViewer):
    """ Visualize the output of protocol volume strain """
    _label = 'viewer normalize strain'
    _targets = [XmippProtNormalizeStrain]
    _environments = [DESKTOP_TKINTER]
    
    def __init__(self, **args):
        ProtocolViewer.__init__(self, **args)

    def _defineParams(self, form):
        form.addSection(label='Visualization')
        # Select the level to show
        form.addParam('show', LabelParam, 
                      label="Strain calculation to visualize")     
        form.addParam('protToShow', PointerParam, pointerClass='XmippProtVolumeStrain',
                      label="Select protocol")     

    def _getVisualizeDict(self):
        return {'show': self._viewStrain}        

    def _viewStrain(self, e=None):
        import os
        protToShow =  self.protToShow.get()
        protId = protToShow.getObjId()
        
        views=[]
        fnCmd = self.protocol._getPath('%d_result_strain_chimera.cmd'%protId)
        if os.path.exists(fnCmd):
            views.append(ChimeraView(fnCmd))
        fnCmd = self.protocol._getPath('%d_result_localrot_chimera.cmd'%protId)
        if os.path.exists(fnCmd):
            views.append(ChimeraView(fnCmd))
        fnCmd = protToShow._getPath('result_morph_chimera.cmd')
        if os.path.exists(fnCmd):
            views.append(ChimeraView(fnCmd))
        return views

