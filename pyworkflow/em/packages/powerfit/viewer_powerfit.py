# **************************************************************************
# *
# * Authors:     Mohsen Kazemi  (mkazemi@cnb.csic.es)
# *              C.O.S. Sorzano (coss@cnb.csic.es)
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os
from pyworkflow.viewer import DESKTOP_TKINTER, WEB_DJANGO, ProtocolViewer
from pyworkflow.em.viewer import ChimeraView, ObjectView
from pyworkflow.em.showj import MODE, MODE_MD, ORDER, VISIBLE, SORT_BY
from protocol_powerfit import PowerfitProtRigidFit
import pyworkflow.protocol.params as params

class PowerfitProtRigidFitViewer(ProtocolViewer):
    """ Wrapper to visualize powerfit results
    """
    
    _label = 'viewer validate_overfitting'
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [PowerfitProtRigidFit]
        
    def _defineParams(self, form):
        form.addSection(label='Show Powerfit')
        form.addParam('doShowFitTable', params.LabelParam, label="Display fitting quality")
        form.addParam('modelNumber', params.IntParam, default=1, label="Model to visualize")
        form.addParam('doShowFit', params.LabelParam, label="Display the fitting")
    
    def _getVisualizeDict(self):
        return {'doShowFit': self._visualizeFit, 'doShowFitTable': self._visualizeFitTable}
    
    def _visualizeFit(self, e=None):
        import os
        fnCmd = self.protocol._getExtraPath('chimera_%d.cmd'%self.modelNumber)
        if os.path.exists(fnCmd):
            return [ChimeraView(fnCmd)]

    def _visualizeFitTable(self, e=None):
        views = []
        if hasattr(self.protocol, "outputPDBs"):
            labels = 'id _filename _powerfit_cc _powerfit_Fish_z _powerfit_rel_z'
            views.append(ObjectView(self._project,
                                    self.protocol.outputPDBs.strId(), 
                                    self.protocol.outputPDBs.getFileName(),
                                    viewParams={MODE: MODE_MD, ORDER: labels, VISIBLE: labels, SORT_BY:"_powerfit_cc desc"}))
        return views 
    