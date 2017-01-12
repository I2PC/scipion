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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
from pyworkflow.em.packages.xmipp3.convert import readSetOfClasses
"""
This module implement the wrappers aroung Xmipp CL2D protocol
visualization program.
"""
from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO
from pyworkflow.em import *
from protocol_split_volume import XmippProtSplitvolume
from pyworkflow.gui.text import *
from pyworkflow.gui.dialog import showError, showWarning
from pyworkflow.protocol.params import LabelParam, LEVEL_ADVANCED
from pyworkflow.em.packages.xmipp3.viewer import XmippViewer
import glob

class XmippViewerSplitVolume(XmippViewer):
    """ Wrapper to visualize different type of data objects
    with the Xmipp program xmipp_showj
    """
    _label = 'viewer split volume'
    _targets = [XmippProtSplitvolume]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    
    def __init__(self, **args):
        XmippViewer.__init__(self, **args)

    def _visualize(self, obj, **args):
        if hasattr(obj, 'outputVolumes'):
            XmippViewer._visualize(self,self.protocol.outputVolumes)
        fnBasis=self.protocol._getExtraPath('split_pc1.vol')
        if exists(fnBasis):
            self._views.append(DataView(fnBasis))
