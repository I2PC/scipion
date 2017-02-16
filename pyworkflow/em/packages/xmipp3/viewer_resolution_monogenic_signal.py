# **************************************************************************
# *
# * Authors:     J.L. Vilas (jlvilas@cnb.csic.es)
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

from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER
from protocol_resolution_monogenic_signal import (XmippProtMonoRes,
                                                  OUTPUT_RESOLUTION_FILE,
                                                  OUTPUT_RESOLUTION_FILE_CHIMERA) 
from pyworkflow.em.viewer_local_resolution import localResolutionViewer

class XmippMonoResViewer(localResolutionViewer):
    """
    Visualization tools for MonoRes results.
    
    MonoRes is a Xmipp packagefor computing the local resolution of 3D
    density maps studied in structural biology, primarily by cryo-electron
    microscopy (cryo-EM).
    """
    _label = 'viewer MonoRes'
    _targets = [XmippProtMonoRes]      
    _environments = [DESKTOP_TKINTER]

    def __init__(self, **kwargs):
        ProtocolViewer.__init__(self, **kwargs)
        localResolutionViewer.OUTPUT_RESOLUTION_FILE = OUTPUT_RESOLUTION_FILE
        localResolutionViewer.OUTPUT_RESOLUTION_FILE_CHIMERA = OUTPUT_RESOLUTION_FILE_CHIMERA
        localResolutionViewer.halves = self.protocol.halfVolumes.get()
        localResolutionViewer.backgroundValue = 0

