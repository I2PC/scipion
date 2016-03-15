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
This module implements viewers for Spider protocols.
"""

from pyworkflow.viewer import DESKTOP_TKINTER, WEB_DJANGO
from pyworkflow.em.viewer import DataView, ObjectView
from pyworkflow.em.packages.xmipp3.viewer import XmippViewer

from spider import PcaFile
from protocol.protocol_custommask import SpiderProtCustomMask


    
class SpiderViewer(XmippViewer):
    """ Wrapper to visualize different type of objects
    with the Xmipp program xmipp_showj. """
    
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [PcaFile, SpiderProtCustomMask]
    _label = 'viewer'

    def _visualize(self, obj, **args):
        
        if isinstance(obj, PcaFile):
            self._views.append(self.textView([obj.getFileName()], "PCA file"))
        
        elif isinstance(obj, SpiderProtCustomMask):
            mask = obj.outputMask
            self._visualize(mask)
            # Visualize the whole stack with all filenames
            self._views.append(DataView(mask.getFileName()))
            
        else:
            XmippViewer._visualize(self, obj)
           
        return self._views
