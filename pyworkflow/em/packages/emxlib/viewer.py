# **************************************************************************
# *
# * Authors:     J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Roberto Marabini       (roberto@cnb.csic.es)
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
"""
This module implement the wrappers around xmipp_showj
visualization program.
"""
from pyworkflow.viewer import TextView
from pyworkflow.em.viewer import DataView
from pyworkflow.viewer import Viewer, DESKTOP_TKINTER, WEB_DJANGO
from dataimport import ProtEmxExport


class EMXViewer(Viewer):
    """ Class to visualize Relion protocols """
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [ProtEmxExport]
    _label = 'viewer emx'
    
    def __init__(self, **args):
        Viewer.__init__(self, **args)

    def visualize(self, obj, **args):
        data = self.protocol._getPath('emxData/data.mrc')
        DataView(data).show()
        
        # text view doesn't work now.
        
#         fn = self.protocol._getPath('emxData/data.emx')
#         self.textView([fn])    
