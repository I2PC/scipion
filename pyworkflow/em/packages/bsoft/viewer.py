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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************


from pyworkflow.em.viewer import CommandView, Viewer, DESKTOP_TKINTER
from pyworkflow.viewer import ProtocolViewer
from pyworkflow.em.data import Volume
from protocol_blocres import BsoftProtBlocres
from convert import getEnviron
from pyworkflow.em.viewer_local_resolution import localResolutionViewer


#------------------------ Some views and  viewers ------------------------
        

class BsoftVolumeView(CommandView):
    def __init__(self, inputFile, **kwargs):

        # Small trick to handle .vol Spider volumes
        if inputFile.endswith('.vol'):
            inputFile += ":spi"

        CommandView.__init__(self, 'bshow "%s" &' % inputFile,
                             env=getEnviron(), **kwargs)

             
class BsoftViewer(Viewer):
    _environments = [DESKTOP_TKINTER]
    _targets = [Volume]
    
    def __init__(self, **kwargs):
        Viewer.__init__(self, **kwargs)

    def visualize(self, obj, **kwargs):
        cls = type(obj)

        fn = obj.getFileName()
        BsoftVolumeView(fn).show()


class BsoftViewerBlocres(localResolutionViewer):
    """
    Visualization tools for blocres results.

    blocres is a bsoft package for computing the local resolution of 3D
    density maps studied in structural biology, primarily by cryo-electron
    microscopy (cryo-EM).
    """
    _label = 'viewer blocres'
    _targets = [BsoftProtBlocres]
    _environments = [DESKTOP_TKINTER]
    

    def __init__(self, **kwargs):
        ProtocolViewer.__init__(self, **kwargs)
        from protocol_blocres import OUTPUT_RESOLUTION_FILE
        localResolutionViewer.OUTPUT_RESOLUTION_FILE = OUTPUT_RESOLUTION_FILE
        localResolutionViewer.OUTPUT_RESOLUTION_FILE_CHIMERA = OUTPUT_RESOLUTION_FILE
        localResolutionViewer.halves = True
        localResolutionViewer.backgroundValue = self.protocol.fill.get()




