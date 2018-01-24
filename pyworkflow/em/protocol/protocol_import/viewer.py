# **************************************************************************
# *
# * Authors:     Marta Martinez (mmmtnez@cnb.csic.es)
# *              Roberto Marabini (roberto@cnb.csic.es)
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
This module implements visualization program
for input volumes.
"""

import os
from distutils.spawn import find_executable
from tkMessageBox import showerror

import pyworkflow.protocol.params as params
from pyworkflow.em.convert import ImageHandler
from pyworkflow.em.protocol.protocol_import.volumes import ProtImportVolumes
from pyworkflow.viewer import DESKTOP_TKINTER, WEB_DJANGO, ProtocolViewer
from pyworkflow.em.utils.chimera_utilities.convert import \
    createCoordinateAxisFile,  adaptOriginFromCCP4ToChimera, getProgram
VOLUME_SLICES = 1
VOLUME_CHIMERA = 0

class viewerProtImportVolumes(ProtocolViewer):
    """ Wrapper to visualize different type of objects
    with the Xmipp program xmipp_showj. """

    _label = 'viewer input volume'
    _targets = [ProtImportVolumes]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]

    def _defineParams(self, form):
        form.addSection(label='Visualization of input volumes')
        form.addParam('displayVol', params.EnumParam,
                      choices=['chimera', 'slices'],
                      default=VOLUME_CHIMERA,
                      display=params.EnumParam.DISPLAY_HLIST,
                      label='Display volume with',
                      help='*chimera*: display volumes as surface with '
                           'Chimera.\n *slices*: display volumes as 2D slices '
                           'along z axis.\n If number of volumes == 1, '
                           'a system of coordinates is shown'
                      )

    def _getVisualizeDict(self):
        return {
            'displayVol': self._showVolumes,
        }

    def _validate(self):
        if find_executable(getProgram()) is None:
            return ["chimera is not available. "
                    "Either install it or choose option 'slices'. "]
        return []

    # =========================================================================
    # ShowVolumes
    # =========================================================================

    def _showVolumes(self, paramName=None):
        if self.displayVol == VOLUME_CHIMERA:
            return self._showVolumesChimera()

        elif self.displayVol == VOLUME_SLICES:
            return self._showVolumesSlices()

    def _createSetOfVolumes(self):
        try:
            setOfVolumes = self.protocol.outputVolumes
            sampling = self.protocol.outputVolumes.getSamplingRate()
        except:
            setOfVolumes = self.protocol._createSetOfVolumes()
            setOfVolumes.append(self.protocol.outputVolume)
            sampling = self.protocol.outputVolume.getSamplingRate()

        return sampling, setOfVolumes

    def _showVolumesChimera(self):
        """ Create a chimera script to visualize selected volumes. """
        tmpFileNameCMD = self.protocol._getTmpPath("chimera.cmd")
        f = open(tmpFileNameCMD, "w")
        sampling, _setOfVolumes = self._createSetOfVolumes()
        count = 0  # first model in chimera is a volume

        if len(_setOfVolumes) == 1:
            count = 1  # first model in chimera is the bild file
            # if we have a single volume then create axis
            # as bild file. Chimera must read the bild file first
            # otherwise system of coordinates will not
            # be in the center of the window

            dim = self.protocol.outputVolume.getDim()[0]
            tmpFileNameBILD = os.path.abspath(self.protocol._getTmpPath(
                "axis.bild"))
            createCoordinateAxisFile(dim,
                                     bildFileName=tmpFileNameBILD,
                                     sampling=sampling)
            f.write("open %s\n" % tmpFileNameBILD)
            count = 1  # skip first model because is not a 3D map

        for vol in _setOfVolumes:
            localVol = os.path.abspath(ImageHandler.removeFileType(
                vol.getFileName()))
            if localVol.endswith("stk"):
                errorWindow(None, "Extension .stk is not supported")
            f.write("open %s\n" % localVol)
            f.write("volume#%d style surface voxelSize %f\n" %
                    (count, sampling))
            count += 1

        if len(_setOfVolumes) > 1:
            f.write('tile\n')
        else:
            x, y, z = adaptOriginFromCCP4ToChimera(vol.getVolOriginAsTuple())
            f.write("volume#1 origin %0.2f,%0.2f,%0.2f\n" % (x, y, z))
        f.close()
        import pyworkflow.em as em
        return [em.ChimeraView(tmpFileNameCMD)]

    def _showVolumesSlices(self):
        # Write an sqlite with all volumes selected for visualization.
        sampling, setOfVolumes = self._createSetOfVolumes()

        if len(setOfVolumes) == 1:
            return [self.objectView(self.protocol.outputVolume)]

        return [self.objectView(setOfVolumes)]


def errorWindow(tkParent, msg):
    try:
        showerror("Error",  # bar title
                  msg,  # message
                  parent=tkParent)
    except:
        print("Error:", msg)

