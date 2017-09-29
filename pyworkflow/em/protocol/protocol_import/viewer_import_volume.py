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
from pyworkflow.em.protocol.protocol_import.volumes import ProtImportVolumes
import pyworkflow.protocol.params as params
from pyworkflow.viewer import DESKTOP_TKINTER, WEB_DJANGO, ProtocolViewer
from distutils.spawn import find_executable
from os.path import exists


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
                       default=VOLUME_CHIMERA, display=params.EnumParam.DISPLAY_HLIST,
                       label='Display volume with',
                       help='*chimera*: display volumes as surface with Chimera.\n '
                            '*slices*: display volumes as 2D slices along z axis.\n'
                            'If number of volumes == 1, a system of coordinates is shown'
                      )


    def _getVisualizeDict(self):
        return {
            'displayVol': self._showVolumes,
        }

    def _validate(self):
        if find_executable("chimera") is None:
            return ["chimera is not available. Either install it or choose option 'slices'. "]
        return []

    # ===============================================================================
    # ShowVolumes
    # ===============================================================================

    def _showVolumes(self, paramName=None):
        if self.displayVol == VOLUME_CHIMERA:
            return self._showVolumesChimera()

        elif self.displayVol == VOLUME_SLICES:
            return self._showVolumesXmipp()

    def _getVolumeName(self, vol):
        volName = vol.getFileName().replace(':mrc', '')
        if not exists(volName):
            print "Volume %s does not exist" % volName
        else:
            return volName

    def _getVolOrigin(self, vol):
        """chimera only handles integer shifts"""
        origin = vol.getOrigin(returnInitIfNone=True).getShifts()
        x = int(origin[0])  # * sampling
        y = int(origin[1])  # * sampling
        z = int(origin[2])  # * sampling
        return x, y, z

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

        if len(_setOfVolumes) ==1:
            count = 1 # first model in chimera is the bild file
            dim = self.protocol.outputVolume.getDim()[0]
            # if we have a single volume then create axis
            # as bild file. Chimera must read the bild file first
            # otherwise system of coordinates will not
            # be in the center of the window
            tmpFileNameBILD = self.protocol._getTmpPath("axis.bild")
            tmpFileNameBILD = os.path.abspath(tmpFileNameBILD)
            f.write("open %s\n" % tmpFileNameBILD)
            ff = open(tmpFileNameBILD, "w+")
            arrowDict = {}
            arrowDict["x"] = arrowDict["y"] = arrowDict["z"] = \
                sampling * dim * 3. / 4.
            arrowDict["r1"] = 0.1 #sampling * dim / 600.
            arrowDict["r2"] = 4 * arrowDict["r1"]
            arrowDict["rho"] = 0.75 #sampling * dim / 150.
            print arrowDict
            ff.write(".color 1 0 0\n"
                    ".arrow 0 0 0 %(x)d 0 0 %(r1)f %(r2)f %(rho)f\n"
                    ".color 1 1 0\n"
                    ".arrow 0 0 0 0 %(y)d 0 %(r1)f %(r2)f %(rho)f\n"
                    ".color 0 0 1\n"
                    ".arrow 0 0 0 0 0 %(z)d %(r1)f %(r2)f %(rho)f\n" %
                     arrowDict)
            ff.close()
            count = 1  # skip first model because is not a 3D map

        for vol in _setOfVolumes:
            localVol = os.path.abspath(self._getVolumeName(vol))
            x, y, z = self._getVolOrigin(vol)
            f.write("open %s\n" % localVol)
            f.write("volume#%d style surface voxelSize %f\n" %
                    (count,sampling))
            count += 1
        if len(_setOfVolumes) > 1:
            f.write('tile\n')
        else:
            f.write("volume#1 originIndex %d,%d,%d\n" % (x, y, z))
        f.close()
        import pyworkflow.em as em
        return [em.ChimeraView(tmpFileNameCMD)]

    def _showVolumesXmipp(self):
        #Write an sqlite with all volumes selected for visualization.

        sampling, setOfVolumes = self._createSetOfVolumes()

        if len(setOfVolumes) == 1:
            return [self.objectView(self.protocol.outputVolume)]

        return [self.objectView(setOfVolumes)]
