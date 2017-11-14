# **************************************************************************
# *
# * Authors:  Roberto Marabini (roberto@cnb.csic.es), May 2013
# *           Marta Martinez (mmmtnez@cnb.csic.es)
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

from pyworkflow.viewer import DESKTOP_TKINTER, WEB_DJANGO, ProtocolViewer
import pyworkflow.protocol.params as params
from pyworkflow.em.packages.chimera.convert import symMapperScipionchimera
from protocol_extract_unit_cell import XmippProtExtractUnit
from pyworkflow.em.constants import SYM_I222
from distutils.spawn import find_executable
import pyworkflow.em as em
import os
from os.path import exists
from pyworkflow.em.data import (SetOfVolumes)

VOLUME_SLICES = 1
VOLUME_CHIMERA = 0


class viewerXmippProtExtractUnit(ProtocolViewer):
    """ Visualize the input and output volumes of protocol XmippProtExtractUnit
        by choosing Chimera (3D) or Xmipp visualizer (2D).
        The axes of coordinates x, y, z will be shown by choosing Chimera"""
    _label = 'viewer extract unit cell'
    _targets = [XmippProtExtractUnit]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]

    # ROB: I know that there is a nice chimera interface but it does not work
    # in this case since I am interested in reading the MRC header. So I will
    # use chimera as an external program

    def _defineParams(self, form):
        form.addSection(label='Visualization of input volume and extracted '
                              'unit cell')
        form.addParam('displayVol', params.EnumParam,
                      choices=['chimera', 'slices'], default=VOLUME_CHIMERA,
                      display=params.EnumParam.DISPLAY_HLIST,
                      label='Display volume with',
                      help='*chimera*: display volumes as surface with '
                           'Chimera.\n*slices*: display volumes as 2D slices '
                           'along z axis.\n')

    def _getVisualizeDict(self):
        return{
            'displayVol': self._showVolumes,
        }

    def _validate(self):
        if find_executable('chimera') is None:
            return ["chimera is not available. Either install it or choose"
                    " option 'slices'. "]
        return []

    # =========================================================================
    # Show Volumes
    # =========================================================================

    def _showVolumes(self, paramName=None):
        if self.displayVol == VOLUME_CHIMERA:
            return self._showVolumesChimera()
        elif self.displayVol == VOLUME_SLICES:
            return self._showVolumesXmipp()

    def _getVolumeName(self, vol):
        volName = vol.getFileName().replace(':mrc', '')
        return volName

    def _getOrigin(self, vol):
        origin = vol.getOrigin(returnInitIfNone=True).getShifts()
        x = int(origin[0])
        y = int(origin[1])
        z = int(origin[2])
        return x, y, z

    def _createSetOfVolumes(self):
        if not exists(self.protocol._getTmpPath() + '/tmpVolumes.sqlite'):
            tmpFileName = self.protocol._getTmpPath("tmpVolumes.sqlite")
            _inputVol = self.protocol.inputVolumes.get()
            _outputVol = self.protocol.outputVolume
            setOfVolumes = SetOfVolumes(filename=tmpFileName)
            setOfVolumes.append(_inputVol)
            setOfVolumes.append(_outputVol)
            setOfVolumes.write()
        else:
            tmpFileName = self.protocol._getTmpPath()+'/tmpVolumes.sqlite'
            setOfVolumes = SetOfVolumes(filename=tmpFileName)
        return setOfVolumes

    def _showVolumesChimera(self):
        tmpFileNameCMD = self.protocol._getTmpPath("chimera.cmd")
        f = open(tmpFileNameCMD, "w")
        dim = self.protocol.inputVolumes.get().getDim()[0]
        sampling = self.protocol.inputVolumes.get().getSamplingRate()
        tmpFileNameBILD = self.protocol._getTmpPath("axis.bild")
        tmpFileNameBILD = os.path.abspath(tmpFileNameBILD)
        f.write("open %s\n" % tmpFileNameBILD)

        ff = open(tmpFileNameBILD, "w+")
        arrowDict = {}
        arrowDict["x"] = arrowDict["y"] = arrowDict["z"] = \
            sampling * dim * 3. / 4.
        arrowDict["r1"] = 0.1  # sampling * dim / 600.
        arrowDict["r2"] = 4 * arrowDict["r1"]
        arrowDict["rho"] = 0.75  # sampling * dim / 150.

        ff.write(".color 1 0 0\n"
                 ".arrow 0 0 0 %(x)d 0 0 %(r1)f %(r2)f %(rho)f\n"
                 ".color 1 1 0\n"
                 ".arrow 0 0 0 0 %(y)d 0 %(r1)f %(r2)f %(rho)f\n"
                 ".color 0 0 1\n"
                 ".arrow 0 0 0 0 0 %(z)d %(r1)f %(r2)f %(rho)f\n" %
                 arrowDict)
        ff.close()
        _inputVol = self.protocol.inputVolumes.get()
        _outputVol = self.protocol.outputVolume
        inputVolFileName = os.path.abspath(self._getVolumeName(_inputVol))
        # input vol origin coordinates
        x_input, y_input, z_input = self. _getOrigin(_inputVol)
        f.write("open %s\n" % inputVolFileName)
        f.write("volume #1 style mesh level 0.001 voxelSize %f originIndex "
                "%d,%d,%d\n"
                % (_inputVol.getSamplingRate(), x_input, y_input, z_input))

        outputVolFileName = os.path.abspath(self._getVolumeName(_outputVol))
        # output vol origin coordinates
        x_output, y_output, z_output = self. _getOrigin(_outputVol)
        f.write("open %s\n" % outputVolFileName)
        f.write("volume #2 style surface level 0.001 voxelSize %f originIndex "
                "%d,%d,%d\n"
                % (_outputVol.getSamplingRate(), x_output, y_output, z_output))

        cMap = ['red', 'yellow', 'green', 'cyan', 'blue']
        d = {}
        d['outerRadius'] = self.protocol.outerRadius.get() * sampling
        d['innerRadius'] = self.protocol.innerRadius.get() * sampling
        d['symmetry'] = \
            symMapperScipionchimera[self.protocol.symmetryGroup.get()]

        if self.protocol.symmetryGroup >= SYM_I222:
            f.write("shape icosahedron mesh true radius %(outerRadius)d "
                    "orientation %(symmetry)s\n" % d)
        step = (d['outerRadius'] - d['innerRadius']) / float(len(cMap) - 1)
        f.write("scolor #2  geom radial center 0,0,0 cmap ")
        counter = 0
        s = ""
        for color in cMap:
            s += "%d,%s:" % (d['innerRadius'] + counter * step, color)
            counter += 1
        f.write(s[:-1] + '\n')

        f.close()

        return [em.ChimeraView(tmpFileNameCMD)]

    def _showVolumesXmipp(self):

        setOfVolumes = self._createSetOfVolumes()

        return [self.objectView(setOfVolumes)]
