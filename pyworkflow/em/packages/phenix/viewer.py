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

import os
from distutils.spawn import find_executable
from os.path import exists

import pyworkflow.em as em
import pyworkflow.protocol.params as params
from pyworkflow.em.convert import ImageHandler
from pyworkflow.em.data import (SetOfVolumes, Volume)
from pyworkflow.em.viewers.chimera_utils import \
    createCoordinateAxisFile, getProgram
from pyworkflow.viewer import DESKTOP_TKINTER, WEB_DJANGO, ProtocolViewer
from protocol_automated_sharpening import PhenixProtAutomatedSharpening

VOLUME_SLICES = 1
VOLUME_CHIMERA = 0


class viewerPhenixAutoSharpen(ProtocolViewer):
    """ Visualize the input and output volumes of protocol XmippProtExtractUnit
        by choosing Chimera (3D) or Xmipp visualizer (2D).
        The axes of coordinates x, y, z will be shown by choosing Chimera"""
    _label = 'viewer extract unit cell'
    _targets = [PhenixProtAutomatedSharpening]
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
        if find_executable(getProgram()) is None:
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



    def _createSetOfVolumes(self):
        if not exists(self.protocol._getTmpPath('tmpVolumes.sqlite')):
            tmpFileName = self.protocol._getTmpPath("tmpVolumes.sqlite")
            _outputVol = self.protocol.outputVol
            setOfVolumes = SetOfVolumes(filename=tmpFileName)
            setOfVolumes.append(_outputVol)
            if _outputVol.hasHalfMaps():
                origin = _outputVol.getOrigin(force=True)
                sampling = _outputVol.getSamplingRate()

                for name in _outputVol.getHalfMaps():                
                    vol = Volume()
                    vol.setFileName(name)
                    vol.setSamplingRate(sampling) 
                    vol.setOrigin(origin)
                    setOfVolumes.append(vol)
            setOfVolumes.write()
        else:
            tmpFileName = self.protocol._getTmpPath('tmpVolumes.sqlite')
            setOfVolumes = SetOfVolumes(filename=tmpFileName)
        return setOfVolumes

    def _showVolumesChimera(self):
        _setOfVolumes = self._createSetOfVolumes()
        tmpFileNameCMD = self.protocol._getTmpPath("chimera.cmd")
        f = open(tmpFileNameCMD, "w")
        dim = self.protocol.inputMap.get().getDim()[0]
        sampling = self.protocol.inputMap.get().getSamplingRate()
        if len(_setOfVolumes)==1:        
            tmpFileName = os.path.abspath(self.protocol._getTmpPath("axis.bild"))
            createCoordinateAxisFile(dim,
                                     bildFileName=tmpFileName,
                                     sampling=sampling)
            f.write("open %s\n" % tmpFileName)
    
            _outputVol = self.protocol.outputVol
            outputVolFileName = os.path.abspath(ImageHandler.removeFileType(
                _outputVol.getFileName()))
    
            # output vol origin coordinates
            x_output, y_output, z_output = _outputVol.getVolOriginAsTuple()
            f.write("open %s\n" % outputVolFileName)
            f.write("volume #1 style surface voxelSize %f\nvolume #1 origin "
                    "%0.2f,%0.2f,%0.2f\n"
                    % (_outputVol.getSamplingRate(), x_output, y_output, z_output))
        else:
            count = 0
            for vol in _setOfVolumes:
                localVol = os.path.abspath(ImageHandler.removeFileType(
                    vol.getFileName()))
                f.write("open %s\n" % localVol)
                f.write("volume #%d style surface voxelSize %f\n" %
                        (count, sampling))
                count += 1
            f.write("tile\n")
        f.close()

        return [em.ChimeraView(tmpFileNameCMD)]

    def _showVolumesXmipp(self):

        setOfVolumes = self._createSetOfVolumes()

        return [self.objectView(setOfVolumes)]
