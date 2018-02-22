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

from pyworkflow.em.convert import ImageHandler
from protocol_fit import ChimeraProtRigidFit
from protocol_operate import ChimeraProtOperate
from pyworkflow.em.viewers.chimera_utils import \
    createCoordinateAxisFile, \
    adaptOriginFromCCP4ToChimera, runChimeraProgram, \
    getProgram, chimeraPdbTemplateFileName
from pyworkflow.viewer import DESKTOP_TKINTER, Viewer


class ChimeraProtRigidFitViewer(Viewer):
    """ Visualize the output of protocols protocol_fit and protocol_operate """
    _label = 'viewer fit'
    _targets = [ChimeraProtRigidFit, ChimeraProtOperate]
    _environments = [DESKTOP_TKINTER]

    def _visualize(self, obj, **args):
        # Construct the coordinate file and visualization
        _inputVol = self.protocol.inputVolume.get()
        if _inputVol is None:
            _inputVol = self.protocol.pdbFileToBeRefined.get().getVolume()

        dim = _inputVol.getDim()[0]
        sampling = _inputVol.getSamplingRate()

        bildFileName = os.path.abspath(self.protocol._getTmpPath(
            "axis_output.bild"))
        createCoordinateAxisFile(dim,
                                 bildFileName=bildFileName,
                                 sampling=sampling)
        fnCmd = self.protocol._getTmpPath("chimera_output.cmd")
        f = open(fnCmd, 'w')
        f.write("open %s\n" % bildFileName)

        try:
            outputVol = self.protocol.output3Dmap
            outputVolFileName = os.path.abspath(outputVol.getFileName())
        except:
            outputVol =  _inputVol
            outputVolFileName = os.path.abspath(
                    ImageHandler.removeFileType(outputVol.getFileName()))
        f.write("open %s\n" % outputVolFileName)
        x, y, z = outputVol.getOrigin().getShifts()
        f.write("volume #1 style surface voxelSize %f origin "
                "%0.2f,%0.2f,%0.2f\n"
                % (outputVol.getSamplingRate(), x, y, z))

        directory = self.protocol._getExtraPath()
        for filename in os.listdir(directory):
            if filename.endswith(".pdb"):
                path = os.path.join(directory, filename)
                f.write("open %s\n" % os.path.abspath(path))

        f.close()

        # run in the background
        runChimeraProgram(getProgram(), fnCmd+"&")
        return []
