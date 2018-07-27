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
from protocol_superpose_pdbs import PhenixProtRunSuperposePDBs
from pyworkflow.em.viewers.chimera_utils import \
    createCoordinateAxisFile, runChimeraProgram, getProgram
from pyworkflow.viewer import DESKTOP_TKINTER, Viewer


class PhenixProtRunSuperposePDBsViewer(Viewer):
    """ Visualize the output of protocols protocol_fit and protocol_operate """
    _environments = [DESKTOP_TKINTER]
    _label = 'Superpose PDBs viewer'
    _targets = [PhenixProtRunSuperposePDBs]

    def _visualize(self, obj, **args):
        # THe input map or pdb may be a parameter from the protocol
        # or from the parent protocol.
        fnCmd = self.protocol._getTmpPath("chimera_output.cmd")
        vol = self.protocol.inputStructureFixed.get().getVolume()
        if vol is not None:
            dim = vol.getDim()[0]
            sampling = vol.getSamplingRate()
            self.volFileName = vol.getFileName()
            if vol.hasOrigin():
                self.x, self.y, self.z = vol.getOrigin().getShifts()
            else:
                self.x, self.y, self.z = vol.getOrigin(force=True).getShifts()
        else:
            # To show pdbs only
            dim = 150.
            sampling = 1.

        bildFileName = os.path.abspath(self.protocol._getTmpPath(
            "axis_output.bild"))
        createCoordinateAxisFile(dim,
                                 bildFileName=bildFileName,
                                 sampling=sampling)

        pdbList = []
        outputPdb = self.protocol.outputPdb.getFileName()
        pdbList.append(outputPdb)
        inputPdbFixed = self.protocol.inputStructureFixed.get().getFileName()
        pdbList.append(inputPdbFixed)
        inputPdbMoving = self.protocol.inputStructureMoving.get().getFileName()
        pdbList.append(inputPdbMoving)

        with open(fnCmd, 'w') as f:
            f.write("open %s\n" % bildFileName)
            if vol is not None:
                f.write("open %s\n" % self.volFileName)
                f.write("volume #1 style surface voxelSize %f\n"
                        "volume #1 origin %0.2f,%0.2f,%0.2f\n"
                        % (vol.getSamplingRate(), self.x, self.y, self.z))
            for filename in pdbList:
                f.write("open %s\n" % os.path.abspath(filename))

        # run in the background
        runChimeraProgram(getProgram(), fnCmd+"&")
        return []


