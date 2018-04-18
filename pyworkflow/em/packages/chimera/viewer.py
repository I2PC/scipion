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
from protocol_restore import ChimeraProtRestore
from pyworkflow.em.viewers.chimera_utils import \
    createCoordinateAxisFile, runChimeraProgram, getProgram, sessionFile
from pyworkflow.viewer import DESKTOP_TKINTER, Viewer


class ChimeraViewerBase(Viewer):
    """ Visualize the output of protocols protocol_fit and protocol_operate """
    _environments = [DESKTOP_TKINTER]

    def _visualize(self, obj, **args):
        # THe input map or pdb may be a parameter from the protocol
        # or from the parent protocol.
        try:
            _inputVol = self.protocol.inputVolume.get()
            directory = self.protocol._getExtraPath()
        except:
            _inputVol = self.protocol.inputProtocol.get().inputVolume.get()
            directory = self.protocol.inputProtocol.get()._getExtraPath()

        if _inputVol is None:
            try:
                _inputVol = self.protocol.pdbFileToBeRefined.get().getVolume()
            except:
                _inputVol = self.protocol.inputProtocol.get().\
                    pdbFileToBeRefined.get().getVolume()

        if _inputVol is not None:
            dim = _inputVol.getDim()[0]
            sampling = _inputVol.getSamplingRate()
            _showVol = _inputVol
        else:
            try:
                outputVol = self.protocol.output3Dmap
                dim = outputVol.getDim()[0]
                sampling = outputVol.getSamplingRate()
                _showVol = outputVol
            except:
                # To show pdbs only
                dim = 150.
                sampling = 1.
                _showVol = None

        bildFileName = os.path.abspath(self.protocol._getTmpPath(
            "axis_output.bild"))
        createCoordinateAxisFile(dim,
                                 bildFileName=bildFileName,
                                 sampling=sampling)
        fnCmd = self.protocol._getTmpPath("chimera_output.cmd")
        f = open(fnCmd, 'w')
        f.write("open %s\n" % bildFileName)

        if _showVol is not None:
        # In case we have PDBs only, the _inputVol is None:
            showVolFileName = os.path.abspath(
                        ImageHandler.removeFileType(_showVol.getFileName()))
            f.write("open %s\n" % showVolFileName)
            if _showVol.hasOrigin():
                x, y, z = _showVol.getOrigin().getShifts()
            else:
                x, y, z = _showVol.getOrigin(force=True).getShifts()

            f.write("volume #1 style surface voxelSize %f origin "
                    "%0.2f,%0.2f,%0.2f\n"
                    % (_showVol.getSamplingRate(), x, y, z))

        for filename in os.listdir(directory):
            if filename.endswith(".pdb") or filename.endswith(".cif"):
                path = os.path.join(directory, filename)
                f.write("open %s\n" % os.path.abspath(path))

        f.close()

        # run in the background
        runChimeraProgram(getProgram(), fnCmd+"&")
        return []

class ChimeraRestoreViewer(Viewer):
    """ Visualize the output of protocols protocol_fit and protocol_operate """
    _label = 'viewer restore'
    _targets = [ChimeraProtRestore]

    def _visualize(self, obj, **args):
        path1 = os.path.join(self.protocol._getExtraPath(), sessionFile)
        if os.path.exists(path1):
            #restored SESSION
            path = os.path.abspath(path1)
        else:
            # SESSION from inputProtocol
            path2 = os.path.join(
                self.protocol.inputProtocol.get()._getExtraPath(), sessionFile)
            path = os.path.abspath(path2)

        runChimeraProgram(getProgram(), path  + "&")
        return []


class ChimeraProtRigidFitViewer(ChimeraViewerBase):
    _label = 'viewer fit'
    _targets = [ChimeraProtRigidFit]

class ChimeraProtOperateViewer(ChimeraViewerBase):
    _label = 'viewer operate'
    _targets = [ChimeraProtOperate]
