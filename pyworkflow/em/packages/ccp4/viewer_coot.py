# **************************************************************************
# *
# * Authors:  Roberto Marabini (roberto@cnb.csic.es), May 2013
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
from pyworkflow.em.packages.xmipp3.viewer import XmippViewer, PdbFile
from pyworkflow.em.utils.chimera_utilities.convert import \
    createCoordinateAxisFile, \
    adaptOriginFromCCP4ToChimera, runChimeraProgram, \
    getProgram, chimeraPdbTemplateFileName
from pyworkflow.em.utils.ccp4_utilities.convert import cootPdbTemplateFileName
from pyworkflow.viewer import DESKTOP_TKINTER, Viewer
from protocol_coot import CCP4ProtCoot

# TODO: very likely this should inherit from ProtocolViewer
# not from XmippViewer. But then I get an empty form :-(


class CootRefineViewer(Viewer):
    """ Visualize the output of protocol volume strain """
    _label = 'coot viewer'
    _targets = [CCP4ProtCoot]
    _environments = [DESKTOP_TKINTER]

    def _visualize(self, obj, **args):
        print "11111111111111111111111111:"
            # TODO if input volume is not mrc this will not work.
        # Construct the coordinate file and visualization
        bildFileName = os.path.abspath(self.protocol._getTmpPath(
            "axis.bild"))
        dims = []
        samplings = []
        if self.protocol.inputVolumes.get() is None:
            dim = self.protocol.pdbFileToBeRefined.get().getVolume().getDim()[0]
            sampling = self.protocol.pdbFileToBeRefined.get().getVolume().getSamplingRate()
            dims.append(dim)
            samplings.append(sampling)
            print "AAAAAAAAAAAAAAAA:" \
                  "dim and sampling from volume associated " \
                  "to pdb"
        else:
            dim = self.protocol.inputVolumes.get().getDim()[0]
            sampling = self.protocol.inputVolumes.get().getSamplingRate()
            dims.append(dim)
            samplings.append(sampling)
            print "BBBBBBBBBBBBBB: dim and sampling from the input volume"
        createCoordinateAxisFile(max(dims),
                                 bildFileName=bildFileName,
                                 sampling=max(samplings))
        fnCmd = self.protocol._getTmpPath("chimera.cmd")
        f = open(fnCmd, 'w')
        f.write("open %s\n" % bildFileName)

        try:
            outputsVol = self.protocol.outputs
            count = 1
            for outputVol in outputsVol:
                outputVolFileName = os.path.abspath(outputVol.getFileName())
                f.write("open %s\n" % outputVolFileName)
                f.write("volume #%d style surface\n" % count)
                count =+ 1
                print "CCCCCCCCCCC: Saved the volume generated"
        except:
            outputsVol = []
            if self.protocol.inputVolumes.get() is None:
                outputVol = self.protocol.pdbFileToBeRefined.get().getVolume()
                outputsVol.append(outputVol)
                print "DDDDDDDDDDD: Saved the starting volume"

            else:
                for outputVol in self.protocol.inputVolumes.get():
                    outputsVol.append(outputVol)
                    print "EEEEEEEEEEEEEEEEEEE: Saved the starting volume"
            count = 1
            for outputVol in outputsVol:
                outputVolFileName = os.path.abspath(
                        ImageHandler.removeFileType(outputVol.getFileName()))
                x, y, z = adaptOriginFromCCP4ToChimera(
                    outputVol.getOrigin().getShifts())
                f.write("open %s\n" % outputVolFileName)
                f.write("volume #%d  style surface voxelSize %f origin "
                        "%0.2f,%0.2f,%0.2f\n"
                        % (count, outputVol.getSamplingRate(),x, y, z))
                count =+ 1

        outputsPDB = self.protocol.outputs
        counter = 1
        for outputPDB in outputsPDB:
            template = self.protocol._getExtraPath(cootPdbTemplateFileName)
            outputPDB = os.path.abspath(template%counter)
            if os.path.exists(outputPDB):
                f.write("open %s\n" % outputPDB)
                counter =+ 1

        f.close()
        #return [em.ChimeraView(fnCmd)]

        #args = ""
        #args += fnCmd + " "
        #self.protocol._log.info('Launching: ' + getProgram() + ' ' + args)

        # run in the background
        runChimeraProgram(getProgram(), fnCmd+"&")
        return []



