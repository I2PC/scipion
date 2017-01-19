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

import os
from os.path import join

import pyworkflow.protocol.params as params
import pyworkflow.utils as pwutils
from pyworkflow import VERSION_1_1
from pyworkflow.em.protocol import EMProtocol
from pyworkflow.em.data import SetOfMicrographs

from convert import writeSetOfMicrographs


class ProtRelionExportCtf(EMProtocol):
    """
    Export a SetOfCTF to the expected Relion STAR file.
    """
    _label = 'export ctf'
    _version = VERSION_1_1
    CTF_STAR_FILE = 'micrographs_ctf.star'
    
    #--------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        
        form.addSection(label='Input')

        form.addParam('inputCTF', params.PointerParam,
                      pointerClass="SetOfCTF",
                      label='Input CTF',
                      help='Select set of CTF that you want to export.')
        
        form.addParallelSection(threads=0, mpi=0)
            
    #--------------------------- INSERT steps functions ------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('writeCtfStarStep')
        
    def writeCtfStarStep(self):
        inputCTF = self.inputCTF.get()
        print "Writing SetOfCtf: %s" % inputCTF
        print " to: %s" % self._getPath(self.CTF_STAR_FILE)

        ctfMicSet = inputCTF.getMicrographs()
        micSet = SetOfMicrographs(filename=':memory:')

        psd = inputCTF.getFirstItem().getPsdFile()
        hasPsd = psd and os.path.exists(psd)

        if hasPsd:
            psdPath = self._getPath('PSD')
            pwutils.makePath(psdPath)
            print "Writing PSD files to %s" % psdPath

        for ctf in inputCTF:
            # Get the corresponding micrograph
            mic = ctfMicSet[ctf.getObjId()]
            mic2 = mic.clone()
            mic2.setCTF(ctf)
            if hasPsd:
                psdFile = ctf.getPsdFile()
                newPsdFile = join('PSD', '%s_psd.mrc' % mic.getMicName())
                pwutils.copyFile(psdFile, self._getPath(newPsdFile))
                ctf.setPsdFile(newPsdFile)
            micSet.append(mic2)

        writeSetOfMicrographs(micSet, self._getPath(self.CTF_STAR_FILE),
                              preprocessImageRow=self.preprocessMicrograph)

    #--------------------------- INFO functions --------------------------------

    def _summary(self):
        summary = []
        ctfStarFn = self._getPath(self.CTF_STAR_FILE)

        if os.path.exists(ctfStarFn):
            summary.append("Output CTF STAR file written to: \n%s" % ctfStarFn)
        else:
            summary.append("No output generated yet.")

        return summary
    
    #--------------------------- UTILS functions -------------------------------
    def preprocessMicrograph(self, mic, micRow):
        micRow.setValue('rlnCtfImage', mic.getCTF().getPsdFile())
