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

from pyworkflow.protocol.params import PathParam, PointerParam
from pyworkflow.em.data import CTFModel
from pyworkflow.em.protocol import ProtImportFiles, ProtCTFMicrographs



class ProtImportCTF(ProtImportFiles, ProtCTFMicrographs):
    """ Protocol to import results from CTFfind.
    Select the micrographs to associate the computed CTFs and
    the path pattern where to find the CTF files, that should
    contains the micrograph id in the name. """
    _label = 'import ctf'
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection('Input')
        form.addParam('inputMicrographs', PointerParam, pointerClass='SetOfMicrographs', 
                          label='Input micrographs',
                          help='Select the particles that you want to update the CTF parameters.')
        form.addParam('pattern', PathParam, 
                      label='CTF pattern',
                      help='Select files containing the CTF estimation.\n'
                           'You should use #### characters in the pattern\n'
                           'to mark where the micrograph id will be taken from. ')
    
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('importCtfStep', 
                                 self.inputMicrographs.get().getObjId(),
                                 self.getPattern())
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def importCtfStep(self, micsId, pattern):
        """ Copy movies matching the filename pattern
        Register other parameters.
        """
        inputMics = self.inputMicrographs.get()
        ctfSet = self._createSetOfCTF()
        ctfSet.setMicrographs(inputMics)
        
        ctfFiles = self._getFilePaths(pattern)
        ctfDict = {}
        for fn in ctfFiles:
            # Try to match the micrograph id from filename
            # this is set by the user by using #### format in the pattern
            match = self._idRegex.match(fn)
            if match is None:
                raise Exception("File '%s' doesn't match the pattern '%s'" % (fn, self.pattern.get()))
            ctfId = int(match.group(1))
            ctfDict[ctfId] = fn
            
        from pyworkflow.em.packages.grigoriefflab.convert import parseCtffindOutput

        for mic in inputMics:
            if mic.getObjId() in ctfDict:
                defocusU, defocusV, defocusAngle = parseCtffindOutput(ctfDict[mic.getObjId()])
            else:
                self.warning("CTF for micrograph id %d was not found." % mic.getObjId())
                defocusU, defocusV, defocusAngle = -999, -1, -999
            
            # save the values of defocus for each micrograph in a list
            ctf = CTFModel()
            ctf.copyObjId(mic)
            ctf.setStandardDefocus(defocusU, defocusV, defocusAngle)
            ctf.setMicrograph(mic)
            ctfSet.append(ctf)
        
        self._defineOutputs(outputCTF=ctfSet)
        self._defineCtfRelation(inputMics, ctfSet)            
    
    def _summary(self):
        summary = []
        return summary    
    
    def _methods(self):
        return []
    
    def _stepsCheck(self):
        # Just to avoid the stream checking inherited from ProtCTFMicrographs
        pass
