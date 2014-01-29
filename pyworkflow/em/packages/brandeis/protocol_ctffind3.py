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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************
"""
This module contains the protocol for CTF estimation with ctffind3
"""


from pyworkflow.em import *  
from pyworkflow.utils.which import which
from pyworkflow.utils.path import makePath, replaceBaseExt, join, basename


class ProtCTFFind(ProtCTFMicrographs):
    """Protocol to perform CTF estimation on a set of micrographs
    using the ctffind3 program"""
    _label = 'ctffind3'
    _references = ['[[http://www.ncbi.nlm.nih.gov/pubmed/12781660][Mindell JA & Grigorieff N., J Struct Biol. (2003)]]']

    def _prepareCommand(self):
        self._params['step_focus'] = 1000.0
        # Convert digital frequencies to spatial frequencies
        sampling = self.inputMics.getSamplingRate()
        self._params['lowRes'] = sampling / self._params['lowRes']
        self._params['highRes'] = sampling / self._params['highRes']        
        
        if which('ctffind3.exe') is '':
            raise Exception('Missing ctffind3.exe')
         
        self._program = 'export NATIVEMTZ=kk ; ' + which('ctffind3.exe')
        self._args = """   << eof > %(ctffindOut)s
%(micFn)s
%(ctffindPSD)s
%(sphericalAberration)f,%(voltage)f,%(ampContrast)f,%(magnification)f,%(scannedPixelSize)f
%(windowSize)d,%(lowRes)f,%(highRes)f,%(minDefocus)f,%(maxDefocus)f,%(step_focus)f
"""

    def _getPsdPath(self, micDir):
        return join(micDir, 'ctffind_psd.mrc')
    
    def _estimateCTF(self, micFn, micDir):
        """ Run ctffind3 with required parameters """
        # Create micrograph dir 
        makePath(micDir)
        # Update _params dictionary
        self._params['micFn'] = micFn
        self._params['micDir'] = micDir
        self._params['ctffindOut'] = join(micDir, 'ctffind.out')
        self._params['ctffindPSD'] = self._getPsdPath(micDir)
                
        self.runJob(None, self._program, self._args % self._params)
        #print "command: ", self._program, self._args % self._params    

    def _parseOutput(self, filename):
        """ Try to find the output estimation parameters
        from filename. It search for a line containing: Final Values.
        """
        f = open(filename)
        result = None
        for line in f:
            if 'Final Values' in line:
                # Take DefocusU, DefocusV and Angle as a tuple
                # that are the first three values in the line
                result = tuple(map(float, line.split()[:3]))
                break
        f.close()
        return result
            
    def _getCTFModel(self, defocusU, defocusV, defocusAngle, psdFile):
        ctf = CTFModel()
#        ctf.copyAttributes(self.inputMics, 'samplingRate')
#        ctf.copyAttributes(self.inputMics.microscope, 'voltage', 'sphericalAberration')
        ctf.defocusU.set(defocusU)
        ctf.defocusV.set(defocusV)
        ctf.defocusAngle.set(defocusAngle)
        ctf.psdFile.set(psdFile)
        
        return ctf
        
    def createOutput(self):
        ctfSet = self._createSetOfCTF()
        
        for fn, micDir, mic in self._iterMicrographs():
            out = join(micDir, 'ctffind.out')
            result = self._parseOutput(out)
            defocusU, defocusV, defocusAngle = result
            #micOut = Micrograph()
            #micOut.setFileName(mic.getFileName())
            #micOut.setCTF(self._getCTFModel(defocusU, defocusV, defocusAngle, 
            #                                    self._getPsdPath(micDir)))
            #micSet.append(micOut)
            ctfModel = self._getCTFModel(defocusU, defocusV, defocusAngle, 
                                                self._getPsdPath(micDir))
            ctfModel.micFile.set(mic.getFileName())
            ctfSet.append(ctfModel)

        # This property should only be set by CTF estimation protocols
        #micSet.setHasCTF(True)     
            
        #micSet.write()
        ctfSet.write() 
        
        #self._defineOutputs(outputMicrographs=micSet)
        self._defineOutputs(outputCTF=ctfSet)
        #self._defineSourceRelation(self.inputMics, micSet)
        self._defineRelation(RELATION_CTF, ctfSet, self.inputMics)

    def _validate(self):
        errors = []
        if which('ctffind3.exe') is '':
            errors.append('Missing ctffind3.exe')
        return errors
            
