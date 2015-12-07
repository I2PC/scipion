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
In this module are protocol base classes related to EM imports of Micrographs, Particles, Volumes...
"""

from pyworkflow.protocol.params import PointerParam
from pyworkflow.em.data import SetOfMicrographs
from pyworkflow.utils.path import removeBaseExt

from base import ProtImportFiles



class ProtImportCTF(ProtImportFiles):
    """Common protocol to import a set of ctfs into the project"""
    # This label should be set in subclasses
    _label = 'import ctf'

    _outputClassName = "SetOfCTF"

    IMPORT_FROM_AUTO = 0
    IMPORT_FROM_XMIPP3 = 1
    IMPORT_FROM_GRIGORIEFF = 2

    #--------------------------- DEFINE param functions --------------------------------------------

    def _defineImportParams(self, form):
        """ Just redefine to put some import parameters.
        """
        form.addParam('inputMicrographs', PointerParam, pointerClass='SetOfMicrographs',
                          label='Input micrographs',
                          help='Select the particles that you want to update the CTF parameters.')

    def _getImportChoices(self):
        """ Return a list of possible choices
        from which the import can be done.
        (usually packages formats such as: xmipp3, eman2, relion...etc.
        """
        return ['auto', 'xmipp','grigorieff']

    def _getDefaultChoice(self):
        return  self.IMPORT_FROM_AUTO

    def _getFilesCondition(self):
        """ Return an string representing the condition
        when to display the files path and pattern to grab
        files.
        """
        return True

    #--------------------------- INSERT functions ---------------------------------------------------
    def _insertAllSteps(self):
        importFrom = self.importFrom.get()
        self._insertFunctionStep('importCTFStep', importFrom)


    def getImportClass(self):
        """ Return the class in charge of importing the files. """
        filesPath = self.filesPath.get()
        importFrom = self.importFrom.get()
        if importFrom == self.IMPORT_FROM_AUTO:
            importFrom = self.getFormat()

        if importFrom == self.IMPORT_FROM_XMIPP3:
            from pyworkflow.em.packages.xmipp3.dataimport import XmippImport
            return XmippImport(self, filesPath)
        elif importFrom == self.IMPORT_FROM_GRIGORIEFF:
            from pyworkflow.em.packages.grigoriefflab.dataimport import GrigorieffLabImportCTF
            return GrigorieffLabImportCTF(self)
        else:
            return None
        
    #--------------------------- STEPS functions ---------------------------------------------------
    def importCTFStep(self, importFrom):
        """ Copy ctfs matching the filename pattern. """
        ci = self.getImportClass()

        inputMics = self.inputMicrographs.get()
        ctfSet = self._createSetOfCTF()
        ctfSet.setMicrographs(inputMics)

        outputMics = self._createSetOfMicrographs()
        SetOfMicrographs.copyInfo(outputMics, inputMics)

        createOutputMics = False
        
        files = [f for f, _ in self.iterFiles()]
        n = len(files)
        if n == 0:
            raise Exception("No files where found in path: '%s'\n"
                            "matching the pattern: '%s'" % (self.filesPath, self.filesPattern))
        print "Matching files: ", len(files)

        for mic in inputMics:
            micName = mic.getMicName()
            micBase = removeBaseExt(mic.getFileName())
            for fileName, fileId in self.iterFiles():
                if (fileId == mic.getObjId() or
                    micBase in fileName or micName in fileName):
                    ctf = ci.importCTF(mic, fileName)
                    ctfSet.append(ctf)
                    outputMics.append(mic)
                    break
            else:
                # If not CTF is found for a micrograph remove it from output mics
                self.warning("CTF for micrograph id %d was not found. Removed from set of micrographs." % mic.getObjId())
                createOutputMics = True

        self._defineOutputs(outputCTF=ctfSet)
        # If some of the micrographs had not ctf a subset of micrographs have been produced
        if createOutputMics:
            self._defineOutputs(outputMicrographs=outputMics)
            self._defineCtfRelation(outputMics, ctfSet)
        else:
            self._defineCtfRelation(inputMics, ctfSet)

    #--------------------------- INFO functions ----------------------------------------------------
    
    def _summary(self):
        summary = []
        if not hasattr(self, 'outputCTF'):
            summary.append("Output " + self._outputClassName + " not ready yet.")
            if self.copyFiles:
                summary.append("*Warning*: You select to copy files into your project.\n"
                               "This will make another copy of your data and may take \n"
                               "more time to import. ")
        else:
            summary.append("Imported *%d* CTFs " % self.outputCTF.getSize())
            summary.append("from %s" % self.getPattern())

        return summary
    
    def _methods(self):
        methods = []

        return methods
    
    def _validate(self):
        errors = []
        files = [f for f, _ in self.iterFiles()]
        if not files:
            errors.append("No files where found in path: '%s'\n"
                          "matching the pattern: '%s'" % (self.filesPath, self.filesPattern))
        return errors
    
    #--------------------------- UTILS functions ---------------------------------------------------
    def getFormat(self):
        for fileName, _ in self.iterFiles():
            if fileName.endswith('.log') or fileName.endswith('.txt'):
                return self.IMPORT_FROM_GRIGORIEFF
            if fileName.endswith('.ctfparam'):
                return self.IMPORT_FROM_XMIPP3
        return -1



