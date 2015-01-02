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
from pyworkflow.em.data import CTFModel

from base import ProtImportFiles



class ProtImportCTF(ProtImportFiles):
    """Common protocol to import a set of ctfs into the project"""
    # This label should be set in subclasses
    _label = 'import ctf'

    _outputClassName = "SetOfCTF"

    IMPORT_FROM_AUTO = 0
    IMPORT_FROM_XMIPP3 = 1
    IMPORT_FROM_BRANDEIS = 2
    IMPORT_FROM_EMX = 3

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
        return ['auto', 'xmipp','brandeis', 'emx']

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
        ci = self.getImportClass()
        self._insertFunctionStep('importCTFStep', self.getPattern())


    def getImportClass(self):
        """ Return the class in charge of importing the files. """

        if self.importFrom == self.IMPORT_FROM_EMX:
            from pyworkflow.em.packages.emxlib import EmxImport
            self.importFilePath = self.filesPath.get().strip()
            return EmxImport(self, self.importFilePath)
        elif self.importFrom == self.IMPORT_FROM_XMIPP3:
            from pyworkflow.em.packages.xmipp3.dataimport import XmippImport
            self.importFilePath = self.mdFile.get('').strip()
            return XmippImport(self, self.filesPath.get())
        elif self.importFrom == self.IMPORT_FROM_BRANDEIS:
            from pyworkflow.em.packages.brandeis.dataimport import BrandeisImport
            return BrandeisImport(self)
        else:
            self.importFilePath = ''
            return None
        
    #--------------------------- STEPS functions ---------------------------------------------------
    def importCTFStep(self, pattern):
        """ Copy ctfs matching the filename pattern.
        """
        ci = self.getImportClass()

        inputMics = self.inputMicrographs.get()
        ctfSet = self._createSetOfCTF()
        ctfSet.setMicrographs(inputMics)

        ctfFiles = self.getMatchFiles(self.getPattern())

        ctfDict = {}

        for i, (fileName, fileId) in enumerate(self.iterFiles()):
            for mic in inputMics:
                if removeBaseExt(mic.getFileName()) == removeBaseExt(fileName):
                    ci.importCoordinates(mic, fileName, coordSet)

        for fn in ctfFiles:
            # Try to match the micrograph id from filename
            # this is set by the user by using #### format in the pattern
            match = self._idRegex.match(fn)
            if match is None:
                raise Exception("File '%s' doesn't match the pattern '%s'" % (fn, self.pattern.get()))
            ctfId = int(match.group(1))
            ctfDict[ctfId] = fn

        for mic in inputMics:
            if mic.getObjId() in ctfDict:
                defocusU, defocusV, defocusAngle = ci.getCTFParams(ctfDict[mic.getObjId()])

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

    
    #--------------------------- INFO functions ----------------------------------------------------
    
    def _summary(self):
        summary = []
        if self.ctfSet is None:
            summary.append("Output " + self._outputClassName + " not ready yet.")
            if self.copyFiles:
                summary.append("*Warning*: You select to copy files into your project.\n"
                               "This will make another copy of your data and may take \n"
                               "more time to import. ")
        else:
            summary.append("Import of *%d* %s from %s" % (self.ctfSet.getSize(), self._outputClassName, self.getPattern()))

        return summary
    
    def _methods(self):
        methods = []

        return methods
    
    #--------------------------- UTILS functions ---------------------------------------------------



