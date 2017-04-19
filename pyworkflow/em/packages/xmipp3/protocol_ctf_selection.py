# **************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
# *              Amaya Jimenez    (ajimenez@cnb.csic.es)
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


import pyworkflow.protocol.params as params
import pyworkflow.em as em
from datetime import datetime
import os
import pyworkflow.utils as pwutils
from pyworkflow.em.data import SetOfCTF


class XmippProtCTFSelection(em.ProtCTFMicrographs):
    """
    Protocol to make a selection of meaningful CTFs in basis of the defocus values,
    the astigmatism, and the resolution.
    """
    _label = 'ctf selection'
    
    def __init__(self, **args):
        em.ProtCTFMicrographs.__init__(self, **args)
        self._freqResol = {}

    #--------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        # Read a ctf estimation
        form.addParam('inputCTFs', params.PointerParam, pointerClass='SetOfCTF',
                      important=True,
                      label="Set of CTFs",
                      help='Estimated CTF to evaluate.')
        line = form.addLine('Defocus (A)',
                            help='Minimum and maximum values for defocus in Angstroms')
        line.addParam('minDefocus', params.FloatParam, default=0, label='Min')
        line.addParam('maxDefocus', params.FloatParam, default=40000, label='Max')

        form.addParam('astigmatism', params.FloatParam, default=1000,
                      label='Astigmatism (A)',
                      help='Maximum value allowed for astigmatism in Angstroms. '
                           'If the evaluated CTF does not fulfill this requirement, it will be discarded.')
        form.addParam('resolution', params.FloatParam, default=5,
                  label='Resolution (A)',
                  help='Minimum value for resolution in Angstroms. '
                       'If the evaluated CTF does not fulfill this requirement, it will be discarded.')

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        """for each ctf insert the steps to compare it
        """
        self._insertFunctionStep('readCTF', wait=True)


    # --------------------------- STEPS functions -------------------------------
    def readCTF(self):
        inputCTFs = self.inputCTFs.get()
        outputCTFs = self._createSetOfCTF()

        #self.minDefocus=20000
        #self.maxDefocus=30000

        for ctf in inputCTFs:
            defocusU = ctf.getDefocusU()
            defocusV = ctf.getDefocusV()
            astigm = defocusU - defocusV
            resol = ctf._ctffind4_ctfResolution.get()
            if defocusU>self.minDefocus and defocusU<self.maxDefocus and \
                defocusV>self.minDefocus and defocusV<self.maxDefocus and \
                astigm <self.astigmatism and resol<self.resolution:

                outputCTFs.append(ctf)

        self._defineOutputs(outputCTF=outputCTFs)
        self._defineTransformRelation(self.inputCTFs.get(), outputCTFs)

    def _stepsCheck(self):
        # Input movie set can be loaded or None when checked for new inputs
        # If None, we load it
        self._checkNewInput()
        self._checkNewOutput()

    def _loadInputList(self):
        """ Load the input set of ctf and create a list. """
        ctfFile = self.inputCTFs.get().getFileName()
        self.debug("Loading input db: %s" % ctfFile)
        ctfSet = SetOfCTF(filename=ctfFile)
        ctfSet.loadAllProperties()
        self.listOfMovies = [m.clone() for m in ctfSet]
        self.streamClosed = ctfSet.isStreamClosed()
        ctfSet.close()
        self.debug("Closed db.")

    def _checkNewInput(self):
        # Check if there are new ctfs to process from the input set
        localFile = self.inputCTFs.get().getFileName()
        now = datetime.now()
        self.lastCheck = getattr(self, 'lastCheck', now)
        mTime = datetime.fromtimestamp(os.path.getmtime(localFile))
        self.debug('Last check: %s, modification: %s'
                  % (pwutils.prettyTime(self.lastCheck),
                     pwutils.prettyTime(mTime)))
        # If the input ctfs.sqlite have not changed since our last check,
        # it does not make sense to check for new input data
        if self.lastCheck > mTime and hasattr(self, 'listOfCTFs'):
            return None

        self.lastCheck = now
        # Open input movies.sqlite and close it as soon as possible
        self._loadInputList()
        newCTFs = any(m.getObjId() not in self.insertedDict
                        for m in self.listOfCTFs)
        outputStep = self._getFirstJoinStep()#####################################

        if newCTFs:
            fDeps = self._insertNewCTFSteps(self.insertedDict,
                                               self.listOfCTFs)
            if outputStep is not None:
                outputStep.addPrerequisites(*fDeps)
            self.updateSteps()

    #--------------------------- INFO functions --------------------------------
    def _summary(self):
        message = []
        return message
    
    def _methods(self):
        pass#nothing here
    
    def _validate(self):
        """ The function of this hook is to add some validation before the protocol
        is launched to be executed. It should return a list of errors. If the list is
        empty the protocol can be executed.
        """
        message = [ ]
        fileCTF = self.inputCTFs.get()
        if fileCTF == None:
            message.append("You must specify a set of CTFs.")
        return message

    def _stepsCheck(self):
        # Just to avoid the stream checking inherited from ProtCTFMicrographs
        pass
    
        
