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

from pyworkflow.em.data import SetOfCTF, SetOfMicrographs
from pyworkflow.protocol.constants import (STATUS_NEW)


class XmippProtCTFSelection(em.ProtCTFMicrographs):
    """
    Protocol to make a selection of meaningful CTFs in basis of the defocus values,
    the astigmatism, and the resolution.
    """
    _label = 'ctf selection'
    
    def __init__(self, **args):
        em.ProtCTFMicrographs.__init__(self, **args)

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
        form.addParam('resolution', params.FloatParam, default=7,
                  label='Resolution (A)',
                  help='Minimum value for resolution in Angstroms. '
                       'If the evaluated CTF does not fulfill this requirement, it will be discarded.')

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        """ This is a very easy program, the createOutputStep function is a dummy function
            There is a loop that calls the function stepsCheck that will be broken when the
            status associated to createOutputStep change from wait=True to wait=False
        """

        # lastly function to be executed when streaming is done
        # it will be executed when waitCondition is set to True
        self._insertFunctionStep('createOutputStep', wait=True)


    # --------------------------- STEPS functions -------------------------------

    def createOutputStep(self):
        # Do nothing now, the output should be ready.
        pass


    def _getFirstJoinStepName(self):
        # This function will be used for streaming, to check which is
        # the first function that need to wait for all micrographs
        # to have completed, this can be overriden in subclasses
        # (e.g., in Xmipp 'sortPSDStep')
        return 'createOutputStep'


    def _getFirstJoinStep(self):
        for s in self._steps:
            if s.funcName == self._getFirstJoinStepName():
                return s
        return None


    def _checkNewCTFs(self, ctfInSet):
        """ Check for already computed CTF and update the output set. """
        ctfDict = {}
        newCTF = False
        micSet = SetOfMicrographs(filename=self._getPath('micrographs.sqlite'))
        ctfSet = SetOfCTF(filename=self._getPath('ctfs.sqlite')) #reading this protocol output
        ctfSet.setMicrographs(ctfInSet.getMicrographs())

        #Create a dict with the processed CTFs
        for ctf in ctfSet:
            ctfDict[ctf.getObjId()] = True

        if ctfDict: # it means there are previous ctfs computed
            ctfSet.loadAllProperties()
            ctfSet.enableAppend()
        else:
            ctfSet.setStreamState(ctfSet.STREAM_OPEN)

        print "EN _checkNewCTFs"

        for ctf in ctfInSet:

            if not ctf.getObjId() in ctfDict:
                defocusU = ctf.getDefocusU()
                defocusV = ctf.getDefocusV()
                astigm = defocusU - defocusV
                resol = ctf._ctffind4_ctfResolution.get()  # PENDIENTEEEEEEEEEEEEEEEEEE
                if defocusU > self.minDefocus and defocusU < self.maxDefocus and \
                defocusV > self.minDefocus and defocusV < self.maxDefocus and \
                astigm < self.astigmatism and resol < self.resolution:
                    pass
                else:
                    ctf.setEnabled(False)
                    #ctf.getMicrograph().setEnabled(False)

                print ctf.getMicrograph().printAll()
                ctfSet.append(ctf)
                micSet.append(ctf.getMicrograph())
                print "Added"
                newCTF = True

        return newCTF, ctfSet, micSet



    def _stepsCheck(self):

        #check if there are new CTFs
        ctfFn = self.inputCTFs.get().getFileName()
        self.micSet = self.inputCTFs.get().getMicrographs()
        ctfSet = SetOfCTF(filename=ctfFn)
        ctfSet.loadAllProperties()
        ctfSet.setMicrographs(self.micSet)
        streamClosed = ctfSet.isStreamClosed()
        outputStep = self._getFirstJoinStep()

        #check if there are new ctfs and process them
        #return number of processed ctfs
        newCTF, ctfOutSet, microOutSet = self._checkNewCTFs(ctfSet)
        ctfOutSet.setMicrographs(self.micSet)

        if newCTF:
            # Check if it is the first time we are registering CTF
            firstTime = not self.hasAttribute('outputCTF')

            if streamClosed:
                streamMode = ctfSet.STREAM_CLOSED
            else:
                streamMode = ctfSet.STREAM_OPEN
            self._updateOutputSet('outputCTF', ctfOutSet, streamMode)
            self._updateOutputSet('outputMicrograph', microOutSet, streamMode)
            if firstTime:  # define relation just once
                self._defineCtfRelation(self.inputCTFs, ctfOutSet)

        else:
            ctfSet.close()

        if outputStep and outputStep.isWaiting() and streamClosed:
            outputStep.setStatus(STATUS_NEW)

        ctfSet.close()


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

