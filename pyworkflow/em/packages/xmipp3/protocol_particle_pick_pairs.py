# **************************************************************************
# *
# * Authors:     Jose Gutierrez Tabuenca (jose.gutierrez@cnb.csic.es)
# *              J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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

from pyworkflow.object import String
import pyworkflow.protocol.params as params
import pyworkflow.utils as pwutils
from pyworkflow.em.protocol import ProtParticlePicking
from pyworkflow.em.data_tiltpairs import CoordinatesTiltPair
from pyworkflow.em.showj import launchTiltPairPickerGUI

from xmipp3 import XmippProtocol
import convert



class XmippProtParticlePickingPairs(ProtParticlePicking, XmippProtocol):
    """ Picks particles in a set of untilted-tilted pairs of micrographs. """
    _label = 'tilt pairs particle picking'

    def __init__(self, **args):
        ProtParticlePicking.__init__(self, **args)
        # The following attribute is only for testing
        self.importFolder = String(args.get('importFolder', None))

    #--------------- DEFINE param functions ---------------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputMicrographsTiltedPair', params.PointerParam,
                      pointerClass='MicrographsTiltPair',
                      label="Micrographs tilt pair",
                      help='Select the MicrographsTiltPair ')
        form.addParam('memory', params.FloatParam, default=2,
                      label='Memory to use (In Gb)', expertLevel=2)

        #----------- INSERT steps functions ----------------------------------
    def _insertAllSteps(self):
        """ The Particle Picking process is realized for a pair
        of set of micrographs
        """
        self.micsFn = self._getPath('input_micrographs.xmd')
        # Convert input into xmipp Metadata format
        self._insertFunctionStep('convertInputStep')

        # Launch Particle Picking GUI
        if not self.importFolder.hasValue():
            self._insertFunctionStep('launchParticlePickGUIStep', interactive=True)
        else: # This is only used for test purposes
            self._insertFunctionStep('_importFromFolderStep')

    #------------------- STEPS functions -----------------------------------
    def convertInputStep(self):
        micTiltPairs = self.inputMicrographsTiltedPair.get()
        # Get the converted input micrographs in Xmipp format
        convert.writeSetOfMicrographsPairs(micTiltPairs.getUntilted(),
                                           micTiltPairs.getTilted(),
                                           self.micsFn)

    def launchParticlePickGUIStep(self):
        process = launchTiltPairPickerGUI(self.micsFn, self._getExtraPath(),
                                          self, memory='%dg' % self.memory.get())
        process.wait()

    def _importFromFolderStep(self):
        """ This function will copy Xmipp .pos files for
        simulating a particle picking run...this is only
        for testing purposes.
        """
        extraDir = self._getExtraPath()

        for f in pwutils.getFiles(self.importFolder.get()):
            pwutils.copyFile(f, extraDir)

        self.registerCoords(extraDir)

    #--------------------------- INFO functions --------------------------------------------
    def _citations(self):
        return []

    #--------------------------- UTILS functions -------------------------------------------
    def __str__(self):
        """ String representation of a Particle Picking Tilt run """
        outputs = self.getOutputsSize()

        if outputs == 0:
            msg = "No particles picked yet."
        elif outputs == 1:
            picked = self.getCoords().getSize()
            mics = self.inputMicrographsTiltedPair.get().getTilted().getSize()
            msg = "Number of particles picked: %d " % picked
            msg += "(from %d micrographs)" % mics
        else:
            msg = 'Number of outputs: %d' % outputs

        return msg

    def getInputMicrographs(self):
        return self.inputMicrographsTiltedPair.get().getTilted()

    def getCoords(self):
        return self.getCoordsTiltPair()

    def _summary(self):
        summary = []
        if self.getInputMicrographs() is  not None:
            summary.append("Number of input micrographs: %d"
                           % self.getInputMicrographs().getSize())

        if self.getOutputsSize() >= 1:
            for key, output in self.iterOutputAttributes(CoordinatesTiltPair):
                summary.append("*%s:*" % key)
                summary.append("  Particles pairs picked: %d" % output.getSize())
                summary.append("  Particle size: %d \n" % output.getBoxSize())
        else:
            summary.append("Output tilpairs not ready yet.")

        return summary

    def __getOutputSuffix(self):
        maxCounter = -1
        for attrName, _ in self.iterOutputAttributes(CoordinatesTiltPair):
            suffix = attrName.replace('outputCoordinatesTiltPair', '')
            try:
                counter = int(suffix)
            except:
                counter = 1 # when there is not number assume 1
            maxCounter = max(counter, maxCounter)

        return str(maxCounter+1) if maxCounter > 0 else '' # empty if not outputs

    def _getBoxSize(self):
        """ Redefine this function to set a specific box size to the output
        coordinates untilted and tilted.
        """
        return None

    def _readCoordinates(self, coordsDir, suffix=''):
        micTiltPairs = self.inputMicrographsTiltedPair.get()
        uSuffix = 'Untilted' + suffix
        tSuffix = 'Tilted' + suffix
        uSet = micTiltPairs.getUntilted()
        tSet = micTiltPairs.getTilted()
        # Create Untilted and Tilted SetOfCoordinates
        uCoordSet = self._createSetOfCoordinates(uSet, suffix=uSuffix)
        convert.readSetOfCoordinates(coordsDir, uSet, uCoordSet)
        uCoordSet.write()
        tCoordSet = self._createSetOfCoordinates(tSet, suffix=tSuffix)
        convert.readSetOfCoordinates(coordsDir, tSet, tCoordSet)
        tCoordSet.write()
        boxSize = self._getBoxSize()
        if boxSize:
            uCoordSet.setBoxSize(boxSize)
            tCoordSet.setBoxSize(boxSize)

        return uCoordSet, tCoordSet

    def _readAngles(self, micsFn, suffix=''):
        # Read Angles from input micrographs
        anglesSet = self._createSetOfAngles(suffix=suffix)
        convert.readAnglesFromMicrographs(micsFn, anglesSet)
        anglesSet.write()
        return anglesSet

    def registerCoords(self, coordsDir, store=True):
        micTiltPairs = self.inputMicrographsTiltedPair.get()
        suffix = self.__getOutputSuffix()

        uCoordSet, tCoordSet = self._readCoordinates(coordsDir, suffix)
        micsFn = self._getPath('input_micrographs.xmd')
        anglesSet = self._readAngles(micsFn, suffix)
        # Create CoordinatesTiltPair object
        outputset = self._createCoordinatesTiltPair(micTiltPairs,
                                                    uCoordSet, tCoordSet,
                                                    anglesSet, suffix)
        summary = self.getSummary(outputset)
        outputset.setObjComment(summary)
        outputName = 'outputCoordinatesTiltPair' + suffix
        outputs = {outputName: outputset}
        self._defineOutputs(**outputs)
        self._defineSourceRelation(self.inputMicrographsTiltedPair, outputset)
        if store:
            self._store()

