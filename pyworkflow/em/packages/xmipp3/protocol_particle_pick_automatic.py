# **************************************************************************
# *
# * Authors:     Jose Gutierrez Tabuenca (jose.gutierrez@cnb.csic.es)
# *              Laura del Cano (laura.cano@cnb.csic.es)
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

from pyworkflow.em import *  
from pyworkflow.utils.path import *  
import xmipp
from xmipp3 import XmippProtocol

from convert import readSetOfCoordinates


MICS_SAMEASPICKING = 0
MICS_OTHER = 1


class XmippParticlePickingAutomatic(ProtParticlePickingAuto, XmippProtocol):
    """Protocol to pick particles automatically in a set of
    micrographs using previous training """
    _label = 'auto-picking (step 2)'  
    
    filesToCopy = ['model_training.txt', 'model_svm.txt',
                   'model_pca_model.stk', 'model_rotpca_model.stk',
                   'model_particle_avg.xmp', 'config.xmd', 'templates.stk']
    
    def __init__(self, **kwargs):
        ProtParticlePickingAuto.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL
    
    # --------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
    
        form.addSection(label='Input')

        form.addParam('xmippParticlePicking', PointerParam,
                      label="Xmipp particle picking run",
                      pointerClass='XmippProtParticlePicking',
                      #pointerCondition='isFinished',
                      help='Select the previous xmipp particle picking run.')

        form.addParam('micsToPick', EnumParam,
                      choices=['Same as supervised', 'Other'],
                      default=0, label='Micrographs to pick',
                      display=EnumParam.DISPLAY_LIST,
                      help="Select from which set of micrographs to pick using "
                           "the training from supervised run."
                           "If you use Same as supervised, the same set of "
                           "micrographs used for training the picker will be "
                           "used at this point. If you select Other, you can "
                           "select another set of micrograph (normally from "
                           "the same specimen) and pick them completely "
                           "automatic using the trained picker.")

        form.addParam('inputMicrographs', PointerParam,
                      pointerClass='SetOfMicrographs',
                      condition='micsToPick==%d' % MICS_OTHER,
                      label="Micrographs",
                      help="Select other set of micrographs to pick using the "
                           "trained picker.")

        form.addParam('memory', FloatParam, default=2,
                      label='Memory to use (In Gb)', expertLevel=2)

        self._defineStreamingParams(form)

        form.addParallelSection(threads=1, mpi=1)
        
    # --------------------------- INSERT steps functions -----------------------
    def _insertInitialSteps(self):
        # Get pointer to input micrographs
        self.particlePickingRun = self.xmippParticlePicking.get()

        copyId = self._insertFunctionStep('copyInputFilesStep')

        return [copyId]

    # --------------------------- STEPS functions ------------------------------
    def copyInputFilesStep(self):
        # Copy training model files to current run
        for f in self.filesToCopy:
            copyFile(self.particlePickingRun._getExtraPath(f),
                     self._getExtraPath(f))

    def _pickMicrograph(self, mic, *args):
        micPath = mic.getFileName()
        # Get particle picking boxsize from the previous run
        boxSize = self.particlePickingRun.outputCoordinates.getBoxSize()
        modelRoot = self._getExtraPath('model')

        micName = removeBaseExt(micPath)
        proceed = True
        if self.micsToPick == MICS_SAMEASPICKING:
            basePos = replaceBaseExt(micPath, "pos")
            fnPos = self.particlePickingRun._getExtraPath(basePos)
            if exists(fnPos):
                blocks = xmipp.getBlocksInMetaDataFile(fnPos)
                copy = True
                if 'header' in blocks:
                    mdheader = xmipp.MetaData("header@" + fnPos)
                    state = mdheader.getValue(xmipp.MDL_PICKING_MICROGRAPH_STATE,
                                              mdheader.firstObject())
                    if state == "Available":
                        copy = False
                if copy:
                    # Copy manual .pos file of this micrograph
                    copyFile(fnPos, self._getExtraPath(basename(fnPos)))
                    proceed = False            

        if proceed:
            args = "-i %s " % micPath
            args += "--particleSize %d " % boxSize
            args += "--model %s " % modelRoot
            args += "--outputRoot %s " % self._getExtraPath(micName)
            args += "--mode autoselect --thr %d" % self.numberOfThreads

            self.runJob("xmipp_micrograph_automatic_picking", args)

    def readSetOfCoordinates(self, workingDir, coordSet):
        readSetOfCoordinates(workingDir, self.getInputMicrographs(), coordSet)

    def readCoordsFromMics(self, workingDir, micList, coordSet):
        readSetOfCoordinates(workingDir, micList, coordSet)
        
    # --------------------------- INFO functions -------------------------------
    def _validate(self):
        validateMsgs = []
        
        if not hasattr(self.xmippParticlePicking.get(),"outputCoordinates"):
            validateMsgs.append("You need to generate coordinates for the "
                                "supervised picking")
   
        srcPaths = [self.xmippParticlePicking.get()._getExtraPath(k)
                    for k in self.filesToCopy]
        # Check that all needed files exist
        if missingPaths(*srcPaths):
            validateMsgs.append('Input picking run has not been trained, '
                                'use *Autopick* for at least one micrograph')
            
        # If other set of micrographs is provided they should have same
        # sampling rate and acquisition
        if self.micsToPick.get() == MICS_OTHER:
            inputMics = self.inputMicrographs.get()
            manualMics = self.xmippParticlePicking.get().inputMicrographs.get()
            pixsizeInput = inputMics.getSamplingRate()

            pixsizeMics = manualMics.getSamplingRate()
            acq = manualMics.getAcquisition()

            if pixsizeInput != pixsizeMics:
                validateMsgs.append('New micrographs should have same sampling '
                                    'rate as the ones already picked.')

            if not inputMics.getAcquisition().equalAttributes(acq):
                validateMsgs.append('New micrographs should have same '
                                    'acquisition parameters as the ones '
                                    'already picked.')

        return validateMsgs
    
    def getSummary(self, coordSet):
        summary = []
        summary.append("Previous run: %s" %
                       self.xmippParticlePicking.get().getNameId())
        configfile = join(self._getExtraPath(), 'config.xmd')
        existsConfig = exists(configfile)
        if existsConfig:
            md = xmipp.MetaData('properties@' + configfile)
            configobj = md.firstObject()
            def _get(label):
                return md.getValue(label, configobj)
            pickingState = _get(xmipp.MDL_PICKING_STATE)
            particleSize = _get(xmipp.MDL_PICKING_PARTICLE_SIZE)
            activeMic = _get(xmipp.MDL_MICROGRAPH)
            isAutopick = pickingState != "Manual"
            manualParticlesSize = _get(xmipp.MDL_PICKING_MANUALPARTICLES_SIZE)
            autoParticlesSize = _get(xmipp.MDL_PICKING_AUTOPARTICLES_SIZE)
            
            summary.append("Manual particles picked: %s" % manualParticlesSize)
            summary.append("Particle size:%d" %(particleSize))
            autopick = "Yes" if isAutopick else "No"
            summary.append("Autopick: " + autopick)
            if isAutopick:
                summary.append("Automatic particles picked: %s"
                               % autoParticlesSize)
            summary.append("Last micrograph: " + activeMic)
        return "\n".join(summary)

    def getMethods(self, output):
        manualPickName = self.xmippParticlePicking.get().getNameId()
        msg = 'Program picked %d particles ' % output.getSize()
        msg += 'of size %d ' % output.getBoxSize()
        msg += 'using training from %s. ' % manualPickName
        msg += 'For more detail see [Abrishami2013]'
        return msg

    def _citations(self):
        return ['Abrishami2013']
    
    # --------------------------- UTILS functions ------------------------------
    def getCoordsDir(self):
        return self._getExtraPath()
    
    def getInputMicrographsPointer(self):
        # Get micrographs to pick
        if self.micsToPick == MICS_SAMEASPICKING:
            inputPicking = self.xmippParticlePicking.get()
            if inputPicking is None:
                return None
            else:
                return inputPicking.inputMicrographs
        else:
            return self.inputMicrographs
        
    def getInputMicrographs(self):
        """ Return the input micrographs that can be the same of the supervised
        picking or other ones selected by the user. (This can be used to pick
        a new set of micrographs with the same properties than a previous
        trained ones. )
        """ 
        if self.getInputMicrographsPointer() is not None:
            return self.getInputMicrographsPointer().get()
        else:
            return None
