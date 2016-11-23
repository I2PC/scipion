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


class XmippParticlePickingAutomatic(ProtParticlePicking, XmippProtocol):
    """Protocol to pick particles automatically in a set of
    micrographs using previous training """
    _label = 'auto-picking (step 2)'  
    
    filesToCopy = ['model_training.txt', 'model_svm.txt',
                   'model_pca_model.stk', 'model_rotpca_model.stk',
                   'model_particle_avg.xmp', 'config.xmd', 'templates.stk']
    
    def __init__(self, **args):        
        ProtParticlePicking.__init__(self, **args)
        self.stepsExecutionMode = STEPS_PARALLEL
    
    # --------------------------- DEFINE param functions -------------------------------------
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
                      help='Select from which set of micrographs to pick using '
                           'the training from supervised run.'
                           'If you use Same as supervised, the same set of micrographs '
                           'used for training the picker will be used at this point. '
                           'If you select Other, you can select another set of micrograph'
                           '(normally from the same specimen) and pick them completely '
                           'automatic using the trained picker.')
        form.addParam('inputMicrographs', PointerParam, label="Micrographs",
                      pointerClass='SetOfMicrographs',
                      condition='micsToPick==%d' % MICS_OTHER,
                      help='Select other set of micrographs to pick using the trained picker.')
        form.addParam('memory', FloatParam, default=2,
                      label='Memory to use (In Gb)', expertLevel=2)
        form.addParallelSection(threads=1, mpi=1)
        
    # --------------------------- INSERT steps functions -------------------------------------
    def _insertAllSteps(self):
        """The Particle Picking process is realized for a set of micrographs"""
        
        # Get pointer to input micrographs 
        self.particlePickingRun = self.xmippParticlePicking.get()
        
        copyId = self._insertFunctionStep('copyInputFilesStep')
        # Get micrographs to pick
        #self.inputMicrographs.set(self.getInputMicrographs())
            
        deps = []
        for mic in self.getInputMicrographs():
            stepId = self._insertFunctionStep('autopickMicrographStep',
                                              mic.getFileName(),
                                              prerequisites=[copyId])
            deps.append(stepId)
                    
        self._insertFunctionStep('_createOutput', self._getExtraPath(),
                                 prerequisites=deps)
        
    # --------------------------- STEPS functions --------------------------------------------
    def copyInputFilesStep(self):
        # Copy training model files to current run
        for f in self.filesToCopy:
            copyFile(self.particlePickingRun._getExtraPath(f),
                     self._getExtraPath(f))
        
    def autopickMicrographStep(self, micPath):
        # Get particle picking boxsize from the previous run
        boxSize = self.particlePickingRun.outputCoordinates.getBoxSize()
        modelRoot = self._getExtraPath('model')

        micName = removeBaseExt(micPath)
        proceed = True
        if self.micsToPick == MICS_SAMEASPICKING:
            fnPos = self.particlePickingRun._getExtraPath(replaceBaseExt(micPath, "pos"))
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
            oroot = self._getExtraPath(micName)
            thr = self.numberOfThreads.get()
            args = "-i %(micPath)s --particleSize %(boxSize)d --model %(modelRoot)s " % locals()
            args += " --outputRoot %(oroot)s --mode autoselect --thr %(thr)d" % locals()
            # TODO: What is this?
#           if self.Fast:
#               args += " --fast "
            self.runJob("xmipp_micrograph_automatic_picking", args)

    def readSetOfCoordinates(self, workingDir, coordSet):
        readSetOfCoordinates(workingDir, self.getInputMicrographs(), coordSet)
        
    # --------------------------- INFO functions ---------------------------------------------
    def _validate(self):
        validateMsgs = []
        
        if not hasattr(self.xmippParticlePicking.get(),"outputCoordinates"):
            validateMsgs.append("You need to generate coordinates for the supervised picking")
   
        srcPaths = [self.xmippParticlePicking.get()._getExtraPath(k) for k in self.filesToCopy]
        # Check that all needed files exist
        if missingPaths(*srcPaths):
            validateMsgs.append('Input picking run has not been trained, '
                                'use *Autopick* for at least one micrograph')
            
        # If other set of micrographs is provided they should have same
        # sampling rate and acquisition
        if self.micsToPick.get() == MICS_OTHER:
            pixsizeInput = self.inputMicrographs.get().getSamplingRate()
            pixsizeMics = self.xmippParticlePicking.get().inputMicrographs.get().getSamplingRate()
            acq = self.xmippParticlePicking.get().inputMicrographs.get().getAcquisition()

            if pixsizeInput != pixsizeMics:
                validateMsgs.append('New micrographs should have same sampling '
                                    'rate as the ones already picked.')

            if not self.inputMicrographs.get().getAcquisition().equalAttributes(acq):
                validateMsgs.append('New micrographs should have same acquisition '
                                    'parameters as the ones already picked.')

        return validateMsgs
    
    def getSummary(self, coordSet):
        summary = []
        summary.append("Previous run: " + self.xmippParticlePicking.get().getNameId())
        configfile = join(self._getExtraPath(), 'config.xmd')
        existsConfig = exists(configfile)
        if existsConfig:
            md = xmipp.MetaData('properties@' + configfile)
            configobj = md.firstObject()
            pickingState = md.getValue(xmipp.MDL_PICKING_STATE, configobj)
            particleSize = md.getValue(xmipp.MDL_PICKING_PARTICLE_SIZE, configobj)
            activeMic = md.getValue(xmipp.MDL_MICROGRAPH, configobj)
            isAutopick = pickingState != "Manual"
            manualParticlesSize = md.getValue(xmipp.MDL_PICKING_MANUALPARTICLES_SIZE, configobj)
            autoParticlesSize = md.getValue(xmipp.MDL_PICKING_AUTOPARTICLES_SIZE, configobj)
            
            summary.append("Manual particles picked: %s" % manualParticlesSize)
            summary.append("Particle size:%d" %(particleSize))
            autopick = "Yes" if isAutopick else "No"
            summary.append("Autopick: " + autopick)
            if isAutopick:
                summary.append("Automatic particles picked: %s" % autoParticlesSize)
            summary.append("Last micrograph: " + activeMic)
        return "\n".join(summary)

    def getMethods(self, output):
        msg = 'Program picked %d particles of size %d ' % (output.getSize(), output.getBoxSize())
        msg += 'using training from %s. ' % self.xmippParticlePicking.get().getNameId()
        msg += 'For more detail see [Abrishami2013]'
        return msg

    def _citations(self):
        return ['Abrishami2013']
    
    # --------------------------- UTILS functions --------------------------------------------
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
        a new set of micrographs with the same properties than a previous trained
        ones. )
        """ 
        if self.getInputMicrographsPointer() is not None:
            return self.getInputMicrographsPointer().get()
        else:
            return None
