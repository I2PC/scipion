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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************
"""
This sub-package contains the XmippParticlePickingAutomatic protocol
"""

from pyworkflow.em import *  
from pyworkflow.utils.path import *  
import xmipp
from xmipp3 import XmippProtocol

from convert import createXmippInputMicrographs, readSetOfCoordinates

MICS_SAMEASPICKING = 0
MICS_OTHER = 1

class XmippParticlePickingAutomatic(ProtParticlePicking, XmippProtocol):
    """Protocol to pick particles automatically in the project"""
    _label = 'automatic picking'
    
    filesToCopy = ['model_training.txt', 'model_svm.txt', 'model_pca_model.stk', 'model_rotpca_model.stk', 
               'model_particle_avg.xmp', 'config.xmd', 'templates.stk']
    
    
    def __init__(self, **args):        
        ProtParticlePicking.__init__(self, **args)
        
    def _defineParams(self, form):
    
        form.addSection(label='Input')
        form.addParam('xmippParticlePicking', PointerParam, label="Xmipp particle picking run",
                      pointerClass='XmippProtParticlePicking', pointerCondition='isFinished',
                      help='Select the previous xmipp particle picking run.')
        form.addParam('micsToPick', EnumParam, choices=['Same as supervised', 'Other'], 
                      default=0, label='Micrographs to pick', display=EnumParam.DISPLAY_LIST,
                      help='Select from which set of micrographs to pick using the training from supervised run.' 
                      'If you use Same as supervised, the same set of micrographs used for training the picker '
                      'will be used at this point. If you select Other, you can select another set of micrograph'
                      '(normally from the same specimen) and pick them completely automatic using the trainned picker.')
        form.addParam('inputMicrographs', PointerParam, label="Micrographs",
                      pointerClass='SetOfMicrographs', condition='micsToPick == 1',
                      help='Select other set of micrographs to pick using the trained picker.')
        form.addParam('memory', FloatParam, default=2,
                   label='Memory to use (In Gb)', expertLevel=2)        
        
    def _defineSteps(self):
        """The Particle Picking proccess is realized for a set of micrographs"""
        
        # Get pointer to input micrographs 
        self.particlePickingRun = self.xmippParticlePicking.get()

        # Get particle picking boxsize from the previous run
        boxSize = self.particlePickingRun.outputCoordinates.getBoxSize()
        
        modelRoot = self._getExtraPath('model')
        
        # Copy training model files to current run
        for file in self.filesToCopy:
            copyFile(self.particlePickingRun._getExtraPath(file), self._getExtraPath(file))
        
        # Get micrographs to pick
        if self.micsToPick.get() == MICS_SAMEASPICKING:
            self.micrographs = self.particlePickingRun.inputMicrographs.get()
            # Store in DB inputMicrographs pointer to be used on summary
            self.inputMicrographs.set(self.particlePickingRun.inputMicrographs.get())
        else:
            self.micrographs = self.inputMicrographs.get()
            
        for mic in self.micrographs:
            micPath = mic.getFileName()
            micName = removeBaseExt(micPath)
            proceed = True
            if self.micsToPick.get() == MICS_SAMEASPICKING:
                fnPos = join(self.particlePickingRun._getExtraPath(), replaceBaseExt(micPath, "pos"))
                if exists(fnPos):
                    blocks = xmipp.getBlocksInMetaDataFile(fnPos)
                    copy = True
                    if 'header' in blocks:
                        mdheader = xmipp.MetaData("header@" + fnPos)
                        state = mdheader.getValue(xmipp.MDL_PICKING_MICROGRAPH_STATE, mdheader.firstObject())
                        if state == "Available":
                            copy = False
                    if copy:
                        # Copy manual .pos file of this micrograph
                        copyFile(fnPos, join(self._getExtraPath(), basename(fnPos)))
                        proceed = False            
    
            if proceed:
                oroot = self._getExtraPath(micName)
                cmd = "-i %(micPath)s --particleSize %(boxSize)d --model %(modelRoot)s --outputRoot %(oroot)s --mode autoselect" % locals()
                #TODO: What is this?
#                if self.Fast:
#                    cmd += " --fast "
                self._insertRunJobStep("xmipp_micrograph_automatic_picking", cmd)
                    
        # Insert step to create output objects       
        self._insertFunctionStep('createOutput')
            
    def createOutput(self):
        posDir = self._getExtraPath()
        coordSet = self._createSetOfCoordinates()
        coordSet.setMicrographs(self.micrographs)
        readSetOfCoordinates(posDir, self.micrographs, coordSet)
        coordSet.write()
        self._defineOutputs(outputCoordinates=coordSet)
        self._defineDataSource(self.micrographs, coordSet)
        
    def _validate(self):
        validateMsgs = []
        srcPaths = [self.xmippParticlePicking.get()._getExtraPath(k) for k in self.filesToCopy]
        print srcPaths
        # Check that all needed files exist
        if missingPaths(*srcPaths):
            validateMsgs.append('Input particle picking run is not valid.')
            
        # If other set of micrographs is provided they should have same sampling rate and acquisition
        if self.micsToPick.get() == MICS_OTHER:
            if self.inputMicrographs.get().getSamplingRate() != self.xmippParticlePicking.get().inputMicrographs.get().getSamplingRate():
                validateMsgs.append('New micrographs should have same sampling rate as the ones already picked.')
            if not self.inputMicrographs.get().getAcquisition().equalAttributes(self.xmippParticlePicking.get().inputMicrographs.get().getAcquisition()):
                validateMsgs.append('New micrographs should have same acquisition parameters as the ones already picked.')
        return validateMsgs
    
    def _summary(self):
        summary = []
        if not hasattr(self, 'outputCoordinates'):
            summary.append("Output coordinates not ready yet.") 
        else:
            summary.append("Previous run: " + self.xmippParticlePicking.get().getNameId())
            summary.append("Number of particles picked: %d (from %d micrographs)" % (self.outputCoordinates.getSize(), self.inputMicrographs.get().getSize()))
        return summary