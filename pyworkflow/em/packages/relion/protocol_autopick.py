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

import os
from os.path import relpath
from pyworkflow.protocol.params import (PointerParam, FloatParam, RelationParam,
                                        IntParam, BooleanParam, LEVEL_ADVANCED, 
                                        LabelParam)
from pyworkflow.protocol.constants import STEPS_PARALLEL
from pyworkflow.em.protocol.protocol_particles import ProtParticlePicking
from pyworkflow.em.constants import RELATION_CTF

from protocol_base import ProtRelionBase
from convert import writeSetOfMicrographs, writeReferences, readSetOfCoordinates
from pyworkflow.em.convert import ImageHandler
import pyworkflow.em.metadata as md
import pyworkflow.utils as pwutils



class ProtRelionAutopickBase(ProtParticlePicking, ProtRelionBase):

    #--------------------------- INSERT steps functions -------------------------------------------- 
    def _insertAllSteps(self): 
        # Convert the input micrographs and references to 
        # the required Relion star files
        convertId = self._insertFunctionStep('convertInputStep', 
                                 self.getInputMicrographs().strId(),
                                 self.getInputReferences().strId())
        # Insert the picking steps for each micrograph separately 
        # (Relion requires in that way if micrographs have different dimensions)
        allMicSteps = []
        for mic in self.getInputMicrographs():
            autopickId = self._insertAutopickStep(mic, convertId)
            allMicSteps.append(autopickId)
        # Register final coordinates as output
        self._insertFunctionStep('createOutputStep', 
                                 prerequisites=allMicSteps)
        
    #--------------------------- STEPS functions --------------------------------------------
    def _preprocessMicrographRow(self, img, imgRow):
        # Temporarly convert the few micrographs to tmp and make sure
        # they are in 'mrc' format
        # Get basename and replace extension by 'mrc'
        newName = pwutils.replaceBaseExt(img.getFileName(), 'mrc')
        newPath = self._getExtraPath(newName)
        
        # If the micrographs are in 'mrc' format just create a link
        # if not, convert to 'mrc'
        if img.getFileName().endswith('mrc'):
            pwutils.createLink(img.getFileName(), newPath)
        else:
            self._ih.convert(img, newPath)
        # The command will be launched from the working dir
        # so, let's make the micrograph path relative to that
        img.setFileName(os.path.join('extra', newName))
        img.setCTF(self.ctfDict[img.getMicName()])
        
    def _getMicStarFile(self, mic):
        return self._getExtraPath(pwutils.replaceBaseExt(mic.getFileName(), 'star'))            
        
    def _postprocessMicrographRow(self, img, imgRow):
        imgRow.writeToFile(self._getMicStarFile(img))
            
    def convertInputStep(self, micsId, refsId):
        self._ih = ImageHandler() # used to convert micrographs
        # Match ctf information against the micrographs
        self.ctfDict = {}
        for ctf in self.ctfRelations.get():
            self.ctfDict[ctf.getMicrograph().getMicName()] = ctf.clone()
        
        micStar = self._getPath('input_micrographs.star')
        writeSetOfMicrographs(self.getInputMicrographs(), micStar, 
                              preprocessImageRow=self._preprocessMicrographRow,
                              postprocessImageRow=self._postprocessMicrographRow)
        
        writeReferences(self.getInputReferences(), self._getPath('input_references'))  
        
    def autopickMicrographStep(self, micStarFile, params, threshold, minDistance, fom):
        """ Launch the 'relion_autopick' for a micrograph with the given parameters. """
        # Call relion_autopick to allow picking of micrographs with different size
        params += ' --i %s' % relpath(micStarFile, self.getWorkingDir())
        params += ' --threshold %0.3f --min_distance %0.3f %s' % (threshold, minDistance, fom)
        
        self.runJob(self._getProgram('relion_autopick'), params, cwd=self.getWorkingDir())      

    #--------------------------- UTILS functions --------------------------------------------
    def getInputReferences(self):
        return self.inputReferences.get()
    
    def getInputMicrographs(self):
        return self.inputMicrographs.get()
        
    def getCoordsDir(self):
        return self._getTmpPath('xmipp_coordinates')
    
    def _writeXmippCoords(self, coordSet):
        micSet = self.getInputMicrographs()
        coordPath = self._getTmpPath('xmipp_coordinates')
        pwutils.cleanPath(coordPath)
        pwutils.makePath(coordPath)
        import pyworkflow.em.packages.xmipp3 as xmipp3
        micPath = os.path.join(coordPath, 'micrographs.xmd')
        xmipp3.writeSetOfMicrographs(micSet, micPath)
        xmipp3.writeSetOfCoordinates(coordPath, coordSet)
        return micPath, coordPath
        
    def writeXmippOutputCoords(self):
        return self._writeXmippCoords(self.outputCoordinates)
        
    def writeXmippCoords(self):
        """ Write the SetOfCoordinates as expected by Xmipp
        to be display with its GUI. 
        """
        micSet = self.getInputMicrographs()
        coordSet = self._createSetOfCoordinates(micSet)
        coordSet.setBoxSize(self.getInputReferences().getDim()[0])
        starFiles = [self._getExtraPath(pwutils.removeBaseExt(mic.getFileName()) + '_autopick.star')
                     for mic in micSet]
        readSetOfCoordinates(coordSet, starFiles)
        return self._writeXmippCoords(coordSet)
        
    
    
class ProtRelionAutopickFom(ProtRelionAutopickBase):
    """    
    This Relion protocol uses 2D class averages as templates to run the auto-picking 
    job-type. In this first stage, the auto-picking will be run just in few micrographs 
    to optimise two of its main parameters ( _Picking threshold_ and _Minimum inter-particle distance_).
    
    In order to save time, only 2 or 3 micrographs should be used with their CTF 
    information. One should use representative micrographs for the entire data set, 
    e.g. a high and a low-defocus one, and/or with thin or thick ice. 

    The expensive part of this calculation is to calculate a probability-based figure-of-merit 
    (related to the cross-correlation coefficient between each rotated reference and all positions 
    in the micrographs. This calculation is followed by a much cheaper peak-detection algorithm that 
    uses the threshold and minimum distance parameters mentioned above. Because these parameters 
    need to be optimised, this first stage of the auto-picking will write out so-called FOM maps.
    These are two large (micrograph-sized) files per reference. To avoid running into hard disc I/O 
    problems, the autopicking program can only be run sequentially (hence there is no option to use 
    more than one single MPI processor).
    """
    _label = 'auto-picking (1)'
    
    #--------------------------- DEFINE param functions --------------------------------------------   
    def _defineParams(self, form):
        form.addSection(label='Input')
        
        form.addParam('inputMicrographs', PointerParam, pointerClass='SetOfMicrographs',
                      label='Input micrographs (a few)', important=True,
                      help='Select a set with just a few micrographs to be used\n'
                           'in the auto-picking job. A few should be used in order\n'
                           'to perform the expensive calculation of the figure-of-merits.\n'
                           'writing out the so-called FOM maps.')
        form.addParam('ctfRelations', RelationParam, allowsNull=True,
                      relationName=RELATION_CTF, attributeName='getInputMicrographs',
                      label='CTF estimation',
                      help='Choose some CTF estimation related to input micrographs. \n')
        
        form.addSection('References')
        form.addParam('inputReferences', PointerParam, pointerClass='SetOfAverages',
                      label='Input references', important=True, 
                      help='Input references (SetOfAverages) for auto-pick. \n\n'
                           'Note that the absolute greyscale needs to be correct, \n'
                           'so only use images with proper normalization.')
        form.addParam('particleDiameter', IntParam, default=200,
                      label='Particle diameter (A)',
                      help='') 
        form.addParam('lowpassFilterRefs', IntParam, default=20,
                      label='Lowpass filter references (A)',
                      help='Lowpass filter that will be applied to the references \n'
                           'before template matching. \n'
                           'Do NOT use very high-resolution templates to search your micrographs. \n'
                           'The signal will be too weak at high resolution anyway,\n' 
                           'and you may find Einstein from noise....')
        form.addParam('angularSampling', IntParam, default=5,
                      label='Angular sampling (deg)',
                      help='Angular sampling in degrees for exhaustive searches \n'
                           'of the in-plane rotations for all references.')
        form.addParam('refsHaveInvertedContrast', BooleanParam, default=False,
                      label='References have inverted contrast',
                      help='Set to Yes to indicate that the reference have inverted \n'
                           'contrast with respect to the particles in the micrographs.')
        form.addParam('refsCtfCorrected', BooleanParam, default=True,
                      label='Are References CTF corrected?',
                      help='Set to Yes if the references were created with CTF-correction inside RELION.\n'
                           'If set to Yes, the input micrographs should contain the CTF information.')
        form.addParam('ignoreCTFUntilFirstPeak', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label='Ignore CTFs until first peak?',
                      help='Set this to Yes, only if this option was also used to generate the references.')
        
        form.addSection('Autopicking')
   
        form.addParam('fomLabel', LabelParam,  
                      important=True,
                      label='Auto-picking parameters will be optimized later interactively.')
        
    #--------------------------- INSERT steps functions --------------------------------------------  
        
    def getAutopickParams(self):
        params = ' --o autopick'
        params += ' --particle_diameter %d' % self.particleDiameter
        params += ' --angpix %0.3f' % self.getInputMicrographs().getSamplingRate()
        params += ' --ref input_references.star'
        
        if self.refsHaveInvertedContrast:
            params += ' --invert'
        
        if self.refsCtfCorrected:
            params += ' --ctf'
            
        params += ' --ang %d' % self.angularSampling
        params += ' --lowpass %d' % self.lowpassFilterRefs
        
        return params

    def autopickStep(self, threshold, minDistance, fom):
        """ This method is used from the wizard to optimize the parameters. """
        params = self.getAutopickParams()
        for mic in self.getInputMicrographs():
            self.autopickMicrographStep(self._getMicStarFile(mic), 
                                        params, threshold, minDistance, fom)  
                  
    def _insertAutopickStep(self, mic, convertId):
        """ Prepare the command line for calling 'relion_autopick' program """
        params = self.getAutopickParams()
        # Use by default the references size as --min_distance
        return self._insertFunctionStep('autopickMicrographStep', self._getMicStarFile(mic), 
                                        params, 0.25, 
                                        self.getInputDimA(), ' --write_fom_maps',
                                        prerequisites=[convertId])
        
    #--------------------------- STEPS functions --------------------------------------------

    def createOutputStep(self):
        self.summaryVar.set('This protocol does not generate any output.\n'
                            'The FOM maps were written to be used later to optimize \n'
                            'the _Threshold_ and _Inter-particle distance_ \n'
                            'parameters.')
    
    #--------------------------- INFO functions -------------------------------------------- 
    def _validate(self):
        """ Should be overriden in subclasses to 
        return summary message for NORMAL EXECUTION. 
        """
        errors = []
        if self.particleDiameter > self.getInputDimA():
            errors.append('Particle diameter (%d) can not be greater than size (%d)' % 
                          (self.particleDiameter, self.getInputDimA()))
        return errors
    
    def _summary(self):
        """ Should be overriden in subclasses to 
        return summary message for NORMAL EXECUTION. 
        """
        return [self.summaryVar.get('')]
    
    #--------------------------- UTILS functions --------------------------------------------
    def getInputDimA(self):
        """ Return the dimension of input references in A. """
        inputRefs = self.getInputReferences()
        if inputRefs is None:
            return None
        else:
            return inputRefs.getDim()[0] * inputRefs.getSamplingRate()
                
        
        
class ProtRelionAutopick(ProtRelionAutopickBase):
    """    
    This Relion protocol uses 2D class averages as templates to run the auto-picking 
    job-type. In this second stage, the protocol reads the FOM maps to optimize
    the 'Threshold' and 'Inter-particle distance'.
    """
    _label = 'auto-picking (2)'
    
    
    def __init__(self, **kwargs):        
        ProtRelionAutopickBase.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL
        
    #--------------------------- DEFINE param functions --------------------------------------------   
    def _defineParams(self, form):
        
        form.addSection(label='Input')
        
        form.addParam('inputAutopickFom', PointerParam, pointerClass='ProtRelionAutopickFom',
                      label='Input auto-pick FOM', important=True,
                      help='Select a previous auto-picking job that generated the FOM maps.')        
        form.addParam('inputMicrographs', PointerParam, pointerClass='SetOfMicrographs',
                      label='Input micrographs (full set)', important=True,
                      help='Select the full set of micrographs to be picked automatically\n'
                           'with the adjusted parameters.')        
        form.addParam('ctfRelations', RelationParam, allowsNull=True,
                      relationName=RELATION_CTF, attributeName='getInputMicrographs',
                      label='CTF estimation',
                      help='Choose some CTF estimation related to input micrographs. \n')
        
        group = form.addGroup('Autopicking')
           
        group.addParam('fomLabel', LabelParam,  
                      important=True,
                      label="Modify parameters and click the 'wizard' to see output coordinates")        
        group.addParam('pickingThreshold', FloatParam, default=0.25,
                      label='Picking threshold',
                      help='Use lower thresholds to pick more particles (and more junk probably)')
        group.addParam('interParticleDistance', IntParam, default=100,
                      label='Minimum inter-particle distance (A)',
                      help='Particles closer together than this distance \n'
                           'will be consider to be a single cluster. \n'
                           'From each cluster, only one particle will be picked.')
        
        form.addParallelSection(threads=4, mpi=1)
        
    #--------------------------- INSERT steps functions --------------------------------------------  
        
    def _insertAutopickStep(self, mic, convertId):
        """ Prepare the command line for calling 'relion_autopick' program """
        params = self.getAutopickParams()
        # Use by default the references size as --min_distance
        return self._insertFunctionStep('autopickMicrographStep', self._getMicStarFile(mic), 
                                        params,
                                        self.pickingThreshold.get(),
                                        self.interParticleDistance.get(),
                                        '', # neither read or write fom
                                        prerequisites=[convertId])
                
    #--------------------------- STEPS functions --------------------------------------------
    
    def createOutputStep(self):
        micSet = self.getInputMicrographs()
        coordSet = self._createSetOfCoordinates(micSet)
        coordSet.setBoxSize(self.getInputReferences().getDim()[0])
        starFiles = [self._getExtraPath(pwutils.removeBaseExt(mic.getFileName()) + '_autopick.star')
                     for mic in micSet]
        readSetOfCoordinates(coordSet, starFiles)
        
        self._defineOutputs(outputCoordinates=coordSet)
    
    #--------------------------- INFO functions -------------------------------------------- 
    def _validate(self):
        """ Should be overriden in subclasses to 
        return summary message for NORMAL EXECUTION. 
        """
        errors = []
        return errors
    
    def _summary(self):
        """ Should be overriden in subclasses to 
        return summary message for NORMAL EXECUTION. 
        """
        return [self.summaryVar.get('')]
    
    #--------------------------- UTILS functions --------------------------------------------
    def getInputAutopick(self):
        return self.inputAutopickFom.get()
    
    def getInputReferences(self):
        return self.getInputAutopick().inputReferences.get()
    
    def getAutopickParams(self):
        return self.getInputAutopick().getAutopickParams()
