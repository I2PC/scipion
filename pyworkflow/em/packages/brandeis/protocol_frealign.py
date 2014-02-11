# **************************************************************************
# *
# * Authors:     Josue Gomez Blanco (jgomez@cnb.csic.es)
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
# *  e-mail address 'jgomez@cnb.csic.es'
# *
# **************************************************************************
"""
This module contains the protocol to obtain a refined 3D recontruction from a set of particles using Frealign
"""

from pyworkflow.em import *
from pyworkflow.utils import *
import brandeis
from data import *
from constants import *
import os
import time


class ProtFrealign(ProtRefine3D):
    """Protocol to perform a volume from a SetOfParticles
    using the frealign program"""
    _label = 'frealign'
    _references = ['[[http://dx.doi.org/10.1016/j.jsb.2006.05.004][Grigorieff N,  JSB (2007)]]',
                   '[[http://www.ncbi.nlm.nih.gov/pubmed/16384646][Wolf M, et.al, Ultramicroscopy (2006)]]',
                   '[[http://www.ncbi.nlm.nih.gov/pubmed/15556702][Stewart A & Grigorieff N, Ultramicroscopy (2004)]]',
                   '[[http://www.ncbi.nlm.nih.gov/pubmed/9571020][Grigorieff N, JMB (1998)]]',
                   '[[http://www.sciencedirect.com/science/article/pii/S104784771200144X][Sindelar CV & Grigorieff N, Ultramicroscopy (2012)]]',
                   '[[http://www.sciencedirect.com/science/article/pii/S1047847713001858][Lyumkis D, et. al, JSB (2013)]]'
                    ]


    def __init__(self, **args):
        ProtRefine3D.__init__(self, **args)
        self.stepsExecutionMode = STEPS_PARALLEL

    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputParticles', PointerParam, label="Input particles", important=True,
                      pointerClass='SetOfParticles',pointerCondition='hasCTF',
                      help='Select the input particles.\n')  

        form.addParam('input3DReferences', PointerParam,
                      pointerClass='SetOfVolumes',
                      label='Initial 3D reference volume:', 
                      help='Input 3D reference reconstruction.\n')

        form.addParam('numberOfIterations', IntParam, default=10,
                      label='Number of iterations:',
                      help='Number of iterations to perform.')

        form.addParam('useInitialAngles', BooleanParam, default=False,
                      label="Use initial angles/shifts ? ", 
                      help='Set to *Yes* if you want to use the projection assignment (angles/shifts) \n '
                      'associated with the input particles (hasProjectionAssigment=True)')

        form.addSection(label='Flow Control')

        form.addParam('Firstmode', EnumParam, condition='not useInitialAngles', choices=['Simple search & Refine', 'Search, Refine, Randomise'],
                      label="Operation mode for iteration 1:", default=brandeis.MOD2_SIMPLE_SEARCH_REFINEMENT,
                      display=EnumParam.DISPLAY_COMBO,
                      help='Parameter *IFLAG* in FREALIGN\n\n'
                           'This option is only for the iteration 1.\n'
                           '_Mode -3_: Simple search, Refine and create a parameter file for the set\n'
                           '_Mode -4_: Search, Refine, Randomize and create a parameterfile for the set')

        form.addParam('mode', EnumParam, choices=['Recontruction only', 'Refinement', 'Random Search & Refine',
                                                   'Simple search & Refine', 'Search, Refine, Randomise'],
                      label="Operation mode:", default=brandeis.MOD_REFINEMENT,
                      display=EnumParam.DISPLAY_COMBO,
                      help='Parameter *IFLAG* in FREALIGN\n\n'
                           '_Mode 0_: Reconstruction only parameters as read in\n'
                           '_Mode 1_: Refinement & Reconstruction\n'
                           '_Mode 2_: Random Search & Refinement\n'
                           '_Mode 3_: Simple search & Refine\n'
                           '_Mode 4_: Search, Refine, Randomise & extend to *RREC*\n')

        form.addParam('doMagRefinement', BooleanParam, default=False,
                      label="Refine magnification?", expertLevel=LEVEL_EXPERT,
                      help='Parameter *FMAG* in FREALIGN\n\n'
                           'Set _True or False_ to enable/disable magnification refinement')

        form.addParam('doDefRefinement', BooleanParam, default=False,
                      label="Refine defocus?", expertLevel=LEVEL_ADVANCED,
                      help='Parameter *FDEF* in FREALIGN\n\n'
                           'Set _True of False_ to enable/disable defocus parameter refinement')

        form.addParam('doAstigRefinement', BooleanParam, default=False,
                      label="Refine astigmatism?", expertLevel=LEVEL_EXPERT,
                      help='Parameter *FASTIG* in FREALIGN\n\n'
                           'Set _True or False_ to enable/disable astigmatic angle refinement')

        form.addParam('doDefPartRefinement', BooleanParam, default=False,
                      label="Refine defocus for individual particles?", expertLevel=LEVEL_EXPERT,
                      help='Parameter *FPART* in FREALIGN\n\n'
                           'Set _True or False_ to enable/disable defocus parameter refinement\n'
                           'for individual particles. Otherwise defocus change is constrained\n'
                           'to be the same for all particles in one image\n')

        form.addParam('methodEwaldSphere', EnumParam, choices=['Disable', 'Simple', 'Reference-based', 'Simple with reversed handedness',
                                                               'Reference-based with reversed handedness'],
                      default=brandeis.EWA_DISABLE, expertLevel=LEVEL_EXPERT,
                      label="Ewald sphere correction", display=EnumParam.DISPLAY_COMBO,
                      help='Parameter *IEWALD* in FREALIGN\n\n'
                           'The options available to Ewald correction are:\n'
                           '_None_: No correction. *Option 0* for parameter *IEWALD* in FREALIGN\n\n'
                           '_Simple_: Do correction, simple insertion method. *Option 1*.\n'
                           '_Reference-based_: Do correction, reference-based method. *Option 2*.\n'
                           '_Simple with reversed handedness_: \n'
                           '   Do correction, simple insertion method with inverted handedness. *Option -1*\n'
                           '_Reference-based with reversed handedness_: \n'
                           '   Do correction, reference-based method with inverted handedness. *Option -2*')

        form.addParam('doExtraRealSpaceSym', BooleanParam, default=False,
                      label="Apply extra real space symmetry?",
                      help='Parameter *FBEAUT* in FREALIGN\n\n'
                           'Apply extra real space symmetry averaging \n'
                           'and masking to beautify final map just prior to output.')

        form.addParam('doWienerFilter', BooleanParam, default=False,
                      label="Apply Wiener filter?", expertLevel=LEVEL_EXPERT,
                      help='Parameter *FFILT* in FREALIGN\n\n'
                           'Apply single particle Wiener filter to final reconstruction.')

        form.addParam('doBfactor', BooleanParam, default=False,
                      label="Calculate and apply Bfactor?", expertLevel=LEVEL_EXPERT,
                      help='Parameter *FBFACT* in FREALIGN\n\n'
                           'Determine and apply B-factor to final reconstruction.')

        form.addParam('writeMatchProjections', BooleanParam, default=True,
                      label="Write matching projections?",
                      help='Parameter *FMATCH* in FREALIGN\n\n'
                           'Set _True or False_ to enable/disable output \n'
                           'of matching projections (for diagnostic purposes).')

        form.addParam('methodCalcFsc', EnumParam, choices=['calculate FSC', 'Calculate one 3DR with odd particles', 
                                                     'Calculate one 3DR with even particles',
                                                     'Calculate one 3DR with all particles'],
                      default=brandeis.FSC_CALC,
                      label="Calculation of FSC", display=EnumParam.DISPLAY_COMBO,
                      help='parameter *IFSC* in FREALIGN\n\n'
                           'Calculation of FSC table:\n'
                           '_Calculate FSC_: Internally calculate two reconstructions with odd and even \n'
                           '   numbered particles and generate FSC table at the end of the run.\n'
                           '   *Option 0* for parameter *IFSC* in FREALIGN\n\n'
                           'The following options reduce memory usage:\n'
                           '_Calculate one 3DR with odd particles_:\n'
                           '   Only calculate one reconstruction using odd particles. *Option 1*.\n'
                           '_Calculate one 3DR with even particles_:\n'
                           '   Only calculate one reconstruction using even particles. *Option 2*.\n'
                           '_Calculate one 3DR with all particles_:\n'
                           '   Only calculate one reconstruction using all particles. *Option 3*.')

        form.addParam('doAditionalStatisFSC', BooleanParam, default=True,
                      label="Calculate aditional statistics in FSC?",
                      help='Parameter *FSTAT* in FREALIGN\n\n'
                           'Calculate additional statistics in resolution table at the end \n'
                           '(*QFACT, SSNR, CC* and related columns). Setting *FSTAT* False saves over 50% of memory!.')

        form.addParam('paddingFactor', EnumParam, choices=['1', '2', '4'], default=brandeis.PAD_4,
                      label='Padding Factor', display=EnumParam.DISPLAY_COMBO,
                      help='Parameter *IBLOW* in FREALIGN\n\n'
                           'Padding factor for reference structure.\n'
                           'Padding factor 4 requires the most memory but results\n'
                           'in the fastest search & refinement.\n')

        form.addSection(label='General Parameters')
        
        form.addParam('innerRadius', FloatParam, default='0.0', 
                      label='Inner radius of reconstruction (A):', 
                      help='Parameter *RI* in FREALIGN\n\n'
                           'In Angstroms from centre of particle.\n'
                           'Enter the inner radius of the volume to be reconstructed.\n' 
                           'This is useful for reconstructions of viruses and other\n' 
                           'particles that might be hollow or have a disordered core.')
              
        form.addParam('outerRadius', FloatParam, default='108.0', 
                      label='Outer radius of reconstruction (A):', 
                      help='Parameter *RO* in FREALIGN\n\n'
                           'In Angstroms from centre of particle.\n'
                           'Enter the outer radius of the volume to be reconstructed.\n'
                           'The program will also apply a mask with a cosine edge to the particle image\n'
                           'before processing (done inside *CTFAPPLY* using  *HALFW=6* pixels for cosine bell).')
        
        form.addParam('ThresholdMask', FloatParam, default='0.0', 
                      label='Threshold to for masking the input 3D structure:', expertLevel=LEVEL_ADVANCED,
                      help='Parameter *XSTD* in FREALIGN\n\n'
                           'Filtered 3D model - note this 3D masking does not use RI.\n'
                           '- if positive, calculates mask with subroutine *D3MASK*, equiv to\n'
                           '  solvent flattening with 5-pixel-cosine-bell smoothed mask\n'
                           '  boundary.  The mask is then used to multiply the input 3D map,\n'
                           '  which is then used for all parameter refinement and subsequent\n'
                           '  calculations.\n'
                           '- if negative, calculates mask with subroutine *D2MASK* resulting\n'
                           '  in a sharp binary (0/1) mask boundary for which is used for\n'
                           '  both parameter refinement and reconstruction, and to mask and\n'
                           '  output the matching projections.  Each matching particle image\n'
                           '  is also always masked with a cosine bell edged function of\n'
                           '  radius RI.\n'
                           'If set 0, disables this function.')
        
        form.addParam('pseudoBFactor', FloatParam, default='100.0', 
                      label='Resol-Dependent weighting of particles for 3D reconstruction:',
                      help='Parameter *PBC* in FREALIGN\n\n'
                           'Automatic weighting is applied to each particle: a pseudo-temperature (B)\n'
                           'factor is applied to each particle according to its relative phase\n'
                           'residual against the reference. The weight is calculated as\n'
                           '          W = exp (-DELTAP/PBC * R^2)\n'
                           'with DELTAP = relative phase residual (actual phase residual minus BOFF),\n'
                           'PBC = conversion constant (5.0 in the example),\n'
                           'and R^2 the squared resolution in Fourier units (R = 0.0 ... 0.5).\n'
                           'A large value for PBC (e.g. 100.0) gives equal weighting to each particle.')

        form.addParam('avePhaseResidual', FloatParam, default='65.0', 
                      label='Average phase residual:',
                      help='Parameter *BOFF* in FREALIGN\n\n'
                           'Approximate average phase residual of all particles,\n'
                           ' used in calculating weights for contributions of different\n'
                           'particles to 3D map (see Grigorieff, 1998).')

        form.addParam('angStepSize', FloatParam, default='15.0', 
                      label='Angular step size for the angular search:',
                      help='Parameter *DANG* in FREALIGN\n\n'
                           'Angular step size for the angular search used in modes *IFLAG*=3,4.\n'
                           'There is a program default value calculated taking resolution into\n'
                           'account, but if this input value is non-zero, the program value is\n'
                           'overridden.')

        form.addParam('numberRandomSearch', FloatParam, default='10.0', 
                      label='Number of randomised search/refinement trials:',
                      help='Parameter *ITMAX* in FREALIGN\n\n'
                           'number of cycles of randomised search/refinement used in modes IFLAG=2,4\n'
                           'There is a program default value (10 cycles), but if this input value is\n'
                           'non-zero, the program value is overridden.\n')

        form.addParam('numberPotentialMatches', FloatParam, default='10.0', 
                      label='number of potential matches:',
                      help='Parameter *IPMAX* in FREALIGN\n\n'
                           'number of potential matches in a search that should be tested further in\n'
                           'a subsequent local refinement.\n')

        form.addParam('paramRefine', EnumParam, choices=['Refine all', 'Refine only euler angles', 'Refine only X & Y shifts'],
                      default=brandeis.REF_ALL,
                      label="Parameters to refine", display=EnumParam.DISPLAY_COMBO,
                      help='Parameter *MASK* in FREALIGN\n\n'
                           'Parameters to include in refinement')

        form.addParam('symmetry', TextParam, default='c1',
                      label='Point group symmetry:',
                      help='Parameter *ASYM* in FREALIGN\n\n'
                           'Specify the symmetry.Choices are: C(n),D(n),T,O,I,I1,I2 or N (can be zero)\n'
                           'n  = rotational symmetry required in pointgroup C(n) or D(n)\n'
                           'N  = number of symmetry matrices to read in.\n'
                           'T  = tetrahedral pointgroup 23\n'
                           'O  = octahedral pointgroup 432\n'
                           'I  = icosahedral 532 symmetry in setting 1 (5fold is on X)\n'
                           'I1 = also in setting 1 (X) - as used by Imagic\n'
                           'I2 = in setting 2 (Y) - as used by Crowther et. al')        

        form.addSection(label='3D Reconstruction')
        
        form.addParam('relMagnification', FloatParam, default='1.0', 
                      label='Relative magnification:',
                      help='Parameter *RELMAG* in FREALIGN\n\n'
                           'Relative magnification of data set. The magnification feature allows\n'
                           'for manual variations of magnification between data sets.')

        form.addParam('targetPhaseResidual', FloatParam, default='10.0',
                      label='Target phase residual:', expertLevel=LEVEL_EXPERT,
                      help='Parameter *TARGET* in FREALIGN\n\n'
                           'Target phase residual (for resolution between RMAX1 and RMAX2)\n'
                           'during parameter saerch and refinement, below which the search and/or\n'
                           'refinement is terminated.  This is normally set low (e.g. THRESH=10.0)\n'
                           'This will give excellent determination of particle orientations.')

        form.addParam('PhaseResidual', FloatParam, default='90.0', 
                      label='Phase residual cut-off:',
                      help='Parameter *THRESH* in FREALIGN\n\n'
                           'Any particles with a higher overall phase residual will not be\n'
                           'included in the reconstruction when IFLAG=0,1,2,3. This variable\n'
                           'is often used with IFLAG=0 in separate runs to calculate maps\n'
                           'using various values of THRESH to find the optimum value to produce\n'
                           'the best map as judged from the statistics.')

        form.addParam('beamTiltX', FloatParam, default='0.0',
                      label='Beam tilt in X direction (in mrad):', expertLevel=LEVEL_EXPERT,
                      help='Parameter *TX* in FREALIGN.')

        form.addParam('beamTiltY', FloatParam, default='0.0',
                      label='Beam tilt in Y direction (in mrad):', expertLevel=LEVEL_EXPERT,
                      help='Parameter *TY* in FREALIGN.')

        form.addParam('resolution', FloatParam, default='10.0', 
                      label='Resol. of reconstruction (A):',
                      help='Parameter *RREC* in FREALIGN\n\n'
                           'Resolution to which the reconstruction is calculated.\n'
                           'If several datasets have different values, the data is individually\n'
                           'limited in the summation to the RREC resolution but symmetry is\n'
                           'applied, statistics output and the final map calculated to the\n'
                           'maximum resolution requested for any dataset.')

        form.addParam('lowResolRefine', FloatParam, default='200.0', 
                      label='Low resolution in refinement (A):',
                      help='Parameter *RMAX1* in FREALIGN\n\n'
                           'Resolution of the data included in the search/refinement. These\n'
                           'two parameters (RMAX1,RMAX2) are very important.  The successful\n'
                           'alignment of particles depends critically on the signal-to-noise\n'
                           'ratio of thecross-correlation or phase residual calculation, and\n'
                           'exclusion of weak data at high resolution or spurious, very strong\n'
                           'artefacts at low resolution can make a big difference.  Success can\n'
                           'be judged by whether the X,Y coordinates of the particle centres are\n'
                           'reasonable.')

        form.addParam('highResolRefine', FloatParam, default='25.0', 
                      label='High resolution in refinement (A):',
                      help='Parameter *RMAX2* in FREALIGN\n\n'
                           'Resolution of the data included in the search/refinement. These\n'
                           'two parameters (RMAX1,RMAX2) are very important.  The successful\n'
                           'alignment of particles depends critically on the signal-to-noise\n'
                           'ratio of thecross-correlation or phase residual calculation, and\n'
                           'exclusion of weak data at high resolution or spurious, very strong\n'
                           'artefacts at low resolution can make a big difference.  Success can\n'
                           'be judged by whether the X,Y coordinates of the particle centres are\n'
                           'reasonable.')

        form.addParam('defocusUncertainty', FloatParam, default='200.0', 
                      label='Defocus uncertainty (A):', expertLevel=LEVEL_EXPERT,
                      help='Parameter *DFSIG* in FREALIGN\n\n'
                           'This will restrain the change in defocus when refining defocus values\n'
                           'for individual particles.')

        form.addParam('Bfactor', FloatParam, default='0.0', 
                      label='B-factor to apply to particle:', expertLevel=LEVEL_EXPERT,
                      help='Parameter *RBFACT* in FREALIGN\n\n'
                           'B-factor to apply to particle image projections before orientation\n'
                           'determination or refinement.  This allows inclusion of high resolution\n'
                           'data but with a reduced (or increased if negative) weight.  *WGH* and\n'
                           '*RBFAC* can be manipulated in particle parameter refinement as if they\n'
                           'were low pass and high pass filters.  *WGH* and the CTF are used to\n'
                           'correct the density in the final map, whereas *RBFAC* is not.\n\n'
                           '_NOTE_: This option should be used with great care as amplification of\n'
                           'noisy high-resolution terms can lead to over-fitting and artificially\n'
                           'high values in the FSC curve (se publication #2 above). FREALIGN uses an\n'
                           'automatic weighting scheme and RBFACT should normally be set to 0.0.')

        form.addParallelSection(threads=1, mpi=1)        
    
    def _insertAllSteps(self):
        """Insert the steps to refine orientations and shifts of the SetOfParticles
        """
       
        numberOfBlocks = self.numberOfThreads.get()
        depsRecons = []
        
        for iter in range(1, self.numberOfIterations.get() + 1):
            initId = self._insertFunctionStep('initIterStep', iter, numberOfBlocks, prerequisites=depsRecons)
            depsRefine = self._insertRefineIterStep(numberOfBlocks, iter, [initId])
            reconsId = self._insertFunctionStep("reconstructVolumeStep", iter, numberOfBlocks, prerequisites=depsRefine)
            depsRecons = [reconsId]
        self._insertFunctionStep("createOutput", prerequisites=depsRecons)
    
    def initIterStep(self, iter, numberOfBlocks):
        """ Prepare files and directories for the current iteration """
        
        self._createIterWorkingDir(iter) # create the working directory for the current iteration.
        iterDir = self._iterWorkingDir(iter)
        prevIter = iter - 1 
        prevIterDir = self._iterWorkingDir(prevIter)
        
        imgSet = self.inputParticles.get()
        imgFn = self._getTmpPath('particles.mrc')
        vol = self.input3DReferences.get()
        volFn = self._getTmpPath('volume.mrc')
        refVol = join(iterDir, 'reference_volume_iter_%03d.mrc' % iter) # reference volume of the step.
        iterVol = join(iterDir, 'volume_iter_%03d.mrc' % iter) # refined volume of the step
        prevIterVol = join(prevIterDir, 'volume_iter_%03d.mrc' % prevIter) # volume of the previous iteration
        
        if iter==1:
            imgSet.writeStack(imgFn) # convert the SetOfParticles into a mrc stack.
            vol.writeStack(volFn) # convert the reference volume into a mrc volume
            copyFile(volFn, refVol)  #Copy the initial volume in the current directory.
        else:
            self._splitParFile(iter, numberOfBlocks)
            copyFile(prevIterVol, refVol)   #Copy the reference volume as refined volume.
        copyFile(refVol, iterVol)   #Copy the reference volume as refined volume.
        
    def _insertRefineIterStep(self, numberOfBlocks, iter, depsInitId):
        """ execute the refinement for the current iteration """
        
        depsRefine = []
        iterDir = self._iterWorkingDir(iter)
        
        if iter == 1:
            if not self.useInitialAngles.get():
                stepConstructId = self._insertFunctionStep("constructParamFilesStep", numberOfBlocks, iterDir, iter, prerequisites=depsInitId)
                depsConstruct = [stepConstructId]
                for block in range(1, numberOfBlocks + 1):
                    refineId = self._insertFunctionStep("refineBlockStep", iterDir, block, prerequisites=depsConstruct)
                    depsRefine.append(refineId)
            else:
                # ToDo: Construct the function to extract the coordinates and euler angles for each particle from SetOfParticles.
                pass
        else:
            for block in range(1, numberOfBlocks + 1):
                refineId = self._insertFunctionStep("refineParticlesStep", iter, block, prerequisites=depsInitId)
                depsRefine.append(refineId)
        return depsRefine
    
    def _particlesPerBlock(self, numberOfBlocks, numberOfParticles):
        """ Return a list with numberOfBlocks values, each value will be
        the number of particles assigned to each block. The number of particles
        will be distributed equally between each block as possible.
        """
        restBlock = numberOfParticles % numberOfBlocks
        colBlock = numberOfParticles / numberOfBlocks
        # Create a list with the number of particles assigned
        # to each block, initially equally distributed
        blockParticles = [colBlock] * numberOfBlocks
        # Now assign the particles in the rest
        for i, v in enumerate(blockParticles):
            if i < restBlock:
                blockParticles[i] += 1
                
        return blockParticles
        
    def _particlesInBlock(self, block):
        """calculate the initial and final particles that belongs to this block"""
        
        imgSet = self.inputParticles.get()
        numberOfBlocks = self.numberOfThreads.get()
        
        blockParticles = self._particlesPerBlock(numberOfBlocks, imgSet.getSize())
        initPart = 0
        lastPart = 0
        for i in range(block):
            initPart = lastPart + 1
            lastPart = lastPart + blockParticles[i]
        particlesInilast = [initPart, lastPart]
        return particlesInilast
    
    def _getParamsIteration(self, imgSet, iter):
        """ Defining the current iteration
        """
        
        iterDir = self._iterWorkingDir(iter)
        samplingRate3DR = imgSet.getSamplingRate()
        paramDic = {}
        #Prepare arguments to call program fralign_v8.exe
        paramsDic = {'mode': self.mode.get(),
                        'useInitialAngles': self.useInitialAngles.get(),
                        'innerRadius': self.innerRadius.get(),
                        'outerRadius': self.outerRadius.get(),
                        'ThresholdMask': self.ThresholdMask.get(),
                        'pseudoBFactor': self.pseudoBFactor.get(),
                        'avePhaseResidual': self.avePhaseResidual.get(),
                        'angStepSize': self.angStepSize.get(),
                        'numberRandomSearch': self.numberRandomSearch.get(),
                        'numberPotentialMatches': self.numberPotentialMatches.get(),
                        'sym': self.symmetry.get(),
                        'relMagnification': self.relMagnification.get(),
                        'targetPhaseResidual': self.targetPhaseResidual.get(),
                        'PhaseResidual': self.PhaseResidual.get(),
                        'beamTiltX': self.beamTiltX.get(),
                        'beamTiltY': self.beamTiltY.get(),
                        'resol': self.resolution.get(),
                        'lowRes': self.lowResolRefine.get(),
                        'highRes': self.highResolRefine.get(),
                        'defocusUncertainty': self.defocusUncertainty.get(),
                        'Bfactor': self.Bfactor.get(),
                        'sampling3DR': samplingRate3DR
                       }
        
        # Get the particles stack
        imgsFn = os.path.relpath(self._getTmpPath('particles.mrc'), iterDir)
        paramsDic['imageFn'] = imgsFn
        acquisition = imgSet.getAcquisition()
        # Get the amplitude Contrast of the micrographs
        paramsDic['ampContrast'] = acquisition.getAmplitudeContrast()
        # Get the scanned pixel size of the micrographs
        paramsDic['scannedPixelSize'] = acquisition.getMagnification() * imgSet.getSamplingRate() / 10000
        # Get the voltage and spherical aberration of the microscope
        paramsDic['voltage'] = acquisition.getVoltage()
        paramsDic['sphericalAberration'] = acquisition.getSphericalAberration()

        # Defining the operation mode
        if self.mode == MOD_RECONSTRUCTION:
            paramsDic['mode'] = 0
        elif self.mode == MOD_REFINEMENT:
            paramsDic['mode'] = 1
        elif self.mode == MOD_RANDOM_SEARCH_REFINEMENT:
            paramsDic['mode'] = 2
        elif self.mode == MOD_SIMPLE_SEARCH_REFINEMENT:
            paramsDic['mode'] = 3
        else:
            paramsDic['mode'] = 4
        
        # Defining the operation mode for iteration 1.
        if self.Firstmode == MOD2_SIMPLE_SEARCH_REFINEMENT:
            paramsDic['mode2'] = -3
        else:
            paramsDic['mode2'] = -4

        # Defining if magnification refinement is going to do
        if self.doMagRefinement:
            paramsDic['doMagRefinement'] = 'T'
        else:
            paramsDic['doMagRefinement'] = 'F'
            
        # Defining if defocus refinement is going to do
        if self.doDefRefinement:
            paramsDic['doDefocusRef'] = 'T'
        else:
            paramsDic['doDefocusRef'] = 'F'

        # Defining if astigmatism refinement is going to do
        if self.doAstigRefinement:
            paramsDic['doAstigRef'] = 'T'
        else:
            paramsDic['doAstigRef'] = 'F'
        
        # Defining if defocus refinement for individual particles is going to do
        if self.doDefPartRefinement:
            paramsDic['doDefocusPartRef'] = 'T'
        else:
            paramsDic['doDefocusPartRef'] = 'F'
        
        if self.methodEwaldSphere == EWA_DISABLE:
            paramsDic['metEwaldSphere'] = 0
        elif self.methodEwaldSphere == EWA_SIMPLE:
            paramsDic['metEwaldSphere'] = 1
        elif self.methodEwaldSphere == EWA_REFERENCE:
            paramsDic['metEwaldSphere'] = 2
        elif self.methodEwaldSphere == EWA_SIMPLE_HAND:
            paramsDic['metEwaldSphere'] = -1
        else:
            paramsDic['metEwaldSphere'] = -2
        
        # Defining if apply extra real space symmetry
        if self.doExtraRealSpaceSym:
            paramsDic['doExtraRealSpaceSym'] = 'T'
        else:
            paramsDic['doExtraRealSpaceSym'] = 'F'
        
        # Defining if wiener filter is going to apply
        if self.doWienerFilter:
            paramsDic['doWienerFilter'] = 'T'
        else:
            paramsDic['doWienerFilter'] = 'F'
        
        # Defining if wiener filter is going to calculate and apply
        if self.doBfactor:
            paramsDic['doBfactor'] = 'T'
        else:
            paramsDic['doBfactor'] = 'F'
        
        # Defining if matching projections is going to write
        if self.writeMatchProjections:
            paramsDic['writeMatchProj'] = 'T'
        else:
            paramsDic['writeMatchProj'] = 'F'
        
        # Defining the method to FSC calcutalion
        if self.methodCalcFsc == FSC_CALC:
            paramsDic['metFsc'] = 0
        elif self.methodCalcFsc == FSC_3DR_ODD:
            paramsDic['metFsc'] = 1
        elif self.methodCalcFsc == FSC_3DR_EVEN:
            paramsDic['metFsc'] = 2
        elif self.methodCalcFsc == FSC_3DR_ALL:
            paramsDic['metFsc'] = 3
        
        
        if self.doAditionalStatisFSC:
            paramsDic['doAditionalStatisFSC'] = 'T'
        else:
            paramsDic['doAditionalStatisFSC'] = 'F'
            
        if self.paddingFactor == PAD_1:
            paramsDic['padFactor'] = 1
        elif self.paddingFactor == PAD_2:
            paramsDic['padFactor'] = 2
        else:
            paramsDic['padFactor'] = 4
            
        if self.paramRefine == REF_ALL:
            paramsDic['paramRefine'] = '1, 1, 1, 1, 1'
        elif self.paramRefine == REF_ANGLES:
            paramsDic['paramRefine'] = '1, 1, 1, 0, 0'
        else:
            paramsDic['paramRefine'] = '0, 0, 0, 1, 1'
        
        return paramsDic
    
    def _createIterWorkingDir(self, iter):
        """create a new directory for the iterarion and change to this directory.
        """
        workDir = self._iterWorkingDir(iter)
        makePath(workDir)   # Create a directory for a current iteration
    
    def _iterWorkingDir(self, iter):
        """ Define which is the directory for the current iteration"""
        iterDir = 'iter_%03d' % iter
        workDir = self._getExtraPath(iterDir)
        return workDir

    def reconstructVolumeStep(self, iter, numberOfBlocks):
        """Reconstruct a volume from a SetOfParticles with its current parameters refined
        """

        imgSet = self.inputParticles.get()
        self._mergeAllParFiles(iter, numberOfBlocks)  # merge all parameter files generated in a refineIterStep function.
        
        initParticle = 1
        finalParticle = imgSet.getSize()
        params = self._getParamsIteration(imgSet, iter)
        
        params['outputParFn'] = 'output_param_file_%06d' % initParticle + '_%06d_' % finalParticle + 'iter_%03d.par' % iter
        params['initParticle'] = initParticle
        params['finalParticle'] = finalParticle

        params2 = self._setParams3DR(iter)
        
        params3DR = dict(params.items() + params2.items())
        
        args = self._prepareCommand()
        self.runJob(self._program, args % params3DR)
        self._leaveDir()
    
    def refineParticlesStep(self, iter, block):
        """Only refine the parameters of the SetOfParticles
        """
        
        param = {}
        imgSet = self.inputParticles.get()
        
        iterDir = self._iterWorkingDir(iter)
        if block==1:
            self._enterDir(iterDir) # enter to the working directory for the current iteration.
        
        iniPart, lastPart = self._particlesInBlock(block)
        prevIter = iter - 1
        param['inputParFn'] = 'particles_%02d_' % block + 'iter_%03d.par' % prevIter
        param['initParticle'] = iniPart
        param['finalParticle'] = lastPart
        param['frealignOut'] = 'frealign_Output_%02d.log' % block

        paramDic = self._setParamsRefineParticles(iter, block)
        initParamsDict = self._getParamsIteration(imgSet, iter)
        
        paramsRefine = dict(initParamsDict.items() + paramDic.items() + param.items())
        args = self._prepareCommand()
        
        self.runJob(self._program, args % paramsRefine)
        
    def _setParamsRefineParticles(self, iter, block):
        paramDics = {}
        paramDics['stopParam'] = -100
        paramDics['volume'] = 'reference_volume_iter_%03d.mrc' % iter
        paramDics['outputParFn'] = 'particles_%02d_' % block + 'iter_%03d.par' % iter
        paramDics['inputParFn'] = paramDics['outputParFn']
        paramDics['imgFnMatch'] = 'particles_match_%02d_' % block + 'iter_%03d.mrc' % iter
        paramDics['outputShiftFn'] = 'particles_shifts_%02d_' % block + 'iter_%03d.shft' % iter
        paramDics['3Dweigh'] = 'volume_weights_iter_%02d_' % block + 'iter_%03d' % iter
        paramDics['FSC3DR1'] = 'volume_1_%02d_' % block + 'iter_%03d' % iter
        paramDics['FSC3DR2'] = 'volume_2_%02d_' % block + 'iter_%03d' % iter
        paramDics['VolPhResidual'] = 'volume_phasediffs_%02d_' % block + 'iter_%03d' % iter
        paramDics['VolpointSpread'] = 'volume_pointspread_%02d_' % block + 'iter_%03d' % iter
        return paramDics
        
    def _setParams3DR(self, iter):
        """ Setting the parameters to reconstruct a new 3DR"""
        paramDic = {}
        paramDic['mode'] = 0
        paramDic['stopParam'] = 0   #The stopParam must be 0 if you want obtain a 3D reconstruction.
        paramDic['volume'] = 'volume_iter_%03d.mrc' % iter
        paramDic['inputParFn'] = 'particles_iter_%03d.par' % iter
        paramDic['imgFnMatch'] = 'particles_match_iter_%03d.mrc' % iter
        paramDic['outputShiftFn'] = 'particles_shifts_iter_%03d.shft' % iter
        paramDic['3Dweigh'] = 'volume_weights_iter_%03d' % iter
        paramDic['FSC3DR1'] = 'volume_1_iter_%03d.mrc' % iter
        paramDic['FSC3DR2'] = 'volume_2_iter_%03d.mrc' % iter
        paramDic['VolPhResidual'] = 'volume_phasediffs_iter_%03d' % iter
        paramDic['VolpointSpread'] = 'volume_pointspread_iter_%03d' % iter
        paramDic['frealignOut'] = 'fralign_volume_reconstruct_iter_%03d' % iter
        return paramDic
        
    def __openParamFile(self, blockNumber, paramsDict):
        """ Open the file and write the first part of the block param file. """
        if which('frealign_v8.exe') is '':
            raise Exception('Missing frealign_v8.exe')
        initaLines = """frealign_v8.exe << eot > %(frealignOut)s
M,%(mode2)s,%(doMagRefinement)s,%(doDefocusRef)s,%(doAstigRef)s,%(doDefocusPartRef)s,%(metEwaldSphere)s,%(doExtraRealSpaceSym)s,%(doWienerFilter)s,%(doBfactor)s,%(writeMatchProj)s,%(metFsc)s,%(doAditionalStatisFSC)s,%(padFactor)s
%(outerRadius)s,%(innerRadius)s,%(sampling3DR)s,%(ampContrast)s,%(ThresholdMask)s,%(pseudoBFactor)s,%(avePhaseResidual)s,%(angStepSize)s,%(numberRandomSearch)s,%(numberPotentialMatches)s
%(paramRefine)s
%(initParticle)s,%(finalParticle)s
%(sym)s
%(relMagnification)s,%(scannedPixelSize)s,%(targetPhaseResidual)s,%(PhaseResidual)s,%(sphericalAberration)s,%(voltage)s,%(beamTiltX)s,%(beamTiltY)s
%(resol)s,%(lowRes)s,%(highRes)s,%(defocusUncertainty)s,%(Bfactor)s
%(imageFn)s
%(imgFnMatch)s
"""
        paramFile = 'block%03d.sh' % blockNumber
        f = open(paramFile, 'w+')
        f.write(initaLines % paramsDict)
        return f
    
    def __writeParamParticle(self, f, particleLine):
        """ Write a particle line to the param file """
        f.write(particleLine)
    
    def __closeParamFile(self, f, paramsDict):
        """ Close the param file for a block. """
        finaLines = """%(outputParFn)s
%(outputShiftFn)s
%(stopParam)s, 0., 0., 0., 0., 0., 0., 0.
%(volume)s
%(3Dweigh)s
%(FSC3DR1)s
%(FSC3DR2)s
%(VolPhResidual)s
%(VolpointSpread)s
eot
"""
        f.write(finaLines % paramsDict)
        f.close()
        
    def constructParamFilesStep(self, numberOfBlocks, iterDir, iter):
        """ This function construct a parameter file (.par) with the information of the SetOfParticle.
        This function will execute only in iteration 1.
        """
        
        self._enterDir(iterDir)
        
        imgSet = self.inputParticles.get()
        magnification = imgSet._acquisition.magnification.get()
        blockParticles = self._particlesPerBlock(numberOfBlocks, imgSet.getSize())
        initParamsDict = self._getParamsIteration(imgSet, iter)
        params = {}
        lastPart = 0
        for block in range(numberOfBlocks):
            more = 1
            params['initParticle'] = 1 + lastPart
            lastPart = lastPart + blockParticles[block]
            params['finalParticle'] = lastPart
            numberOfBlock = block + 1
            params['frealignOut'] = 'frealign_Output_%02d.log' % numberOfBlock
            paramDic = self._setParamsRefineParticles(iter, numberOfBlock)
            paramsRefine = dict(initParamsDict.items() + params.items() + paramDic.items())
            f = self.__openParamFile(block + 1, paramsRefine)
            
            # ToDo: Implement a better method to get the info particles. Now, you iterate several times over the SetOfParticles (as many threads as you have)
            for i, img in enumerate(imgSet):
                film = img.getObjId()
                ctf = img.getCTF()
                defocusU, defocusV, astig = ctf.getDefocusU(), ctf.getDefocusV(), ctf.getDefocusAngle()
                partCounter = i + 1
                
                if partCounter == lastPart: # The last particle in the block
                    more = 0
                particleLine = '1, %05d,' % magnification + ' %05d,'% film + ' %05f,'% defocusU + ' %05f,'% defocusV + ' %02f,'% astig + ' %01d\n' % more
                self.__writeParamParticle(f, particleLine)
                
                if more == 0: # close the block.
                    self.__closeParamFile(f, paramsRefine)
                    break
        self._leaveDir()
                
    def refineBlockStep(self, iterDir, block):
        """ This function execute the bash script for refine a subset(block) of images.
        It will enter in the iteration dir and execute the script there. 
        """
        if block == 1:
            self._enterDir(iterDir)
            
        # time.sleep(block)
        program = "./block%03d.sh" % block
        os.chmod(program, 0775)
        self.runJob(program, "")
        
        
        #self._leaveDir()
        
    def _prepareCommand(self):
        """ prepare the command to execute"""

        if which('frealign_v8.exe') is '':
            raise Exception('Missing frealign_v8.exe')
        self._program = which('frealign_v8.exe')
        args = """  << eot >> %(frealignOut)s
M,%(mode)s,%(doMagRefinement)s,%(doDefocusRef)s,%(doAstigRef)s,%(doDefocusPartRef)s,%(metEwaldSphere)s,%(doExtraRealSpaceSym)s,%(doWienerFilter)s,%(doBfactor)s,%(writeMatchProj)s,%(metFsc)s,%(doAditionalStatisFSC)s,%(padFactor)s
%(outerRadius)s,%(innerRadius)s,%(sampling3DR)s,%(ampContrast)s,%(ThresholdMask)s,%(pseudoBFactor)s,%(avePhaseResidual)s,%(angStepSize)s,%(numberRandomSearch)s,%(numberPotentialMatches)s
%(paramRefine)s
%(initParticle)s,%(finalParticle)s
%(sym)s
%(relMagnification)s,%(scannedPixelSize)s,%(targetPhaseResidual)s,%(PhaseResidual)s,%(sphericalAberration)s,%(voltage)s,%(beamTiltX)s,%(beamTiltY)s
%(resol)s,%(lowRes)s,%(highRes)s,%(defocusUncertainty)s,%(Bfactor)s
%(imageFn)s
%(imgFnMatch)s
%(inputParFn)s
%(outputParFn)s
%(outputShiftFn)s
%(stopParam)s, 0., 0., 0., 0., 0., 0., 0.
%(volume)s
%(3Dweigh)s
%(FSC3DR1)s
%(FSC3DR2)s
%(VolPhResidual)s
%(VolpointSpread)s
eot
"""
        return args
    
    def _mergeAllParFiles(self, iter, numberOfBlocks):
        """ This method merge all parameters files that has been created in a refineIterStep """
        
        if numberOfBlocks != 1:
            #iterDir = self._iterWorkingDir(iter)
            file2 = 'particles_iter_%03d.par' % iter
            f2 = open(file2, 'w+')
            
            for block in range(1, numberOfBlocks + 1):
                file1 = 'particles_%02d_' % block + 'iter_%03d.par' % iter
                f1 = open(file1)
                
                if block == 1:
                    lines = f1.readlines()
                    f2.writelines(lines[:-1])
                else:
                    for l in f1:
                        if not l.startswith('C'):
                            f2.write(l)
                f1.close()
            f2.close()
    
    def _splitParFile(self, iter, numberOfBlocks):
        """ This method split the parameter files that has been previosuly merged """
        
        if numberOfBlocks != 1:
            prevIter = iter -1
            prevIterDir = self._iterWorkingDir(prevIter)
            iterDir = self._iterWorkingDir(iter)
            file1 = join(prevIterDir, 'particles_iter_%03d.par' % prevIter)
            
            for block in range(1, numberOfBlocks + 1):
                f1 = open(file1)
                file2 = join(iterDir, 'particles_%02d_' % block + 'iter_%03d.par' % prevIter)
                f2 = open(file2, 'w+')
                initpart, finalPart = self._particlesInBlock(block)
                
                for l in f1:
                    
                    if l.startswith('C'):
                        f2.write(l)
                    else:
                        split = l.split()
                        numPart = int(''.join(split[:1]))
                        
                        if numPart <= finalPart:
                            f2.write(l)
                        else:
                            break
                f2.close()
                f1.close()

    
    def createOutput(self):
        
        lastIter = self.numberOfIterations.get()
        lastIterDir = self._iterWorkingDir(lastIter)
        volFn = join(lastIterDir, 'volume_iter_%03d.mrc' % lastIter)
        vol = Volume()
        vol.setSamplingRate(self.inputParticles.get().getSamplingRate())
        vol.setFileName(volFn)
        self._defineOutputs(outputVolume=vol)
        
    def _validate(self):
        errors = []
        if which('frealign_v8.exe') is '':
            errors.append('Missing frealign_v8.exe')
        return errors
    
    def _summary(self):
        summary = []
        summary.append("Input particles:  %s" % self.inputParticles.get().getNameId())
        summary.append("Input volumes:  %s" % self.input3DReferences.get().getNameId())
        
        if not hasattr(self, 'outputVolumes'):
            summary.append("Output volumes not ready yet.")
        else:
            summary.append("Output volumes: %s" % self.vol.get())
        
        return summary
