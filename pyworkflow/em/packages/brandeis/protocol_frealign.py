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
This module contains the protocol for CTF estimation with ctffind3
"""

from pyworkflow.em import *
from pyworkflow.utils import *
import brandeis
from data import *
from constants import *


class BrandeisDefFrealign(Form):
    """Create the definition of parameters for
    the Frealign protocol"""
    def __init__(self):
        Form.__init__(self)
        
        self.addSection(label='Input')

        self.addParam('inputParticles', PointerParam, label="Input particles", important=True,
                      pointerClass='SetOfParticles',
                      help='Select the input particles.\n')  

        self.addParam('input3DReference', PointerParam,
                      pointerClass='Volume',
                      label='Initial 3D reference volume:', 
                      help='Input 3D reference reconstruction.\n')
        
        self.addParam('Firstmode', EnumParam, choices=['Simple search & Refine', 'Search, Refine, Randomise'],
                      label="Operation mode for iteration 1:", default=brandeis.MOD2_SIMPLE_SEARCH_REFINEMENT,
                      display=EnumParam.DISPLAY_COMBO,
                      help='Parameter <IFLAG> in FREALIGN\n'
                           'This option is only for the iteration 1'
                           'Mode -3: Simple search, Refine and create a parameter file for the set\n'
                           'Mode -4: Search, Refine, Randomise and create a parameterfile for the set')

        self.addParam('mode', EnumParam, choices=['Recontruction only', 'Refinement', 'Random Search & Refine',
                                                   'Simple search & Refine', 'Search, Refine, Randomise'],
                      label="Operation mode:", default=brandeis.MOD_REFINEMENT,
                      display=EnumParam.DISPLAY_COMBO,
                      help='Parameter <IFLAG> in FREALIGN'
                           'Mode 0: Reconstruction only parameters as read in\n'
                           'Mode 1: Refinement & Reconstruction\n'
                           'Mode 2: Random Search & Refinement\n'
                           'Mode 3: Simple search & Refine\n'
                           'Mode 4: Search, Refine, Randomise & extend to RREC\n')

        self.addParam('useInitialAngles', BooleanParam, default=False,
                      label="Use initial angles/shifts ? ", 
                      help='Set to <Yes> if you want to use the projection assignment (angles/shifts)\n'
                      'associated with the input particles (hasProjectionAssigment=True)')

        self.addParam('numberOfIterations', IntParam, default=4,
                      label='Number of iterations',
                      help='Number of iterations to perform.')

        self.addParam('doMagRefinement', BooleanParam, default=False,
                      label="Refine magnification?", expertLevel=LEVEL_EXPERT,
                      help='Set True or False to enable/disable magnification refinement\n'
                           'Parameter <FMAG> in FREALIGN')

        self.addParam('doDefRefinement', BooleanParam, default=False,
                      label="Refine defocus?", expertLevel=LEVEL_ADVANCED,
                      help='Set True of False to enable/disable defocus parameter refinement\n'
                           'Parameter <FDEF> in FREALIGN')

        self.addParam('doAstigRefinement', BooleanParam, default=False,
                      label="Refine astigmatism?", expertLevel=LEVEL_EXPERT,
                      help='Set True or False to enable/disable astigmatic angle refinement\n'
                           'Parameter <FASTIG> in FREALIGN')

        self.addParam('doDefPartRefinement', BooleanParam, default=False,
                      label="Refine defocus for individual particles?", expertLevel=LEVEL_EXPERT,
                      help='Set True of False to enable/disable defocus parameter refinement \n'
                           'for individual particles. Otherwise defocus change is constrained \n'
                           'to be the same for all particles in one image\n'
                           'Parameter <FPART> in FREALIGN')

        self.addParam('methodEwaldSphere', EnumParam, choices=['Disable', 'Simple', 'Reference-based', 'Simple with reversed handedness',
                                                               'Reference-based with reversed handedness'],
                      default=brandeis.EWA_DISABLE, expertLevel=LEVEL_EXPERT,
                      label="Ewald sphere correction", display=EnumParam.DISPLAY_COMBO,
                      help='Ewald correction are (parameter <IEWALD> in FREALIGN):\n'
                           'Disable: No correction. Option <0> for parameter <IEWALD> in FREALIGN.\n'
                           'Simple: Do correction, simple insertion method. Option <1>.\n'
                           'Reference-based: Do correction, reference-based method. Option <2>.\n'
                           'Simple with reversed handedness: \n'
                           '   Do correction, simple insertion method with inverted handedness. Option -1\n'
                           'Reference-based with reversed handedness: \n'
                           '   Do correction, reference-based method with inverted handedness. Option <-2>')

        self.addParam('doExtraRealSpaceSym', BooleanParam, default=False,
                      label="Apply extra real space symmetry?",
                      help='Apply extra real space symmetry averaging \n'
                           'and masking to beautify final map just prior to output.\n'
                           'Parameter <FBEAUT> in FREALIGN')

        self.addParam('doWienerFilter', BooleanParam, default=False,
                      label="Apply Wiener filter?", expertLevel=LEVEL_EXPERT,
                      help='Apply single particle Wiener filter to final reconstruction.\n'
                           'Parameter <FFILT> in FREALIGN')

        self.addParam('doBfactor', BooleanParam, default=False,
                      label="Calculate and apply Bfactor?", expertLevel=LEVEL_EXPERT,
                      help='Determine and apply B-factor to final reconstruction.\n'
                           'Parameter <FBFACT> in FREALIGN')

        self.addParam('writeMatchProjections', BooleanParam, default=True,
                      label="Write matching projections?",
                      help='Set True or False to enable/disable output \n'
                           'of matching projections (for diagnostic purposes)\n'
                           'Parameter <FMATCH> in FREALIGN')

        self.addParam('methodCalcFsc', EnumParam, choices=['calculate FSC', 'Calculate one 3DR with odd particles', 
                                                     'Calculate one 3DR with even particles',
                                                     'Calculate one 3DR with all particles'],
                      default=brandeis.FSC_CALC,
                      label="Calculation of FSC", display=EnumParam.DISPLAY_COMBO,
                      help='Calculation of FSC table (parameter <IFSC> in FREALIGN):\n'
                           'calculate FSC: Internally calculate two reconstructions with odd and even \n'
                           '   numbered particles and generate FSC table at the end of the run.\n'
                           '   Option <0> for parameter <IFSC> in FREALIGN.\n'
                           'The following options reduce memory usage:'
                           'Calculate one 3DR with odd particles: \n'
                           '   Only calculate one reconstruction using odd particles. Option <1>.\n'
                           'Calculate one 3DR with even particles: \n'
                           '   Only calculate one reconstruction using even particles. Option <2>.\n'
                           'Calculate one 3DR with all particles: \n'
                           '   Only calculate one reconstruction using all particles. Option <3>.')

        self.addParam('doAditionalStatisFSC', BooleanParam, default=True,
                      label="Calculate aditional statistics in FSC?",
                      help='Parameter <FSTAT> in FREALIGN\n'
                           'Calculate additional statistics in resolution table at the end \n'
                           '(QFACT, SSNR, CC and related columns). Setting FSTAT=F saves over 50% of memory!.')

        self.addParam('paddingFactor', EnumParam, choices=['1', '2', '4'], default=brandeis.PAD_4,
                      label='Padding Factor', display=EnumParam.DISPLAY_COMBO,
                      help='Parameter <IBLOW> in FREALIGN\n'
                           'Padding factor for reference structure.\n'
                           'Padding factor 4 requires the most memory but results\n'
                           'in the fastest search & refinement.\n')

        self.addSection(label='Projection Matching')
        
        self.addParam('innerRadius', FloatParam, default='0.0', 
                      label='Inner radius of reconstruction (in Amgs):', 
                      help='Parameter <RI> in FREALIGN\n'
                           'In Angstroms from centre of particle.\n'
                           'Enter the inner radius of the volume to be reconstructed.\n' 
                           'This is useful for reconstructions of viruses and other\n' 
                           'particles that might be hollow or have a disordered core.')
              
        self.addParam('outerRadius', FloatParam, default='108.0', 
                      label='Outer radius of reconstruction (in Amgs):', 
                      help='Parameter <RO> in FREALIGN\n'
                           'In Angstroms from centre of particle.\n'
                           'Enter the outer radius of the volume to be reconstructed.\n'
                           'he program will also apply a mask with a cosine edge to the particle image\n'
                           'before processing (done inside CTFAPPLY using  HALFW=6 pixels for cosine bell).')
        
        self.addParam('ThresholdMask', FloatParam, default='0.0', 
                      label='Threshold to for masking the input 3D structure:', expertLevel=LEVEL_ADVANCED,
                      help='Parameter <XSTD> in FREALIGN.\n'
                           'filtered 3D model - note this 3D masking does not use RI.\n'
                           '- if positive, calculates mask with subroutine D3MASK, equiv to\n'
                           '  solvent flattening with 5-pixel-cosine-bell smoothed mask\n'
                           '  boundary.  The mask is then used to multiply the input 3D map,\n'
                           '  which is then used for all parameter refinement and subsequent\n'
                           '  calculations.\n'
                           '- if negative, calculates mask with subroutine D2MASK resulting\n'
                           '  in a sharp binary (0/1) mask boundary for which is used for\n'
                           '  both parameter refinement and reconstruction, and to mask and\n'
                           '  output the matching projections.  Each matching particle image\n'
                           '  is also always masked with a cosine bell edged function of\n'
                           '  radius RI.\n'
                           'If set 0, disables this function.')
        
        self.addParam('pseudoBFactor', FloatParam, default='5.0', 
                      label='Resol-Dependent weighting of particles for 3D reconstruction:',
                      help='Parameter <PBC> in FREALIGN.\n'
                           'Automatic weighting is applied to each particle: a pseudo-temperature (B)\n'
                           'factor is applied to each particle according to its relative phase\n'
                           'residual against the reference. The weight is calculated as\n'
                           '          W = exp (-DELTAP/PBC * R^2)\n'
                           'with DELTAP = relative phase residual (actual phase residual minus BOFF),\n'
                           'PBC = conversion constant (5.0 in the example),\n'
                           'and R^2 the squared resolution in Fourier units (R = 0.0 ... 0.5).\n'
                           'A large value for PBC (e.g. 100.0) gives equal weighting to each particle.')

        self.addParam('avePhaseResidual', FloatParam, default='60.0', 
                      label='Average phase residual:',
                      help='Parameter <BOFF> in FREALIGN.\n'
                           'Approximate average phase residual of all particles,\n'
                           ' used in calculating weights for contributions of different\n'
                           'particles to 3D map (see Grigorieff, 1998).')

        self.addParam('angStepSize', FloatParam, default='10.0', 
                      label='Angular step size for the angular search:',
                      help='Parameter <DANG> in FREALIGN.\n'
                           'Angular step size for the angular search used in modes IFLAG=3,4.\n'
                           'There is a program default value calculated taking resolution into\n'
                           'account, but if this input value is non-zero, the program value is\n'
                           'overridden.')

        self.addParam('numberRandomSearch', FloatParam, default='10.0', 
                      label='Number of randomised search/refinement trials:',
                      help='Parameter <ITMAX> in FREALIGN.\n'
                           'number of cycles of randomised search/refinement used in modes IFLAG=2,4\n'
                           'There is a program default value (10 cycles), but if this input value is\n'
                           'non-zero, the program value is overridden.\n')

        self.addParam('numberPotentialMatches', FloatParam, default='10.0', 
                      label='number of potential matches:',
                      help='Parameter <IPMAX> in FREALIGN.\n'
                           'number of potential matches in a search that should be tested further in\n'
                           'a subsequent local refinement.\n')

        self.addParam('paramRefine', EnumParam, choices=['Refine All', 'Refine only Euler angles', 'Refine Only X & Y shifts'],
                      default=brandeis.REF_ALL,
                      label="Parameters to refine", display=EnumParam.DISPLAY_COMBO,
                      help='Parameter <MASK> in FREALIGN.\n'
                           'Parameters to include in refinemet')

        self.addParam('symmetry', TextParam, default='c1',
                      label='Point group symmetry:',
                      help='Parameter <ASYM> in FREALIGN.\n'
                           'Specify the symmetry.Choices are: C(n),D(n),T,O,I,I1,I2 or N (can be zero)\n'
                           'n  = rotational symmetry required in pointgroup C(n) or D(n)\n'
                           'N  = number of symmetry matrices to read in.\n'
                           'T  = tetrahedral pointgroup 23\n'
                           'O  = octahedral pointgroup 432\n'
                           'I  = icosahedral 532 symmetry in setting 1 (5fold is on X)\n'
                           'I1 = also in setting 1 (X) - as used by Imagic\n'
                           'I2 = in setting 2 (Y) - as used by Crowther et. al')        

        self.addParam('relMagnification', FloatParam, default='1.0', 
                      label='Relative magnification:',
                      help='Parameter <RELMAG> in FREALIGN.\n'
                           'Relative magnification of data set. The magnification feature allows\n'
                           'for manual variations of magnification between data sets.')

        self.addParam('targetPhaseResidual', FloatParam, default='10.0',
                      label='Target phase residual:', expertLevel=LEVEL_EXPERT,
                      help='Parameter <TARGET> in FREALIGN.\n'
                           'Target phase residual (for resolution between RMAX1 and RMAX2)\n'
                           'during parameter saerch and refinement, below which the search and/or\n'
                           'refinement is terminated.  This is normally set low (e.g. THRESH=10.0)\n'
                           'This will give excellent determination of particle orientations.')

        self.addParam('PhaseResidual', FloatParam, default='90.0', 
                      label='Phase residual cut-off:',
                      help='Parameter <THRESH> in FREALIGN.\n'
                           'Any particles with a higher overall phase residual will not be\n'
                           'included in the reconstruction when IFLAG=0,1,2,3. This variable\n'
                           'is often used with IFLAG=0 in separate runs to calculate maps\n'
                           'using various values of THRESH to find the optimum value to produce\n'
                           'the best map as judged from the statistics.')

        self.addParam('beamTiltX', FloatParam, default='0.0',
                      label='Beam tilt in X direction (in mrad):', expertLevel=LEVEL_EXPERT,
                      help='Parameter <TX> in FREALIGN.')

        self.addParam('beamTiltY', FloatParam, default='0.0',
                      label='Beam tilt in Y direction (in mrad):', expertLevel=LEVEL_EXPERT,
                      help='Parameter <TY> in FREALIGN.')

        self.addParam('resolution', FloatParam, default='10.0', 
                      label='Resol. of reconstruction (in Amgs):',
                      help='Parameter <RREC> in FREALIGN.\n'
                           'Resolution to which the reconstruction is calculated.\n'
                           'If several datasets have different values, the data is individually\n'
                           'limited in the summation to the RREC resolution but symmetry is\n'
                           'applied, statistics output and the final map calculated to the\n'
                           'maximum resolution requested for any dataset.')

        self.addParam('lowResolRefine', FloatParam, default='200.0', 
                      label='Low resolution in refinement (in Amgs):',
                      help='Parameter <RMAX1> in FREALIGN.\n'
                           'Resolution of the data included in the search/refinement. These\n'
                           'two parameters (RMAX1,RMAX2) are very important.  The successful\n'
                           'alignment of particles depends critically on the signal-to-noise\n'
                           'ratio of thecross-correlation or phase residual calculation, and\n'
                           'exclusion of weak data at high resolution or spurious, very strong\n'
                           'artefacts at low resolution can make a big difference.  Success can\n'
                           'be judged by whether the X,Y coordinates of the particle centres are\n'
                           'reasonable.')

        self.addParam('highResolRefine', FloatParam, default='25.0', 
                      label='High resolution in refinement (in Amgs):',
                      help='Parameter <RMAX2> in FREALIGN.\n'
                           'Resolution of the data included in the search/refinement. These\n'
                           'two parameters (RMAX1,RMAX2) are very important.  The successful\n'
                           'alignment of particles depends critically on the signal-to-noise\n'
                           'ratio of thecross-correlation or phase residual calculation, and\n'
                           'exclusion of weak data at high resolution or spurious, very strong\n'
                           'artefacts at low resolution can make a big difference.  Success can\n'
                           'be judged by whether the X,Y coordinates of the particle centres are\n'
                           'reasonable.')

        self.addParam('defocusUncertainty', FloatParam, default='200.0', 
                      label='Defocus uncertainty (in Amgs):', expertLevel=LEVEL_EXPERT,
                      help='Parameter <DFSIG> in FREALIGN.\n'
                           'This will restrain the change in defocus when refining defocus values\n'
                           'for individual particles.')

        self.addParam('Bfactor', FloatParam, default='0.0', 
                      label='B-factor to apply to particle:', expertLevel=LEVEL_EXPERT,
                      help='Parameter <RBFACT> in FREALIGN.\n'
                           'B-factor to apply to particle image projections before orientation\n'
                           'determination or refinement.  This allows inclusion of high resolution\n'
                           'data but with a reduced (or increased if negative) weight.  WGH and\n'
                           'RBFAC can be manipulated in particle parameter refinement as if they\n'
                           'were low pass and high pass filters.  WGH and the CTF are used to\n'
                           'correct the density in the final map, whereas RBFAC is not.\n'
                           'NOTE: This option should be used with great care as amplification of\n'
                           'noisy high-resolution terms can lead to over-fitting and artificially\n'
                           'high values in the FSC curve (se publication #2 above). FREALIGN uses an\n'
                           'automatic weighting scheme and RBFACT should normally be set to 0.0.')

        self.addParallelSection(threads=1, mpi=1)

class ProtFrealign(ProtRefine3D):
    """Protocol to perform a volume from a SetOfParticles
    using the frealign program"""
    _definition = BrandeisDefFrealign()
    _label = 'Frealign Protocol'
        
    def _defineSteps(self):
        # for entire set of particles
        maxIter = self.numberOfIterations.get()
        initVol = self.input3DReference.get()
        imgSet = self.inputParticles.get()
        coordSet = imgSet.getCoordinates()
        coordSet.loadIfEmpty()
        micSet = coordSet.getMicrographs()
        micSet.loadIfEmpty()
        
        while iter < maxIter:
            prevIter = iter
            iter = iter + 1
            self._createWorkingDir(iter)
            vol = 'volume_iter_%03d.mrc' % iter
            self._defCurrentIteration()
            # -------------------------------------
            # for each sub-set of particles
            if iter == 1:
                self._constructParamFile()
            else:
                self._params['stringConstructor'] = '../iter_%03d/particles_%06d_%06d_iter_%03d.par' % prevIter, iniParticle, finalParticle, prevIter
                self._refineParticles()
            #--------------------------------------
            # for entire set of particles
            if self._param['mode'] != 0:
                self._params['stringConstructor'] = '../iter_%03d/particles_%06d_%06d_iter_%03d.par' % prevIter, iniParticle, finalParticle, prevIter
                self._generate3DR()
            self._leaveWorkingDir()
        
    def _defCurrentIteration(self):
        """ Defining the current iteration
        """        
        #Prepare arguments to call program fralign_v8.exe
        self._params = {'imgsFn': imgsFn,
                        'mode': self.mode.get(),
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
                        'sampling3DR': self.input3DRef.getSamplingRate()
                       }

        # Get reference map for the iteration.
        self._params['initVol'] = self.input3DRef.get()
        
        # Get the amplitude Contrast of the micrographs
        # We need the amplitude Contrast for CTF correction
        self._params['ampContrast'] = imgSet.ampContrast.get()

        # Defining the operation mode
        if self.mode == MOD_RECONSTRUCTION:
            self._params['mode'] = 0
        elif self.mode == MOD_REFINEMENT:
            self._params['mode'] = 1
        elif self.mode == MOD_RANDOM_SEARCH_REFINEMENT:
            self._params['mode'] = 2
        elif self.mode == MOD_SIMPLE_SEARCH_REFINEMENT:
            self._params['mode'] = 3
        else:
            self._params['mode'] = 4
        
        # Defining the operation mode for iteration 1.
        if self.Firstmode == MOD2_SIMPLE_SEARCH_REFINEMENT:
            self._params['mode2'] = -3
        else:
            self._params['mode2'] = -4

        # Defining if magnification refinement is going to do
        if self.doMagRefinement:
            self._params['doMagRefinement'] = T
        else:
            self._params['doMagRefinement'] = F
            
        # Defining if defocus refinement is going to do
        if self.doDefRefinement:
            self._params['doDefocusRef'] = T
        else:
            self._params['doDefocusRef'] = F

        # Defining if astigmatism refinement is going to do
        if self.doAstigRefinement:
            self._params['doAstigRef'] = T
        else:
            self._params['doAstigRef'] = F
        
        # Defining if defocus refinement for individual particles is going to do
        if self.doDefPartRefinement:
            self._params['doDefocusPartRef'] = T
        else:
            self._params['doDefocusPartRef'] = F
        
        if self.methodEwaldSphere == EWA_DISABLE:
            self._params['metEwaldSphere'] = 0
        elif self.methodEwaldSphere == EWA_SIMPLE:
            self._params['metEwaldSphere'] = 1
        elif self.methodEwaldSphere == EWA_REFERENCE:
            self._params['metEwaldSphere'] = 2
        elif self.methodEwaldSphere == EWA_SIMPLE_HAND:
            self._params['metEwaldSphere'] = -1
        else:
            self._params['metEwaldSphere'] = -2
        
        # Defining if apply extra real space symmetry
        if self.doExtraRealSpaceSym:
            self._params['doExtraRealSpaceSym'] = T
        else:
            self._params['doExtraRealSpaceSym'] = F
        
        # Defining if wiener filter is going to apply
        if self.doWienerFilter:
            self._params['doWienerFilter'] = T
        else:
            self._params['doWienerFilter'] = F
        
        # Defining if wiener filter is going to calculate and apply
        if self.doBfactor:
            self._params['doBfactor'] = T
        else:
            self._params['doBfactor'] = F
        
        # Defining if matching projections is going to write
        if self.writeMatchProjections:
            self._params['writeMatchProj'] = T
        else:
            self._params['writeMatchProj'] = F
            
        if self.methodCalcFsc == FSC_CALC:
            self._params['metFsc'] = 0
        elif self.methodCalcFsc == FSC_3DR_ODD:
            self._params['metFsc'] = 1
        elif self.methodCalcFsc == FSC_3DR_EVEN:
            self._params['metFsc'] = 2
        elif self.methodCalcFsc == FSC_3DR_ALL:
            self._params['metFsc'] = 3
        
        
        if self.doAditionalStatisFSC:
            self._params['doAditionalStatisFSC'] = T
        else:
            self._params['doAditionalStatisFSC'] = F
            
        if self.paddingFactor == PAD_1:
            self._params['padFactor'] = 1
        elif self.paddingFactor == PAD_2:
            self._params['padFactor'] = 2
        else:
            self._params['padFactor'] = 4
            
        if self.paramRefine == REF_ALL:
            self._params['paramRefine'] = '1, 1, 1, 1, 1'
        elif self.paramRefine == REF_ANGLES:
            self._params['paramRefine'] = '1, 1, 1, 0, 0'
        else:
            self._params['paramRefine'] = '0, 0, 0, 1, 1'
        
        self._params['scannedPixelSize'] = micSet.scannedPixelSize.get()
        self._params['voltage'] = micSet.voltage.get()
        self._params['sphericalAberration'] = micSet.sphericalAberration.get()
        
    def _createWorkingDir(self, iter):        
        iterDir = 'iter_%03d' % iter
        workDir = self._getPath(iterDir)    
        # Create dir for current iteration
        makePath(workDir)
        os.chdir(workDir)
    
    def _generate3DR(self):        
        mode = self._params['mode']
        self._params['mode'] = 0
        # The stopParam must be 0 if you want to obtain a 3D reconstruction.
        self._params['stopParam'] = 0
        
        self._prepareCommand(iniParticle, finalParticle)
        self._params['mode'] = mode
    
    def _refineParticles(self):
        self._params['stopParam'] = -100
        self._prepareCommand(iniParticle, finalParticle)
        
    def _constructParamFile(self):
        """ This function construct a unique parameter file (.par) with the information of the SetOfParticle.
        This function will execute it only in iteration 1.
        """
        mode = self._params['mode']
        self._params['mode'] = selself._params['Firstmode']
        self._params['stopParam'] = 0
        part = iniParticle
        for part,img in enumerate(imgSet):
            #film = coordSet[img.getId()].getmicId()
            #mag = micSet.magnification.get()
           # defocus1 = img.
            if part == finalParticle:
                break
            
            
        self._params['stringConstructor'] = 1, magnification, film, def1, def2, astig, more
        self._prepareCommand(iniParticle, finalParticle)
        self._params['mode'] = mode

        

    def _prepareCommand(self, iniParticle, finalParticle):
#        TODO: change this for writeSetOfParticles and writeVolume functions.
#        TODO: Uncomment the following line when writeSetOfParticles function is implemented
#        imgsFn = writeSetOfParticles(self.inputParticles.get())

        self._params['initParticle'] = iniParticle
        self._params['finalParticle'] = finalParticle
        self._params['imgFnMatch'] = 'particles_match_iter_%03d.mrc' % iter
        self._params['outputParFn'] = 'particles_%06d_%06d_iter_%03d.par' % iniParticle, finalParticle, iter
        self._params['outputShiftFn'] = 'particles_shifts_iter_%03d.shft' % iter
        self._params['volume'] = vol
        self._params['3Dweigh'] = 'volume_weights_iter_%03d' % iter
        self._params['FSC3DR1'] = 'volume_1_iter_%03d' % iter
        self._params['FSC3DR2'] = 'volume_2_iter_%03d' % iter
        self._params['VolPhResidual'] = 'volume_phasediffs_iter_%03d' % iter
        self._params['VolpointSpread'] = 'volume_pointspread_iter_%03d' % iter

        if which('frealign_v8.exe') is '':
            raise Exception('Missing frealign_v8.exe')
        self._program = which('frealign_v8.exe')
        self._args = """   << eot > %(frealignOut)s
M,%(mode)i,%(doMagRefinement)s,%(doDefocusRef)s,%(doAstigRef)s,%(doDefocusPartRef)s,%(metEwaldSphere)d,%(doExtraRealSpaceSym)s,%(doWienerFilter)s,%(doBfactor)s,%(writeMatchProj)s,%(metFsc)d,%(doAditionalStatisFSC)s,%(padFactor)d
%(outerRadius)f,%(innerRadius)f,%(sampling3DR)f,%(ampContrast)f,%(ThresholdMask)f,%(pseudoBFactor)f,%(avePhaseResidual)f,%(angStepSize)f,%(numberRandomSearch)f,%(numberPotentialMatches)f
%(paramRefine)s
%(initParticle)d,%(finalParticle)d
%(sym)s
%(relMagnification)f,%(scannedPixelSize)f,%(targetPhaseResidual)f,%(PhaseResidual)f,%(sphericalAberration)f,%(voltage)f,%(beamTiltX)f,%(beamTiltY)f
%(resol)f,%(lowRes)f,%(highRes)f,%(defocusUncertainty)f,%(Bfactor)f
%(imgsFn)s
%(stringConstructor)s
%(inputParFn)s
%(outputParFn)s
%(outputShiftFn)s
%(stopParam)i, 0., 0., 0., 0., 0., 0., 0.
%(volume)s
%(3Dweigh)s
%(FSC3DR1)s
%(FSC3DR2)s
%(VolPhResidual)s
%(VolpointSpread)s
eot
"""
# 
#     def _getPsdPath(self, micDir):
#         return join(micDir, 'ctffind_psd.mrc')
#     
#     def _estimateCTF(self, micFn, micDir):
#         """ Run ctffind3 with required parameters """
#         # Create micrograph dir 
#         makePath(micDir)
#         # Update _params dictionary
#         self._params['micFn'] = micFn
#         self._params['micDir'] = micDir
#         self._params['frealignOut'] = join(micDir, 'frealign_protocol.log')
#         self._params['ctffindPSD'] = self._getPsdPath(micDir)
#                 
#         self.runJob(None, self._program, self._args % self._params)
#         # print "command: ", self._program, self._args % self._params    
# 
#     def _parseOutput(self, filename):
#         """ Try to find the output estimation parameters
#         from filename. It search for a line containing: Final Values.
#         """
#         f = open(filename)
#         result = None
#         for line in f:
#             if 'Final Values' in line:
#                 # Take DefocusU, DefocusV and Angle as a tuple
#                 # that are the first three values in the line
#                 result = tuple(map(float, line.split()[:3]))
#                 break
#         f.close()
#         return result
#             
#     def _getCTFModel(self, defocusU, defocusV, defocusAngle, psdFile):
#         ctf = CTFModel()
#         ctf.copyAttributes(self.inputMics, 'samplingRate')
#         ctf.copyAttributes(self.inputMics.microscope, 'voltage', 'sphericalAberration')
#         
#         ctf.defocusU.set(defocusU)
#         ctf.defocusV.set(defocusV)
#         ctf.defocusAngle.set(defocusAngle)
#         ctf.psdFile.set(psdFile)
#         
#         return ctf
#         
#     def createOutput(self):
#         path = self._getPath('micrographs.sqlite')
#         micSet = SetOfMicrographs(path)
#         micSet.copyInfo(self.inputMics)
#         
#         for fn, micDir, mic in self._iterMicrographs():
#             out = join(micDir, 'ctffind.out')
#             result = self._parseOutput(out)
#             defocusU, defocusV, defocusAngle = result
#             micOut = Micrograph(fn)
#             micOut.setCTF(self._getCTFModel(defocusU, defocusV, defocusAngle,
#                                                 self._getPsdPath(micDir)))
#             micSet.append(micOut)
#             
#         micSet.write()
#         # This property should only be set by CTF estimation protocols
#         micSet.setCTF(True)     
#         self._defineOutputs(outputMicrographs=micSet)
# 	
#     def _validate(self):
#         errors = []
#         if which('ctffind3.exe') is '':
#             errors.append('Missing ctffind3.exe')
# 	return errors
#             
