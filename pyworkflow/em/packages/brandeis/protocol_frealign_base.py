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
"""This module contains the protocol base class for frealign protocols"""

import os
from os.path import join, exists, basename

from pyworkflow.object import Integer
from pyworkflow.utils.path import copyFile, createLink, makePath
from pyworkflow.protocol.constants import STEPS_PARALLEL, LEVEL_ADVANCED, LEVEL_EXPERT
from pyworkflow.protocol.params import (StringParam, BooleanParam, IntParam, PointerParam,
                                        EnumParam, FloatParam, TextParam)
from pyworkflow.em.protocol import EMProtocol
from pyworkflow.em.convert import ImageHandler

from constants import (MOD2_SIMPLE_SEARCH_REFINEMENT, MOD_REFINEMENT, EWA_DISABLE, FSC_CALC, MEM_0,
                       INTERPOLATION_1, REF_ALL, MOD_RECONSTRUCTION, MOD_RANDOM_SEARCH_REFINEMENT,
                       MOD_SIMPLE_SEARCH_REFINEMENT, EWA_REFERENCE, EWA_SIMPLE_HAND, EWA_SIMPLE,
                       FSC_3DR_ODD, FSC_3DR_EVEN, FSC_3DR_ALL, MEM_1, MEM_2, INTERPOLATION_0, REF_ANGLES, REF_SHIFTS)
from brandeis import FREALIGN, FREALIGN_PATH, FREALIGNMP_PATH


class ProtFrealignBase(EMProtocol):
    """ This class cointains the common functionalities for all Frealign protocols.
    In subclasses there should be little changes about the steps to execute and several parameters
    """
    IS_REFINE = True

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)
        self.stepsExecutionMode = STEPS_PARALLEL
        self._lastIter = Integer(0)
    
    def _createFilenameTemplates(self):
        """ Centralize how files are called for iterations and references. """
        
        def iterFile(suffix=''):
            iterDn = 'iter_%(iter)03d'
            iterDir = self._getExtraPath(iterDn)
            return join(iterDir, suffix)
        
        myDict = {
                  'particles': self._getTmpPath('particles.mrc'),
                  'init_vol': self._getTmpPath('volume.mrc'),
                  # Volumes for the iteration
                  'ref_vol': iterFile('reference_volume_iter_%(iter)03d.mrc'),
                  'iter_vol': iterFile('volume_iter_%(iter)03d.mrc'),
                  'output_vol_par': 'output_vol_iter_%(iter)03d.par',
                  # dictionary for all set
                  'output_par': iterFile('particles_iter_%(iter)03d.par'),
                  'shift' : 'particles_shifts_iter_%(iter)03d.shft',
                  'match' : 'particles_match_iter_%(iter)03d.mrc',
                  'weight' : 'volume_weights_iter_%(iter)03d.mrc',
                  'vol1' : 'volume_1_iter_%(iter)03d.mrc',
                  'vol2' : 'volume_2_iter_%(iter)03d.mrc',
                  'phase' : 'volume_phasediffs_iter_%(iter)03d.mrc',
                  'spread' : 'volume_pointspread_iter_%(iter)03d.mrc',
                  # dictionary for each processing block
                  'input_par_block': iterFile('particles_iter_%(prevIter)03d_%(block)02d.par'),
                  'output_par_block': iterFile('particles_iter_%(iter)03d_%(block)02d.par'),
                  'shift_block' : 'particles_shifts_iter_%(iter)03d_%(block)02d.shft',
                  'match_block' : 'particles_match_iter_%(iter)03d_%(block)02d.mrc', 
                  'weight_block' : 'volume_weights_iter_%(iter)03d_%(block)02d',
                  'vol1_block' : 'volume_1_%(iter)03d_%(block)02d_iter',
                  'vol2_block' : 'volume_2_iter_%(iter)03d_%(block)02d',
                  'phase_block' : 'volume_phasediffs_iter_%(iter)03d_%(block)02d',
                  'spread_block' : 'volume_pointspread_iter_%(iter)03d_%(block)02d',
                  # each class volumes for the iteration
                  'ref_vol_class': iterFile('reference_volume_iter_%(iter)03d_class_%(ref)02d.mrc'),
                  'iter_vol_class': iterFile('volume_iter_%(iter)03d_class_%(ref)02d.mrc'),
                  'output_vol_par_class': 'output_vol_iter_%(iter)03d_class_%(ref)02d.par',
                  # dictionary for each class
                  'output_par_class': iterFile('particles_iter_%(iter)03d_class_%(ref)02d.par'),
                  'output_par_class_tmp': iterFile('particles_iter_%(iter)03d_class_0.par'),
                  'shift_class' : 'particles_shifts_iter_%(iter)03d_class_%(ref)02d.shft',
                  'match_class' : 'particles_match_iter_%(iter)03d_class_%(ref)02d.mrc',
                  'weight_class' : 'volume_weights_iter_%(iter)03d_class_%(ref)02d.mrc',
                  'vol1_class' : 'volume_1_iter_%(iter)03d_class_%(ref)02d.mrc',
                  'vol2_class' : 'volume_2_iter_%(iter)03d_class_%(ref)02d.mrc',
                  'phase_class' : 'volume_phasediffs_iter_%(iter)03d_class_%(ref)02d.mrc',
                  'spread_class' : 'volume_pointspread_iter_%(iter)03d_class_%(ref)02d.mrc',
                  'logFileRecons' : 'logRecons_iter_%(iter)03d_class_%(ref)02d.log',
                  # dictionary for each processing block and class
                  'input_par_block_class': iterFile('particles_iter_%(prevIter)03d_class_%(ref)02d_%(block)02d.par'),
                  'output_par_block_class': iterFile('particles_iter_%(iter)03d_class_%(ref)02d_%(block)02d.par'),
                  'shift_block_class' : 'particles_shifts_iter_%(iter)03d_class_%(ref)02d_%(block)02d.shft',
                  'match_block_class' : 'particles_match_iter_%(iter)03d_class_%(ref)02d_%(block)02d.mrc', 
                  'weight_block_class' : 'volume_weights_iter_%(iter)03d_class_%(ref)02d_%(block)02d',
                  'vol1_block_class' : 'volume_1_%(iter)03d_class_%(ref)02d_%(block)02d_iter',
                  'vol2_block_class' : 'volume_2_iter_%(iter)03d_class_%(ref)02d_%(block)02d',
                  'phase_block_class' : 'volume_phasediffs_iter_%(iter)03d_class_%(ref)02d_%(block)02d',
                  'spread_block_class' : 'volume_pointspread_iter_%(iter)03d_class_%(ref)02d_%(block)02d',
                  'logFileRefine' : 'logRefine_iter_%(iter)03d_class_%(ref)02d_%(block)02d.log'
                  }
        
        self._updateFilenamesDict(myDict)
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        
        form.addParam('doContinue', BooleanParam, default=False,
              label='Continue from a previous run?',
              help='If you set to *Yes*, you should select a previous'
              'run of type *%s* class and most of the input parameters'
              'will be taken from it.' % self.getClassName())
        form.addParam('continueRun', PointerParam, pointerClass=self.getClassName(),
                      condition='doContinue', allowsNull=True,
                      label='Select previous run',
                      help='Select a previous run to continue from.')
        form.addParam('continueIter', StringParam, default='last',
                      condition='doContinue', 
                      label='Continue from iteration',
                      help='Select from which iteration do you want to continue.'
                           'if you use *last*, then the last iteration will be used.'
                           'otherwise, a valid iteration number should be provided.')        
        form.addParam('inputParticles', PointerParam, label="Input particles", important=True,
                      pointerClass='SetOfParticles',pointerCondition='hasCTF',
                      condition='not doContinue',
                      help='Select the input particles.\n')  
        form.addParam('input3DReference', PointerParam,
                      pointerClass='Volume',
                      label='Initial 3D reference volume:',
                      condition='not doContinue', 
                      help='Input 3D reference reconstruction.\n')
        form.addParam('numberOfIterations', IntParam, default=10,
                      label='Number of iterations:',
                      help='Number of iterations to perform. If continue option is True,'
                           'you going to do this number of new iterations (e.g. if'
                           '*Continue from iteration* is set 3 and this param is set 10,'
                           ' the final iteration of the protocol will be the 13th.')
        
        if not self.IS_REFINE:
            form.addParam('numberOfClasses', IntParam, default=3, 
                          label='Number of classes:',
                          condition='not doContinue',
                          help='The number of classes (K) for a multi-reference refinement.'
                               'These classes will be made in a random manner from a single'
                               'reference by division of the data into random subsets during the'
                               'first iteration.')
            form.addParam('itRefineAngles', IntParam, default=5, 
                          label='Every how many iterations refine the angles?',
                          help='The number of classes (K) for a multi-reference refinement.'
                               'These classes will be made in a random manner from a single'
                               'reference by division of the data into random subsets during the'
                               'first iteration.')
            form.addParam('itRefineShifts', IntParam, default=10, 
                          label='Every how many iterations refine the shifts?',
                          help='The number of classes (K) for a multi-reference refinement.'
                               'These classes will be made in a random manner from a single'
                               'reference by division of the data into random subsets during the'
                               'first iteration.')
        form.addParam('useInitialAngles', BooleanParam, default=False,
                      label="Use initial angles/shifts ? ", 
                      help='Set to *Yes* if you want to use the projection assignment (angles/shifts) \n '
                      'associated with the input particles (hasProjectionAssigment=True)')
        
        form.addSection(label='Flow Control')
        
        form.addParam('Firstmode', EnumParam, condition='not useInitialAngles and not doContinue',
                      choices=['Simple search & Refine', 'Search, Refine, Randomise'],
                      label="Operation mode for iteration 1:", default=MOD2_SIMPLE_SEARCH_REFINEMENT,
                      display=EnumParam.DISPLAY_COMBO,
                      help='Parameter *IFLAG* in FREALIGN\n\n'
                           'This option is only for the iteration 1.\n'
                           '_Mode -3_: Simple search, Refine and create a parameter file for the set\n'
                           '_Mode -4_: Search, Refine, Randomize and create a parameterfile for the set')
        form.addParam('mode', EnumParam, choices=['Recontruction only', 'Refinement', 'Random Search & Refine',
                                                   'Simple search & Refine', 'Search, Refine, Randomise'],
                      label="Operation mode:", default=MOD_REFINEMENT,
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
                      label="Refine defocus for individual particles?", expertLevel=LEVEL_ADVANCED,
                      condition='doDefRefinement',
                      help='Parameter *FPART* in FREALIGN\n\n'
                           'Set _True or False_ to enable/disable defocus parameter refinement\n'
                           'for individual particles if Refine defocus is True, otherwise defocus\n'
                           'change is constrained to be the same for all particles in one image\n')
        form.addParam('methodEwaldSphere', EnumParam, choices=['Disable', 'Simple', 'Reference-based', 'Simple with reversed handedness',
                                                               'Reference-based with reversed handedness'],
                      default=EWA_DISABLE, expertLevel=LEVEL_EXPERT,
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
                      label="Apply Wiener filter?",
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
                      default=FSC_CALC,
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
        form.addParam('memory', EnumParam, choices=['NO pad - NO multi-vol', 'pad - NO multi-vol',
                                                         'NO pad- multi-vol', 'pad - multi-vol'],
                      default=MEM_0,
                      label='Memory usage', display=EnumParam.DISPLAY_COMBO,
                      help='Parameter *IMEM* in FREALIGN\n\n'
                            '_NO pad - NO multi-vol_: no padding of reference during refinement,\n'
                            '  no multi-volume parallelization during reconstruction (least memory usage).\n'
                            '  *Option 0* for parameter *IMEM* in FREALIGN\n'
                           '_pad - NO multi-vol_: padding of reference during refinement,\n'
                           '   no multi-volume parallelization during reconstruction. *Option 1*.\n'
                           '_NO pad - multi-vol_:no padding of reference during refinement,\n'
                           '   multi-volume parallelization during reconstruction. *Option 2*.\n'
                           '_pad - multi-vol_: padding of reference during refinement,\n'
                           '   multi-volume parallelization during reconstruction (most memory usage). *Option 3*.')
        form.addParam('interpolationScheme', EnumParam, choices=['Nearest neighbor', 'Trilinear'],
                      default=INTERPOLATION_1,
                      label='Interpolation Scheme', display=EnumParam.DISPLAY_COMBO,
                      help='Parameter *INTERP* in FREALIGN\n\n'
                            'The options are:\n'
                            '_Nearest neighbor_. *Option 0* for parameter *INTERP* in FREALIGN\n'
                           '_Trilinear_. More time-consuming. *Option 1* for parameter *INTERP* in FREALIGN\n')
        
        form.addSection(label='General Parameters')
        
        line = form.addLine('Reconstruction radius (A):',
                      help='Parameters *RI* and *RO* in FREALIGN\n\n'
                           'In Angstroms from centre of particle.\n'
                           'Enter the inner and outer radius of the volume to be reconstructed.\n' 
                           'This is useful for reconstructions of viruses and other\n' 
                           'particles that might be hollow or have a disordered core.\n'
                           'The program will also apply a mask with a cosine edge to \n'
                           'the particle image before processing \n'
                           '(done inside *CTFAPPLY* using  *HALFW=6* pixels for cosine bell).')
        line.addParam('innerRadius', FloatParam, default='0.0', 
                      label='Inner')
        line.addParam('outerRadius', FloatParam, default='108.0', 
                      label='Outer')
                       
        form.addParam('molMass', FloatParam, default='500.0', 
                      label='Molecular mass of the specimen (kDa):',
                      condition="doWienerFilter",
                      help='Parameter *MW* in FREALIGN\n\n'
                           'Approximate molecular mass of the particle, in kDa. This is used to\n'
                           'calculate the optimal filter for the 3D reconstruction.\n')
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
        form.addParam('pseudoBFactor', FloatParam, default='20.0', 
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
        form.addParam('avePhaseResidual', FloatParam, default='35.0', 
                      label='Average phase residual:',
                      help='Parameter *BOFF* in FREALIGN\n\n'
                           'Approximate average phase residual of all particles,\n'
                           ' used in calculating weights for contributions of different\n'
                           'particles to 3D map (see Grigorieff, 1998).')
        form.addParam('angStepSize', FloatParam, default='0.0',
                      condition="mode==3 or mode==4",
                      label='Angular step size for the angular search:',
                      help='Parameter *DANG* in FREALIGN\n\n'
                           'Angular step size for the angular search used in modes *IFLAG*=3,4.\n'
                           'There is a program default value calculated taking resolution into\n'
                           'account, but if this input value is non-zero, the program value is\n'
                           'overridden.')
        form.addParam('numberRandomSearch', FloatParam, default='10.0', 
                      label='Number of randomised search/refinement trials:',
                      condition="mode==2 or mode==4",
                      help='Parameter *ITMAX* in FREALIGN\n\n'
                           'number of cycles of randomised search/refinement used in modes IFLAG=2,4\n'
                           'There is a program default value (10 cycles), but if this input value is\n'
                           'non-zero, the program value is overridden.\n')
        form.addParam('numberPotentialMatches', FloatParam, default='20.0', 
                      label='number of potential matches:',
                      help='Parameter *IPMAX* in FREALIGN\n\n'
                           'number of potential matches in a search that should be tested further in\n'
                           'a subsequent local refinement.\n')
        form.addParam('paramRefine', EnumParam, choices=['Refine all', 'Refine only euler angles',
                                                         'Refine only X & Y shifts', 'None'],
                      default=REF_ALL,
                      label="Parameters to refine", display=EnumParam.DISPLAY_COMBO,
                      help='Parameter *MASK* in FREALIGN\n\n'
                           'Parameters to include in refinement')
        form.addParam('symmetry', TextParam, default='C1',
                      label='Point group symmetry:',
                      condition='not doContinue',
                      help='Parameter *ASYM* in FREALIGN\n\n'
                           'Specify the symmetry.Choices are: Cn,Dn,T,O,I,I1,I2,N or H (can be zero)\n'
                           'n  = rotational symmetry required in pointgroup C(n) or D(n)\n'
                           'N  = number of symmetry matrices to read in.\n'
                           'T  = tetrahedral pointgroup 23\n'
                           'O  = octahedral pointgroup 432\n'
                           'I  = icosahedral 532 symmetry in setting 1 (5fold is on X)\n'
                           'I1 = also in setting 1 (X) - as used by Imagic\n'
                           'I2 = in setting 2 (Y) - as used by Crowther et. al\n'
                           'H  = helical symmetry')

        form.addSection(label='3D Reconstruction')
        
        form.addParam('relMagnification', FloatParam, default='1.0', 
                      label='Relative magnification:',
                      help='Parameter *RELMAG* in FREALIGN\n\n'
                           'Relative magnification of data set. The magnification feature allows\n'
                           'for manual variations of magnification between data sets.')
        form.addParam('targetScore', FloatParam, default='90.0',
                      label='Target score:', expertLevel=LEVEL_EXPERT,
                      help='Parameter *TARGET* in FREALIGN\n\n'
                           'Target score (for resolution between RMAX1 and RMAX2)\n'
                           'during parameter search and refinement, above which the search and/or\n'
                           'refinement is terminated.  This is normally set high (e.g. THRESH=90.0)\n'
                           'This will give excellent determination of particle orientations.')
        form.addParam('score', FloatParam, default='10.0', 
                      label='Score cut-off:',
                      help='Parameter *THRESH* in FREALIGN\n\n'
                           'Any particles with a lower overall score will not be included\n'
                           'in the reconstruction when IFLAG=0,1,2,3. This variable\n'
                           'is often used with IFLAG=0 in separate runs to calculate maps\n'
                           'using various values of THRESH to find the optimum value to produce\n'
                           'the best map as judged from the statistics.')

        line = form.addLine('Beam tilt in direction: ',
                      help='Parameters *TX* and *TY* in FREALIGN in mrad')
        line.addParam('beamTiltX', FloatParam, default='0.0', label='X ')
        line.addParam('beamTiltY', FloatParam, default='0.0', label='Y ')
        
        form.addParam('resolution', FloatParam, default='10.0', 
                      label='Resolution of reconstruction (A):',
                      help='Parameter *RREC* in FREALIGN\n\n'
                           'Resolution to which the reconstruction is calculated.\n'
                           'If several datasets have different values, the data is individually\n'
                           'limited in the summation to the RREC resolution but symmetry is\n'
                           'applied, statistics output and the final map calculated to the\n'
                           'maximum resolution requested for any dataset.')
        
        line = form.addLine('Resolution in refinement (A)',
                      help='Parameters *RMIN* and *RMAX* in FREALIGN\n\n'
                           'Resolution of the data included in the search/refinement. These\n'
                           'two parameters (RMIN,RMAX) are very important.  The successful\n'
                           'alignment of particles depends critically on the signal-to-noise\n'
                           'ratio of thecross-correlation or phase residual calculation, and\n'
                           'exclusion of weak data at high resolution or spurious, very strong\n'
                           'artefacts at low resolution can make a big difference.  Success can\n'
                           'be judged by whether the X,Y coordinates of the particle centres are\n'
                           'reasonable.')
        line.addParam('lowResolRefine', FloatParam, default='200.0', label='Low')
        line.addParam('highResolRefine', FloatParam, default='25.0', label='High')

        form.addParam('resolClass', FloatParam, default='25.0', 
                      label='High resolution for classification (A):',
                      help='Parameter *RCLAS* in FREALIGN\n\n'
                            'High-resolution limit used for classification.\n'
                            'It should typically be set to the same resolution\n'
                            'limit used also for the refinement, or a bit lower.\n'
                            'Resolution of the data included in the search/refine')
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

        form.addParallelSection(threads=1, mpi=2)
    
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        """Insert the steps to refine orientations and shifts of the SetOfParticles
        """
        self.numberOfBlocks = max(self.numberOfMpi.get() - 1, self.numberOfThreads.get() - 1, 1)
        
        self._createFilenameTemplates()
        self._insertContinueStep()
        self._insertItersSteps()
        self._insertFunctionStep("createOutputStep")
    
    def _insertContinueStep(self):
        if self.doContinue:
            continueRun = self.continueRun.get()
            self.inputParticles.set(continueRun.inputParticles.get())
            self.symmetry.set(continueRun.symmetry.get())
            self.input3DReference.set(None)
            if self.continueIter.get() == 'last':
                self.initIter = continueRun._getCurrIter()
            else:
                self.initIter = int(self.continueIter.get()) + 1
            self._insertFunctionStep('continueStep', self.initIter)
        else:
            self.initIter = 1
            
        self.finalIter = self.initIter + self.numberOfIterations.get()
    
    def _insertItersSteps(self):
        """ Insert the steps for all iters """
        
        for iterN in self._allItersN():
            initId = self._insertFunctionStep('initIterStep', iterN)
            paramsDic = self._getParamsIteration(iterN)
            depsRefine = self._insertRefineIterStep(iterN, paramsDic, [initId])
            self._insertFunctionStep("reconstructVolumeStep", iterN, paramsDic, prerequisites=depsRefine)
        
    def _insertRefineIterStep(self, iterN, paramsDic, depsInitId):
        """ execute the refinement for the current iteration """
        
        depsRefine = []
        if iterN == 1:
            if not self.useInitialAngles.get():
                stepConstructId = self._insertFunctionStep("constructParamFilesStep", paramsDic, prerequisites=depsInitId)
                depsConstruct = [stepConstructId]
                for block in self._allBlocks():
                    refineId = self._insertFunctionStep("refineBlockStep", block, prerequisites=depsConstruct)
                    depsRefine.append(refineId)
            else:
                initAngStepId = self._insertFunctionStep("writeInitialAnglesStep", prerequisites=depsInitId)
                paramsDic['paramRefine'] = '0, 0, 0, 1, 1'
                for block in self._allBlocks():
                    refineId = self._insertFunctionStep("refineParticlesStep", iterN, block, paramsDic, prerequisites=[initAngStepId])
                    depsRefine.append(refineId)
        else:
            for block in self._allBlocks():
                refineId = self._insertFunctionStep("refineParticlesStep", iterN, block, paramsDic, prerequisites=depsInitId)
                depsRefine.append(refineId)
        return depsRefine
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def continueStep(self, iterN):
        """Create a symbolic link of a previous iteration from a previous run."""
        iterN = iterN - 1
        self._setLastIter(iterN)
        continueRun = self.continueRun.get()
        prevDir = continueRun._iterWorkingDir(iterN)
        currDir = self._iterWorkingDir(iterN)
        createLink(prevDir, currDir)
        
        imgSet = self.inputParticles.get()
        imgFn = self._getFileName('particles')
        imgSet.writeStack(imgFn)
        
    def initIterStep(self, iterN):
        """ Prepare files and directories for the current iteration """
        
        self._createIterWorkingDir(iterN) # create the working directory for the current iteration.
        prevIter = iterN - 1
        refVol = self._getFileName('ref_vol', iter=iterN) # reference volume of the step.
        iterVol =  self._getFileName('iter_vol', iter=iterN) # refined volume of the step
        prevIterVol = self._getFileName('iter_vol', iter=prevIter) # volume of the previous iteration
        
        if iterN == 1:
            imgSet = self.inputParticles.get()
            vol = self.input3DReference.get()
            
            imgFn = self._getFileName('particles')
            volFn = self._getFileName('init_vol')
            #TODO check if the input is already a single mrc stack
            imgSet.writeStack(imgFn) # convert the SetOfParticles into a mrc stack.
            ImageHandler().convert(vol.getLocation(), volFn) # convert the reference volume into a mrc volume
            copyFile(volFn, refVol)  #Copy the initial volume in the current directory.
        else:
            self._splitParFile(iterN, self.numberOfBlocks)
            copyFile(prevIterVol, refVol)   #Copy the reference volume as refined volume.
        copyFile(refVol, iterVol)   #Copy the reference volume as refined volume.
    
    def constructParamFilesStep(self, paramsDic):
        """ Construct a parameter file (.par) with the information of the SetOfParticles. """

        #  This function will be called only in iteration 1.

        iterN = 1
        iterDir = self._iterWorkingDir(iterN)
        self._enterDir(iterDir)
        
        imgSet = self.inputParticles.get()
        magnification = imgSet.getAcquisition().getMagnification()
        params = {}

        #frealign need to have a numeric micId not longer than 5 digits
        micIdList = imgSet.aggregate(['count'],'_micId',['_micId'])
        micIdMap={}
        counter = 0;
        for mic in micIdList:
            micIdMap[mic['_micId']]=counter
            counter = counter +1

        for block in self._allBlocks():
            more = 1
            initPart, lastPart = self._initFinalBlockPaticles(block)
            params['initParticle'] = initPart
            params['finalParticle'] = lastPart
            paramDic = self._setParamsRefineParticles(iterN, block)
            paramsRefine = dict(paramsDic.items() + params.items() + paramDic.items())
            f = self.__openParamFile(block, paramsRefine)
            
            # ToDo: Implement a better method to get the info particles.
            #  Now, you iterate several times over the SetOfParticles
            # (as many threads as you have)

            for i, img in enumerate(imgSet):
                film = micIdMap[img.getMicId()]
                ctf = img.getCTF()
                defocusU, defocusV, astig = ctf.getDefocusU(), ctf.getDefocusV(), ctf.getDefocusAngle()
                partCounter = i + 1
                
                if partCounter == lastPart: # The last particle in the block
                    more = 0
                particleLine = ('1, %05d, %05d, %05f, %05f, %02f, %01d\n' %
                                (magnification, film, defocusU, defocusV, astig, more))
                self.__writeParamParticle(f, particleLine)
                
                if more == 0: # close the block.
                    self.__closeParamFile(f, paramsRefine)
                    break
        self._leaveDir()
    
    def refineBlockStep(self, block):
        """ This function execute the bash script for refine a subset(block) of images.
        It will enter in the iteration dir and execute the script there. 
        """
        iterDir = self._iterWorkingDir(1)
        program = "./block%03d.sh" % block
        os.chmod(join(iterDir, program), 0775)
        self.runJob(program, "", cwd=iterDir)
    
    def writeInitialAnglesStep(self):
        """This function write a .par file with all necessary information for a refinement"""
        
        imgSet = self.inputParticles.get()
        #frealign need to have a numeric micId not longer than 5 digits
        micIdList = imgSet.aggregate(['count'],'_micId',['_micId'])
        self.micIdMap={}
        counter = 0;
        for mic in micIdList:
            self.micIdMap[mic['_micId']]=counter
            counter = counter +1

        for block in self._allBlocks():
            more = 1
            _, lastPart = self._initFinalBlockPaticles(block)
            parFn = self._getFileName('input_par_block', block= block, iter=1, prevIter=0)
            f = open(parFn, 'w')
            f.write("C           PSI   THETA     PHI       SHX       SHY     MAG  FILM      DF1"
                    "      DF2  ANGAST     OCC     -LogP      SIGMA   SCORE  CHANGE\n")
            
            # ToDo: Implement a better method to get the info particles. Now, you iterate several times over the SetOfParticles (as many threads as you have)
            for i, img in enumerate(imgSet):
                partCounter = i + 1

                if partCounter == lastPart: # The last particle in the block
                    more = 0
                
                self.writeAnglesLines(partCounter, img, f)
            
                if more == 0: # close the block.
                    f.close()
                    break
    
    def refineParticlesStep(self, iterN, block, paramsDic):
        """Only refine the parameters of the SetOfParticles
        """
        param = {}
        
        iterDir = self._iterWorkingDir(iterN)
        
        iniPart, lastPart = self._initFinalBlockPaticles(block)
        prevIter = iterN - 1
        param['inputParFn'] = self._getBaseName('input_par_block', block= block, iter=iterN, prevIter=prevIter)
        param['initParticle'] = iniPart 
        param['finalParticle'] = lastPart
        
        paramDic = self._setParamsRefineParticles(iterN, block)
        
        paramsRefine = dict(paramsDic.items() + paramDic.items() + param.items())
        args = self._prepareCommand()
        
        if self.mode.get() != 0:
            # frealign program is already in the args script, that's why runJob('')
            self.runJob('', args % paramsRefine, cwd=iterDir)
        else:
            pass
            ##ugly hack when for reconstruction only, just copy the input files
            #inFile  = self._getFileName('input_par_block', block= block, iter=iterN, prevIter=prevIter)
            #outFile = self._getFileName('output_par_block', block=block, iter=iterN)
            #print "I am in dir: ", os.getcwd()
            #print "copying params files", inFile, outFile
            #copyFile(inFile, outFile)
    
    def reconstructVolumeStep(self, iterN, paramsDic):
        """Reconstruct a volume from a SetOfParticles with its current parameters refined
        """
        
        imgSet = self.inputParticles.get()
        self._mergeAllParFiles(iterN, self.numberOfBlocks)  # merge all parameter files generated in a refineIterStep function.
        
        initParticle = 1
        finalParticle = imgSet.getSize()
        
        os.environ['NCPUS'] = str(self.numberOfBlocks)
        paramsDic['frealign'] = FREALIGNMP_PATH
        paramsDic['outputParFn'] = self._getFileName('output_vol_par', iter=iterN)
        paramsDic['initParticle'] = initParticle
        paramsDic['finalParticle'] = finalParticle
#         paramsDic['paramRefine'] = '0, 0, 0, 0, 0'
        
        params2 = self._setParams3DR(iterN)
        
        params3DR = dict(paramsDic.items() + params2.items())

        args = self._prepareCommand()
        iterDir = self._iterWorkingDir(iterN)
        # frealign program is already in the args script, that's why runJob('')
        self.runJob('', args % params3DR, cwd=iterDir)
        self._setLastIter(iterN)
    
    def createOutputStep(self):
        pass # should be implemented in subclasses
    
    #--------------------------- INFO functions ----------------------------------------------------
    def _validate(self):
        errors = []
        if self.doContinue:
            continueRun = self.continueRun.get()
            imgSet = continueRun.inputParticles.get()
        else:
            imgSet = self.inputParticles.get()
        
        if not exists(FREALIGN_PATH):
            errors.append('Missing ' + FREALIGN_PATH)
        
        partSizeX, _, _ = imgSet.getDim()
        if not self.doContinue:
            volSizeX, _, _ = self.input3DReference.get().getDim()
            if partSizeX != volSizeX:
                errors.append('Volume and particles dimensions must be equal!!!')
        
        halfX = partSizeX % 2
        if halfX != 0:
            errors.append('Particle dimensions must be even!!!')
        if not imgSet.hasAlignmentProj() and self.useInitialAngles.get():
            errors.append("Particles has not initial angles !!!")
        return errors
    
    def _summary(self):
        summary = []
        if self.inputParticles.hasValue():
            summary.append("Number of particles:  %d" % self.inputParticles.get().getSize())
        if self.input3DReference.hasValue():
            summary.append("Input volume:  %s" % self.input3DReference.get().getFileName())
        
        if not hasattr(self, 'outputVolume'):
            summary.append("Output volumes not ready yet.")
        else:
            summary.append("Number of iterations: %d" % self.numberOfIterations.get())
#             summary.append("Angular step size: %f" % self.angStepSize.get())
            summary.append("symmetry: %s" % self.symmetry.get())
            summary.append("Final volume: %s" % self.outputVolume.getFileName())
        
        return summary
    
    def _methods(self):
        # ToDo: implement this method
        return self._summary()
    
    #--------------------------- UTILS functions ---------------------------------------------------
    def _getParamsIteration(self, iterN):
        """ Defining the current iteration
        """
        imgSet = self.inputParticles.get()
        
        #Prepare arguments to call program fralign_v9.exe
        paramsDic = {'frealign': FREALIGN_PATH,
                     'mode': self.mode.get(),
                     'innerRadius': self.innerRadius.get(),
                     'outerRadius': self.outerRadius.get(),
                     'molMass': self.molMass.get(),
                     'ThresholdMask': self.ThresholdMask.get(),
                     'pseudoBFactor': self.pseudoBFactor.get(),
                     'avePhaseResidual': self.avePhaseResidual.get(),
                     'angStepSize': self.angStepSize.get(),
                     'numberRandomSearch': self.numberRandomSearch.get(),
                     'numberPotentialMatches': self.numberPotentialMatches.get(),
                     'sym': self.symmetry.get(),
                     'relMagnification': self.relMagnification.get(),
                     'targetScore': self.targetScore.get(),
                     'score': self.score.get(),
                     'beamTiltX': self.beamTiltX.get(),
                     'beamTiltY': self.beamTiltY.get(),
                     'resol': self.resolution.get(),
                     'lowRes': self.lowResolRefine.get(),
                     'highRes': self.highResolRefine.get(),
                     'resolClass': self.resolClass.get(),
                     'defocusUncertainty': self.defocusUncertainty.get(),
                     'Bfactor': self.Bfactor.get(),
                     'sampling3DR': imgSet.getSamplingRate()
                    }
        
        # Get the particles stack
        iterDir = self._iterWorkingDir(iterN)
        paramsDic['imageFn'] = os.path.relpath(self._getFileName("particles"), iterDir)
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
        if self.doMagRefinement and iterN != 1:
            paramsDic['doMagRefinement'] = 'T'
        else:
            paramsDic['doMagRefinement'] = 'F'
            
        # Defining if defocus refinement is going to do
        if self.doDefRefinement and iterN != 1:
            paramsDic['doDefocusRef'] = 'T'
        else:
            paramsDic['doDefocusRef'] = 'F'

        # Defining if astigmatism refinement is going to do
        if self.doAstigRefinement and iterN != 1:
            paramsDic['doAstigRef'] = 'T'
        else:
            paramsDic['doAstigRef'] = 'F'
        
        # Defining if defocus refinement for individual particles is going to do
        if self.doDefPartRefinement and iterN != 1:
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
            
        if self.memory == MEM_0:
            paramsDic['memory'] = 0
        elif self.memory == MEM_1:
            paramsDic['memory'] = 1
        elif self.memory == MEM_2:
            paramsDic['memory'] = 2
        else:
            paramsDic['memory'] = 3
        
        if self.interpolationScheme == INTERPOLATION_0:
            paramsDic['interpolation'] = 0
        else:
            paramsDic['interpolation'] = 1
            
        if self.paramRefine == REF_ALL:
            paramsDic['paramRefine'] = '1, 1, 1, 1, 1'
        elif self.paramRefine == REF_ANGLES:
            paramsDic['paramRefine'] = '1, 1, 1, 0, 0'
        elif self.paramRefine == REF_SHIFTS:
            paramsDic['paramRefine'] = '0, 0, 0, 1, 1'
        else:
            paramsDic['paramRefine'] = '0, 0, 0, 0, 0'
        
        return paramsDic
    
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
        particles_iter_001.par
    def _createIterWorkingDir(self, iterN):
        """create a new directory for the iterarion and change to this directory.
        """
        workDir = self._iterWorkingDir(iterN)
        makePath(workDir)   # Create a directory for a current iteration
    
    def _iterWorkingDir(self, iterN, *paths):
        """ Define which is the directory for the current iteration"""
        iterDir = 'iter_%03d' % iterN
        workDir = self._getExtraPath(iterDir, *paths)
        return workDir

    def _getBaseName(self, key, **args):
        """ Remove the folders and return the file from the filename. """
        return basename(self._getFileName(key, **args))
    
    def _setParamsRefineParticles(self, iterN, block):
        paramDics = {}
        paramDics['stopParam'] = -100
        paramDics['volume'] = self._getBaseName('ref_vol', iter=iterN)
        paramDics['outputParFn'] = self._getBaseName('output_par_block', block=block, iter=iterN)
        paramDics['inputParFn'] = paramDics['outputParFn']
        paramDics['imgFnMatch'] = self._getFileName('match_block', block=block, iter=iterN)
        paramDics['outputShiftFn'] = self._getFileName('shift_block', block=block, iter=iterN)
        paramDics['3Dweigh'] = self._getFileName('weight_block', block=block, iter=iterN)
        paramDics['FSC3DR1'] = self._getFileName('vol1_block', block=block, iter=iterN)
        paramDics['FSC3DR2'] = self._getFileName('vol2_block', block=block, iter=iterN)
        paramDics['VolPhResidual'] = self._getFileName('phase_block', block=block, iter=iterN)
        paramDics['VolpointSpread'] = self._getFileName('spread_block', block=block, iter=iterN)
        paramDics['logFile'] = self._getFileName('logFileRefine', block=block, iter=iterN, ref=1)
        return paramDics
    
    def _setParams3DR(self, iterN):
        """ Setting the parameters to reconstruct a new 3DR"""
        paramDics = {}
        paramDics['mode'] = 0
        paramDics['stopParam'] = 0   #The stopParam must be 0 if you want obtain a 3D reconstruction.
        paramDics['volume'] = self._getBaseName('iter_vol', iter=iterN)
        paramDics['inputParFn'] = self._getBaseName('output_par', iter=iterN)
        paramDics['imgFnMatch'] = self._getFileName('match', iter=iterN)
        paramDics['outputShiftFn'] = self._getFileName('shift', iter=iterN)
        paramDics['3Dweigh'] = self._getFileName('weight', iter=iterN)
        paramDics['FSC3DR1'] = self._getFileName('vol1', iter=iterN)
        paramDics['FSC3DR2'] = self._getFileName('vol2', iter=iterN)
        paramDics['VolPhResidual'] = self._getFileName('phase', iter=iterN)
        paramDics['VolpointSpread'] = self._getFileName('spread', iter=iterN)
        paramDics['logFile'] = self._getFileName('logFileRecons', iter=iterN, ref=1)
        return paramDics
        
    def __openParamFile(self, blockNumber, paramsDict):
        """ Open the file and write the first part of the block param file. """
        if not exists(FREALIGN_PATH):
            raise Exception('Missing ' + FREALIGN)
        initaLines = """%(frealign)s << eot > %(logFile)s
M,%(mode2)s,%(doMagRefinement)s,%(doDefocusRef)s,%(doAstigRef)s,%(doDefocusPartRef)s,%(metEwaldSphere)s,%(doExtraRealSpaceSym)s,%(doWienerFilter)s,%(doBfactor)s,%(writeMatchProj)s,%(metFsc)s,%(doAditionalStatisFSC)s,%(memory)s,%(interpolation)s
%(outerRadius)s,%(innerRadius)s,%(sampling3DR)s,%(molMass)s,%(ampContrast)s,%(ThresholdMask)s,%(pseudoBFactor)s,%(avePhaseResidual)s,%(angStepSize)s,%(numberRandomSearch)s,%(numberPotentialMatches)s
%(paramRefine)s
%(initParticle)s,%(finalParticle)s
%(sym)s
%(relMagnification)s,%(scannedPixelSize)s,%(targetScore)s,%(score)s,%(sphericalAberration)s,%(voltage)s,%(beamTiltX)s,%(beamTiltY)s
%(resol)s,%(lowRes)s,%(highRes)s,%(resolClass)s,%(defocusUncertainty)s,%(Bfactor)s
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
        
    def _prepareCommand(self):
        """ prepare the command to execute"""

        if not exists(FREALIGN_PATH):
            raise Exception('Missing ' + FREALIGN)
        args = """%(frealign)s  << eot > %(logFile)s
M,%(mode)s,%(doMagRefinement)s,%(doDefocusRef)s,%(doAstigRef)s,%(doDefocusPartRef)s,%(metEwaldSphere)s,%(doExtraRealSpaceSym)s,%(doWienerFilter)s,%(doBfactor)s,%(writeMatchProj)s,%(metFsc)s,%(doAditionalStatisFSC)s,%(memory)s,%(interpolation)s
%(outerRadius)s,%(innerRadius)s,%(sampling3DR)s,%(molMass)s,%(ampContrast)s,%(ThresholdMask)s,%(pseudoBFactor)s,%(avePhaseResidual)s,%(angStepSize)s,%(numberRandomSearch)s,%(numberPotentialMatches)s
%(paramRefine)s
%(initParticle)s,%(finalParticle)s
%(sym)s
%(relMagnification)s,%(scannedPixelSize)s,%(targetScore)s,%(score)s,%(sphericalAberration)s,%(voltage)s,%(beamTiltX)s,%(beamTiltY)s
%(resol)s,%(lowRes)s,%(highRes)s,%(resolClass)s,%(defocusUncertainty)s,%(Bfactor)s
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
    
    def _mergeAllParFiles(self, iterN, numberOfBlocks):
        """ This method merge all parameters files that has been created in a refineIterStep """

        #if we only want to reconstruct then use the initial par file
        #instead of the output one since they are empty
        file2 = self._getFileName('output_par', iter=iterN)
        if (self.mode.get()==0):
            inFile = self._getFileName('input_par_block', block= numberOfBlocks, iter=1, prevIter=0)
            print inFile, file2
            copyFile(inFile,file2)
        else:
            if numberOfBlocks != 1:
                f2 = open(file2, 'w+')
                f2.write("C           PSI   THETA     PHI       SHX       SHY     MAG  FILM      DF1"
                         "      DF2  ANGAST     OCC     -LogP      SIGMA   SCORE  CHANGE\n")
                for block in range(1, numberOfBlocks + 1):
                    file1 = self._getFileName('output_par_block', block=block, iter=iterN)
                    if not os.path.exists(file1):
                         raise Exception ("Error: file %s does not exists" % file1)
                    f1 = open(file1)

    #                 if block == 1:
    #                     lines = f1.readlines()
    #                     f2.writelines(lines[:-2])
    #                 else:
                    for l in f1:
                        if not l.startswith('C'):
                            f2.write(l)
                    f1.close()
                f2.close()
            else:
                file1 = self._getFileName('output_par_block', block=1, iter=iterN)
                copyFile(file1, file2)
    
    def _splitParFile(self, iterN, numberOfBlocks):
        """ This method split the parameter files that has been previosuly merged """
        
        prevIter = iterN -1
        file1 = self._getFileName('output_par', iter=prevIter)
        if numberOfBlocks != 1:
            for block in range(1, numberOfBlocks + 1):
                f1 = open(file1)
                file2 = self._getFileName('input_par_block', block=block, iter=iterN, prevIter=prevIter)
                f2 = open(file2, 'w+')
                f2.write("C           PSI   THETA     PHI       SHX       SHY     MAG  FILM      DF1"
                         "      DF2  ANGAST     OCC     -LogP      SIGMA   SCORE  CHANGE\n")

                _, finalPart = self._initFinalBlockPaticles(block)
                
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
        else:
            file2 = self._getFileName('input_par_block', block=1, iter=iter, prevIter=prevIter)
            copyFile(file1, file2)
    
    def _setLastIter(self, iterN):
        self._lastIter.set(iterN)
    
    def _getLastIter(self):
        return self._lastIter.get()
    
    def _getCurrIter(self):
        return self._getLastIter() + 1
    
    def _allItersN(self):
        """ Iterate over iterations steps """
        for iterN in range(self.initIter, self.finalIter):
            yield iterN
    
    def _allBlocks(self):
        """ Iterate over all numberOfCPUs. """
        for i in range(1, self.numberOfBlocks+1):
            yield i
    
    def _initFinalBlockPaticles(self, block):
        """ return initial and final particle number for a determined block """
        finalPart = 0
        for i in range(block):
            particlesPerBlock = self._particlesPerBlock(self.numberOfBlocks, self.inputParticles.get().getSize())
            initPart = 1 + finalPart
            finalPart = finalPart + particlesPerBlock[i]
        return initPart, finalPart
    
    def writeAnglesLines(self, counter, img, filePar):
        
        objId = self.micIdMap[img.getMicId()]
        
        # get alignment parameters for each particle
        from convert import geometryFromMatrix
        shifts, angles = geometryFromMatrix(img.getTransform().getMatrix())
        #TODO: check if can use shiftZ
        shiftX, shiftY, _ = shifts * img.getSamplingRate()
        psi, theta, phi = angles
#        shiftX = float(str(align._xmipp_shiftX)) * img.getSamplingRate()
#        shiftY = float(str(align._xmipp_shiftY)) * img.getSamplingRate()
#        psi   = float(str(align._xmipp_anglePsi))
#        theta = float(str(align._xmipp_angleTilt))
#        phi   = float(str(align._xmipp_angleRot))
                    
        # get ctfModel for each particle
        ctfModel = img.getCTF()
        defU     = ctfModel.getDefocusU()
        defV     = ctfModel.getDefocusV()
        defAngle = ctfModel.getDefocusAngle()
        
        # get the adquisition info
        acquisition = img.getAcquisition()
        mag = acquisition.getMagnification()
        
        filePar.write("%(counter)7d %(psi)7.2f %(theta)7.2f %(phi)7.2f %(shiftX)9.2f %(shiftY)9.2f"
                   " %(mag)7.0f %(objId)5d %(defU)8.1f %(defV)8.1f %(defAngle)7.2f  100.00      0000     0.5000   00.00   00.00\n" % locals())
