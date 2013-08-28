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
        
        self.addParam('mode', EnumParam, choices=['Recontruction only', 'Refinement', 'Random Search & Refine',
                                                   'Simple search & Refine', 'Search, Refine, Randomise', 'Create new param file, Simple search & Refine', 
                                                   'Create new param file, Search, Refine & Randomise'],
                      label="Operation mode", default=brandeis.MOD_REFINEMENT,
                      display=EnumParam.DISPLAY_COMBO,
                      help='Parameter <IFLAG> in FREALIGN'
                           'Mode  0: Reconstruction only parameters as read in\n'
                           'Mode  1: Refinement & Reconstruction\n'
                           'Mode  2: Random Search & Refinement\n'
                           'Mode  3: Simple search & Refine\n'
                           'Mode  4: Search,Refine,Randomise & extend to RREC\n'
                           'Mode -3: bootstrap parameter file, then Mode 3\n'
                           'Mode -4: bootstrap parameter file, then Mode 4\n'
                           '\n'
                           '\n'
                           '\n'
                           '\n'
                           '')

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






class ProtFrealign():
    """
    """

    def _prepareCommand(self):
        self._params['step_focus'] = 1000.0
        # Convert digital frequencies to spatial frequencies
        sampling = self.inputMics.samplingRate.get()
        self._params['lowRes'] = sampling / self._params['lowRes']
        self._params['highRes'] = sampling / self._params['highRes']        
        
        if which('ctffind3.exe') is '':
            raise Exception('Missing ctffind3.exe')
         
        self._program = 'export NATIVEMTZ=kk ; ' + which('ctffind3.exe')
        self._args = """   << eof > %(ctffindOut)s
%(micFn)s
%(ctffindPSD)s
%(sphericalAberration)f,%(voltage)f,%(ampContrast)f,%(magnification)f,%(scannedPixelSize)f
%(windowSize)d,%(lowRes)f,%(highRes)f,%(minDefocus)f,%(maxDefocus)f,%(step_focus)f
"""

    def _getPsdPath(self, micDir):
        return join(micDir, 'ctffind_psd.mrc')
    
    def _estimateCTF(self, micFn, micDir):
        """ Run ctffind3 with required parameters """
        # Create micrograph dir 
        makePath(micDir)
        # Update _params dictionary
        self._params['micFn'] = micFn
        self._params['micDir'] = micDir
        self._params['ctffindOut'] = join(micDir, 'ctffind.out')
        self._params['ctffindPSD'] = self._getPsdPath(micDir)
                
        self.runJob(None, self._program, self._args % self._params)
        # print "command: ", self._program, self._args % self._params    

    def _parseOutput(self, filename):
        """ Try to find the output estimation parameters
        from filename. It search for a line containing: Final Values.
        """
        f = open(filename)
        result = None
        for line in f:
            if 'Final Values' in line:
                # Take DefocusU, DefocusV and Angle as a tuple
                # that are the first three values in the line
                result = tuple(map(float, line.split()[:3]))
                break
        f.close()
        return result
            
    def _getCTFModel(self, defocusU, defocusV, defocusAngle, psdFile):
        ctf = CTFModel()
        ctf.copyAttributes(self.inputMics, 'samplingRate')
        ctf.copyAttributes(self.inputMics.microscope, 'voltage', 'sphericalAberration')
        
        ctf.defocusU.set(defocusU)
        ctf.defocusV.set(defocusV)
        ctf.defocusAngle.set(defocusAngle)
        ctf.psdFile.set(psdFile)
        
        return ctf
        
    def createOutput(self):
        path = self._getPath('micrographs.sqlite')
        micSet = SetOfMicrographs(path)
        micSet.copyInfo(self.inputMics)
        
        for fn, micDir, mic in self._iterMicrographs():
            out = join(micDir, 'ctffind.out')
            result = self._parseOutput(out)
            defocusU, defocusV, defocusAngle = result
            micOut = Micrograph(fn)
            micOut.setCTF(self._getCTFModel(defocusU, defocusV, defocusAngle,
                                                self._getPsdPath(micDir)))
            micSet.append(micOut)
            
        micSet.write()
        # This property should only be set by CTF estimation protocols
        micSet.setCTF(True)     
        self._defineOutputs(outputMicrographs=micSet)
	
    def _validate(self):
        errors = []
        if which('ctffind3.exe') is '':
            errors.append('Missing ctffind3.exe')
	return errors
            
