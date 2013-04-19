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
"""
Protocols related to EM
"""
import os
import shutil
from pyworkflow.object import String, Float
from pyworkflow.protocol import Protocol
from pyworkflow.protocol.params import *
from pyworkflow.em import SetOfMicrographs
from pyworkflow.utils.path import removeBaseExt


class DefImportMicrographs(Form):
    """Create the definition of parameters for
    the ImportMicrographs protocol"""
    def __init__(self):
        Form.__init__(self)
    
        self.addSection(label='Input')
        self.addParam('pattern', StringParam, label="Pattern")
        self.addParam('tiltPairs', BooleanParam, default=False, important=True,
                   label='Are micrographs tilt pairs?')
        
        self.addSection(label='Microscope description')
        self.addParam('voltage', FloatParam, default=200,
                   label='Microscope voltage (in kV)')
        self.addParam('sphericalAberration', FloatParam, default=2.26,
                   label='Spherical aberration (in mm)')
        self.addParam('samplingRateMode', EnumParam, default=0,
                   label='Sampling rate mode',
                   choices=['From image', 'From scanner'])
        self.addParam('samplingRate', FloatParam,
                   label='Sampling rate (A/px)', condition='samplingRateMode==0')
        self.addParam('magnification', IntParam, default=60000,
                   label='Magnification rate', condition='samplingRateMode==1')
        self.addParam('scannedPixelSize', FloatParam,
                   label='Scanned pixel size', condition='samplingRateMode==1')
        

class ProtImportMicrographs(Protocol):
    """Protocol to import a set of micrographs in the project"""
    _definition = DefImportMicrographs()
    
    def __init__(self, **args):
        Protocol.__init__(self, **args)         
        
    def _defineSteps(self):
        self._insertFunctionStep('importMicrographs', self.pattern.get(), self.tiltPairs.get(),
                                self.voltage.get(), self.sphericalAberration.get(),
                                self.samplingRate.get())
        
    def importMicrographs(self, pattern, tiltPairs, voltage, sphericalAberration, samplingRate):
        """Copy micrographs matching the filename pattern
        Register other parameters"""
        from glob import glob
        files = glob(pattern)
        if len(files) == 0:
            raise Exception('importMicrographs:There is not files matching pattern')
        path = self._getPath('micrographs.txt')
        micSet = SetOfMicrographs(filename=path, tiltPairs=tiltPairs)
        micSet.microscope.voltage.set(voltage)
        micSet.microscope.sphericalAberration.set(sphericalAberration)
        micSet.samplingRate.set(samplingRate)

        for f in files:
            dst = self._getPath(os.path.basename(f))            
            shutil.copyfile(f, dst)
            micSet.append(dst)
        
        micSet.writeToFile(path)
        self._defineOutputs(outputMicrographs=micSet)
        outFiles = micSet._files + [path]
        
        return outFiles


class DefCTFMicrographs(Form):
    """Create the definition of parameters for
    the XmippCtfMicrographs protocol"""
    def __init__(self):
        Form.__init__(self)
    
        self.addSection(label='Input')
        self.addParam('inputMicrographs', PointerParam, label="Micrographs", pointerClass='SetOfMicrographs')
        
        self.addSection(label='CTF Estimation')
        self.addParam('ampContrast', FloatParam, default=0.1,
                      label='Amplitude Contrast',
                      help='It should be a positive number, typically between 0.05 and 0.3.')
        self.addParam('lowRes', FloatParam, default=0.05,
                      label='Lowest resolution',
                      help='Give a value in digital frequency (i.e. between 0.0 and 0.5). '
                           'This cut-off prevents the typically peak at the center of the PSD '
                           'to interfere with CTF estimation. The default value is 0.05, but for '
                           'micrographs with a very fine sampling this may be lowered towards 0.0')
        self.addParam('highRes', FloatParam, default=0.35,
                      label='Highest resolution', 
                      help='Give a value in digital frequency (i.e. between 0.0 and 0.5). '
                           'This cut-off prevents high-resolution terms where only noise exists '
                           'to interfere with CTF estimation. The default value is 0.35, but it should '
                           'be increased for micrographs with signals extending beyond this value. '
                           'However, if your micrographs extend further than 0.35, you should consider '
                           'sampling them at a finer rate.')
        self.addParam('minDefocus', FloatParam, default=0.5,
                      label='Minimum defocus to search (in microns)',
                      help=' Minimum defocus value (in microns) to include in defocus search. ' 
                      'Underfocus is represented by a positive number.',
                      expertLevel=LEVEL_ADVANCED)
        self.addParam('maxDefocus', FloatParam, default=10.,
                      label='Maximum defocus to search (in microns)',
                      help='Maximum defocus value (in microns) to include in defocus search. '
                           'Underfocus is represented by a positive number.',
                      expertLevel=LEVEL_ADVANCED)
        self.addParam('windowSize', IntParam, default=256,
                      label='Window size',
                      help='The PSD is estimated from small patches of this size. Bigger patches '
                           'allow identifying more details. However, since there are fewer windows, '
                           'estimations are noisier',
                      expertLevel=LEVEL_ADVANCED)        


class ProtCTFMicrographs(Protocol):
    
    _definition = DefCTFMicrographs()
        
#    def __init__(self, **args):
#        
#        Protocol.__init__(self, **args)
            
    
    def _iterMicrographs(self):
        """Iterate over micrographs and yield
        micrograph name and a directory to process """
        for mic in self.inputMics:
            fn = mic.getFileName()
            micrographDir = self._getExtraPath(removeBaseExt(fn)) 
            yield (fn, micrographDir)  
        
    def _defineSteps(self):
        ''' insert the steps to perform ctf estimation on a set of micrographs
        '''
        # Get pointer to input micrographs 
        self.inputMics = self.inputMicrographs.get() 
        
        self.params = {'kV': self.inputMics.microscope.voltage.get(),
                       'Cs': self.inputMics.microscope.sphericalAberration.get(),
                       'sampling_rate': self.inputMics.samplingRate.get(),
                       'ctfmodelSize': self.windowSize.get(),
                       'Q0': self.ampContrast.get(),
                       'min_freq': self.lowRes.get(),
                       'max_freq': self.highRes.get(),
                       'pieceDim': self.windowSize.get(),
                       'defocus_range': (self.maxDefocus.get()-self.minDefocus.get())*10000/2,
                       'defocusU': (self.maxDefocus.get()+self.minDefocus.get())*10000/2
                       }
        
        # For each micrograph insert the steps to process it
        for fn, micrographDir in self._iterMicrographs():
            
            # CTF estimation with Xmipp
            self._insertFunctionStep('estimateCTF', fn, micrographDir)
                    
        # Insert step to create output objects       
        self._insertFunctionStep('createOutput')


class ProtPreprocessMicrographs(Protocol):
    pass


class ProtExtractParticles(Protocol):
    pass


class ProtParticlePicking(Protocol):
    pass


class ProtAlign(Protocol):
    pass


class ProtClassify(Protocol):
    pass


class ProtAlignClassify(Protocol):
    pass
