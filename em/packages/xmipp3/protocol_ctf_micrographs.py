# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Laura del Cano (ldelcano@cnb.csic.es)
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
This sub-package contains the XmippCtfMicrographs protocol
"""


from pyworkflow.em import *  
from pyworkflow.utils import *  
from xmipp import MetaData, MDL_MICROGRAPH, MDL_MICROGRAPH_ORIGINAL, MDL_MICROGRAPH_TILTED, MDL_MICROGRAPH_TILTED_ORIGINAL
from pyworkflow.em.packages.xmipp3.data import *
from os.path import join


class XmippDefCTFMicrograph(Form):
    """Create the definition of parameters for
    the XmippCtfMicrographs protocol"""
    def __init__(self):
        Form.__init__(self)
    
        self.addSection(label='Input')
        self.addParam('inputMicrographs', PointerParam, label="Micrographs", pointerClass='SetOfMicrographs')
        
        self.addSection(label='CTF Estimation')
        self.addParam('ctfDownFactor', FloatParam, default=2.,
                      label='CTF Downsampling factor',
                      help='Set to 1 for no downsampling. Non-integer downsample factors are possible. '
                           'This downsampling is only used for estimating the CTF and it does not affect '
                           'any further calculation. Ideally the estimation of the CTF is optimal when the '
                           'Thon rings are not too concentrated at the origin (too small to be seen) and not '
                           'occupying the whole power spectrum (since this downsampling might entail aliasing).')
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
        self.addParam('fastDefocusEst', BooleanParam, default=True,
                      label='Fast defocus estimate')
        self.addParam('minDefocus', FloatParam, default=0.5,
                      label='Minimum defocus to search (in microns)',
                      help=' Minimum defocus value (in microns) to include in defocus search. ' 
                      'Underfocus is represented by a positive number.')
        self.addParam('maxDefocus', FloatParam, default=10.,
                      label='Maximum defocus to search (in microns)',
                      help='Maximum defocus value (in microns) to include in defocus search. '
                           'Underfocus is represented by a positive number.')
        self.addParam('windowSize', IntParam, default=256,
                      label='Window size',
                      help='The PSD is estimated from small patches of this size. Bigger patches '
                           'allow identifying more details. However, since there are fewer windows, '
                           'estimations are noisier')

class XmippProtCTFMicrographs(ProtCTFMicrographs):
    """Protocol to perform CTF estimation on a set of micrographs in the project"""
    
    __prefix = join('%(micrographDir)s','xmipp_ctf')
    __templateDict = {
        # This templates are relative to a micrographDir
        'prefix': __prefix,
        'ctfparam': __prefix +  '.ctfparam',
        'psd': __prefix + '.psd',
        'enhanced_psd': __prefix + '_enhanced_psd.xmp',
        'ctfmodel_quadrant': __prefix + '_ctfmodel_quadrant.xmp',
        'ctfmodel_halfplane': __prefix + '_ctfmodel_halfplane.xmp',
        'ctffind_ctfparam': join('%(micrographDir)s', 'ctffind.ctfparam'),
        'ctffind_spectrum': join('%(micrographDir)s', 'ctffind_spectrum.mrc')
        }

    _definition = XmippDefCTFMicrograph()
        
    def __init__(self, **args):
        
        Protocol.__init__(self, **args)
        
    def _defineSteps(self):
        '''for each micrograph insert the steps to preprocess it
        '''
        # Get pointer to input micrographs 
        self.inputMics = self.inputMicrographs.get() 
        extraDir = os.path.join(self.workingDir,'extra')
        
        
        # For each micrograph insert the steps to preprocess it
        for mic in self.inputMics:
            fn = mic.getFileName()
            shortname = os.path.basename(fn)
            micrographDir = os.path.join(extraDir,shortname) 
            # Downsample if necessary
            if self.ctfDownFactor.get() != 1:
                finalname = self.tmpPath(os.path.splitext(fn)[0] + "_tmp.mrc")
                self._insertRunJobStep("xmipp_transform_downsample",
                                      "-i %s -o %s --step %f --method fourier" % (fn,finalname,self.ctfDownFactor.get())) 
            else :
                finalname = fn
            
            # CTF estimation with Xmipp
            args="--micrograph "+finalname+\
                 " --oroot " + self._getFilename('prefix', micrographDir=micrographDir)+\
                 " --kV "+str(self.inputMics.microscope.voltage.get())+\
                 " --Cs "+str(self.inputMics.microscope.sphericalAberration.get())+\
                 " --sampling_rate "+str(self.inputMics.samplingRate.get()*self.ctfDownFactor.get())+\
                 " --downSamplingPerformed "+str(self.ctfDownFactor.get())+\
                 " --ctfmodelSize 256"+\
                 " --Q0 "+str(self.ampContrast.get())+\
                 " --min_freq "+str(self.lowRes.get())+\
                 " --max_freq "+str(self.highRes.get())+\
                 " --pieceDim "+str(self.windowSize.get())+\
                 " --defocus_range "+str((self.maxDefocus.get()-self.minDefocus.get())*10000/2)+\
                 " --defocusU "+str((self.maxDefocus.get()+self.minDefocus.get())*10000/2)
            if (self.fastDefocusEst.get()):
                args+=" --fastDefocus"
            self._insertRunJobStep('xmipp_ctf_estimate_from_micrograph', args)
        
        # Insert step to create output objects       
        self._insertFunctionStep('createOutput')


    def createOutput(self, IOTable):

        #TODO: Define what kind of output produce
        
        # self.defineOutputs(micrograph=self.outputMicrographs)
