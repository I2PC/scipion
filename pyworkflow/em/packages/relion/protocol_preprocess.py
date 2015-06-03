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

from glob import glob

from pyworkflow.utils.path import moveFile
from pyworkflow.em.protocol.protocol_particles import ProtProcessParticles
from pyworkflow.protocol.params import PointerParam, BooleanParam, FloatParam, IntParam, Positive

from pyworkflow.em.packages.relion.convert import writeSetOfParticles
from pyworkflow.em.packages.relion.protocol_base import ProtRelionBase
import pyworkflow.em as em

class ProtRelionPreprocessParticles(ProtProcessParticles, ProtRelionBase):
    """ Wrapper to Relion preprocess program.
    This protocol provides an easy way to execute *relion_preprocess* program
    to perform operations such us: normalize, filtering or scaling on the particles.
    """
    _label = 'preprocess particles'
    
    def __init__(self, **args):
        ProtProcessParticles.__init__(self, **args)
                   
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', PointerParam, pointerClass='SetOfParticles',
                      label="Input particles", important=True, 
                      help='Select the input images from the project.')   

        form.addParam('doNormalize', BooleanParam, default=True,
                      label='Normalize', important=True,
                      help='If set to True, particles will be normalized in the way RELION prefers it.\n'
                           'It is recommended to *always normalize your particles*, and use a reasonable\n'
                           'radius for the circle around your particles outside of which the standard\n'
                           'deviation and average values for the noise are calculated.\n\n'
                           '*Note*: if the particles are re-scaled, the radius for normalize will be\n'
                           'taken over the new dimensions.')
        form.addParam('backRadius', IntParam, default=-1,
                      condition='doNormalize',
                      label='Background radius (px)',
                      help='Pixels outside this circle are assumed to be noise and their stddev '
                      'is set to 1. Radius for background circle definition (in pixel).')        

        form.addParam('doRemoveDust', BooleanParam, default=False,
                      label='Remove dust from particles', 
                      help='If there are white or black artefacts on the micrographs \n'
                           '(e.g. caused by dust or hot/dead pixels), these may be removed \n'
                           'by using a positive value for the dust removal options.'
                           'All black/white pixels with values above the given parameter times\n'
                           'the standard deviation of the noise are replaced by random values from'
                           'a Gaussian distribution. For cryo-EM data, values around 3.5-5 are '
                           'often useful. Make sure you do not erase part of the true signal. ')
        line = form.addLine('Dust sigmas', condition='doRemoveDust',
                            help='Sigma-values above which white/black dust will be removed\n'
                                 '(negative value means no dust removal)')
        line.addParam('whiteDust', FloatParam, default=-1., label='White')
        line.addParam('blackDust', FloatParam, default=-1., label='Black')

        form.addParam('doInvert', BooleanParam, default=False,
                      label='Invert contrast', 
                      help='Invert the contrast if your particles are black over a white background.')

        form.addSection('Scale and window')
        form.addParam('doScale', BooleanParam, default=False,
                      label='Scale particles?',
                      help='Re-scale the particles to this size (in pixels).')
        form.addParam('scaleSize', IntParam, default=0, validators=[Positive],
                      condition='doScale',
                      label='Scale size (px)',
                      help='New particle size in pixels.')  
        form.addParam('doWindow', BooleanParam, default=False,
                      label='Window particles?', 
                      help='Re-window the particles to this size (in pixels).')
        form.addParam('windowSize', IntParam, default=0, validators=[Positive],
                      condition='doWindow',
                      label='Window size (px)',
                      help='New particles windows size (in pixels).')

    #--------------------------- INSERT steps functions --------------------------------------------
    
    def _insertAllSteps(self):
        self._insertFunctionStep("convertInputStep", self.inputParticles.get().getObjId())
        self._insertFunctionStep('processStep')
        self._insertFunctionStep('createOutputStep')
    
    #--------------------------- STEPS functions --------------------------------------------
    def convertInputStep(self, particlesId):
        """ Create the input file in STAR format as expected by Relion.
        If the input particles comes from Relion, just link the file. 
        """
        imgSet = self.inputParticles.get()
        # Create links to binary files and write the relion .star file
        writeSetOfParticles(imgSet, self._getPath('input_particles.star'), 
                            outputDir=self._getExtraPath(),
                            writeAlignment=False,
                            postprocessImageRow=self._postprocessImageRow)
    
    def processStep(self):
        # Enter here to generate the star file or to preprocess the images
        
        outputRadius = self._getOutputRadius()
        params = ' --operate_on input_particles.star'
        
        if self.doNormalize:
            radius = self.backRadius.get()
            if radius <= 0:
                radius = outputRadius
            params += ' --norm --bg_radius %d' % radius
            
        if self.doRemoveDust:
            wDust = self.whiteDust.get()
            if wDust > 0:
                params += ' --white_dust %f' % wDust
            bDust = self.blackDust.get()
            if bDust > 0:
                params += ' --black_dust %f' % bDust
                
        if self.doInvert:
            params += ' --invert_contrast'            
        
        if self.doScale:
            params += ' --scale %d' % self.scaleSize.get()

        if self.doWindow:
            params += ' --window %d' % self.windowSize.get()

        self.runJob(self._getProgram('relion_preprocess'), params, cwd=self._getPath())
        
        outputMrcs = glob(self._getPath('particles*.mrcs'))[0] # In Relion 1.3 it is produces particles.mrcs.mrcs
        partFn = self._getPath('particles.mrcs')
        if outputMrcs != partFn:
            # Override the initial converted mrcs particles stack
            # It also make easy to use the same .star file as output
            moveFile(outputMrcs, partFn)
    
    def createOutputStep(self):
        inputSet = self.inputParticles.get()
        
        if isinstance(inputSet, em.SetOfAverages):
            imgSet = self._createSetOfAverages()
        else:
            imgSet = self._createSetOfParticles()
            
        imgSet.copyInfo(inputSet)
        if self.doScale:
            oldSampling = inputSet.getSamplingRate()
            scaleFactor = self._getScaleFactor(inputSet)
            newSampling = oldSampling * scaleFactor
            imgSet.setSamplingRate(newSampling)
        
        for i, img in enumerate(inputSet):
            img.setLocation(i+1, self._getPath('particles.mrcs'))
            if self.doScale and inputSet.hasAlignment():
                a = img.getTransform()
                m = a.getMatrix()
                # Check if for scale we only need to scale shifts
                m[0, 3] *= 1/scaleFactor
                m[1, 3] *= 1/scaleFactor
            imgSet.append(img)
        
        self._defineOutputs(outputParticles=imgSet)
        self._defineTransformRelation(inputSet, imgSet)

#--------------------------- INFO functions -------------------------------------------- 
    def _validate(self):
        """ Should be overriden in subclasses to 
        return summary message for NORMAL EXECUTION. 
        """
        validateMsgs = []
        if self.doScale and self.scaleSize.get() % 2 != 0: 
            validateMsgs += ["Only re-scaling to even-sized images is allowed in RELION."]
        if self.doWindow and self.windowSize.get() % 2 != 0:
            validateMsgs += ["Only re-windowing to even-sized images is allowed in RELION."]

        if self.doNormalize:
            outputRadius = self._getOutputRadius()
            if self.backRadius > outputRadius:
                validateMsgs.append('Set a normalization background radius less than '
                                    'the particles output radius (%d pixels). ' % outputRadius)
        return validateMsgs
    
    def _summary(self):
        summary = []
        summary.append('Operations applied:')
        if self.doScale:
            summary.append("- Particles scaled to *%d* pixels" % self.scaleSize.get())
        if self.doWindow:
            summary.append("- Particles windowed to *%d* pixels" % self.windowSize.get())
        if self.doNormalize:
            summary.append("- Normalization, background radius *%d* pixels" % self.backRadius.get())
        if self.doRemoveDust:
            if self.whiteDust > 0:
                summary.append('- Removed white dust (sigma=%0.3f)' % self.whiteDust.get())
            if self.blackDust > 0:
                summary.append('- Removed black dust (sigma=%0.3f)' % self.blackDust.get())
        if self.doInvert:
            summary.append('- Inverted contrast')   

        return summary
    
    def _methods(self):
        summary = 'The input particles were preprocessed using Relion.'
        if self.doScale:
            summary += " Scaled to *%d* pixels." % self.scaleSize.get()
        if self.doWindow:
            summary += " Windowed to *%d* pixels." % self.windowSize.get()
        if self.doNormalize:
            summary += " The particles were normalized using a background radius of *%d* pixels." % self.backRadius.get()
        if self.doRemoveDust:
            if self.whiteDust > 0:
                summary += ' White dust was removed (pixels with a sigma > %0.3f).' % self.whiteDust.get()
            if self.blackDust > 0:
                summary += ' Black dust was removed (pixels with a sigma > %0.3f).' % self.blackDust.get()
        if self.doInvert:
            summary += ' The original constrast was inverted.'
            
        return [summary]                
    #--------------------------- UTILS functions ---------------------------------------------------
    def _getOutputRadius(self):
        """ Get the radius of the output particles"""
        if self.doScale:
            radius = self.scaleSize.get()/2
        else:
            xdim = self.inputParticles.get().getDimensions()[0]
            radius = xdim / 2
        return radius
    
    def _postprocessImageRow(self, img, imgRow):
        """ Since relion_preprocess will runs in its working directory
        we need to modify the default image path (from project dir)
        and make them relative to run working dir.
        """
        from convert import relativeFromFileName
        relativeFromFileName(imgRow, self._getPath())
    
    def _getScaleFactor(self, inputSet):
        xdim = inputSet.getDim()[0]
        scaleFactor = xdim/float(self.scaleSize.get())
        return scaleFactor
        
