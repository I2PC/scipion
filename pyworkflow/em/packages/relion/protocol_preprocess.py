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
This module contains the protocol base class for Relion protocols
"""
from pyworkflow.em import * 
from pyworkflow.utils import environAdd, moveFile, cleanPath
from convert import writeSetOfParticles, readSetOfParticles

from protocol_base import ProtRelionBase

class ProtRelionPreprocessParticles(ProtProcessParticles, ProtRelionBase):
    """ Wrapper to Relion preprocess program.
    """
    _label = 'preprocess particles'
    
    def __init__(self, **args):
        ProtProcessParticles.__init__(self, **args)
                   
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', PointerParam, pointerClass='SetOfParticles',
                      label="Input particles",  
                      help='Select the input images from the project.')   
        form.addParam('doNormalize', BooleanParam, default=False,
                      label='Normalize',
                      help='If set to True, particles will be normalized in the way RELION prefers it.')
        form.addParam('backRadius', IntParam, default=0,
                      condition='doNormalize',
                      label='Background radius',
                      help='Pixels outside this circle are assumed to be noise and their stddev '
                      'is set to 1. Radius for background circle definition (in pix.).')
    
    #--------------------------- INSERT steps functions --------------------------------------------
    
    def _insertAllSteps(self):
        self.imgStar = self._getPath('input_particles.star')
        self.imgFn = self._getPath('input_particles.mrcs')
        self._insertFunctionStep("convertInputStep")
        print "IMAGES STAR=%s" % self.imgStar
        self._insertProcessStep()
        
        self._insertFunctionStep('createOutputStep')
        
    def convertInputStep(self):
        """ Create the input file in STAR format as expected by Relion.
        If the input particles comes from Relion, just link the file. 
        """
        imgSet = self.inputParticles.get()
        writeSetOfParticles(imgSet, self.imgStar, self.imgFn)
        
    def _insertProcessStep(self):
        if self.doNormalize:
            self._insertFunctionStep("normalizeStep")
            
    def normalizeStep(self):

        # Enter here to generate the star file or to normalize the images
        imgSet = self.inputParticles.get()
        Xdim = imgSet.getDimensions()[0]
        
        self._enterDir(self._getPath())
        imgFn = os.path.relpath(self.imgFn, self._getPath())
        radius = self.backRadius.get()
        if radius <= 0:
            radius = Xdim / 2
        params = '--operate_on %(imgFn)s --norm --bg_radius %(radius)s'
        self.runJob(self._getProgram('relion_preprocess'), params % locals())
         
        outputMrcs = glob('particles*.mrcs') # In Relion 1.3 it is produces particles.mrcs.mrcs
        if not outputMrcs:
            raise Exception("Not particles produced by 'relion_preprocess'")
        outFn = outputMrcs[0]
        outStar = 'particles.star'
         
        moveFile(outFn, imgFn)
        cleanPath(outStar)
        self._leaveDir()

              
    def createOutputStep(self):
        imgSet = self._createSetOfParticles()
        imgSet.copyInfo(self.inputParticles.get())
        readSetOfParticles(self.imgStar, imgSet)
        self._defineOutputs(outputParticles=imgSet)
        self._defineTransformRelation(self.inputParticles.get(), self.outputParticles)

#--------------------------- INFO functions -------------------------------------------- 
    def _validate(self):
        """ Should be overriden in subclasses to 
        return summary message for NORMAL EXECUTION. 
        """
        return []
    
    def _summary(self):
        """ Should be overriden in subclasses to 
        return summary message for NORMAL EXECUTION. 
        """
        return []
            