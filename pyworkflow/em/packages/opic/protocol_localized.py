# *****************************************************************************
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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# *****************************************************************************
"""
This module contains the protocol for localized reconstruction.
"""

from pyworkflow.em import ALIGN_PROJ
from pyworkflow.protocol.params import (PointerParam, BooleanParam, 
                                        StringParam, IntParam, Positive,
                                        EnumParam, NumericRangeParam,
                                        PathParam, FloatParam)
from pyworkflow.em.protocol import ProtParticles

from convert import (writeSetOfParticles, convertBinaryVol, getEnviron,
                     getProgram, readSetOfParticles)
# import re
# from glob import glob
# from os.path import exists
# 
# from pyworkflow.utils.path import cleanPath

# import pyworkflow.em.metadata as md
# from pyworkflow.em.data import SetOfClasses3D

# from constants import ANGULAR_SAMPLING_LIST, MASK_FILL_ZERO


CMM = 0
HAND = 1


class ProtLocalizedRecons(ProtParticles):
    """ This class cointains a re-implementation to a method for the
    localized three-dimensional reconstruction of such subunits. 
    After determining the particle orientations, local areas 
    corresponding to the subunits can be extracted and treated as 
    single particles.
    """
    
    _label = 'localized reconstruction'
    
    def __init__(self, **args):        
        ProtParticles.__init__(self, **args)
    
    def _initialize(self):
        self._createFilenameTemplates()
    
    def _createFilenameTemplates(self):
        """ Centralize how files are called. """
        myDict = {
                  'input_star': self._getPath('input_particles.star'),
                  'output': self._getExtraPath('output_particles'),
                  'output_star': self._getExtraPath('output_particles.star')
                  }
        self._updateFilenamesDict(myDict)
    
    #--------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', PointerParam,
                      pointerClass='SetOfParticles',
                      pointerCondition='hasAlignmentProj',
                      important=True,
                      label="Input particles",  
                      help='Select the input images from the project.')

        group = form.addGroup('Symmetry')
        group.addParam('symmetryGroup', StringParam, default='c1',
                      label="Symmetry", 
                      help='If the molecule is asymmetric, set Symmetry group '
                           'to C1. Note their are multiple possibilities for '
                           'icosahedral symmetry: \n'
                           '* I1: No-Crowther 222 (standard in Heymann, '
                           'Chagoyen & Belnap, JSB, 151 (2005) 196-207)\n'
                           '* I2: Crowther 222 \n'
                           '* I3: 52-setting (as used in SPIDER?) \n'
                           '* I4: A different 52 setting \n')

        group.addParam('randomize', BooleanParam, default=False,
                      label='Randomize the order of the symmetry matrices?',
                      help='Useful for preventing preferred orientations.')
        group.addParam('relaxSym', BooleanParam, default=False,
                      label='Relax symmetry?',
                      help='Create one random subparticle for each particle ')

        group = form.addGroup('Vectors')
        group.addParam('defineVector',  EnumParam, default=CMM,
                      label='Is vector defined by?',
                      choices=['cmm file', 'string'],
                      display=EnumParam.DISPLAY_HLIST,
                      help='')
        group.addParam('vector', NumericRangeParam, default='0,0,1',
                      label='Location vectors', condition="defineVector==1",
                      help='Vector defining the location of the '
                           'subparticles. The vector is defined by 3 '
                           'values x,y,z separated by comma. \n'
                           'More than one vector can be specified separated by'
                           'semicolon. For example: \n'
                           '0,0,1            # Defines only one vector.\n'
                           '0,0,1; 1,0,0;    # Defines two vectors.'
                       )
        group.addParam('vectorFile', PathParam, default='',
                      condition="defineVector==0",
                      label='file obtained by Chimera: ',
                      help='CMM file defining the location(s) of the '
                           'sub-particle(s). Use instead of vector. ')
        group.addParam('length', FloatParam, default=-1,
                      label='Alternative length of the vector (A)',
                      help='Use to adjust the sub-particle center. If it '
                           'is <= 0, the length of the given vector. '
                           'Different values must be separated by commas.')

        form.addSection('Sub-particles')
        form.addParam('alignSubparticles', BooleanParam, default=True,
                      label='Align the subparticles?',
                      help='Align sub-particles to the standard orientation. ')
        form.addParam('unique', FloatParam, default=-1,
                      label='Angle to keep unique sub-particles: ',
                      help='Keep only unique subparticles within angular '
                           'distance. It is useful to remove overlapping '
                           'sub-particles on symmetry axis.')
        form.addParam('mindist', FloatParam, default=-1,
                      label='Minimum distance between sub-particles',
                      help='In pixels. Minimum distance between the '
                           'subparticles in the image. All overlapping ones '
                           'will be discarded.')
        form.addParam('side', FloatParam, default=-1,
                      label='Angle to keep sub-particles from side views',
                      help='Keep only particles within specified angular '
                           'distance from side views. All others will be '
                           'discarded. ')
        form.addParam('top', FloatParam, default=-1,
                      label='Angle to keep sub-particles from top views',
                      help='Keep only particles within specified angular '
                           'distance from top views. All others will be '
                           'discarded. ')
        form.addParallelSection(threads=1, mpi=3)
    
    #--------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        partsId = self.inputParticles.get().getObjId()
        self._initialize()
        self._insertFunctionStep('convertInputStep', partsId)
        self._insertFunctionStep('localizedReconsStep')
        self._insertFunctionStep('createOutputStep')
    
    #--------------------------- STEPS functions ------------------------------
    def convertInputStep(self, particlesId):
        """ Create the input file in STAR format as expected by Relion.
        Params:
            particlesId: use this parameters just to force redo of convert if 
                the input particles are changed.
        """
        
        imgSet = self._getInputParticles()
        imgStar = self._getFileName('input_star')

        self.info("Converting set from '%s' into '%s'" %
                           (imgSet.getFileName(), imgStar))
        
        # Pass stack file as None to avoid write the images files
        writeSetOfParticles(imgSet, imgStar, self._getExtraPath())
    
    def localizedReconsStep(self):
        params = {"symmetryGroup" : self.symmetryGroup.get(),
                  "output" : self._getFileName('output'),
                  "vector" : self.vector.get(),
                  "vectorFile" : self.vectorFile.get(),
                  "length" : self.length.get(),
                  "unique" : self.unique.get(),
                  "mindist" : self.mindist.get(),
                  "side" : self.side.get(),
                  "top" : self.top.get(),
                  "pxSize" : self.inputParticles.get().getSamplingRate(),
                  "dim" : self.inputParticles.get().getXDim()
                  }
        
        args = ""

        if self.randomize:
            args += "--randomize "
        
        if self.relaxSym:
            args += "--relax_symmetry "
        
        if self.alignSubparticles:
            args += "--align_subparticles "
        
        args += "--j %d --np %d " % (self.numberOfThreads, self.numberOfMpi)
        
        if self.defineVector == CMM:
            args += "--cmm %s " % self.vectorFile
        else:
            args += "--vector %s " % self.vector
        
        if params["length"] > 0:
            args += "--length %f " % params["length"]

        # We use --subparitcle_size of 1 because we are not going to
        # extract the subparticles at this point, so we really need the
        # subparticle box size later
        args += ("--output %(output)s --angpix %(pxSize)f "
                 "--sym %(symmetryGroup)s --particle_size %(dim)d "
                 "--subparticle_size 1 " ) % params
        
        if params["unique"] > 0:
            args += "--unique %f " % params["unique"]
        
        if params["mindist"] > 0:
            args += "--mindist %f " % params["mindist"]
        
        if params["side"] > 0:
            args += "--side %f " % params["side"]
        
        if params["top"] > 0:
            args += "--top %f " % params["top"]
        
        args += "%s " % self._getFileName('input_star')
        
        self.runJob(getProgram(), args, env=getEnviron(), numberOfMpi=1)
    
    def createOutputStep(self):
        inputSet = self._getInputParticles()
        outputSet = self._createSetOfParticles()
        outputSet.copyInfo(inputSet)
        readSetOfParticles(self._getFileName('output_star'),
                           outputSet, alignType=ALIGN_PROJ)
        
        self._defineOutputs(outputParticles=outputSet)
        self._defineSourceRelation(self.inputParticles, outputSet)
    
    #--------------------------- INFO functions -------------------------------
    def _validate(self):
        errors = []
        return errors
    
    def _citations(self):
        return ['Serban2015']

    def _summary(self):
        summary = []
        self._initialize()
        return summary
    
    def _methods(self):
        return []
    
    #--------------------------- UTILS functions ------------------------------
    def _getInputParticles(self):
        return self.inputParticles.get()
    
