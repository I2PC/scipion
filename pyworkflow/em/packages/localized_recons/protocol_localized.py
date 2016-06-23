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

from pyworkflow.em.protocol import ProtParticles
import pyworkflow.protocol.params as params

from pyworkflow.em.packages.relion.convert import (convertBinaryVol, 
                                                   writeSetOfParticles, 
                                                   getVersion)
from convert import vectorsFromCmm
from data import Vector, SetOfVectors
# import re
# from glob import glob
# from os.path import exists
# 
from pyworkflow.protocol.constants import LEVEL_ADVANCED
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
    label='localized reconstruction'
    
    def __init__(self, **args):        
        ProtParticles.__init__(self, **args)
    
    def _initialize(self):
        self._createFilenameTemplates()
    
    def _createFilenameTemplates(self):
        """ Centralize how files are called. """
        myDict = {
                  'input_star': self._getPath('input_particles.star'),
                  'output': self._getExtraPath('output_particles'),
                  'output_star': self._getExtraPath('output_particles.star'),
                  'volume_masked': self._getTmpPath('volume_masked.mrc'),
                  'dummy_mask': self._getTmpPath('dummy_mask.mrc')
                  }
        self._updateFilenamesDict(myDict)
    
    #--------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', params.PointerParam,
                      pointerClass='SetOfParticles',
                      pointerCondition='hasAlignmentProj',
                      important=True,
                      label="Input particles",  
                      help='Select the input images from the project.')
        form.addParam('splitParticles', params.BooleanParam, default=False,
                      label='Do split particle stacks?',
                      help='Split particle stacks (needs to be done once).')
        form.addParam('inputVolume', params.PointerParam, 
                      pointerClass='Volume',
                      label="Input volume", allowsNull=True,
                      help='Map with density that to be subtracted from '
                           'particle images.')
        form.addParam('refMask', params.PointerParam,
                      pointerClass='VolumeMask',
                      label='Reference mask (optional)', allowsNull=True,
                      help="The volume will be masked applying this mask.")
        form.addParam('symmetryGroup', params.StringParam, default='c1',
                      label="Symmetry", 
                      help='If the molecule is asymmetric, set Symmetry group '
                           'to C1. Note their are multiple possibilities for '
                           'icosahedral symmetry: \n'
                           '* I1: No-Crowther 222 (standard in Heymann, '
                           'Chagoyen & Belnap, JSB, 151 (2005) 196-207)\n'
                           '* I2: Crowther 222 \n'
                           '* I3: 52-setting (as used in SPIDER?) \n'
                           '* I4: A different 52 setting \n')
        form.addParam('boxSize', params.IntParam, default=0,
                      label='Sub-particle box size',
                      validators=[params.Positive],
                      help='In pixels.  Size of the sub-particle box')
        form.addParam('randomize', params.BooleanParam, default=False,
                      label='Randomize the order of the symmetry matrices?',
                      help='Useful for preventing preferred orientations.')
        form.addParam('relaxSym', params.BooleanParam, default=False,
                      label='Relax symmetry?',
                      help='Create one random subparticle for each particle ')
        form.addParam('defineVector',  params.EnumParam, default=CMM,
                      label='Is vector defined by?',
                      choices=['cmm file', 'by hand'], 
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='')
        form.addParam('vector', params.NumericRangeParam, default='0,0,1',
                      label='Location vectors', condition="defineVector==1",
                      help='Vector(s) defining the location of the '
                           'subparticle(s). The vector is defined by 3 '
                           'values(x,y,z). To define different vectors, must '
                           'be separated by dot and comma(;) among them.')
        form.addParam('vectorFile', params.PathParam, default='',
                      condition="defineVector==0",
                      label='file obtained by Chimera: ',
                      help='CMM file defining the location(s) of the '
                           'sub-particle(s). Use instead of vector. ')
        form.addParam('length', params.NumericRangeParam, default=-1,
                      experlevel=LEVEL_ADVANCED,
                      label='Alternative length of the vector (A)',
                      help='Use to adjust the sub-particle center. If it '
                           'is <= 0, the length of the given vector. '
                           'Different values must be separated by commas.')
        form.addParam('alignSubparticles', params.BooleanParam, default=False,
                      label='Align the subparticles?',
                      help='Align sub-particles to the standard orientation. ')
        form.addParam('unique', params.FloatParam, default=1,
                      label='Angle to keep unique sub-particles: ',
                      help='Keep only unique subparticles within angular '
                           'distance. It is useful to remove overlapping '
                           'sub-particles on symmetry axis.')
        form.addParam('mindist', params.FloatParam, default=5,
                      label='Minimum distance between sub-particles',
                      help='In pixels. Minimum distance between the '
                           'subparticles in the image. All overlapping ones '
                           'will be discarded.')
        form.addParam('side', params.FloatParam, default=25,
                      label='Angle to keep sub-particles from side views',
                      help='Keep only particles within specified angular '
                           'distance from side views. All others will be '
                           'discarded. ')
        form.addParam('top', params.FloatParam, default=25,
                      label='Angle to keep sub-particles from top views',
                      help='Keep only particles within specified angular '
                           'distance from top views. All others will be '
                           'discarded. ')
        form.addParallelSection(threads=1, mpi=3)
    
    #--------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        partsId = self.inputParticles.get().getObjId()
        volId = self.inputVolume.get().getObjId()
        
        self._initialize()
        self._insertFunctionStep('convertInputStep', partsId, volId)
        self._insertFunctionStep('splitStep')
        self._insertFunctionStep('localizedReconsStep')
        self._insertFunctionStep('createOutputStep')
    
    #--------------------------- STEPS functions ------------------------------
    def convertInputStep(self, particlesId, volId):
        """ Create the input file in STAR format as expected by Relion.
        Params:
            particlesId: use this parameters just to force redo of convert if 
                the input particles are changed.
            volId: use this parameters just to force redo of convert if 
                the input volume is changed.
        """
        
        imgSet = self._getInputParticles()
        imgStar = self._getFileName('input_star')

        self.info("Converting set from '%s' into '%s'" %
                           (imgSet.getFileName(), imgStar))
        
        # Pass stack file as None to avoid write the images files
        writeSetOfParticles(imgSet, imgStar, self._getExtraPath())
    
    def splitStep(self):
        import pyworkflow.em.packages.bsoft.getEnviron as bsEnviron
        import pyworkflow.em.packages.relion.getEnviron as relEnviron
        if self.splitParticles:
            print "Creating and splitting the particle stack..."
            print "Creating a normal stack from which nothing is subtracted..."
            x = self.inputParticles.get().getXDim()
            args = "-create %d,%d,%d " % (x, x, x)
            args += " -fill 0 %s" % self._getFileName('dummy_mask')
            self.runJob('beditimg', args, env=bsEnviron)
            
            
            
#             
#             run_command("beditimg -create " + str(part_image_size) + "," + str(part_image_size) + "," + str(part_image_size) + " -fill 0 dummy_mask.mrc", "")
#             run_command("relion_project --i dummy_mask.mrc --o " + path + "particles" + " --ang " + input_star_filename + " --subtract_exp --angpix " + str(angpix), "")
#             run_command("rm -f dummy_mask.mrc", "")
#             run_command("bsplit -digits 6 -first 1 " + path + "particles" + ".mrcs:mrc " + path + "particles" + ".mrc", "")
#             print "Finished splitting the particle stack!"
#             print " "
#     
#         if subtract_masked_map:
#             print "Creating and splitting the particle stack..."
#             print " Creating a stack from which the projections of the masked particle have been subtracted..."
#             run_command("relion_project --i " + masked_map + " --o " + path + "particles_subtracted" + " --ang " + input_star_filename + " --subtract_exp --ctf --angpix " + str(angpix))
#             print "Finished subtracting the model projections!"
#             print " "
#     
#             print "Splitting the subtracted particle stack..."
#             run_command("bsplit -digits 6 -first 1 " + path + "particles_subtracted" + ".mrcs:mrc " + path + "particles_subtracted" + ".mrc", "")
#             print "Finished splitting the subtracted particle stack!"
#             print " "
    
    def localizedReconsStep(self):
        vectorSet = SetOfVectors()
        pxSize = self.inputParticles.get().getSamplingRate()
        
        lengthList = self.getListFromValues(self.length.get(), float, ",")
        lengthList = [x/pxSize for x in lengthList]
        
        if self.defineVector.get () == CMM:
            cmmFn = self.vectorFile.get()
            vectorsFromCmm(cmmFn, vectorSet, pxSize)
            
            if lengthList[0] > 0:
                if len(lengthList) == vectorSet.getSize():
                    for i, vector in enumerate(vectorSet):
                        vector.setDistance(lengthList[i])
                else:
                    raise Exception("Error: The number of distances doesn't "
                                    "match the number of vectors!")
        else:
            vector = Vector()
            vectorList = self.getListFromValues(self.vector.get(), spacer=";")
            for i, l in enumerate(vectorList):
                vec = [float(a) for a in l.split(",")]
                vector.setVector(vec)
                
                if lengthList[0] > 0:
                    if len(lengthList) == len(vectorList):
                        vector.setDistance(lengthList[i])
                    elif len(lengthList) == 1:
                        vector.setDistance(lengthList[0])
                    else:
                        raise Exception("Error: The number of distances doesn't "
                                        "match the number of vectors!")
                else:
                    vector.setDistance(vector.length())
                vectorSet.append(vector)
        
        print "Creating subparticles using:"
        
        for vector in vectorSet:
            print "Vector: ", vector.getVector()
            print ""
            print "Length: %.2f pixels" % vector.getDistance()
            




    def createOutputStep(self):
        pass
    
    #--------------------------- INFO functions -------------------------------
    def _validate(self):
        errors = []
        self.validatePackageVersion('RELION_HOME', errors)
        if self._getInputParticles().isOddX():
            errors.append("Relion only works with even values for the image dimensions!")
        
            errors += self._validateNormal()
        return errors
    
    def _validateNormal(self):
        """ Should be overriden in subclasses to 
        return summary message for NORMAL EXECUTION. 
        """
        return []
    
    def _validateContinue(self):
        """ Should be overriden in subclasses to
        return summary messages for CONTINUE EXECUTION.
        """
        return []
    
    def _citations(self):
        cites = []
        return cites
    
    def _summary(self):
        self._initialize()

        lastIter = self._lastIter()

        if lastIter is not None:
            iterMsg = 'Iteration %d' % lastIter
            if self.hasAttribute('numberOfIterations'):
                iterMsg += '/%d' % self._getnumberOfIters()
        else:
            iterMsg = 'No iteration finished yet.'
        summary = [iterMsg]

        flip = '' if self._getInputParticles().isPhaseFlipped() else 'not '
        flipMsg = "Your images have %sbeen ctf-phase corrected" % flip
        summary.append(flipMsg)
        
        if self.doContinue:
            summary += self._summaryContinue()
        summary += self._summaryNormal()
        return summary
    
    def _summaryNormal(self):
        """ Should be overriden in subclasses to 
        return summary message for NORMAL EXECUTION. 
        """
        return []
    
    def _summaryContinue(self):
        """ Should be overriden in subclasses to
        return summary messages for CONTINUE EXECUTION.
        """
        return []
    
    def _methods(self):
        """ Should be overriden in each protocol.
        """
        return []
    
    #--------------------------- UTILS functions ------------------------------
    def getListFromValues(self, valuesStr, valType=str, spacer=" "):
        listValues = []
        for values in valuesStr.split(spacer):
            listValues.append(valType(values))
        return listValues

#     def _getProgram(self, program='relion_refine'):
#         """ Get the program name depending on the MPI use or not. """
#         if self.numberOfMpi > 1:
#             program += '_mpi'
#         return program
#     
#     def _getInputParticles(self):
#         return self.inputParticles.get()
#     
#     def _getIterNumber(self, index):
#         """ Return the list of iteration files, give the iterTemplate. """
#         result = None
#         files = sorted(glob(self._iterTemplate))
#         if files:
#             f = files[index]
#             s = self._iterRegex.search(f)
#             if s:
#                 result = int(s.group(1)) # group 1 is 3 digits iteration number
#         return result
#         
#     def _lastIter(self):
#         return self._getIterNumber(-1)
# 
#     def _firstIter(self):
#         return self._getIterNumber(0) or 1
#     
#     def _getIterClasses(self, it, clean=False):
#         """ Return a classes .sqlite file for this iteration.
#         If the file doesn't exists, it will be created by 
#         converting from this iteration data.star file.
#         """
#         data_classes = self._getFileName('classes_scipion', iter=it)
#         
#         if clean:
#             cleanPath(data_classes)
#         
#         if not exists(data_classes):
#             clsSet = self.OUTPUT_TYPE(filename=data_classes)
#             clsSet.setImages(self.inputParticles.get())
#             self._fillClassesFromIter(clsSet, it)
#             clsSet.write()
#             clsSet.close()
# 
#         return data_classes
#     
#     def _getIterData(self, it, **kwargs):
#         """ Sort the it??.data.star file by the maximum likelihood. """
#         data_sqlite = self._getFileName('data_scipion', iter=it)
#         
#         if not exists(data_sqlite):
#             iterImgSet = em.SetOfParticles(filename=data_sqlite)
#             iterImgSet.copyInfo(self._getInputParticles())
#             self._fillDataFromIter(iterImgSet, it)
#             iterImgSet.write()
#             iterImgSet.close()
#         
#         return data_sqlite
#     
#     def _splitInCTFGroups(self, imgStar):
#         """ Add a new colunm in the image star to separate the particles into ctf groups """
#         from convert import splitInCTFGroups
#         splitInCTFGroups(imgStar
#                          , self.defocusRange.get()
#                          , self.numParticles.get())
#     
#     def _getContinueIter(self):
#         continueRun = self.continueRun.get()
#         
#         if continueRun is not None:
#             continueRun._initialize()
#         
#         if self.doContinue:
#             if self.continueIter.get() == 'last':
#                 continueIter = continueRun._lastIter()
#             else:
#                 continueIter = int(self.continueIter.get())
#         else:
#             continueIter = 0
#             
#         return continueIter
#     
#     def _getnumberOfIters(self):
#         return self._getContinueIter() + self.numberOfIterations.get()
#     
#     def _postprocessImageRow(self, img, imgRow):
# #         from convert import locationToRelion
#         partId = img.getParticleId()
#         magnification = img.getAcquisition().getMagnification()
#         imgRow.setValue(md.RLN_PARTICLE_ID, long(partId))
#         imgRow.setValue(md.RLN_CTF_MAGNIFICATION, magnification)
#         imgRow.setValue(md.RLN_MICROGRAPH_NAME, "%06d@fake_movie_%06d.mrcs" %(img.getFrameId() + 1, img.getMicId()))
        
