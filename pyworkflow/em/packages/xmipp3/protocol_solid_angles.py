# **************************************************************************
# *
# * Authors:     Javier Vargas (jvargas@cnb.csic.es)
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

from os.path import join, exists

from pyworkflow.protocol.params import PointerParam, FloatParam, StringParam, LEVEL_ADVANCED
from pyworkflow.utils.path import makePath
from pyworkflow.em.convert import ImageHandler
from pyworkflow.em.data import Image
from pyworkflow.em.protocol import ProtAnalysis3D
from pyworkflow.em.packages.xmipp3.convert import writeSetOfParticles, xmippToLocation
from pyworkflow.em.metadata.utils import getSize, keepColumns

import xmipp
from convert import readSetOfClasses2D, setXmippAttributes
from constants import SYM_URL


class XmippProtSolidAngles(ProtAnalysis3D):
    """    
    Construct image groups based on the angular assignment. All images assigned within a solid angle
    are assigned to a class. Classes are not exclusive and an image may be assigned to multiple classes
    """

    _label = 'solid angles'
    
    def __init__(self, *args, **kwargs):
        ProtAnalysis3D.__init__(self, *args, **kwargs)
        
    #--------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        
        form.addParam('inputVolume', PointerParam, pointerClass='Volume',
                      label="Input volume",  
                      help='Select the input volume.')     
        form.addParam('inputParticles', PointerParam,
                      pointerClass='SetOfParticles',
                      pointerCondition='hasAlignmentProj',
                      label="Input particles",  
                      help='Select the input experimental images with an '
                           'angular assignment.')
        form.addParam('symmetryGroup', StringParam, default='c1',
                      label="Symmetry group", 
                      help='See %s page for a description of the symmetries '
                           'accepted by Xmipp' % SYM_URL)

        form.addParam('angularSampling', FloatParam, default=5,
                      label='Angular sampling',
                      expertLevel=LEVEL_ADVANCED, help="In degrees")
        form.addParam('angularDistance', FloatParam, default=10,
                      label='Angular distance', expertLevel=LEVEL_ADVANCED,
                      help="In degrees. An image belongs to a group if its "
                           "distance is smaller than this value")
        form.addParam('maxShift', FloatParam, default=15,
                      label='Maximum shift', expertLevel=LEVEL_ADVANCED,
                      help="In pixels")
        
        form.addParallelSection(threads=0, mpi=8)
    
    #--------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):        
        self._insertFunctionStep('convertInputStep',
                                 self.inputParticles.get().getObjId(),
                                 self.inputVolume.get().getObjId())
        self._insertFunctionStep('constructGroupsStep',
                                 self.inputParticles.get().getObjId(),
                                 self.angularSampling.get(),
                                 self.angularDistance.get(),
                                 self.symmetryGroup.get())
        self._insertFunctionStep('classifyGroupsStep')
        self._insertFunctionStep('createOutputStep')
    
    #--------------------------- STEPS functions -------------------------------
    def convertInputStep(self, particlesId, volId):
        """ Write the input images as a Xmipp metadata file. 
        particlesId: is only need to detect changes in
        input particles and cause restart from here.
        """
        writeSetOfParticles(self.inputParticles.get(),
                            self._getExpParticlesFn())

        img = ImageHandler()
        img.convert(self.inputVolume.get(), self._getTmpPath("volume.vol"))
   
    def constructGroupsStep(self, particlesId, angularSampling,
                            angularDistance, symmetryGroup):
        args = '-i %s ' % self._getTmpPath("volume.vol")
        args += '-o %s ' % self._getExtraPath("gallery.stk")
        args += '--sampling_rate %f ' % self.angularSampling
        args += '--sym %s ' % self.symmetryGroup
        args += '--method fourier 1 0.25 bspline --compute_neighbors '
        args += '--angular_distance %f ' % self.angularDistance
        args += '--experimental_images %s ' % self._getExpParticlesFn()
        args += '--max_tilt_angle 90 '

        # Create a gallery of projections of the input volume
        # with the given angular sampling
        self.runJob("xmipp_angular_project_library", args)

        args = '--i1 %s ' % self._getExpParticlesFn()
        args += '--i2 %s ' % self._getExtraPath("gallery.doc")
        args += '-o %s ' % self._getExtraPath("neighbours.xmd")
        args += '--dist %f ' % self.angularDistance
        args += '--sym %s ' % self.symmetryGroup
        args += '--check_mirrors '

        # Compute several groups of the experimental images into
        # different angular neighbourhoods
        self.runJob("xmipp_angular_neighbourhood", args, numberOfMpi=1)

    def classifyOneGroup(self, projNumber, projMdBlock, projRef, mdOut):
        """ Classify one of the neighbourhood groups if not empty.
         Class information will be stored in output metadata: mdOut
        """
        blockSize = getSize(projMdBlock)

        if blockSize == 0: # Skip empty blocks, i.e., directions with no images
            return

        fnDir = self._getExtraPath("direction_%s" % projNumber)
        makePath(fnDir)
        fnOut = join(fnDir, "level_00/class_classes.stk")

        if not exists(fnOut):
            # Run CL2D classification for the images assigned to one direction
            args = "-i %s " % projMdBlock
            args += "--odir %s " % fnDir
            args += "--ref0 %s --iter 1 --nref 1 " % projRef
            args += "--distance correlation --classicalMultiref "
            args += "--maxShift %f " % self.maxShift
            self.runJob("xmipp_classify_CL2D", args)

            # Run align of the class and the projection
            fnAlignRoot = join(fnDir, "classes")
            args = "-i %s " % fnOut
            args += "--ref %s " % projRef
            args += " --oroot %s --iter 1" % fnAlignRoot
            self.runJob("xmipp_image_align", args, numberOfMpi=1)

            # Apply alignment
            args = "-i %s_alignment.xmd --apply_transform" % fnAlignRoot
            self.runJob("xmipp_transform_geometry", args, numberOfMpi=1)

        # Add a new entry in the output metadata for this class
        imagesMd = xmipp.MetaData("class000001_images@%s"
                                  % join(fnDir, "level_00/class_classes.xmd"))
        fnDirectional = self._getPath("directionalClasses.xmd")
        imagesMd.write("class%06d_images@%s" % (projNumber, fnDirectional),
                      xmipp.MD_APPEND)

        objId = mdOut.addObject()
        mdOut.setValue(xmipp.MDL_REF, projNumber, objId)
        mdOut.setValue(xmipp.MDL_IMAGE, projRef, objId)
        mdOut.setValue(xmipp.MDL_IMAGE1, "1@%s" % fnOut, objId)
        mdOut.setValue(xmipp.MDL_CLASS_COUNT, imagesMd.size(), objId)

    def classifyGroupsStep(self):
        # For each proje
        mdOut = xmipp.MetaData()

        fnNeighbours = self._getExtraPath("neighbours.xmd")
        fnGallery = self._getExtraPath("gallery.stk")
        fnDirectional = self._getPath("directionalClasses.xmd")
        # Write an empty block to put classes as first block
        mdOut.write("classes@" + fnDirectional)

        for block in xmipp.getBlocksInMetaDataFile(fnNeighbours):
            # Figure out the projection number from the block name
            projNumber = int(block.split("_")[1])

            self.classifyOneGroup(projNumber,
                                  projMdBlock="%s@%s" % (block, fnNeighbours),
                                  projRef="%06d@%s" % (projNumber, fnGallery),
                                  mdOut=mdOut)

        galleryMd = xmipp.MetaData(self._getExtraPath("gallery.doc"))
        # Increment the reference number to starts from 1
        galleryMd.operate("ref=ref+1")
        mdJoined = xmipp.MetaData()
        # Add extra information from the gallery metadata
        mdJoined.join1(mdOut, galleryMd, xmipp.MDL_REF)
        # Remove unnecessary columns
        keepColumns(mdJoined, "ref", "image", "image1",
                    "classCount", "angleRot", "angleTilt")
        outputBlock = "classes@" + fnDirectional
        self.info("Writting classes info to: %s" % outputBlock)
        mdJoined.write(outputBlock, xmipp.MD_APPEND)

    def createOutputStep(self):
        fnDirectional=self._getPath("directionalClasses.xmd")
        if exists(fnDirectional):
            inputParticles = self.inputParticles.get()
            self._sampling = inputParticles.getSamplingRate()
            classes2DSet = self._createSetOfClasses2D(inputParticles)
            readSetOfClasses2D(classes2DSet, fnDirectional, 'classes',
                               preprocessClass=self._readRow)
            self._defineOutputs(outputClasses=classes2DSet)
            self._defineSourceRelation(self.inputParticles, classes2DSet)
    
    def _readRow(self, classItem, row):
        img = Image()
        img.setSamplingRate(self._sampling)
        img.setLocation(xmippToLocation(row.getValue(xmipp.MDL_IMAGE1)))
        classItem.reprojection = img
        
        setXmippAttributes(classItem, row, xmipp.MDL_ANGLE_ROT)
        setXmippAttributes(classItem, row, xmipp.MDL_ANGLE_TILT)

    # --------------------------- INFO functions -------------------------------
    def _validate(self):
        validateMsgs = []
        # if there are Volume references, it cannot be empty.
        if self.inputVolume.get() and not self.inputVolume.hasValue():
            validateMsgs.append('Please provide an input reference volume.')
        if self.inputParticles.get() and not self.inputParticles.hasValue():
            validateMsgs.append('Please provide input particles.')            
        return validateMsgs
    
    def _summary(self):
        summary = []
        return summary


    #----------------------- UTILITY FUNCTIONS ---------------------------------
    def _getExpParticlesFn(self):
        return self._getPath('input_particles.xmd')

    def _getVolFn(self):
        return self._getExtraPath('volume.vol')
    