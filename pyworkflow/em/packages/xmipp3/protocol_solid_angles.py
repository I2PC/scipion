# **************************************************************************
# *
# * Authors:     C.O.S. Sorzano (coss@cnb.csic.es)
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

from os.path import join, exists
import math

import pyworkflow.protocol.params as params
from pyworkflow import VERSION_1_1
from pyworkflow.utils.path import makePath
from pyworkflow.em.convert import ImageHandler, ALIGN_PROJ
from pyworkflow.em.data import Image
from pyworkflow.em.protocol import ProtAnalysis3D
from pyworkflow.em.packages.xmipp3.convert import writeSetOfParticles
import pyworkflow.em.metadata as md

import xmipp
from convert import rowToAlignment, setXmippAttributes, xmippToLocation
from xmipp3 import findRow
from constants import SYM_URL


class XmippProtSolidAngles(ProtAnalysis3D):
    """    
    Construct image groups based on the angular assignment. All images assigned within a solid angle
    are assigned to a class. Classes are not exclusive and an image may be assigned to multiple classes
    """

    _label = 'solid angles'
    _lastUpdateVersion = VERSION_1_1
    
    def __init__(self, *args, **kwargs):
        ProtAnalysis3D.__init__(self, *args, **kwargs)
        
    #--------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        
        form.addParam('inputVolume', params.PointerParam, pointerClass='Volume',
                      label="Input volume",  
                      help='Select the input volume.')

        form.addParam('inputParticles', params.PointerParam,
                      pointerClass='SetOfParticles',
                      pointerCondition='hasAlignmentProj',
                      label="Input particles",  
                      help='Select the input experimental images with an '
                           'angular assignment.')

        form.addParam('symmetryGroup', params.StringParam, default='c1',
                      label="Symmetry group", 
                      help='See %s page for a description of the symmetries '
                           'accepted by Xmipp' % SYM_URL)

        form.addParam('angularSampling', params.FloatParam, default=5,
                      label='Angular sampling',
                      expertLevel=params.LEVEL_ADVANCED, help="In degrees")

        form.addParam('angularDistance', params.FloatParam, default=10,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Angular distance',
                      help="In degrees. An image belongs to a group if its "
                           "distance is smaller than this value")

        form.addParam('maxShift', params.FloatParam, default=15,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Maximum shift',
                      help="In pixels")

        form.addSection("Directional Classes")

        form.addParam('directionalClasses', params.IntParam, default=1,
                      label='Number of directional classes',
                      help="By default only one class will be computed for "
                           "each projection direction. More classes could be"
                           "computed and this is needed for protocols "
                           "split-volume. ")

        form.addParam('targetResolution', params.FloatParam, default=10,
                      condition="directionalClasses > 1",
                      label='Target resolution (A)')

        form.addParam('refineAngles', params.BooleanParam, default=True,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Refine angles',
                      help="Refine the angles of the classes using a "
                           "continuous angular assignment")

        form.addParam('cl2dIterations', params.IntParam, default=5,
                      expertLevel=params.LEVEL_ADVANCED,
                      condition="directionalClasses > 1",
                      label='Number of CL2D iterations')

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

        if self.refineAngles:
            self._insertFunctionStep('refineAnglesStep')

        self._insertFunctionStep('createOutputStep')
    
    #--------------------------- STEPS functions -------------------------------

    def convertInputStep(self, particlesId, volId):
        """ Write the input images as a Xmipp metadata file. 
        particlesId: is only need to detect changes in
        input particles and cause restart from here.
        """
        inputParticles = self.inputParticles.get()
        inputVolume = self.inputVolume.get()

        writeSetOfParticles(inputParticles, self._getExpParticlesFn())

        img = ImageHandler()
        img.convert(inputVolume, self._getInputVolFn())

        if self._useSeveralClasses():
            # Scale particles
            Xdim = inputParticles.getXDim()
            Ts = inputParticles.getSamplingRate()
            newTs = self.targetResolution.get() * 0.4
            newTs = max(Ts, newTs)
            newXdim = Xdim * Ts / newTs
            self.runJob("xmipp_image_resize",
                        "-i %s -o %s --save_metadata_stack %s --dim %d" %
                        (self._getExpParticlesFn(),
                         self._getTmpPath('scaled_particles.stk'),
                         self._getTmpPath('scaled_particles.xmd'),
                         newXdim))
            # Scale volume
            Xdim = inputVolume.getXDim()
            if Xdim != newXdim:
                self.runJob("xmipp_image_resize", "-i %s --dim %d"
                            % (self._getInputVolFn(), newXdim), numberOfMpi=1)

    def constructGroupsStep(self, particlesId, angularSampling,
                            angularDistance, symmetryGroup):
        args = '-i %s ' % self._getInputVolFn()
        args += '-o %s ' % self._getExtraPath("gallery.stk")
        args += '--sampling_rate %f ' % self.angularSampling
        args += '--sym %s ' % self.symmetryGroup
        args += '--method fourier 1 0.25 bspline --compute_neighbors '
        args += '--angular_distance %f ' % self.angularDistance
        args += '--experimental_images %s ' % self._getInputParticlesFn()
        args += '--max_tilt_angle 90 '

        # Create a gallery of projections of the input volume
        # with the given angular sampling
        self.runJob("xmipp_angular_project_library", args)

        args = '--i1 %s ' % self._getInputParticlesFn()
        args += '--i2 %s ' % self._getExtraPath("gallery.doc")
        args += '-o %s ' % self._getExtraPath("neighbours.xmd")
        args += '--dist %f ' % self.angularDistance
        args += '--sym %s ' % self.symmetryGroup
        args += '--check_mirrors '

        # Compute several groups of the experimental images into
        # different angular neighbourhoods
        self.runJob("xmipp_angular_neighbourhood", args, numberOfMpi=1)

    def classifyOneGroup(self, projNumber, projMdBlock, projRef,
                         mdClasses, mdImages):
        """ Classify one of the neighbourhood groups if not empty.
         Class information will be stored in output metadata: mdOut
        """
        blockSize = md.getSize(projMdBlock)
        Nclasses = self.directionalClasses.get()
        Nlevels = int(math.ceil(math.log(Nclasses) / math.log(2)))

        # Skip projection directions with not enough images to
        # create a given number of classes
        if blockSize / Nclasses < 10:
            return

        fnDir = self._getExtraPath("direction_%s" % projNumber)
        makePath(fnDir)

        # Run CL2D classification for the images assigned to one direction
        args = "-i %s " % projMdBlock
        args += "--odir %s " % fnDir
        args += "--ref0 %s --iter 1 --nref %d " % (projRef, Nclasses)
        args += "--distance correlation --classicalMultiref "
        args += "--maxShift %f " % self.maxShift
	try:
            self.runJob("xmipp_classify_CL2D", args)
	except:
	    return 

        # After CL2D the stk and xmd files should be produced
        classesXmd = join(fnDir, "level_%02d/class_classes.xmd" % Nlevels)
        classesStk = join(fnDir, "level_%02d/class_classes.stk" % Nlevels)

        # Let's check that the output was produced
        if not exists(classesStk):
            return

        # Run align of the class average and the projection representative
        fnAlignRoot = join(fnDir, "classes")
        args = "-i %s " % classesStk
        args += "--ref %s " % projRef
        args += " --oroot %s --iter 1" % fnAlignRoot
        self.runJob("xmipp_image_align", args, numberOfMpi=1)

        # Apply alignment
        args = "-i %s_alignment.xmd --apply_transform" % fnAlignRoot
        self.runJob("xmipp_transform_geometry", args, numberOfMpi=1)

        for classNo in range(1, Nclasses+1):
            localImagesMd = xmipp.MetaData("class%06d_images@%s"
                                           % (classNo, classesXmd))

            # New class detected
            self.classCount += 1
            # Check which images have not been assigned yet to any class
            # and assign them to this new class
            for objId in localImagesMd:
                imgId = localImagesMd.getValue(xmipp.MDL_ITEM_ID, objId)
                # Add images not classify yet and store their class number
                if imgId not in self.classImages:
                    self.classImages.add(imgId)
                    newObjId = mdImages.addObject()
                    mdImages.setValue(xmipp.MDL_ITEM_ID, imgId, newObjId)
                    mdImages.setValue(xmipp.MDL_REF2, self.classCount, newObjId)

            newClassId = mdClasses.addObject()
            mdClasses.setValue(xmipp.MDL_REF, projNumber, newClassId)
            mdClasses.setValue(xmipp.MDL_REF2, self.classCount, newClassId)
            mdClasses.setValue(xmipp.MDL_IMAGE, "%d@%s" %
                               (classNo, classesStk), newClassId)
            mdClasses.setValue(xmipp.MDL_IMAGE1, projRef, newClassId)
            mdClasses.setValue(xmipp.MDL_CLASS_COUNT,localImagesMd.size(),newClassId)

    def classifyGroupsStep(self):
        # Create two metadatas, one for classes and another one for images
        mdClasses = xmipp.MetaData()
        mdImages = xmipp.MetaData()

        fnNeighbours = self._getExtraPath("neighbours.xmd")
        fnGallery = self._getExtraPath("gallery.stk")

        self.classCount = 0
        self.classImages = set()

        for block in xmipp.getBlocksInMetaDataFile(fnNeighbours):
            # Figure out the projection number from the block name
            projNumber = int(block.split("_")[1])

            self.classifyOneGroup(projNumber,
                                  projMdBlock="%s@%s" % (block, fnNeighbours),
                                  projRef="%06d@%s" % (projNumber, fnGallery),
                                  mdClasses=mdClasses,
                                  mdImages=mdImages)

        galleryMd = xmipp.MetaData(self._getExtraPath("gallery.doc"))
        # Increment the reference number to starts from 1
        galleryMd.operate("ref=ref+1")
        mdJoined = xmipp.MetaData()
        # Add extra information from the gallery metadata
        mdJoined.join1(mdClasses, galleryMd, xmipp.MDL_REF)
        # Remove unnecessary columns
        md.keepColumns(mdJoined, "ref", "ref2", "image", "image1",
                    "classCount", "angleRot", "angleTilt")

        # Write both classes and images
        fnDirectional = self._getDirectionalClassesFn()
        self.info("Writting classes info to: %s" % fnDirectional)
        mdJoined.write(fnDirectional)

        fnDirectionalImages = self._getDirectionalImagesFn()
        self.info("Writing images info to: %s" % fnDirectionalImages)
        mdImages.write(fnDirectionalImages)

    def refineAnglesStep(self):
        # TODO: Implement this step from the current DirectionalClasses protocol
        pass

    def createOutputStep(self):
        self.mdClasses = xmipp.MetaData(self._getDirectionalClassesFn())
        self.mdImages = xmipp.MetaData(self._getDirectionalImagesFn())

        inputParticles = self.inputParticles.get()
        classes2D = self._createSetOfClasses2D(inputParticles)

        self.averageSet = self._createSetOfAverages()
        self.averageSet.copyInfo(inputParticles)
        self.averageSet.setAlignmentProj()

        # Let's use a SetMdIterator because it could be less particles
        # in the metadata produced than in the input set
        iterator = md.SetMdIterator(self.mdImages, sortByLabel=md.MDL_ITEM_ID,
                                    updateItemCallback=self._updateParticle,
                                    skipDisabled=True)

        classes2D.classifyItems(updateItemCallback=iterator.updateItem,
                                updateClassCallback=self._updateClass)

        self._defineOutputs(outputClasses=classes2D)
        self._defineOutputs(outputAverages=self.averageSet)
        self._defineSourceRelation(self.inputParticles, classes2D)
        self._defineSourceRelation(self.inputParticles, self.averageSet)

    def _updateParticle(self, item, row):
        item.setClassId(row.getValue(xmipp.MDL_REF2))

        # item.setTransform(rowToAlignment(row, em.ALIGN_2D))
        #
        # item._rlnNormCorrection = em.Float(row.getValue('rlnNormCorrection'))
        # item._rlnLogLikeliContribution = em.Float(
        #     row.getValue('rlnLogLikeliContribution'))
        # item._rlnMaxValueProbDistribution = em.Float(
        #     row.getValue('rlnMaxValueProbDistribution'))

    def _updateClass(self, item):
        classId = item.getObjId()
        classRow = findRow(self.mdClasses, xmipp.MDL_REF2, classId)

        representative = item.getRepresentative()
        representative.setTransform(rowToAlignment(classRow, ALIGN_PROJ))
        representative.setLocation(xmippToLocation(classRow.getValue(xmipp.MDL_IMAGE)))
        setXmippAttributes(representative, classRow, xmipp.MDL_ANGLE_ROT)
        setXmippAttributes(representative, classRow, xmipp.MDL_ANGLE_TILT)
        setXmippAttributes(representative, classRow, xmipp.MDL_CLASS_COUNT)

        self.averageSet.append(representative)

        reprojection = Image()
        reprojection.setLocation(xmippToLocation(classRow.getValue(xmipp.MDL_IMAGE1)))
        item.reprojection = reprojection

    # def createOutputStep(self):
    #     fnDirectional=self._getPath("directionalClasses.xmd")
    #     if exists(fnDirectional):
    #         inputParticles = self.inputParticles.get()
    #         self._sampling = inputParticles.getSamplingRate()
    #         classes2DSet = self._createSetOfClasses2D(inputParticles)
    #         readSetOfClasses2D(classes2DSet, fnDirectional, 'classes',
    #                            preprocessClass=self._readRow)
    #         self._defineOutputs(outputClasses=classes2DSet)
    #         self._defineSourceRelation(self.inputParticles, classes2DSet)
    #
    # def _readRow(self, classItem, row):
    #     img = Image()
    #     img.setSamplingRate(self._sampling)
    #     img.setLocation(xmippToLocation(row.getValue(xmipp.MDL_IMAGE1)))
    #     classItem.reprojection = img
    #
    #     setXmippAttributes(classItem, row, xmipp.MDL_ANGLE_ROT)
    #     setXmippAttributes(classItem, row, xmipp.MDL_ANGLE_TILT)

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

    def _useSeveralClasses(self):
        return self.directionalClasses > 1

    def _getExpParticlesFn(self):
        return self._getPath('input_particles.xmd')

    def _getInputParticlesFn(self):
        if self._useSeveralClasses():
            return self._getTmpPath('scaled_particles.xmd')
        else:
            return self._getExpParticlesFn()

    def _getInputVolFn(self):
        return self._getTmpPath('volume.vol')

    def _getDirectionalClassesFn(self):
        return self._getPath("directional_classes.xmd")

    def _getDirectionalImagesFn(self):
        return self._getPath("directional_images.xmd")
