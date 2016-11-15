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
from pyworkflow.em.data import Image
from pyworkflow.em.protocol import ProtAnalysis3D
from pyworkflow.em.packages.xmipp3.convert import writeSetOfParticles, xmippToLocation
from pyworkflow.em.metadata.utils import getSize

import xmipp
from convert import readSetOfClasses2D, setXmippAttributes


class XmippProtSolidAngles(ProtAnalysis3D):
    """    
    Construct image groups based on the angular assignment. All images assigned within a solid angle
    are assigned to a class. Classes are not exclusive and an image may be assigned to multiple classes
    """

    _label = 'solid angles'
    
    def __init__(self, *args, **kwargs):
        ProtAnalysis3D.__init__(self, *args, **kwargs)
        
    #--------------------------- DEFINE param functions --------------------------------------------   
    def _defineParams(self, form):
        form.addSection(label='Input')
        
        form.addParam('inputVolume', PointerParam, pointerClass='Volume',
                      label="Input volume",  
                      help='Select the input volume.')     
        form.addParam('inputParticles', PointerParam, pointerClass='SetOfParticles', 
                      label="Input particles",  
                      help='Select the input projection images with an angular assignment.') 
        form.addParam('symmetryGroup', StringParam, default='c1',
                      label="Symmetry group", 
                      help='See [[Xmipp Symmetry][http://www2.mrc-lmb.cam.ac.uk/Xmipp/index.php/Conventions_%26_File_formats#Symmetry]] page '
                           'for a description of the symmetry format accepted by Xmipp') 
        form.addParam('angularSampling', FloatParam, default=5, label='Angular sampling', expertLevel=LEVEL_ADVANCED, help="In degrees")
        form.addParam('angularDistance', FloatParam, default=10, label='Angular distance', expertLevel=LEVEL_ADVANCED,
                      help="In degrees. An image belongs to a group if its distance is smaller than this value")
        form.addParam('maxShift', FloatParam, default=15, label='Maximum shift', expertLevel=LEVEL_ADVANCED,
                      help="In pixels")
        
        form.addParallelSection(threads=0, mpi=8)
    
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):        
        self._insertFunctionStep('convertInputStep', self.inputParticles.get().getObjId(), self.inputVolume.get().getObjId())
        self._insertFunctionStep('constructGroupsStep', self.inputParticles.get().getObjId(),
                                 self.angularSampling.get(), self.angularDistance.get(), self.symmetryGroup.get())
        self._insertFunctionStep('classifyGroupsStep')
        self._insertFunctionStep('createOutputStep')
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def convertInputStep(self, particlesId, volId):
        """ Write the input images as a Xmipp metadata file. 
        particlesId: is only need to detect changes in
        input particles and cause restart from here.
        """
        writeSetOfParticles(self.inputParticles.get(), self._getPath('input_particles.xmd'))
        from pyworkflow.em.convert import ImageHandler
        img = ImageHandler()
        img.convert(self.inputVolume.get(), self._getTmpPath("volume.vol"))
   
    def constructGroupsStep(self, particlesId, angularSampling, angularDistance, symmetryGroup):
        # Generate projections from this reconstruction        
        params = {"inputVol" : self._getTmpPath("volume.vol"),
                  "galleryStk" : self._getExtraPath("gallery.stk"),
                  "galleryXmd" : self._getExtraPath("gallery.doc"),
                  "neighborhoods": self._getExtraPath("neighbours.xmd"),
                  "symmetry" : self.symmetryGroup.get(),
                  "angularSampling" : self.angularSampling.get(),
                  "angularDistance" : self.angularDistance.get(),
                  "expParticles" : self._getPath('input_particles.xmd')
                }
        args = '-i %(inputVol)s -o %(galleryStk)s --sampling_rate %(angularSampling)f --sym %(symmetry)s'
        args += ' --method fourier 1 0.25 bspline --compute_neighbors --angular_distance %(angularSampling)f'
        args += ' --experimental_images %(expParticles)s --max_tilt_angle 90'
        
        self.runJob("xmipp_angular_project_library", args % params)
        
        args = '--i1 %(expParticles)s --i2 %(galleryXmd)s -o %(neighborhoods)s --dist %(angularDistance)f --sym %(symmetry)s --check_mirrors'
        self.runJob("xmipp_angular_neighbourhood", args % params, numberOfMpi=1)
   
    def classifyGroupsStep(self):
        from pyworkflow.em.metadata.utils import getSize
        mdOut = xmipp.MetaData()

        fnNeighbours = self._getExtraPath("neighbours.xmd")
        fnGallery=self._getExtraPath("gallery.stk")
        for block in xmipp.getBlocksInMetaDataFile(fnNeighbours):
            imgNo = block.split("_")[1]
            fnDir = self._getExtraPath("direction_%s"%imgNo)
            makePath(fnDir)
            fnOut = join(fnDir,"level_00/class_classes.stk")
            fnRef = "%s@%s"%(imgNo,fnGallery)
            fnBlock = "%s@%s"%(block,fnNeighbours)
            blockSize = getSize(fnBlock)
            if blockSize!=0:
                if not exists(fnOut):
                    args="-i %s --odir %s --ref0 %s --iter 1 --nref 1 --distance correlation --classicalMultiref --maxShift %d"%\
                        (fnBlock,fnDir,fnRef,self.maxShift.get())
                    self.runJob("xmipp_classify_CL2D", args)
                    fnAlignRoot = join(fnDir,"classes")
                    self.runJob("xmipp_image_align","-i %s --ref %s --oroot %s --iter 1"%(fnOut,fnRef,fnAlignRoot),numberOfMpi=1)
                    self.runJob("xmipp_transform_geometry","-i %s_alignment.xmd --apply_transform"%fnAlignRoot,numberOfMpi=1)
    
    #           # Construct output metadata
                objId = mdOut.addObject()
                mdOut.setValue(xmipp.MDL_REF,int(imgNo),objId)
                mdOut.setValue(xmipp.MDL_IMAGE,"%s"%fnRef,objId)
                mdOut.setValue(xmipp.MDL_IMAGE1,"1@%s"%fnOut,objId)
                mdOut.setValue(xmipp.MDL_CLASS_COUNT,getSize("class000001_images@%s"%join(fnDir,"level_00/class_classes.xmd")),objId)
        fnDirectional=self._getPath("directionalClasses.xmd")
        mdOut.write("classes@"+fnDirectional)
        self.runJob("xmipp_metadata_utilities",'-i %s --operate modify_values "ref=ref+1"'%self._getExtraPath("gallery.doc"), numberOfMpi=1)
        self.runJob("xmipp_metadata_utilities","-i classes@%s --set join %s ref"%(fnDirectional,self._getExtraPath("gallery.doc")), numberOfMpi=1)
        self.runJob("xmipp_metadata_utilities",'-i classes@%s --operate keep_column "ref image image1 classCount angleRot angleTilt"'%fnDirectional, numberOfMpi=1)

        for block in xmipp.getBlocksInMetaDataFile(fnNeighbours):
            imgNo = block.split("_")[1]
            fnDir = self._getExtraPath("direction_%s"%imgNo)
            try:
                mdClass=xmipp.MetaData("class000001_images@%s"%join(fnDir,"level_00/class_classes.xmd"))
                mdClass.write("class%s_images@%s"%(imgNo,fnDirectional),xmipp.MD_APPEND)
            except Exception as e:
                pass

    def createOutputStep(self):
        fnDirectional=self._getPath("directionalClasses.xmd")
        if exists(fnDirectional):
            inputParticles = self.inputParticles.get()
            self._sampling = inputParticles.getSamplingRate()
            classes2DSet = self._createSetOfClasses2D(inputParticles)
            readSetOfClasses2D(classes2DSet, fnDirectional, 'classes', preprocessClass=self._readRow)
            self._defineOutputs(outputClasses=classes2DSet)
            self._defineSourceRelation(self.inputParticles, classes2DSet)
    
    def _readRow(self, classItem, row):
        img = Image()
        img.setSamplingRate(self._sampling)
        img.setLocation(xmippToLocation(row.getValue(xmipp.MDL_IMAGE1)))
        classItem.reprojection = img
        
        setXmippAttributes(classItem, row, xmipp.MDL_ANGLE_ROT)
        setXmippAttributes(classItem, row, xmipp.MDL_ANGLE_TILT)

    #--------------------------- INFO functions -------------------------------------------- 
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

    
    