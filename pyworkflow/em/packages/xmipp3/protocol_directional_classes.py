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
from pyworkflow.object import Float, String
from pyworkflow.protocol.params import (PointerParam, FloatParam, IntParam,
                                        StringParam, LEVEL_ADVANCED, BooleanParam)
from pyworkflow.em.data import Volume
from pyworkflow.em.protocol import ProtAnalysis3D
from pyworkflow.utils.path import moveFile, makePath, cleanPath, cleanPattern
from pyworkflow.em.packages.xmipp3.convert import writeSetOfParticles, readSetOfParticles
from pyworkflow.em.metadata.utils import getSize
import xmipp
import math

class XmippProtDirectionalClasses(ProtAnalysis3D):
    """    
    Analyze 2D classes as assigned to the different directions
    """

    _label = 'directional_classes'
    
    def __init__(self, *args, **kwargs):
        ProtAnalysis3D.__init__(self, *args, **kwargs)
        
    #--------------------------- DEFINE param functions --------------------------------------------   
    def _defineParams(self, form):
        form.addSection(label='Input')
        
        form.addParam('inputVolume', PointerParam, pointerClass='Volume',
                      label="Input volume",  
                      help='Select the input volume.')     
        form.addParam('inputParticles', PointerParam, pointerClass='SetOfParticles', 
                      label="Input particles",   pointerCondition='hasAlignmentProj',
                      help='Select the input projection images with an angular assignment.') 
        form.addParam('symmetryGroup', StringParam, default='c1',
                      label="Symmetry group", 
                      help='See [[Xmipp Symmetry][http://www2.mrc-lmb.cam.ac.uk/Xmipp/index.php/Conventions_%26_File_formats#Symmetry]] page '
                           'for a description of the symmetry format accepted by Xmipp') 
        form.addParam('targetResolution', FloatParam, default=10, label='Target resolution (A)', expertLevel=LEVEL_ADVANCED)
        form.addParam('angularSampling', FloatParam, default=5, label='Angular sampling', expertLevel=LEVEL_ADVANCED, help="In degrees")
        form.addParam('angularDistance', FloatParam, default=10, label='Angular distance', expertLevel=LEVEL_ADVANCED,
                      help="In degrees. An image belongs to a group if its distance is smaller than this value")
        form.addParam('maxShift', FloatParam, default=15, label='Maximum shift', expertLevel=LEVEL_ADVANCED,
                      help="In pixels")
        form.addParam('refineAngles', BooleanParam, default=True, label='Refine angles', expertLevel=LEVEL_ADVANCED,
                      help="Refine the angles of the classes using a continuous angular assignment")
        form.addParam('directionalClasses', IntParam, default=2, label='Number of directional classes', expertLevel=LEVEL_ADVANCED)
        form.addParam('cl2dIterations', IntParam, default=5, label='Number of CL2D iterations', expertLevel=LEVEL_ADVANCED)
        
        form.addParallelSection(threads=0, mpi=8)
    
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):        
        self._insertFunctionStep('convertInputStep', self.inputParticles.get().getObjId(), self.inputVolume.get().getObjId(), 
                                 self.targetResolution.get())
        self._insertFunctionStep('constructGroupsStep', self.inputParticles.get().getObjId(),
                                 self.angularSampling.get(), self.angularDistance.get(), self.symmetryGroup.get())
        self._insertFunctionStep('classifyGroupsStep')
        if self.refineAngles:
            self._insertFunctionStep('refineAnglesStep')
        self._insertFunctionStep('cleanStep')
        self._insertFunctionStep('createOutputStep',1)
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def convertInputStep(self, particlesId, volId, targetResolution):
        """ Write the input images as a Xmipp metadata file. 
        particlesId: is only need to detect changes in
        input particles and cause restart from here.
        """
        writeSetOfParticles(self.inputParticles.get(), self._getPath('input_particles.xmd'))
        Xdim = self.inputParticles.get().getDimensions()[0]
        Ts = self.inputParticles.get().getSamplingRate()
        newTs = self.targetResolution.get()*0.4
        newTs = max(Ts,newTs)
        newXdim = Xdim*Ts/newTs
        self.runJob("xmipp_image_resize","-i %s -o %s --save_metadata_stack %s --dim %d"%\
                    (self._getPath('input_particles.xmd'),
                    self._getExtraPath('scaled_particles.stk'),
                    self._getExtraPath('scaled_particles.xmd'),
                    newXdim))
    
        from pyworkflow.em.convert import ImageHandler
        img = ImageHandler()
        img.convert(self.inputVolume.get(), self._getExtraPath("volume.vol"))
        Xdim = self.inputVolume.get().getDim()[0]
        if Xdim!=newXdim:
            self.runJob("xmipp_image_resize","-i %s --dim %d"%\
                        (self._getExtraPath("volume.vol"),
                        newXdim), numberOfMpi=1)
   
    def constructGroupsStep(self, particlesId, angularSampling, angularDistance, symmetryGroup):
        # Generate projections from this reconstruction        
        params = {"inputVol" : self._getExtraPath("volume.vol"),
                  "galleryStk" : self._getExtraPath("gallery.stk"),
                  "galleryXmd" : self._getExtraPath("gallery.doc"),
                  "neighborhoods": self._getExtraPath("neighbours.xmd"),
                  "symmetry" : self.symmetryGroup.get(),
                  "angularSampling" : self.angularSampling.get(),
                  "angularDistance" : self.angularDistance.get(),
                  "expParticles" : self._getExtraPath('scaled_particles.xmd')
                }
        args = '-i %(inputVol)s -o %(galleryStk)s --sampling_rate %(angularSampling)f --sym %(symmetry)s'
        args += ' --method fourier 1 0.25 bspline --compute_neighbors --angular_distance %(angularSampling)f'
        args += ' --experimental_images %(expParticles)s --max_tilt_angle 90'
        
        self.runJob("xmipp_angular_project_library", args % params)
        
        args = '--i1 %(expParticles)s --i2 %(galleryXmd)s -o %(neighborhoods)s --dist %(angularDistance)f --sym %(symmetry)s --check_mirrors'
        self.runJob("xmipp_angular_neighbourhood", args % params, numberOfMpi=1)               
   
    def classifyGroupsStep(self):
        mdOut = xmipp.MetaData()

        fnNeighbours = self._getExtraPath("neighbours.xmd")
        fnGallery=self._getExtraPath("gallery.stk")
        for block in xmipp.getBlocksInMetaDataFile(fnNeighbours):
            imgNo = block.split("_")[1]
            fnDir = self._getExtraPath("direction_%s"%imgNo)
            makePath(fnDir)
            Nlevels = int(math.ceil(math.log(self.directionalClasses.get())/math.log(2)))
            fnOut = join(fnDir,"level_%02d/class_classes.stk"%Nlevels)
            if not exists(fnOut):
                fnBlock="%s@%s"%(block,fnNeighbours)
                if getSize(fnBlock)>25:
                    try:
                        args="-i %s --odir %s --ref0 %s@%s --iter %d --nref %d --distance correlation --classicalMultiref --maxShift %d"%\
                            (fnBlock,fnDir,imgNo,fnGallery,self.cl2dIterations.get(),self.directionalClasses.get(),self.maxShift.get())
                        self.runJob("xmipp_classify_CL2D", args)
                        fnAlignRoot = join(fnDir,"classes")
                        self.runJob("xmipp_image_align","-i %s --ref %s@%s --oroot %s --iter 1"%(fnOut,imgNo,fnGallery,fnAlignRoot),numberOfMpi=1)
                        self.runJob("xmipp_transform_geometry","-i %s_alignment.xmd --apply_transform"%fnAlignRoot,numberOfMpi=1)
                    
                        # Construct output metadata
                        if exists(fnOut):
                            for i in range(self.directionalClasses.get()):
                                objId = mdOut.addObject()
                                mdOut.setValue(xmipp.MDL_REF,int(imgNo)-1,objId)
                                mdOut.setValue(xmipp.MDL_IMAGE,"%d@%s"%(i+1,fnOut),objId)
                    except:
                        print("The classification failed, probably because of a low number of images.")
                        print("However, this classification does not hinder the protocol to continue")
                        
        fnDirectional=self._getPath("directionalClasses.xmd")
        mdOut.write(fnDirectional)
        self.runJob("xmipp_metadata_utilities","-i %s --set join %s ref"%(fnDirectional,self._getExtraPath("gallery.doc")), numberOfMpi=1)
        
    def refineAnglesStep(self):
        fnDirectional = self._getPath("directionalClasses.xmd")
        newTs = self.targetResolution.get()*0.4
        self.runJob("xmipp_angular_continuous_assign2","-i %s --ref %s --max_resolution %f --sampling %f --optimizeAngles --optimizeShift"%\
                    (fnDirectional,self._getExtraPath("volume.vol"),self.targetResolution.get(),newTs))
    
    def cleanStep(self):
        cleanPath(self._getExtraPath('scaled_particles.stk'))
        cleanPath(self._getExtraPath('scaled_particles.xmd'))
        cleanPath(self._getExtraPath('volume.vol'))
        cleanPattern(self._getExtraPath("direction_*/level_00"))

    def createOutputStep(self, numeroFeo):
        fnDirectional=self._getPath("directionalClasses.xmd")
        if exists(fnDirectional):
            imgSetOut = self._createSetOfParticles()
            imgSetOut.setSamplingRate(self.targetResolution.get()*0.4)
            imgSetOut.setAlignmentProj()
            readSetOfParticles(fnDirectional,imgSetOut)
            self._defineOutputs(outputParticles=imgSetOut)
            self._defineSourceRelation(self.inputParticles, imgSetOut)
            self._defineSourceRelation(self.inputVolume, imgSetOut)

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
    
    #--------------------------- UTILS functions --------------------------------------------

    
    