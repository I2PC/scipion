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

from os.path import join
from pyworkflow.object import Float, String
from pyworkflow.protocol.params import (PointerParam, FloatParam, STEPS_PARALLEL,
                                        StringParam, EnumParam, LEVEL_ADVANCED)
from pyworkflow.em.data import Volume
from pyworkflow.em.protocol import ProtAnalysis3D
from pyworkflow.utils.path import moveFile, makePath
from pyworkflow.em.packages.xmipp3.convert import writeSetOfParticles
import pyworkflow.em.metadata as md


PROJECTION_MATCHING = 0
SIGNIFICANT = 1


class XmippProtValidateNonTilt(ProtAnalysis3D):
    """    
    Ranks a set of volumes according to their alignment reliability obtained from a clusterability test.
    """

    _label = 'validate_nontilt'
    WEB = 0

    
    def __init__(self, *args, **kwargs):
        ProtAnalysis3D.__init__(self, *args, **kwargs)
        
        if (self.WEB == 1):
            self.stepsExecutionMode = STEPS_PARALLEL
        
    #--------------------------- DEFINE param functions --------------------------------------------   
    def _defineParams(self, form):
        form.addSection(label='Input')
        
        form.addParam('inputVolumes', PointerParam, pointerClass='SetOfVolumes, Volume',
                      label="Input volumes",  
                      help='Select the input volumes.')     
        form.addParam('inputParticles', PointerParam, pointerClass='SetOfParticles, SetOfClasses2D', 
                      label="Input particles",  
                      help='Select the input projection images .') 
        form.addParam('symmetryGroup', StringParam, default='c1',
                      label="Symmetry group", 
                      help='See [[Xmipp Symmetry][http://www2.mrc-lmb.cam.ac.uk/Xmipp/index.php/Conventions_%26_File_formats#Symmetry]] page '
                           'for a description of the symmetry format accepted by Xmipp') 
        form.addParam('alignmentMethod', EnumParam, label='Image alignment', choices=['Projection_Matching','Significant'], default=SIGNIFICANT)
        line = form.addLine('Resolution to filter (A)')
        line.addParam('highPassFilter', FloatParam, default=150, label='High')
        line.addParam('lowPassFilter', FloatParam, default=15, label='Low')
        
        form.addParam('angularSampling', FloatParam, default=5, expertLevel=LEVEL_ADVANCED,
                      label="Angular Sampling (degrees)",
                      help='Angular distance (in degrees) between neighboring projection points ')
        form.addParam('numOrientations', FloatParam, default=10, expertLevel=LEVEL_ADVANCED,
                      label="Number of orientations per particle",  
                      help='Number of possible orientations in which a particle can be \n')
        form.addParam('significanceNoise', FloatParam, default=0.95, expertLevel=LEVEL_ADVANCED,
                      label="Significance", help='Significance of the aligniability with respect'
                                                 ' to a a set of uniformly distributed random points \n')
        
        form.addParallelSection(threads=0, mpi=4)
    
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):        
        deps = [] # store volumes steps id to use as dependencies for last step
        self.partSet = self.inputParticles.get()
        
        convertId = self._insertFunctionStep('convertInputStep', self.partSet.getObjId())
        
        for vol in self._iterInputVols():
            filterId=self._insertFunctionStep('filterVolumeStep', vol.getObjId(), vol.getFileName(), prerequisites=[convertId])
            pmStepId = self._insertFunctionStep('projectionLibraryStep', vol.getObjId(), prerequisites=[filterId])
            
            if self.alignmentMethod == SIGNIFICANT:
                sigStepId = self._insertFunctionStep('significantStep', vol.getObjId(), prerequisites=[pmStepId])
            else:
                sigStepId = self._insertFunctionStep('projectionMatchingStep', vol.getObjId(), prerequisites=[pmStepId])
            
            volStepId = self._insertFunctionStep('validationStep', vol.getObjId(), prerequisites=[sigStepId])
            deps.append(volStepId)
        
        self._insertFunctionStep('createOutputStep', prerequisites=deps)
    
    #--------------------------- STEPS functions ---------------------------------------------------
    def convertInputStep(self, particlesId):
        """ Write the input images as a Xmipp metadata file. 
        particlesId: is only need to detect changes in
        input particles and cause restart from here.
        """
        writeSetOfParticles(self.partSet, 
                            self._getMdParticles())
    
    def filterVolumeStep(self, volId, volFn):
        params = {"inputVol" : volFn,
                  "filtVol" : self._getVolFiltered(volId),
                  "highPass" : self.partSet.getSamplingRate()/self.highPassFilter.get(),
                  "lowPass" : self.partSet.getSamplingRate()/self.lowPassFilter.get()
                  }
        
        args = ' -i %(inputVol)s -o %(filtVol)s --fourier band_pass %(highPass)f %(lowPass)f'
        self.runJob('xmipp_transform_filter', args % params, numberOfMpi=1, numberOfThreads=1)
    
    def projectionLibraryStep(self, volId):
        # Generate projections from this reconstruction        
        volDir = self._getVolDir(volId)
        makePath(volDir)
        
        params = {"inputVol" : self._getVolFiltered(volId),
                  "gallery" : self._getGalleryStack(volId),
                  "sampling" : self.partSet.getSamplingRate(),
                  "symmetry" : self.symmetryGroup.get(),
                  "angSampling" : self.angularSampling.get(),
                  "expParticles" : self._getMdParticles()
                  }
        
        args = '-i %(inputVol)s -o %(gallery)s --sampling_rate %(angSampling)f --sym %(symmetry)s'
        args += ' --method fourier 1 0.25 bspline --compute_neighbors --angular_distance -1'
        args += ' --experimental_images %(expParticles)s --max_tilt_angle 90'
        
        self.runJob("xmipp_angular_project_library", args % params)                    

    def significantStep(self, volId):
        params = {"inputParts" : self._getMdParticles(),
                  "symmetry" : self.symmetryGroup.get(),
                  "angSampling" : self.angularSampling.get(),
                  "orientations" : self.numOrientations.get(),
                  "gallery" : self._getGalleryMd(volId),
                  "outDir" : self._getVolDir(volId)
                  }

        args = ' -i %(inputParts)s --sym %(symmetry)s --angularSampling %(angSampling)0.3f --dontReconstruct'
        args += ' --useForValidation %(orientations)0.3f --initgallery  %(gallery)s --odir %(outDir)s --iter 1'
        
        self.runJob('xmipp_reconstruct_significant', args % params)
    
    def projectionMatchingStep(self, volId):
        params = {"inputParts" : self._getMdParticles(),
                  "outerRadius" : self.partSet.getDimensions()[0]/2,
                  "shift" : self.partSet.getDimensions()[0]/10,
                  "search5D" : self.partSet.getDimensions()[0]/10,
                  "gallery" : self._getGalleryStack(volId),
                  "orientations" : self.numOrientations.get(),
                  "output" : self._getAnglesMd(volId)
                  }
        
        args = ' -i %(inputParts)s --Ri 0.0 --Ro %(outerRadius)0.3f --max_shift %(shift)0.3f --append'
        args += ' --search5d_shift %(search5D)0.3f --number_orientations %(orientations)0.3f -o %(output)s --ref %(gallery)s'
        
        self.runJob('xmipp_angular_projection_matching', args % params)
    
    def validationStep(self, volId):
        params = {"inputAngles" : self._getAnglesMd(volId),
                  "filtVol" : self._getVolFiltered(volId),
                  "symmetry" : self.symmetryGroup.get(),
                  "significance" : self.significanceNoise.get(),
                  "outDir" : self._getVolDir(volId)
                  }
        
        args = ' --i %(inputAngles)s --volume %(filtVol)s --odir %(outDir)s'
        args += ' --significance_noise %(significance)0.2f --sym %(symmetry)s'
        
        if (self.alignmentMethod == SIGNIFICANT):
            args += ' --useSignificant '
        
        self.runJob('xmipp_validation_nontilt', args % params)
    
    def createOutputStep(self):
        outputVols = self._createSetOfVolumes()
        
        for vol in self._iterInputVols():
            volume = vol.clone()
            volDir = self._getVolDir(vol.getObjId())
            volPrefix = 'vol%03d_' % (vol.getObjId())
            validationMd = self._getExtraPath(volPrefix + 'validation.xmd')
            moveFile(join(volDir, 'validation.xmd'), 
                     validationMd)
            clusterMd = self._getExtraPath(volPrefix + 'clusteringTendency.xmd')
            moveFile(join(volDir, 'clusteringTendency.xmd'), clusterMd)
            
            mData = md.MetaData(validationMd)
            weight = mData.getValue(md.MDL_WEIGHT, mData.firstObject())
            volume._xmipp_weight = Float(weight)
            volume.clusterMd = String(clusterMd)
            volume.cleanObjId() # clean objects id to assign new ones inside the set
            outputVols.append(volume)
        
        outputVols.setSamplingRate(self.partSet.getSamplingRate())
        self._defineOutputs(outputVolumes=outputVols)
        self._defineTransformRelation(self.inputVolumes, outputVols)
    
    #--------------------------- INFO functions -------------------------------------------- 
    def _validate(self):
        validateMsgs = []
        # if there are Volume references, it cannot be empty.
        if self.inputVolumes.get() and not self.inputVolumes.hasValue():
            validateMsgs.append('Please provide an input reference volume.')
        if self.inputParticles.get() and not self.inputParticles.hasValue():
            validateMsgs.append('Please provide input particles.')            
        return validateMsgs
    
    def _summary(self):
        summary = []

        if  (not hasattr(self,'outputVolumes')):
            summary.append("Output volumes not ready yet.")
        else:
            size = 0
            for i, vol in enumerate(self._iterInputVols()):
                size +=1
            summary.append("Volumes to validate: *%d* " % size)
            summary.append("Angular sampling: %s" % self.angularSampling.get())
            summary.append("Significance value: %s" % self.significanceNoise.get())

        return summary
    
    def _methods(self):
        messages = []
        if (hasattr(self,'outputVolumes')):
            messages.append('The quality parameter(s) has been obtained using the approach [Vargas2014a] with angular sampling of %f and significant value of %f' % (self.angularSampling.get(), self.alpha.get()))
        return messages
    
    def _citations(self):
        return ['Vargas2014a']
    
    #--------------------------- UTILS functions --------------------------------------------
    def _getVolDir(self, volIndex):
        return self._getExtraPath('vol%03d' % volIndex)
    
    def _getVolFiltered(self, volIndex):
        return self._getVolDir(volIndex) + "_filt.vol"
    
    def _iterInputVols(self):
        """ In this function we will encapsulate the logic
        to iterate through the input volumes.
        This give the flexibility of having Volumes, SetOfVolumes or 
        a combination of them as input and the protocol code
        remain the same.
        """
        inputVols = self.inputVolumes.get()
        
        if isinstance(inputVols, Volume):
            yield inputVols
        else:
            for vol in inputVols:
                yield vol
    
    def _defineMetadataRootName(self, mdrootname,volId):
        if mdrootname=='P':
            VolPrefix = 'vol%03d_' % (volId)
            return self._getExtraPath(VolPrefix+'clusteringTendency.xmd')
        if mdrootname=='Volume':

            VolPrefix = 'vol%03d_' % (volId)
            return self._getExtraPath(VolPrefix+'validation.xmd')
            
    def _definePName(self):
        fscFn = self._defineMetadataRootName('P')
        return fscFn
    
    def _defineVolumeName(self,volId):
        fscFn = self._defineMetadataRootName('Volume', volId)
        return fscFn
    
    def _getMdParticles(self):
        return self._getPath('input_particles.xmd')
    
    def _getGalleryStack(self, volIndex):
        return join(self._getVolDir(volIndex), 'gallery.stk')
    
    def _getGalleryMd(self, volIndex):
        return join(self._getVolDir(volIndex), 'gallery.doc')
    
    def _getAnglesMd(self, volIndex):
        return join(self._getVolDir(volIndex), 'angles_iter001_00.xmd')
    
    
