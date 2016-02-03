# **************************************************************************
# *
# * Authors:         Jose Luis Vilas (jvilas@cnb.csic.es)
#                    Javier Vargas (jvargas@cnb.csic.es)
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

from os.path import join, isfile
from shutil import copyfile
from pyworkflow.object import Float, String
from pyworkflow.protocol.params import (PointerParam, FloatParam, STEPS_PARALLEL,
                                        StringParam, BooleanParam, LEVEL_ADVANCED)
from pyworkflow.em.data import Volume
from pyworkflow.em import Viewer
import pyworkflow.em.metadata as md
from pyworkflow.em.protocol import ProtAnalysis3D
from pyworkflow.utils.path import moveFile, makePath
from pyworkflow.em.packages.xmipp3.convert import (writeSetOfParticles,
                                                   writeSetOfVolumes,
                                                   getImageLocation)


class XmippProtMultiRefAlignability(ProtAnalysis3D):
    """    
    Reconstruct a volume using Xmipp_reconstruct_fourier from a given set of particles.
    The alignment parameters will be converted to a Xmipp xmd file
    and used as direction projections to reconstruct.
    """
    _label = 'multireference_aligneability'
    
    def __init__(self, *args, **kwargs):
        ProtAnalysis3D.__init__(self, *args, **kwargs)
        
    #--------------------------- DEFINE param functions --------------------------------------------   
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputVolumes', PointerParam, pointerClass='SetOfVolumes,Volume',
                      label="Input volume(s)",  
                      help='Select the input volume(s).')     
                
        form.addParam('inputParticles', PointerParam, label="Input particles", important=True, 
                      pointerClass='SetOfParticles', pointerCondition='hasAlignment',
                      help='Select the input projection images.')
            
        form.addParam('symmetryGroup', StringParam, default='c1',
                      label="Symmetry group", 
                      help='See [[Xmipp Symmetry][http://www2.mrc-lmb.cam.ac.uk/Xmipp/index.php/Conventions_%26_File_formats#Symmetry]] page '
                           'for a description of the symmetry format accepted by Xmipp') 
        
        form.addParam('angularSampling', FloatParam, default=5, expertLevel=LEVEL_ADVANCED,
                      label="Angular Sampling (degrees)",  
                      help='Angular distance (in degrees) between neighboring projection points ')

        form.addParam('numOrientations', FloatParam, default=10, expertLevel=LEVEL_ADVANCED,
                      label="Number of Orientations for particle",  
                      help='Parameter to define the number of most similar volume \n' 
                      '    projected images for each projection image')

        form.addParam('angDist', FloatParam, default=10, expertLevel=LEVEL_ADVANCED,
                      label="Maximum angular neighbourhood distance (degrees)",  
                      help='Parameter to define the maximum angular neighbourhood distance to determine close particles ')
           
        form.addParam('doNotUseWeights', BooleanParam, default=False, expertLevel=LEVEL_ADVANCED,
                      label="Do not use the weights",
                      help='Do not use the weights in the clustering calculation')
        
        form.addParallelSection(threads=1, mpi=1)

    #--------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):   
             
        convertId = self._insertFunctionStep('convertInputStep', 
                                             self.inputParticles.get().getObjId())
        deps = [] # store volumes steps id to use as dependencies for last step
        commonParams    = self._getCommonParams()
        commonParamsRef = self._getCommonParamsRef()

        sym = self.symmetryGroup.get()
        for i, vol in enumerate(self._iterInputVols()):
            
            volName = getImageLocation(vol)
            volDir = self._getVolDir(i+1)
            
            pmStepId = self._insertFunctionStep('projectionLibraryStep',                                                    
                                                volName, volDir,self.angularSampling.get(),
                                                prerequisites=[convertId])
            
            sigStepId1 = self._insertFunctionStep('significantStep', 
                                                 volName, volDir,
                                                 'exp_particles.xmd',
                                                 commonParams, 
                                                 prerequisites=[pmStepId])
            
            phanProjStepId = self._insertFunctionStep('phantomProject', 
                                                 volName, 
                                                 prerequisites=[sigStepId1])


            sigStepId2 = self._insertFunctionStep('significantStep', 
                                                 volName, volDir,
                                                 'ref_particles.xmd',
                                                 commonParamsRef, 
                                                 prerequisites=[phanProjStepId])


            volStepId = self._insertFunctionStep('angularPrecisionStep', 
                                                 volName, volDir,
                                                 sym,
                                                 prerequisites=[sigStepId2])
            
            
            pmStepId2 = self._insertFunctionStep('projectionLibraryStep',                                                    
                                                volName, volDir,10,
                                                prerequisites=[volStepId])

            
            pmNeigDirId = self._insertFunctionStep('neighbourhoodDirectionStep',                                                    
                                                volName, volDir,sym,
                                                prerequisites=[pmStepId2])
            

            pmAngAccId = self._insertFunctionStep('angularAccuracyStep',                                                    
                                                volName, volDir,i,
                                                prerequisites=[pmNeigDirId])
            
           
            
            deps.append(pmAngAccId)
          
        self._insertFunctionStep('createOutputStep', 
                                 prerequisites=deps)
        
    def convertInputStep(self, particlesId):
        """ Write the input images as a Xmipp metadata file. 
        particlesId: is only need to detect changes in
        input particles and cause restart from here.
        """
        a = self.inputParticles.get() 
        print a.printAll() 
        writeSetOfParticles(self.inputParticles.get(), 
                            self._getPath('input_particles.xmd'))
                    
    def _getCommonParams(self):
        params =  '  -i %s' % self._getPath('input_particles.xmd')        
        params += ' --sym %s' % self.symmetryGroup.get()
        params += ' --dontReconstruct'
        params += ' --useForValidation %0.3f' % (self.numOrientations.get())        
        return params
        
    
    def _getCommonParamsRef(self):
        params =  '  -i %s' % self._getPath('reference_particles.xmd')        
        params += ' --sym %s' % self.symmetryGroup.get()
        params += ' --dontReconstruct'
        params += ' --useForValidation %0.3f' % (self.numOrientations.get())       
        return params
    
    def phantomProject(self,volName):
        
        nproc = self.numberOfMpi.get()
        nT=self.numberOfThreads.get()         
        pathParticles = self._getPath('input_particles.xmd')
        Nx,Ny,Nz = self.inputParticles.get().getDim()
        R = -int(Nx/2)

        f = open(self._getExtraPath('params'),'w')          
        f.write("""# XMIPP_STAR_1 *
#
data_block1
_dimensions2D '%d %d'
_projAngleFile %s
_ctfPhaseFlipped %d
_applyShift 0
_noisePixelLevel   '0 0'""" % (Nx, Ny, pathParticles, self.inputParticles.get().isPhaseFlipped()))
        f.close()
        param =     ' -i %s' % volName
        param +=    ' --params %s' % self._getExtraPath('params')
        param +=    ' -o %s' % self._getPath('reference_particles.xmd')
        param +=    ' --sampling_rate % 0.3f' % self.inputParticles.get().getSamplingRate()
                
        #while (~isfile(self._getExtraPath('params'))):
        #    print 'No created'
        
        self.runJob('xmipp_phantom_project', 
                    param, numberOfMpi=1,numberOfThreads=1)

        param =     ' -i %s' % self._getPath('reference_particles.stk')
        param +=    ' --mask circular %d' % R
        self.runJob('xmipp_transform_mask',param, numberOfMpi=nproc,numberOfThreads=nT)
                
 
    def projectionLibraryStep(self, volName, volDir, angularSampling):
        
        # Generate projections from this reconstruction        
        nproc = self.numberOfMpi.get()
        nT=self.numberOfThreads.get() 
        
        makePath(volDir)
        fnGallery= (volDir+'/gallery.stk')
        params = '-i %s -o %s --sampling_rate %f --sym %s --method fourier 1 0.25 bspline --compute_neighbors --angular_distance %f --experimental_images %s --max_tilt_angle 90'\
                    %(volName,fnGallery, angularSampling, self.symmetryGroup.get(), -1, self._getPath('input_particles.xmd'))
        
        self.runJob("xmipp_angular_project_library", params, numberOfMpi=nproc, numberOfThreads=nT)                    
        
    def significantStep(self, volName, volDir, anglesPath, params):

        nproc = self.numberOfMpi.get()
        nT=self.numberOfThreads.get() 

        fnGallery= (volDir+'/gallery.doc')          
        params += ' --initgallery  %s' % fnGallery
        params += ' --odir %s' % volDir
        params += ' --iter %d' % 1
        self.runJob('xmipp_reconstruct_significant', 
                    params, numberOfMpi=nproc,numberOfThreads=nT)
        copyfile(volDir+'/angles_iter001_00.xmd', self._getExtraPath(anglesPath))
        
    def angularPrecisionStep(self, volName,volDir,sym):
        makePath(volDir)  
        aFile = self._getExtraPath('exp_particles.xmd')
        aFileRef =self._getExtraPath('ref_particles.xmd')
        aFileGallery =(volDir+'/gallery.doc')
        params = '  --volume %s' % volName  
        params += '  --angles_file %s' % aFile
        params += '  --angles_file_ref %s' % aFileRef
        params += '  --gallery %s' % aFileGallery
        params += ' --odir %s' % volDir
        params += ' --sym %s' % sym
        
        if (self.doNotUseWeights.get()) :
            params += ' --dontUseWeights'
            
        self.runJob('xmipp_multireference_aligneability', params,numberOfMpi=1,numberOfThreads=1)


    def neighbourhoodDirectionStep(self, volName,volDir,sym):
          
        aFileGallery =(volDir+'/gallery.doc')
        neighbours = (volDir+'/neighbours.xmd')
        
        params = '  --i1 %s' % self._getPath('input_particles.xmd')   
        params += ' --i2 %s' % aFileGallery
        params += ' -o %s' % neighbours
        params += ' --dist %s' % (self.angDist.get()+1)
        params += ' --sym %s' % sym        
                    
        self.runJob('xmipp_angular_neighbourhood', params,numberOfMpi=1,numberOfThreads=1)
        
        
    def angularAccuracyStep(self, volName,volDir,indx):
         
        nproc = self.numberOfMpi.get()
        nT=self.numberOfThreads.get()
          
        neighbours = (volDir+'/neighbours.xmd')
        
        params =  ' -i %s' %  volName
        params += ' --i2 %s' % neighbours
        params += ' -o %s' %  (volDir+'/pruned_particles_alignability_accuracy.xmd')
                    
        self.runJob('xmipp_angular_accuracy_pca', params,numberOfMpi=nproc,numberOfThreads=nT)
        
    def createOutputStep(self):
        
        for i, vol in enumerate(self._iterInputVols()):        
        
            volDir = self._getVolDir(i+1)
            volume = vol.clone()
            volPrefix = 'vol%03d_' % (i+1)
            
            m1_pruned = md.MetaData()
            m2_pruned = md.MetaData()
            mdJoin_pruned = md.MetaData()            
            m1_pruned.read(volDir+'/pruned_particles_alignability_precision.xmd')
            m2_pruned.read(volDir+'/pruned_particles_alignability_accuracy.xmd')                
            mdJoin_pruned = m1_pruned
            mdJoin_pruned.merge(m2_pruned)
            mdJoin_pruned.write(volDir+'/pruned_particles_alignability.xmd')
            
            prunedMd = self._getExtraPath(volPrefix + 'pruned_particles_alignability.xmd')
            moveFile(join(volDir, 'pruned_particles_alignability.xmd'), 
                     prunedMd)

            m1_volScore = md.MetaData()
            m2_volScore = md.MetaData()
            mdJoin_volScore = md.MetaData()
	                
            m1_volScore.read(volDir+'/validationAlignabilityPrecision.xmd')
            m2_volScore.read(volDir+'/validationAlignabilityAccuracy.xmd') 
	    
            mdJoin_volScore = m1_volScore           
            mdJoin_volScore.merge(m2_volScore)
            mdJoin_volScore.write(volDir+'/validation_alignability.xmd')
            
            validationMd = self._getExtraPath(volPrefix + 'validation_alignability.xmd')
            moveFile(join(volDir, 'validation_alignability.xmd'), validationMd)
                           
            outputVols = self._createSetOfVolumes()
            imgSet = self.inputParticles.get()                  

            outImgSet = self._createSetOfParticles(volPrefix)            
            outImgSet.copyInfo(imgSet)

            outImgSet.copyItems(imgSet,
                                updateItemCallback=self._setWeight,
                                itemDataIterator=md.iterRows(mdJoin_pruned, sortByLabel=md.MDL_ITEM_ID))
                        
            mdValidatoin = md.MetaData(validationMd)
	    
	   
            weight = mdValidatoin.getValue(md.MDL_WEIGHT_PRECISION_ALIGNABILITY, mdValidatoin.firstObject())	    
            volume.weightAlignabilityPrecision  = Float(weight)
	    
            weight = mdValidatoin.getValue(md.MDL_WEIGHT_PRECISION_MIRROR, mdValidatoin.firstObject())	    
            volume.weightMirror  = Float(weight)
	    
    	    weight = mdValidatoin.getValue(md.MDL_SCORE_BY_PCA_RESIDUAL_PROJ, mdValidatoin.firstObject())	    
            volume.AlignabilityPCAProj  = Float(weight)

    	    weight = mdValidatoin.getValue(md.MDL_SCORE_BY_PCA_RESIDUAL_EXP, mdValidatoin.firstObject())	    
            volume.AlignabilityPCAExp  = Float(weight)

            weight = mdValidatoin.getValue(md.MDL_SCORE_BY_PCA_RESIDUAL, mdValidatoin.firstObject())        
            volume.AlignabilityPCA  = Float(weight)
            
    	    weight = mdValidatoin.getValue(md.MDL_SCORE_BY_ZSCORE, mdValidatoin.firstObject())	    
            volume.ZScore  = Float(weight)
            	    
            volume.clusterMd = String(mdJoin_pruned)
            volume.cleanObjId() # clean objects id to assign new ones inside the set            
            outputVols.append(volume)
            self._defineOutputs(outputParticles=outImgSet)
        
        outputVols.setSamplingRate(volume.getSamplingRate())
        self._defineOutputs(outputVolumes=outputVols)
        #self._defineTransformRelation(self.inputVolumes.get(), volume)
        
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

        summary.append("Input particles:  %s" % self.inputParticles.get().getNameId())

        summary.append("-----------------")
        if self.inputVolumes.get():
            for i, vol in enumerate(self._iterInputVols()):
                summary.append("Input volume(s)_%d: [%s]" % (i+1,vol))

        summary.append("-----------------")
        if  (not hasattr(self,'outputVolumes')):
            summary.append("Output volumes not ready yet.")
        else:
            for i, vol in enumerate(self._iterInputVols()):
                            
                VolPrefix = 'vol%03d_' % (i+1)
                mdVal = md.MetaData(self._getExtraPath(VolPrefix+'validation_alignability.xmd'))                
                weight = mdVal.getValue(md.MDL_WEIGHT, mdVal.firstObject())
                summary.append("Output volume(s)_%d : %s" % (i+1,self.outputVolumes.getNameId()))
                summary.append("Quality parameter_%d : %f" % (i+1,weight))
                summary.append("-----------------")        
        return summary
    
    def _methods(self):
        messages = []
        if (hasattr(self,'outputVolumes')):
            messages.append('The quality parameter(s) has been obtained using the approach [Vargas2014a] with angular sampling of %f and number of orientations of %f' % (self.angularSampling.get(), self.numOrientations.get()))
        return messages
    
    def _citations(self):
        return ['Vargas2014a']
    
    #--------------------------- UTILS functions --------------------------------------------
    def _getVolDir(self, volIndex):
        return self._getExtraPath('vol%03d' % volIndex)
    
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
        fscFn = self._defineMetadataRootName('Volume',volId)
        return fscFn
    
    def _setWeight(self, item, row):  
        item._xmipp_scoreAlignabilityPrecision    = Float(row.getValue(md.MDL_SCORE_BY_ALIGNABILITY))
        item._xmipp_scoreAlignabilityAccuracyProj = Float(row.getValue(md.MDL_SCORE_BY_PCA_RESIDUAL_PROJ))
        item._xmipp_scoreAlignabilityAccuracyExp = Float(row.getValue(md.MDL_SCORE_BY_PCA_RESIDUAL_EXP))
        item._xmipp_scoreAlignabilityAccuracy = Float(row.getValue(md.MDL_SCORE_BY_PCA_RESIDUAL))
        item._xmipp_scoreZscore = Float(row.getValue(md.MDL_SCORE_BY_ZSCORE))
        item._xmipp_scoreMirror = Float(row.getValue(md.MDL_SCORE_BY_MIRROR))	
        item._xmipp_weight = Float( float(item._xmipp_scoreZscore)*float(item._xmipp_scoreAlignabilityAccuracy)*float(item._xmipp_scoreAlignabilityPrecision))
