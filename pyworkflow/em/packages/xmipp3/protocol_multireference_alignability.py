# **************************************************************************
# *
# * Authors:         Javier Vargas (jvargas@cnb.csic.es)
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
from pyworkflow.gui.plotter import Plotter


from pyworkflow.em.packages.xmipp3.convert import (writeSetOfParticles,
                                                   writeSetOfVolumes,
                                                   getImageLocation)


class XmippProtMultiRefAlignability(ProtAnalysis3D):
    """    
    Performs soft alignment validation of a set of particles confronting them
    against a given 3DEM map. This protocol produces particle alignment
    precision and accuracy parameters.
    """
    _label = 'multireference aligneability'
    
    def __init__(self, *args, **kwargs):
        ProtAnalysis3D.__init__(self, *args, **kwargs)
        
    #--------------------------- DEFINE param functions --------------------------------------------   
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputVolumes', PointerParam, pointerClass='Volume',
                      label="Input volume(s)",  
                      help='Select the input volume(s).')     
                
        form.addParam('inputParticles', PointerParam,
                      pointerClass='SetOfParticles', pointerCondition='hasAlignment',
                      label="Input particles", important=True,
                      help='Select the input projection images.')
            
        form.addParam('symmetryGroup', StringParam, default='c1',
                      label="Symmetry group", 
                      help='See [[Xmipp Symmetry][http://www2.mrc-lmb.cam.ac.uk/Xmipp/index.php/Conventions_%26_File_formats#Symmetry]] page '
                           'for a description of the symmetry format accepted by Xmipp') 

        form.addParam('isCTFCorrected', BooleanParam, default=False,
                      label="Has been the CTF corrected previously using Wiener filter?",  
                      help='Select true if the CTF has been previously corrected through Wiener filtering')
                
        form.addParam('angularSampling', FloatParam, default=5, expertLevel=LEVEL_ADVANCED,
                      label="Angular Sampling (degrees)",  
                      help='Angular distance (in degrees) between neighboring projection points ')

        form.addParam('numOrientations', FloatParam, default=6, expertLevel=LEVEL_ADVANCED,
                      label="Number of Orientations for particle",  
                      help='Parameter to define the number of most similar volume \n' 
                      '    projected images for each projection image')

        form.addParam('doNotUseWeights', BooleanParam, default=False, expertLevel=LEVEL_ADVANCED,
                      label="Do not use the weights",
                      help='Do not use the weights in the clustering calculation')
        
        form.addParam('pseudoSymmetryGroup', StringParam, default='', expertLevel=LEVEL_ADVANCED,
                      label="Pseudo symmetry group", 
                      help='Add only in case the map is close to a symmetry different and more restrict than the one reported in the parameter Symmetry group.'
                      'See [[Xmipp Symmetry][http://www2.mrc-lmb.cam.ac.uk/Xmipp/index.php/Conventions_%26_File_formats#Symmetry]] page '
                           'for a description of the symmetry format accepted by Xmipp')
        
        form.addParam('minTilt', FloatParam, default=0, expertLevel=LEVEL_ADVANCED,
                      label="Minimum allowed tilt angle",  
                      help='Tilts below this value will not be considered for the alignment')
        
        form.addParam('maxTilt', FloatParam, default=180, expertLevel=LEVEL_ADVANCED,
                      label="Maximum allowed tilt angle without mirror check",  
                      help='Tilts above this value will not be considered for the alignment without mirror check')

                
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

            if ( not (self.pseudoSymmetryGroup.get() == '') ):
                sym = self.pseudoSymmetryGroup.get()
                
            volStepId = self._insertFunctionStep('alignabilityStep', 
                                                 volName, volDir,
                                                 sym,
                                                 prerequisites=[sigStepId2])
            
            deps.append(volStepId)
          

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
        params += ' --useForValidation %0.3f' % (self.numOrientations.get()-1)
        params += ' --dontCheckMirrors'
        return params
        
    
    def _getCommonParamsRef(self):
        params =  '  -i %s' % self._getPath('reference_particles.xmd')        
        params += ' --sym %s' % self.symmetryGroup.get()
        params += ' --dontReconstruct'
        params += ' --useForValidation %0.3f' % (self.numOrientations.get()-1)
        params += ' --dontCheckMirrors'
        return params
    
    def phantomProject(self,volName):
        nproc = self.numberOfMpi.get()
        nT=self.numberOfThreads.get()         
        pathParticles = self._getPath('input_particles.xmd')
        Nx,Ny,Nz = self.inputParticles.get().getDim()
        R = -int(Nx/2)

        f = open(self._getExtraPath('params'),'w')          
        print self.isCTFCorrected.get()
        f.write("""# XMIPP_STAR_1 *
#
data_block1
_dimensions2D '%d %d'
_projAngleFile %s
_ctfPhaseFlipped %d
_ctfCorrected %d
_applyShift 0
_noisePixelLevel   '0 0'""" % (Nx, Ny, pathParticles, self.inputParticles.get().isPhaseFlipped(), self.isCTFCorrected.get()))
        f.close()
        param =  ' -i %s' % volName
        param += ' --params %s' % self._getExtraPath('params')
        param += ' -o %s' % self._getPath('reference_particles.xmd')
        param += ' --sampling_rate % 0.3f' % self.inputParticles.get().getSamplingRate()
                
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
        params = '-i %s -o %s --sampling_rate %f --sym %s --method fourier 1 0.25 bspline --compute_neighbors --angular_distance %f --experimental_images %s --max_tilt_angle %f --min_tilt_angle %f'\
                    %(volName,fnGallery, angularSampling, self.symmetryGroup.get(), -1, self._getPath('input_particles.xmd'), self.maxTilt.get(), self.minTilt.get())
        
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
        copyfile(volDir+'/angles_iter001_00.xmd', self._getTmpPath(anglesPath))
        
    def alignabilityStep(self, volName,volDir,sym):
        
        nproc = self.numberOfMpi.get()
        nT=self.numberOfThreads.get()
         
        makePath(volDir)  
        inputFile = self._getPath('input_particles.xmd') 
        inputFileRef = self._getPath('reference_particles.xmd')
        aFile = self._getTmpPath('exp_particles.xmd')
        aFileRef =self._getTmpPath('ref_particles.xmd')
        aFileGallery =(volDir+'/gallery.doc')

        params = '  -i %s'  % inputFile
        params += ' -i2 %s' % inputFileRef  
        params += '  --volume %s' % volName  
        params += '  --angles_file %s' % aFile
        params += '  --angles_file_ref %s' % aFileRef
        params += '  --gallery %s' % aFileGallery
        params += ' --odir %s' % volDir
        params += ' --sym %s' % sym
        
        if self.doNotUseWeights:
            params += ' --dontUseWeights'
            
        self.runJob('xmipp_multireference_aligneability', params,numberOfMpi=nproc,numberOfThreads=nT)

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
        
        outputVols = self._createSetOfVolumes()

        for i, vol in enumerate(self._iterInputVols()):        
        
            volDir = self._getVolDir(i+1)
            volume = vol.clone()
            volPrefix = 'vol%03d_' % (i+1)

            m_pruned = md.MetaData()
            m_pruned.read(volDir+'/pruned_particles_alignability.xmd')
            prunedMd = self._getExtraPath(volPrefix + 'pruned_particles_alignability.xmd')
            
            moveFile(join(volDir, 'pruned_particles_alignability.xmd'), prunedMd)
            m_volScore = md.MetaData()
            m_volScore.read(volDir+'/validationAlignability.xmd')
            validationMd = self._getExtraPath(volPrefix + 'validation_alignability.xmd')
            moveFile(join(volDir, 'validationAlignability.xmd'), validationMd)
            
            imgSet = self.inputParticles.get()                  

            outImgSet = self._createSetOfParticles(volPrefix)            
            outImgSet.copyInfo(imgSet)

            outImgSet.copyItems(imgSet,
                                updateItemCallback=self._setWeight,
                                itemDataIterator=md.iterRows(prunedMd, sortByLabel=md.MDL_ITEM_ID))
                        
            mdValidatoin = md.getFirstRow(validationMd)        
       
            weight = mdValidatoin.getValue(md.MDL_WEIGHT_PRECISION_ALIGNABILITY)        
            volume.weightAlignabilityPrecision  = Float(weight)
        
            weight = mdValidatoin.getValue(md.MDL_WEIGHT_ACCURACY_ALIGNABILITY)        
            volume.weightAlignabilityAccuracy  = Float(weight)
                    
            weight = mdValidatoin.getValue(md.MDL_WEIGHT_PRECISION_MIRROR)        
            volume.weightMirror  = Float(weight)
                    
            volume.cleanObjId() # clean objects id to assign new ones inside the set            
            outputVols.append(volume)
            self._defineOutputs(outputParticles=outImgSet)
            
            self.createPlot2D(volPrefix,m_pruned)
       
        outputVols.setSamplingRate(volume.getSamplingRate())
        self._defineOutputs(outputVolumes=outputVols)
        
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
        summary = ["Input particles:  %s" % self.inputParticles.get().getNameId()]
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
                weightAccuracy = mdVal.getValue(md.MDL_WEIGHT_ACCURACY_ALIGNABILITY, mdVal.firstObject())
                weightPrecision = mdVal.getValue(md.MDL_WEIGHT_PRECISION_ALIGNABILITY, mdVal.firstObject())
                weightAlignability = mdVal.getValue(md.MDL_WEIGHT_ALIGNABILITY, mdVal.firstObject())
                 
                summary.append("ALIGNABILITY ACCURACY parameter_%d : %f" % (i+1,weightAccuracy))
                summary.append("ALIGNABILITY PRECISION parameter_%d : %f" % (i+1,weightPrecision))
                summary.append("ALIGNABILITY ACCURACY & PRECISION parameter_%d : %f" % (i+1,weightAlignability))
                
                summary.append("-----------------")        
        return summary
    
    def _methods(self):
        messages = []
        if (hasattr(self,'outputVolumes')):
            messages.append('The quality parameter(s) has been obtained using '
                            'the approach [Vargas2014a] with angular sampling '
                            'of %f and number of orientations of %f' % (self.angularSampling.get(),
                                                                        self.numOrientations.get()))
        return messages
    
    def _citations(self):
        return ['Vargas2014a']
    
    #--------------------------- UTILS functions --------------------------------------------
    def _getVolDir(self, volIndex):
        return self._getTmpPath('vol%03d' % volIndex)
    
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
        item._xmipp_scoreAlignabilityPrecision    = Float(row.getValue(md.MDL_SCORE_BY_ALIGNABILITY_PRECISION))
        item._xmipp_scoreAlignabilityAccuracy = Float(row.getValue(md.MDL_SCORE_BY_ALIGNABILITY_ACCURACY))
        item._xmipp_scoreMirror = Float(row.getValue(md.MDL_SCORE_BY_MIRROR))
        item._xmipp_weight = Float( float(item._xmipp_scoreAlignabilityAccuracy)*float(item._xmipp_scoreAlignabilityPrecision))
        
    def createPlot2D(self,volPrefix,md):
        
        import xmipp
        
        figurePath = self._getExtraPath(volPrefix + 'softAlignmentValidation2D.png')
        figureSize = (8, 6)
    
        #alignedMovie = mic.alignMetaData
        plotter = Plotter(*figureSize)
        figure = plotter.getFigure()
    
        ax = figure.add_subplot(111)
        ax.grid()
        ax.set_title('Soft alignment validation map')
        ax.set_xlabel('Angular Precision')
        ax.set_ylabel('Angular Accuracy')

        for objId in md:
            x = md.getValue(xmipp.MDL_SCORE_BY_ALIGNABILITY_PRECISION, objId)
            y = md.getValue(xmipp.MDL_SCORE_BY_ALIGNABILITY_ACCURACY, objId)
            ax.plot(x, y, 'r.',markersize=1)

        ax.grid(True, which='both')
        ax.autoscale_view(True,True,True)
                
        plotter.savefig(figurePath)
        plotter.show()
        return plotter    

