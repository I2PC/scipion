# **************************************************************************
# *
# * Authors:     Javier Vargas and Adrian Quintana (jvargas@cnb.csic.es aquintana@cnb.csic.es)
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
This sub-package contains wrapper around ML2D Xmipp program
"""
from pyworkflow.em import *  
from convert import readSetOfClasses2D, createXmippInputClasses2D, createXmippInputVolumes, readSetOfVolumes
from math import floor
from xmipp import MetaData, MD_APPEND, MDL_MAXCC, MDL_WEIGHT, MDL_IMAGE, \
    MDL_VOLUME_SCORE_SUM, MDL_VOLUME_SCORE_SUM_TH, MDL_VOLUME_SCORE_MEAN, MDL_VOLUME_SCORE_MIN
#    removeFilenamePrefix
from pyworkflow.utils.path import moveFile, cleanPath, copyFile, removeExt
from protlib_xmipp import getMdSize

class XmippProtRansac(ProtInitialVolume):
    """ Protocol to obtain a set of initial volumes. """
    _label = 'ransac'
    
    def __init__(self, **args):
        ProtInitialVolume.__init__(self, **args)
        
#        self.progId = "ransac"
#        self.oroot = ""
        
    def _defineParams(self, form):
        form.addSection(label='Input')
        # JUST for testing purpose        
#        pointerClassList='SetOfClasses2D,SetOfVolumes'
        
        form.addParam('inputClasses', PointerParam, label="Input classes", important=True, 
                      pointerClass='SetOfClasses2D', pointerCondition='hasAverages',
                      help='Select the input classes from the project.'
                           'It should be a SetOfClasses2D class')        
        form.addParam('symmetryGroup', StringParam, default="c1",
                      label='Symmetry group', 
                      help='See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Symmetry for a description of the symmetry groups format'
                        'If no symmetry is present, give c1')
        form.addParam('angularSampling', FloatParam, default=5, expertLevel=LEVEL_EXPERT,
                      label='Angular sampling rate',
                      help='In degrees.'
                      ' This sampling defines how fine the projection gallery from the volume is explored.')
        form.addParam('nRansac', IntParam, default="400", expertLevel=LEVEL_EXPERT,
                      label="Number of RANSAC iterations", 
                      help='Number of initial volumes to test by RANSAC')
        
        form.addParam('dimRed', BooleanParam, default=True, expertLevel=LEVEL_EXPERT,
                      label='Perform dimensionality reduction', 
                      help='The dimensionality reduction is performed using the Local Tangent Space'
                      'Alignment. See http://www.stat.missouri.edu/~ys873/research/LTSA11.pdf')
        form.addParam('numGrids', IntParam, default=3, condition='dimRed', expertLevel=LEVEL_EXPERT,
                      label='Number of grids per dimension',
                      help='Number of squares to sample the classes')
        form.addParam('numSamples', IntParam, default=8, condition='not dimRed', expertLevel=LEVEL_EXPERT,
                      label='Number of random samples',
                      help='Number of squares to sample the classes')
        
        
        form.addParam('corrThresh', FloatParam, default=0.77, expertLevel=LEVEL_EXPERT,
                      label='Inliers threshold',
                      help='Correlation value threshold to determine if an experimental projection is an inlier or outlier.')        
        form.addParam('numVolumes', IntParam, default=10, expertLevel=LEVEL_EXPERT,
                      label='Number of best volumes to refine',
                      help='Number of best volumes to refine using projection matching approach and the input classes')
        form.addParam('numIter', IntParam, default=10, expertLevel=LEVEL_EXPERT,
                      label='Number of iterations to refine the volumes',
                      help='Number of iterations to refine the best volumes using projection matching approach and the input classes')
        form.addParam('initialVolume', PointerParam, label="Initial volume",  
                      pointerClass='SetOfVolumes', allowNull=True,
                      help='You may provide a very rough initial volume as a way to constraint the angular search.'
                            'For instance, when reconstructing a fiber, you may provide a cylinder so that side views'
                            'are assigned to the correct tilt angle, although the rotational angle may be completely wrong')           
                
        form.addParam('maxFreq', IntParam, default=5,
                      label='Max frequency of the initial volume',
                      help=' Max frequency of the initial volume in Angstroms')
        
        form.addParam('useSA', BooleanParam, default=False,
                      label='Combine simulated annealing and RANSAC', 
                      help='This option produces better results at a higher computational cost')
        form.addParam('nIterRandom', IntParam, default=10, condition='useSA',
                      label='Number of simulated annealing iterations',
                      help='During the simulated annealing iterations, all those particles positively'
                        'contributing to the improvement of the volume are considered. In this way,' 
                        'the same image may participate several times from different projection '
                        'directions (but different weights) depending on whether it improves the'
                        'correlation with the volume or not')
        form.addParam('rejection', IntParam, default=50, condition='useSA',
                      label='Percentage of rejected particles',
                       help='At each iteration, the lowest correlated particles are'
                            'removed from the 3D reconstruction, although they may '
                            'participate in the next iteration')
        

        form.addParam('useAll', BooleanParam, default=False, expertLevel=LEVEL_EXPERT,
                      label='Use all images to refine', 
                      help=' When refining a RANSAC volume, use all images to refine it instead of only inliers')
        
        form.addParallelSection(mpi=2)
            
         
        
    def _defineSteps(self):
        # Convert input images if necessary
        self.imgsFn = createXmippInputClasses2D(self, self.inputClasses.get())
        
        """ Filter the """
        initialStepId = self.initialize()
        
        # Generate projection gallery from the initial volume
        if self.initialVolume.hasValue():
            self._insertFunctionStep("projectInitialVolume",self)
        
        deps = [] # Store all steps ids, final step createOutput depends on all of them    
        for n in range(self.nRansac.get()):
            # CTF estimation with Xmipp
            stepId = self._insertFunctionStep('ransacIteration', n,
                                    prerequisites=[initialStepId]) # Make estimation steps indepent between them
            deps.append(stepId)
        
        # Look for threshold, evaluate volumes and get the best
        if self.initialVolume.hasValue():
            self._insertFunctionStep("runJob",programname="rm", params=self._getTmpPath("gallery_InitialVolume*"), NumberOfMpi=1)
            
        self._insertFunctionStep("getCorrThresh",
                                 prerequisites=deps) # Make estimation steps indepent between them)
        self._insertFunctionStep("evaluateVolumes")
        bestVolumesStepId = self._insertFunctionStep("getBestVolumes")        
        
        deps = [] # Store all steps ids, final step createOutput depends on all of them
        # Refine the best volumes
        for n in range(self.numVolumes.get()):
            fnBase='proposedVolume%05d'%n
            fnRoot=self._getPath(fnBase)
                    
            # Simulated annealing
            self._insertFunctionStep('reconstruct',fnRoot,
                                               prerequisites=[bestVolumesStepId]) # Make estimation steps indepent between them
            if self.useSA.get():
                self._insertFunctionStep('simulatedAnnealing',fnRoot)
        
            for it in range(self.numIter.get()):    
                self._insertFunctionStep('reconstruct',fnRoot) 
                self._insertFunctionStep('projMatch',fnBase)
            
            stepId =  self._insertRunJobStep("xmipp_image_resize","-i %s.vol -o %s.vol --dim %d %d" 
                                          %(fnRoot,fnRoot,self.Xdim,self.Xdim))
            
            deps.append(stepId)
        
        # Score each of the final volumes
        self._insertFunctionStep("scoreFinalVolumes",
                                 prerequisites=deps) # Make estimation steps indepent between them
        
        self._insertFunctionStep('createOutput')
        
    def initialize(self):
        self.Xdim=self.inputClasses.get().getDimensions()[0]
        
        fnOutputReducedClass = self._getExtraPath("reducedClasses.xmd")
        fnOutputReducedClassNoExt = os.path.splitext(fnOutputReducedClass)[0]
    
        # Low pass filter and resize        
        maxFreq = self.maxFreq.get()
        ts = self.inputClasses.get().getAverages().getSamplingRate()
        K = 0.25*(maxFreq/ts)
        if K<1:
            K=1
        self.Xdim2 = self.Xdim/K
        if (self.Xdim2 < 32):
            self.Xdim2 = 32
            K = self.Xdim/self.Xdim2
            
        freq = ts/maxFreq
        ts = K*ts

        self._insertRunJobStep("xmipp_transform_filter","-i %s -o %s --fourier low_pass %f --oroot %s"
                                                %(self.imgsFn,fnOutputReducedClass,freq,fnOutputReducedClassNoExt))
        return self._insertRunJobStep("xmipp_image_resize","-i %s --fourier %d -o %s" %(fnOutputReducedClass,self.Xdim2,fnOutputReducedClassNoExt))
    
    def ransacIteration(self, n):
    
        fnOutputReducedClass = self._getExtraPath("reducedClasses.xmd")  
        fnBase="ransac%05d"%n
        fnRoot=self._getTmpPath(fnBase)
        
    
        if self.dimRed:
            # Get a random sample of images
            self.runJob(None,"xmipp_transform_dimred","-i %s --randomSample %s.xmd  %d -m LTSA "%(fnOutputReducedClass,fnRoot,self.numGrids.get()))
        else:        
            self.runJob(None,"xmipp_metadata_utilities","-i %s -o %s.xmd  --operate random_subset %d --mode overwrite "%(fnOutputReducedClass,fnRoot,self.numSamples.get()))
            self.runJob(None,"xmipp_metadata_utilities","-i %s.xmd --fill angleRot rand_uniform -180 180 "%(fnRoot))
            self.runJob(None,"xmipp_metadata_utilities","-i %s.xmd --fill angleTilt rand_uniform 0 180 "%(fnRoot))
            self.runJob(None,"xmipp_metadata_utilities","-i %s.xmd --fill anglePsi  rand_uniform 0 360 "%(fnRoot)) 
    
        # If there is an initial volume, assign angles        
        if self.initialVolume.hasValue():
            fnGallery=self._getTmpPath('gallery_InitialVolume.stk')
            self.runJob(None,"xmipp_angular_projection_matching", "-i %s.xmd -o %s.xmd --ref %s --Ri 0 --Ro %s --max_shift %s --append"\
                   %(fnRoot,fnRoot,fnGallery,str(self.Xdim/2),str(self.Xdim/20)))
    
        # Reconstruct with the small sample
        self.reconstruct(fnRoot)
        
        fnVol = fnRoot+'.vol'
        
        # Simulated annealing
        if self.useSA.get():
            smallIter=int(min(floor(self.nIterRandom.get()/5.0),0));
            self.runJob(None,"xmipp_volume_initial_simulated_annealing","-i %s --initial %s --oroot %s_sa --sym %s --randomIter %d --rejection %f --dontApplyPositive"
                      %(fnRoot+".xmd",fnVol,fnRoot,self.symmetryGroup.get(),smallIter,self.rejection.get()))
            moveFile(fnRoot+"_sa.vol", fnVol)
            cleanPath(fnRoot+"_sa.xmd")
    
        # Generate projections from this reconstruction
        fnGallery=self._getTmpPath('gallery_'+fnBase+'.stk')
        self.runJob(None,"xmipp_angular_project_library", "-i %s -o %s --sampling_rate %f --sym %s --method fourier 1 0.25 bspline --compute_neighbors --angular_distance -1 --experimental_images %s"\
                    %(fnVol,fnGallery,self.angularSampling.get(),self.symmetryGroup.get(),fnOutputReducedClass))
            
        # Assign angles to the rest of images
        fnAngles=self._getTmpPath('angles_'+fnBase+'.xmd')
        self.runJob(None,"xmipp_angular_projection_matching", "-i %s -o %s --ref %s --Ri 0 --Ro %s --max_shift %s --append"\
                              %(fnOutputReducedClass,fnAngles,fnGallery,str(self.Xdim/2),str(self.Xdim/20)))
       
        # Delete intermediate files 
        cleanPath(fnGallery)
        cleanPath(self._getTmpPath('gallery_'+fnBase+'_sampling.xmd'))
        cleanPath(self._getTmpPath('gallery_'+fnBase+'.doc'))
        cleanPath(fnVol)
        cleanPath(self._getTmpPath(fnBase+'.xmd'))
    
    
    def simulatedAnnealing(self, fnRoot):
        self._insertRunJobStep("xmipp_volume_initial_simulated_annealing","-i %s.xmd --initial %s.vol --oroot %s_sa --sym %s --randomIter %d --rejection %f --dontApplyPositive"\
                         %(fnRoot,fnRoot,fnRoot,self.symmetryGroup.get(),self.nIterRandom.get(),self.rejection.get()))
        self._insertRunJobStep('moveFile', source=fnRoot+"_sa.vol", dest=fnRoot+".vol")
        self._insertRunJobStep('deleteFile', filename=fnRoot+"_sa.xmd") 
    
    def reconstruct(self, fnRoot):
        self.runJob(None,"xmipp_reconstruct_fourier","-i %s.xmd -o %s.vol --sym %s " %(fnRoot,fnRoot,self.symmetryGroup.get()))
        self.runJob(None,"xmipp_transform_mask","-i %s.vol --mask circular -%d "%(fnRoot,self.Xdim2/2))
     
    def getCorrThresh(self):
        corrVector = []
        fnCorr=self._getExtraPath("correlations.xmd")               
        mdCorr= MetaData()
    
        for n in range(self.nRansac.get()):
            fnRoot="ransac%05d"%n
            fnAngles=self._getTmpPath("angles_"+fnRoot+".xmd")
            md=MetaData(fnAngles)
            
            for objId in md:
                corr = md.getValue(MDL_MAXCC, objId)
                corrVector.append(corr)
                objIdCorr = mdCorr.addObject()
                mdCorr.setValue(MDL_MAXCC,float(corr),objIdCorr)
    
        mdCorr.write("correlations@"+fnCorr,MD_APPEND)                            
        mdCorr= MetaData()
        sortedCorrVector = sorted(corrVector)
        indx = int(floor(self.corrThresh.get()*(len(sortedCorrVector)-1)))    
        
        #With the line below commented the percentil is not used for the threshold and is used the value introduced in the form
        #CorrThresh = sortedCorrVector[indx]#
            
        objId = mdCorr.addObject()
        mdCorr.setValue(MDL_WEIGHT,self.corrThresh.get(),objId)
        mdCorr.write("corrThreshold@"+fnCorr,MD_APPEND)
        print "Correlation threshold: "+str(self.corrThresh.get())
 
    def evaluateVolumes(self):
        fnCorr=self._getExtraPath("correlations.xmd")
        fnCorr = 'corrThreshold@'+fnCorr
        mdCorr= MetaData(fnCorr)
        objId = mdCorr.firstObject()    
        CorrThresh = mdCorr.getValue(MDL_WEIGHT,objId)
        for n in range(self.nRansac.get()):        
            fnRoot="ransac%05d"%n              
            fnAngles=self._getTmpPath("angles_"+fnRoot+".xmd")    
            md=MetaData(fnAngles)
            numInliers=0
            for objId in md:
                corr = md.getValue(MDL_MAXCC, objId)
               
                if (corr >= CorrThresh) :
                    numInliers = numInliers+corr
    
            md= MetaData()
            objId = md.addObject()
            md.setValue(MDL_WEIGHT,float(numInliers),objId)
            md.write("inliers@"+fnAngles,MD_APPEND)
        
    def getBestVolumes(self):
        volumes = []
        inliers = []
        
        for n in range(self.nRansac.get()):
            fnAngles = self._getTmpPath("angles_ransac%05d"%n+".xmd")
            md=MetaData("inliers@"+fnAngles)
            numInliers=md.getValue(MDL_WEIGHT,md.firstObject())
            volumes.append(fnAngles)
            inliers.append(numInliers)
        
        index = sorted(range(inliers.__len__()), key=lambda k: inliers[k])
        fnBestAngles = ''
        threshold=self.getCCThreshold()
     
        i=self.nRansac.get()-1
        indx = 0
        while i>=0 and indx<self.numVolumes.get():
            fnBestAngles = volumes[index[i]]
            fnBestAnglesOut=self._getPath("proposedVolume%05d"%indx+".xmd")
            copyFile(fnBestAngles,fnBestAnglesOut)
            print("Best volume "+str(indx)+" = "+fnBestAngles)
            if not self.useAll.get():
                self.runJob(None,"xmipp_metadata_utilities","-i %s -o %s --query select \"maxCC>%f \" --mode append" %(fnBestAnglesOut,fnBestAnglesOut,threshold))
                if getMdSize(fnBestAnglesOut) > 0:
                    indx += 1
            else:
                indx += 1
            i -= 1
            
        # Remove unnecessary files
        for n in range(self.nRansac.get()):
            fnAngles = self._getTmpPath("angles_ransac%05d"%n+".xmd")
            cleanPath(fnAngles)
            
            
    def projMatch(self,fnBase):
        fnRoot=self._getPath(fnBase)
        fnGallery=self._getTmpPath('gallery_'+fnBase+'.stk')
        fnOutputReducedClass = self._getExtraPath("reducedClasses.xmd") 
        
        AngularSampling=int(max(floor(self.angularSampling.get()/2.0),2));
        self.runJob(None,"xmipp_angular_project_library", "-i %s.vol -o %s --sampling_rate %f --sym %s --method fourier 1 0.25 bspline --compute_neighbors --angular_distance -1 --experimental_images %s"\
                              %(fnRoot,fnGallery,float(AngularSampling),self.symmetryGroup.get(),fnOutputReducedClass))
    
        self.runJob(None,"xmipp_angular_projection_matching", "-i %s.xmd -o %s.xmd --ref %s --Ri 0 --Ro %s --max_shift %s --append"\
               %(fnRoot,fnRoot,fnGallery,str(self.Xdim/2),str(self.Xdim/20)))
                
        cleanPath(self._getTmpPath('gallery_'+fnBase+'_sampling.xmd'))
        cleanPath(self._getTmpPath('gallery_'+fnBase+'.doc'))
        cleanPath(self._getTmpPath('gallery_'+fnBase+'.stk'))
            
    def getCCThreshold(self):
        fnCorr=self._getExtraPath("correlations.xmd")               
        mdCorr=MetaData("corrThreshold@"+fnCorr)
        return mdCorr.getValue(MDL_WEIGHT, mdCorr.firstObject())
    
    def scoreFinalVolumes(self):
        threshold=self.getCCThreshold()
        mdOut=MetaData()
        for n in range(self.numVolumes.get()):
            fnRoot=self._getPath('proposedVolume%05d'%n)
            fnAssignment=fnRoot+".xmd"
            if exists(fnAssignment):
                self.runJob(None,"xmipp_metadata_utilities","-i %s --fill weight constant 1"%fnAssignment)
                MDassignment=MetaData(fnAssignment)
                sum=0
                thresholdedSum=0
                N=0
                minCC=2
                for id in MDassignment:
                    cc=MDassignment.getValue(MDL_MAXCC,id)
                    sum+=cc
                    thresholdedSum+=cc-threshold
                    if cc<minCC:
                        minCC=cc
                    N+=1
                avg=sum/N
                id=mdOut.addObject()
                mdOut.setValue(MDL_IMAGE,fnRoot+".vol",id)
                mdOut.setValue(MDL_VOLUME_SCORE_SUM,float(sum),id)
                mdOut.setValue(MDL_VOLUME_SCORE_SUM_TH,float(thresholdedSum),id)
                mdOut.setValue(MDL_VOLUME_SCORE_MEAN,float(avg),id)
                mdOut.setValue(MDL_VOLUME_SCORE_MIN,float(minCC),id)
        mdOut.write(self._getPath("proposedVolumes.xmd"))

    def projectInitialVolume(self):
        self.volFn = createXmippInputVolumes(self, self.initialVolume.get())
        
        fnOutputInitVolume=self._getTmpPath("initialVolume.vol")
        self.runJob(None,'xmipp_image_convert',"-i %s -o %s"%(removeExt(self.volFn),fnOutputInitVolume))
        self.runJob(None,"xmipp_image_resize","-i %s --dim %d %d"%(fnOutputInitVolume,self.Xdim2,self.Xdim2))
        fnGallery=self._getTmpPath('gallery_InitialVolume.stk')
        fnOutputReducedClass = self._getExtraPath("reducedClasses.xmd") 
        self.runJob(None,"xmipp_angular_project_library", "-i %s -o %s --sampling_rate %f --sym %s --method fourier 1 0.25 bspline --compute_neighbors --angular_distance -1 --experimental_images %s"\
                              %(fnOutputInitVolume,fnGallery,self.angularSampling.get(),self.symmetryGroup.get(),fnOutputReducedClass))
            
    def createOutput(self):
        fn = self._getPath('proposedVolumes.xmd')
        md = xmipp.MetaData(fn)
        md.addItemId()
        md.write(fn)
        
        volumesSet = self._createSetOfVolumes()
        readSetOfVolumes(fn, volumesSet)
        volumesSet.setSamplingRate(self.inputClasses.get().getAverages().getSamplingRate())
        volumesSet.write()
        
        self._defineOutputs(outputVolumes=volumesSet)

    def _summary(self):
        summary = []
        if not hasattr(self, 'outputVolumes'):
            summary.append("Output volumes not ready yet.")
        else:
            summary.append("RANSAC iterations: %d"%self.nRansac.get())
        
            for n in range(self.numVolumes.get()):
    
                fnBase='proposedVolume%05d'%n
                fnRoot=self._getPath(fnBase+".xmd")
                               
                if os.path.isfile(fnRoot):
                    md=MetaData(fnRoot)
                    if (md.size()< 5) :
                        summary.append("Num of inliers for %s too small and equal to %d"%(fnRoot,md.size()))
                        summary.append("Decrease the value of Inlier Threshold parameter and run again")
                                    
            fnRoot=self._getTmpPath("ransac00000.xmd")    
            
            if os.path.isfile(fnRoot):
                md=MetaData(fnRoot)
            
                if (md.size()< 5) :
                    summary.append("Num of random samples too small and equal to %d"%(md.size()))
                    summary.append("If the option Dimensionality reduction is on, increase the number of grids per dimension")
                    summary.append("If the option Dimensionality reduction is off, increase the number of random samples")
                
            if self.useSA:
                summary.append("Simulated annealing used")
            return summary
            