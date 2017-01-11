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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os
from math import floor

import xmipp
from pyworkflow.object import Float
from pyworkflow.utils.path import cleanPath, copyFile, cleanPattern
from pyworkflow.protocol.params import (PointerParam, FloatParam, BooleanParam,
                                        IntParam, StringParam, 
                                        STEPS_PARALLEL, LEVEL_ADVANCED)
from pyworkflow.em.protocol import ProtInitialVolume
from pyworkflow.em.data import SetOfClasses2D

from convert import writeSetOfClasses2D, readSetOfVolumes, writeSetOfParticles
from utils import isMdEmpty
from pyworkflow.em.convert import ImageHandler


class XmippProtRansac(ProtInitialVolume):
    """ 
    Computes an initial 3d model from a set of projections/classes 
    using RANSAC algorithm.
    
    This method is based on an initial non-lineal dimensionality
    reduction approach which allows to select representative small 
    sets of class average images capturing the most of the structural 
    information of the particle under study. These reduced sets are 
    then used to generate volumes from random orientation assignments. 
    The best volume is determined from these guesses using a random 
    sample consensus (RANSAC) approach.    
     """
    _label = 'ransac'
    
    def __init__(self, **args):
        ProtInitialVolume.__init__(self, **args)
        self.stepsExecutionMode = STEPS_PARALLEL

    #--------------------------- DEFINE param functions --------------------------------------------        
    def _defineParams(self, form):
        form.addSection(label='Input')
         
        form.addParam('inputSet', PointerParam, label="Input averages", important=True, 
                      pointerClass='SetOfClasses2D, SetOfAverages',# pointerCondition='hasRepresentatives',
                      help='Select the input images from the project.'
                           'It should be a SetOfClasses2D object')  
        form.addParam('symmetryGroup', StringParam, default="c1",
                      label='Symmetry group',  
                      help="See http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Symmetry"
                           " for a description of the symmetry groups format in Xmipp.\n"
                           "If no symmetry is present, use _c1_.")
        form.addParam('angularSampling', FloatParam, default=5, expertLevel=LEVEL_ADVANCED,
                      label='Angular sampling rate',
                      help='In degrees.'
                      ' This sampling defines how fine the projection gallery from the volume is explored.')
        form.addParam('nRansac', IntParam, default="400", expertLevel=LEVEL_ADVANCED,
                      label="Number of RANSAC iterations", 
                      help='Number of initial volumes to test by RANSAC')
        
        form.addParam('dimRed', BooleanParam, default=False, expertLevel=LEVEL_ADVANCED,
                      label='Perform dimensionality reduction', 
                      help='The dimensionality reduction is performed using the Local Tangent Space'
                      'Alignment. See http://www.stat.missouri.edu/~ys873/research/LTSA11.pdf')
        form.addParam('numGrids', IntParam, default=3, condition='dimRed', expertLevel=LEVEL_ADVANCED,
                      label='Number of grids per dimension',
                      help='Number of squares to sample the classes')
        form.addParam('numSamples', IntParam, default=8, condition='not dimRed', expertLevel=LEVEL_ADVANCED,
                      label='Number of random samples',
                      help='Number of squares to sample the classes')
        
        form.addParam('corrThresh', FloatParam, default=0.77, expertLevel=LEVEL_ADVANCED,
                      label='Inliers threshold',
                      help='Correlation value threshold to determine if an experimental projection is an inlier or outlier.')        
        form.addParam('numVolumes', IntParam, default=10, expertLevel=LEVEL_ADVANCED,
                      label='Number of best volumes to refine',
                      help='Number of best volumes to refine using projection matching approach and the input classes')
        form.addParam('numIter', IntParam, default=10, expertLevel=LEVEL_ADVANCED,
                      label='Number of iterations to refine the volumes',
                      help='Number of iterations to refine the best volumes using projection matching approach and the input classes')
        form.addParam('initialVolume', PointerParam, label="Initial volume",  expertLevel=LEVEL_ADVANCED,
                      pointerClass='Volume', allowsNull=True,
                      help='You may provide a very rough initial volume as a way to constraint the angular search.'
                            'For instance, when reconstructing a fiber, you may provide a cylinder so that side views'
                            'are assigned to the correct tilt angle, although the rotational angle may be completely wrong')           
                
        form.addParam('maxFreq', IntParam, default=12, expertLevel=LEVEL_ADVANCED,
                      label='Max frequency of the initial volume',
                      help=' Max frequency of the initial volume in Angstroms')
        
        form.addParam('useAll', BooleanParam, default=False, expertLevel=LEVEL_ADVANCED,
                      label='Use all images to refine', 
                      help=' When refining a RANSAC volume, use all images to refine it instead of only inliers')
        
        form.addParallelSection(threads=8, mpi=1)
            
         
    #--------------------------- INSERT steps functions --------------------------------------------    
    def _insertAllSteps(self):
        # Insert some initialization steps
        initialStepId = self._insertInitialSteps()
        
        deps = [] # Store all steps ids, final step createOutput depends on all of them    
        for n in range(self.nRansac.get()):
            # CTF estimation with Xmipp
            stepId = self._insertFunctionStep('ransacIterationStep', n,
                                    prerequisites=[initialStepId]) # Make estimation steps indepent between them
            deps.append(stepId)
        
        # Look for threshold, evaluate volumes and get the best
        self._insertFunctionStep("getCorrThreshStep", prerequisites=deps) # Make estimation steps indepent between them)
        self._insertFunctionStep("evaluateVolumesStep")
        bestVolumesStepId = self._insertFunctionStep("getBestVolumesStep")        
        
        deps = [] # Store all steps ids, final step createOutput depends on all of them
        # Refine the best volumes
        for n in range(self.numVolumes.get()):
            fnBase='proposedVolume%05d' % n
            fnRoot=self._getPath(fnBase)
                    
            for it in range(self.numIter.get()):    
                if it==0:
                    self._insertFunctionStep('reconstructStep',fnRoot, prerequisites=[bestVolumesStepId])
                else:
                    self._insertFunctionStep('reconstructStep',fnRoot)
                self._insertFunctionStep('projMatchStep',fnBase)
            
            stepId =  self._insertFunctionStep("resizeStep",fnRoot,self.Xdim)
            
            deps.append(stepId)
        
        # Score each of the final volumes
        self._insertFunctionStep("scoreFinalVolumes",
                                 prerequisites=deps) # Make estimation steps indepent between them
        
        self._insertFunctionStep('createOutputStep')
    
    def _insertInitialSteps(self):
        # Convert the input classes to a metadata ready for xmipp
        self.imgsFn = self._getExtraPath('input_classes.xmd')
        self._insertFunctionStep('convertInputStep', self.imgsFn)
        
        inputSet = self.inputSet.get()
        self.Xdim = inputSet.getDimensions()[0]
        
        fnOutputReducedClass = self._getExtraPath("reducedClasses.xmd")
        fnOutputReducedClassNoExt = os.path.splitext(fnOutputReducedClass)[0]
    
        # Low pass filter and resize        
        maxFreq = self.maxFreq.get()
        ts = inputSet.getSamplingRate()
        K = 0.25 * (maxFreq / ts)
        if K < 1:
            K = 1
        self.Xdim2 = self.Xdim / K
        if self.Xdim2 < 32:
            self.Xdim2 = 32
            K = self.Xdim / self.Xdim2
            
        freq = ts / maxFreq
        ts = K * ts

        self._insertRunJobStep("xmipp_transform_filter","-i %s -o %s --fourier low_pass %f --oroot %s"
                                                %(self.imgsFn,fnOutputReducedClass,freq,fnOutputReducedClassNoExt))
        lastId = self._insertRunJobStep("xmipp_image_resize","-i %s --fourier %d -o %s" %(fnOutputReducedClass,self.Xdim2,fnOutputReducedClassNoExt))

        # Generate projection gallery from the initial volume
        if self.initialVolume.hasValue():
            lastId = self._insertFunctionStep("projectInitialVolume")
            
        return lastId

    #--------------------------- STEPS functions --------------------------------------------

    def convertInputStep(self, classesFn):
        inputSet = self.inputSet.get()
        
        if isinstance(inputSet, SetOfClasses2D):
            writeSetOfClasses2D(inputSet, classesFn)
        else:
            writeSetOfParticles(inputSet, classesFn)
        
    def ransacIterationStep(self, n):
    
        fnOutputReducedClass = self._getExtraPath("reducedClasses.xmd")  
        fnBase = "ransac%05d"%n
        fnRoot = self._getTmpPath(fnBase)
        
    
        if self.dimRed:
            # Get a random sample of images
            self.runJob("xmipp_transform_dimred","-i %s --randomSample %s.xmd  %d -m LTSA "%(fnOutputReducedClass,fnRoot,self.numGrids.get()))
        else:        
            self.runJob("xmipp_metadata_utilities","-i %s -o %s.xmd  --operate random_subset %d --mode overwrite "%(fnOutputReducedClass,fnRoot,self.numSamples.get()))
            self.runJob("xmipp_metadata_utilities","-i %s.xmd --fill angleRot rand_uniform -180 180 "%(fnRoot))
            self.runJob("xmipp_metadata_utilities","-i %s.xmd --fill angleTilt rand_uniform 0 180 "%(fnRoot))
            self.runJob("xmipp_metadata_utilities","-i %s.xmd --fill anglePsi  rand_uniform 0 360 "%(fnRoot)) 
    
        # If there is an initial volume, assign angles        
        if self.initialVolume.hasValue():
            fnGallery=self._getTmpPath('gallery_InitialVolume.stk')
            self.runJob("xmipp_angular_projection_matching", "-i %s.xmd -o %s.xmd --ref %s --Ri 0 --Ro %s --max_shift %s --append"\
                   %(fnRoot,fnRoot,fnGallery,str(self.Xdim/2),str(self.Xdim/20)))
    
        # Reconstruct with the small sample
        self.reconstructStep(fnRoot)
        
        fnVol = fnRoot+'.vol'
        
        # Generate projections from this reconstruction
        fnGallery=self._getTmpPath('gallery_'+fnBase+'.stk')
        self.runJob("xmipp_angular_project_library", "-i %s -o %s --sampling_rate %f --sym %s --method fourier 1 0.25 bspline --compute_neighbors --angular_distance -1 --experimental_images %s --max_tilt_angle 90"\
                    %(fnVol,fnGallery,self.angularSampling.get(),self.symmetryGroup.get(),fnOutputReducedClass))
            
        # Assign angles to the rest of images
        fnAngles=self._getTmpPath('angles_'+fnBase+'.xmd')
        self.runJob("xmipp_angular_projection_matching", "-i %s -o %s --ref %s --Ri 0 --Ro %s --max_shift %s --append"\
                              %(fnOutputReducedClass,fnAngles,fnGallery,str(self.Xdim/2),str(self.Xdim/20)))
       
        # Delete intermediate files 
        cleanPath(fnGallery)
        cleanPath(self._getTmpPath('gallery_'+fnBase+'_sampling.xmd'))
        cleanPath(self._getTmpPath('gallery_'+fnBase+'.doc'))
        cleanPath(fnVol)
        cleanPath(self._getTmpPath(fnBase+'.xmd'))
    
    def reconstructStep(self, fnRoot):
        from pyworkflow.em.metadata.utils import getSize
        if os.path.exists(fnRoot+".xmd"):
            Nimages=getSize(fnRoot+".xmd")
            if Nimages>0:
                self.runJob("xmipp_reconstruct_fourier","-i %s.xmd -o %s.vol --sym %s " %(fnRoot,fnRoot,self.symmetryGroup.get()))
                self.runJob("xmipp_transform_mask","-i %s.vol --mask circular -%d "%(fnRoot,self.Xdim2/2))
        else:
            print fnRoot+".xmd is empty. The corresponding volume is not generated." 
    
    def resizeStep(self,fnRoot,Xdim):
        if os.path.exists(fnRoot+".vol"):
            self.runJob("xmipp_image_resize","-i %s.vol -o %s.vol --dim %d %d" %(fnRoot,fnRoot,Xdim,Xdim))
        
     
    def getCorrThreshStep(self):
        corrVector = []
        fnCorr=self._getExtraPath("correlations.xmd")               
        mdCorr= xmipp.MetaData()
    
        for n in range(self.nRansac.get()):
            fnRoot="ransac%05d"%n
            fnAngles=self._getTmpPath("angles_"+fnRoot+".xmd")
            md = xmipp.MetaData(fnAngles)
            
            for objId in md:
                corr = md.getValue(xmipp.MDL_MAXCC, objId)
                corrVector.append(corr)
                objIdCorr = mdCorr.addObject()
                mdCorr.setValue(xmipp.MDL_MAXCC,float(corr),objIdCorr)
    
        mdCorr.write("correlations@"+fnCorr,xmipp.MD_APPEND)                            
        mdCorr= xmipp.MetaData()
        sortedCorrVector = sorted(corrVector)
        indx = int(floor(self.corrThresh.get()*(len(sortedCorrVector)-1)))    
        
        #With the line below commented the percentil is not used for the threshold and is used the value introduced in the form
        #CorrThresh = sortedCorrVector[indx]#
            
        objId = mdCorr.addObject()
        mdCorr.setValue(xmipp.MDL_WEIGHT,self.corrThresh.get(),objId)
        mdCorr.write("corrThreshold@"+fnCorr,xmipp.MD_APPEND)
        print "Correlation threshold: "+str(self.corrThresh.get())
    
    
    def evaluateVolumesStep(self):
        fnCorr=self._getExtraPath("correlations.xmd")
        fnCorr = 'corrThreshold@'+fnCorr
        mdCorr= xmipp.MetaData(fnCorr)
        objId = mdCorr.firstObject()    
        CorrThresh = mdCorr.getValue(xmipp.MDL_WEIGHT,objId)
        for n in range(self.nRansac.get()):        
            fnRoot="ransac%05d"%n              
            fnAngles=self._getTmpPath("angles_"+fnRoot+".xmd")    
            md = xmipp.MetaData(fnAngles)
            numInliers=0
            for objId in md:
                corr = md.getValue(xmipp.MDL_MAXCC, objId)
               
                if (corr >= CorrThresh) :
                    numInliers = numInliers+corr
    
            md= xmipp.MetaData()
            objId = md.addObject()
            md.setValue(xmipp.MDL_WEIGHT,float(numInliers),objId)
            md.write("inliers@"+fnAngles,xmipp.MD_APPEND)
    
    def getBestVolumesStep(self):
        volumes = []
        inliers = []
        
        for n in range(self.nRansac.get()):
            fnAngles = self._getTmpPath("angles_ransac%05d"%n+".xmd")
            md=xmipp.MetaData("inliers@"+fnAngles)
            numInliers=md.getValue(xmipp.MDL_WEIGHT,md.firstObject())
            volumes.append(fnAngles)
            inliers.append(numInliers)
        
        index = sorted(range(inliers.__len__()), key=lambda k: inliers[k])
        fnBestAngles = ''
        threshold=self.getCCThreshold()
     
        i = self.nRansac.get()-1
        indx = 0
        while i >= 0 and indx < self.numVolumes:
            fnBestAngles = volumes[index[i]]
            fnBestAnglesOut = self._getPath("proposedVolume%05d"%indx+".xmd")
            copyFile(fnBestAngles, fnBestAnglesOut)
            self._log.info("Best volume %d = %s" % (indx, fnBestAngles))
            if not self.useAll:
                self.runJob("xmipp_metadata_utilities","-i %s -o %s --query select \"maxCC>%f \" --mode append" %(fnBestAnglesOut,fnBestAnglesOut,threshold))
                if not isMdEmpty(fnBestAnglesOut):
                    indx += 1
            else:
                indx += 1
            i -= 1
            
        # Remove unnecessary files
        for n in range(self.nRansac.get()):
            fnAngles = self._getTmpPath("angles_ransac%05d"%n+".xmd")
            cleanPath(fnAngles)
             
    def projMatchStep(self,fnBase):
        fnRoot=self._getPath(fnBase)
        if os.path.exists(fnRoot+".vol"):
            fnGallery=self._getTmpPath('gallery_'+fnBase+'.stk')
            fnOutputReducedClass = self._getExtraPath("reducedClasses.xmd") 
            
            AngularSampling=max(self.angularSampling.get()/2.0,7.5);
            self.runJob("xmipp_angular_project_library", "-i %s.vol -o %s --sampling_rate %f --sym %s --method fourier 1 0.25 bspline --compute_neighbors --angular_distance -1 --experimental_images %s"\
                                  %(fnRoot,fnGallery,float(AngularSampling),self.symmetryGroup.get(),fnOutputReducedClass))
        
            self.runJob("xmipp_angular_projection_matching", "-i %s.xmd -o %s.xmd --ref %s --Ri 0 --Ro %s --max_shift %s --append"\
                   %(fnRoot,fnRoot,fnGallery,str(self.Xdim/2),str(self.Xdim/20)))
                    
            cleanPath(self._getTmpPath('gallery_'+fnBase+'_sampling.xmd'))
            cleanPath(self._getTmpPath('gallery_'+fnBase+'.doc'))
            cleanPath(self._getTmpPath('gallery_'+fnBase+'.stk'))
                 
    def scoreFinalVolumes(self):
        threshold=self.getCCThreshold()
        mdOut=xmipp.MetaData()
        for n in range(self.numVolumes.get()):
            fnRoot=self._getPath('proposedVolume%05d'%n)
            fnAssignment=fnRoot+".xmd"
            if os.path.exists(fnAssignment):
                self.runJob("xmipp_metadata_utilities","-i %s --fill weight constant 1"%fnAssignment)
                MDassignment=xmipp.MetaData(fnAssignment)
                sum=0
                thresholdedSum=0
                N=0
                minCC=2
                for id in MDassignment:
                    cc=MDassignment.getValue(xmipp.MDL_MAXCC,id)
                    sum+=cc
                    thresholdedSum+=cc-threshold
                    if cc<minCC:
                        minCC=cc
                    N+=1
                if N>0:
                    avg=sum/N
                else:
                    avg=0.0
                id=mdOut.addObject()
                mdOut.setValue(xmipp.MDL_IMAGE,fnRoot+".vol",id)
                mdOut.setValue(xmipp.MDL_VOLUME_SCORE_SUM,float(sum),id)
                mdOut.setValue(xmipp.MDL_VOLUME_SCORE_SUM_TH,float(thresholdedSum),id)
                mdOut.setValue(xmipp.MDL_VOLUME_SCORE_MEAN,float(avg),id)
                mdOut.setValue(xmipp.MDL_VOLUME_SCORE_MIN,float(minCC),id)
        mdOut.write(self._getPath("proposedVolumes.xmd"))


    def projectInitialVolume(self):
        fnOutputInitVolume=self._getTmpPath("initialVolume.vol")
        img = ImageHandler()
        img.convert(self.initialVolume.get(), fnOutputInitVolume)
        self.runJob("xmipp_image_resize","-i %s --dim %d %d"%(fnOutputInitVolume,self.Xdim2,self.Xdim2))
        fnGallery=self._getTmpPath('gallery_InitialVolume.stk')
        fnOutputReducedClass = self._getExtraPath("reducedClasses.xmd") 
        self.runJob("xmipp_angular_project_library", "-i %s -o %s --sampling_rate %f --sym %s --method fourier 1 0.25 bspline --compute_neighbors --angular_distance -1 --experimental_images %s"\
                              %(fnOutputInitVolume,fnGallery,self.angularSampling.get(),self.symmetryGroup.get(),fnOutputReducedClass))
   
    def _postprocessVolume(self, vol, row):
        self._counter += 1
        vol.setObjComment('ransac volume %02d' % self._counter)
        vol._xmipp_volScoreSum   = Float(row.getValue(xmipp.MDL_VOLUME_SCORE_SUM))
        vol._xmipp_volScoreSumTh = Float(row.getValue(xmipp.MDL_VOLUME_SCORE_SUM_TH))
        vol._xmipp_volScoreMean = Float(row.getValue(xmipp.MDL_VOLUME_SCORE_MEAN))
        vol._xmipp_volScoreMin = Float(row.getValue(xmipp.MDL_VOLUME_SCORE_MIN))
        
    def createOutputStep(self):
        inputSet = self.inputSet.get()
        fn = self._getPath('proposedVolumes.xmd')
#        md = xmipp.MetaData(fn)
#        md.addItemId()
#        md.write(fn)
        
        volumesSet = self._createSetOfVolumes()
        volumesSet.setSamplingRate(inputSet.getSamplingRate())
        self._counter = 0
        readSetOfVolumes(fn, volumesSet, postprocessImageRow=self._postprocessVolume)
        
        # Set a meanful comment
#         for vol in volumesSet:
#             vol.setObjComment('ransac volume %02d' % vol.getObjId())
#             volumesSet.update(vol)
        
        self._defineOutputs(outputVolumes=volumesSet)
        self._defineSourceRelation(self.inputSet, volumesSet)
        self._storeSummaryInfo(self.numVolumes.get())
    
    #--------------------------- INFO functions --------------------------------------------
    def _validate(self):
        errors = []
        inputSet = self.inputSet.get()
        if isinstance(inputSet, SetOfClasses2D):
            if not self.inputSet.get().hasRepresentatives():
                errors.append("The input classes should have representatives.")
                
                
        if self.dimRed:
            nGrids = self.numGrids.get()
            if (nGrids * nGrids) > inputSet.getSize():
                errors.append('Dimensionaly reduction could not be applied')
                errors.append('if the number of classes is less than the number')
                errors.append('of grids squared. \n')
                errors.append('Consider either provide more classes or')
                errors.append('disable dimensionality reduction')
        return errors
    
    def _summary(self):
        summary = []
        if hasattr(self, 'outputVolumes'):                        
            summary.append("RANSAC iterations: %d" % self.nRansac)
            summary.append("Number of volumes to refine: %d" % self.numVolumes.get())
            summary.append("Number of refinement interations: %d" % self.numIter.get())
            if self.summaryVar.hasValue():
                summary.append(self.summaryVar.get())
        else:
            summary.append("Output volumes not ready yet.")

        return summary
        
    def _methods(self):
        strline=''
        if hasattr(self, 'outputVolumes'):
            strline+='We obtained %d initial volume(s) from %s set of images and symmetry %s.'%\
                           (self.numVolumes.get(),self.getObjectTag('inputSet'),self.symmetryGroup.get())
            strline+='The number of RANSAC volumes generated are %d.From this set, we selected the best %d volume(s) that was refined by %d projection matching iteration(s)'%\
                           (self.nRansac.get(),self.numVolumes.get(),self.numIter.get())
            strline+= ' [Vargas2014]'
            if self.dimRed:
                strline+='We performed dimensionality reduction using a %d x %d grid'%\
                            (self.numGrids.get(),self.numGrids.get())
        return [strline]               
    #--------------------------- UTILS functions --------------------------------------------        
    def getCCThreshold(self):
        fnCorr = self._getExtraPath("correlations.xmd")               
        mdCorr = xmipp.MetaData("corrThreshold@"+fnCorr)
        return mdCorr.getValue(xmipp.MDL_WEIGHT, mdCorr.firstObject())
    
    def _storeSummaryInfo(self, numVolumes):
        """ Store some information when the protocol finishes. """
        msg1 = ''
        msg2 = ''
        
        for n in range(numVolumes):
            fnBase = 'proposedVolume%05d' % n
            fnRoot = self._getPath(fnBase + ".xmd")
                               
            if os.path.isfile(fnRoot):
                md = xmipp.MetaData(fnRoot)
                size = md.size()
                if (size < 5):
                    msg1 = "Num of inliers for model %d too small and equal to %d \n" % (n, size)
                    msg1 += "Decrease the value of Inlier Threshold parameter and run again \n"
                 
        fnRoot = self._getTmpPath("ransac00000.xmd")
        if os.path.isfile(fnRoot):
            md = xmipp.MetaData(fnRoot)
            size = md.size()
            if (size < 5):
                msg2 = "Num of random samples too small and equal to %d.\n" % size
                msg2 += "If the option Dimensionality reduction is on, increase the number of grids per dimension (In this case we recommend to put Dimensionality reduction off).\n"
                msg2 += "If the option Dimensionality reduction is off, increase the number of random samples.\n"
                
        msg = msg1 + msg2
        self.summaryVar.set(msg)

    def _citations(self):
        return ['Vargas2014']
