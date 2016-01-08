# **************************************************************************
# *
# * Authors:     Mohsen Kazemi  (mkazemi@cnb.csic.es)
# *              C.O.S. Sorzano (coss@cnb.csic.es)
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

from pyworkflow.protocol.params import (PointerParam, FloatParam, NumericListParam, IntParam,
                                        StringParam, BooleanParam, LEVEL_ADVANCED)
from pyworkflow.em.data import Volume
from pyworkflow.em.protocol import ProtReconstruct3D
from pyworkflow.em.packages.xmipp3.convert import writeSetOfParticles
from pyworkflow.utils import getFloatListFromValues
from pyworkflow.utils.path import cleanPattern, cleanPath, copyFile
import os
import xmipp
import glob
from pyworkflow.object import Float, String
from math import sqrt
from plotter import XmippPlotter


class XmippProtValidateOverfitting(ProtReconstruct3D):
    """    
    Check how the FSC changes with the number of projections used for 3D reconstruction. This method
    has been proposed by B. Heymann "Validation of 3D EM Reconstructions", 2015
    """
    _label = 'validate overfitting'
    
    #--------------------------- DEFINE param functions --------------------------------------------   
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputParticles', PointerParam, pointerClass='SetOfParticles',pointerCondition='hasAlignmentProj',
                      label="Input particles",  
                      help='Select the input images from the project.')     
        form.addParam('input3DReference', PointerParam,
                 pointerClass='Volume,SetOfVolumes',
                 label='Initial 3D reference volume', 
                 help='Input 3D reference reconstruction to alignment gaussian nise.') 
        form.addParam('symmetryGroup', StringParam, default='c1',
                      label="Symmetry group", 
                      help='See [[Xmipp Symmetry][http://www2.mrc-lmb.cam.ac.uk/Xmipp/index.php/Conventions_%26_File_formats#Symmetry]] page '
                           'for a description of the symmetry format accepted by Xmipp')
        form.addParam('numberOfParticles', NumericListParam, default="100 200 500 1000 2000 5000", expertLevel=LEVEL_ADVANCED,
                      label="Number of particles") 
        form.addParam('numberOfIterations', IntParam, default=10, expertLevel=LEVEL_ADVANCED,
                      label="Number of times the randomization is performed") 
        form.addParam('maxRes', FloatParam, default = 0.5, expertLevel=LEVEL_ADVANCED,
                      label="Maximum resolution (dig.freq)",  
                      help='Nyquist is 0.5') 
        form.addParam('angSampRate', FloatParam, default = 5, expertLevel=LEVEL_ADVANCED,
                      label="Angular sampling rate")  
                      
        form.addParallelSection(threads=4, mpi=1)

    #--------------------------- INSERT steps functions --------------------------------------------

    def _createFilenameTemplates(self):
        """ Centralize how files are called for iterations and references. """
        myDict = {
            'input_xmd': self._getExtraPath('input_particles.xmd')
            }
        self._updateFilenamesDict(myDict)

    def _insertAllSteps(self):
        self._createFilenameTemplates()
        
        #convertInputStep 
        particlesMd = self._getFileName('input_xmd')
        imgSet = self.inputParticles.get()
        writeSetOfParticles(imgSet, particlesMd)
                
        #for debugging purpose
        debugging = False      
                 
        #projections from reference volume
        
        args="-i %s -o %s --sampling_rate %f --sym %s --min_tilt_angle %f --max_tilt_angle %f --compute_neighbors --angular_distance -1 --experimental_images %s"%\
            (self.input3DReference.get().getFileName(),self._getExtraPath("Ref_Projections.stk"),
             self.angSampRate.get(),self.symmetryGroup.get(),0,90, self._getFileName('input_xmd'))
        
        self.runJob("xmipp_angular_project_library",args, numberOfMpi=self.numberOfMpi.get()*self.numberOfThreads.get())
        
        
        numberOfParticles=getFloatListFromValues(self.numberOfParticles.get())
        fractionCounter=0
        for number in numberOfParticles:
            if number<self.inputParticles.get().getSize():
                for iteration in range(0,self.numberOfIterations.get()):
                    self._insertFunctionStep('reconstructionStep',number,fractionCounter,iteration, debugging)
                fractionCounter+=1     
        self._insertFunctionStep('gatherResultsStep', debugging)
        
    #--------------------------- STEPS functions --------------------------------------------
    def reconstructionStep(self, numberOfImages, fractionCounter, iteration, debugging):
        fnRoot = self._getExtraPath("fraction%02d"%fractionCounter)
        Ts = self.inputParticles.get().getSamplingRate()
        
        #for noise
        fnRootN = self._getExtraPath("Nfraction%02d"%fractionCounter)
             
        
        for i in range(0,2):
            fnImgs = fnRoot+"_images_%02d_%02d.xmd"%(i, iteration)
            self.runJob("xmipp_metadata_utilities","-i %s -o %s --operate random_subset %d"%\
                        (self._getFileName('input_xmd'),fnImgs,numberOfImages),numberOfMpi=1)
        
            params =  '  -i %s' % fnImgs
            params += '  -o %s' % fnRoot+"_%02d_%02d.vol"% (i, iteration)
            params += ' --sym %s' % self.symmetryGroup.get()
            params += ' --max_resolution %0.3f' % self.maxRes.get()
            params += ' --padding 2'
            params += ' --thr %d' % self.numberOfThreads.get()
            params += ' --sampling %f' % Ts
            self.runJob('xmipp_reconstruct_fourier', params)
            
            #for noise
            noiseStk = fnRoot+"_noises_%02d.stk"%i
            self.runJob ("xmipp_image_convert", "-i %s -o %s"% (fnImgs, noiseStk))
            self.runJob("xmipp_image_operate", "-i %s --mult 0"% noiseStk)
            self.runJob("xmipp_transform_add_noise", "-i %s --type gaussian 3"% noiseStk)
            fnImgsNL = fnRoot+"_noisesL_%02d.xmd"%i
            noiseStk2 = fnRoot+"_noises2_%02d.stk"%i
            self.runJob ("xmipp_image_convert", "-i %s -o %s --save_metadata_stack %s"% (noiseStk, noiseStk2, fnImgsNL))
            fnImgsNoiseOld = fnRoot+"_noisesOld_%02d.xmd"%i
            fnImgsN = fnRoot+"_noises_%02d_%02d.xmd"%(i, iteration)
            self.runJob("xmipp_metadata_utilities",'-i %s -o %s --operate drop_column "image"'% (fnImgs,fnImgsNoiseOld))
            self.runJob("xmipp_metadata_utilities","-i %s  --set merge %s -o %s"% (fnImgsNL,fnImgsNoiseOld, fnImgsN))
           
            
            
            #alignment gaussian noise
            fnImgsAlign = self._getExtraPath("Nfraction_alignment%02d"%fractionCounter)
            fnImgsAlignN = fnImgsAlign + "_%02d_%02d.xmd"%(i, iteration)
                             
            args="-i %s -o %s -r %s --Ri 0 --Ro -1 --mem 2 --append "%\
                (fnImgsN,fnImgsAlignN,self._getExtraPath("Ref_Projections.stk"))
            
            self.runJob('xmipp_angular_projection_matching',args, numberOfMpi=self.numberOfMpi.get()*self.numberOfThreads.get())
           
           
           
            params =  '  -i %s' % fnImgsAlignN
            params += '  -o %s' % fnRootN+"_%02d_%02d.vol"%(i, iteration)
            params += ' --sym %s' % self.symmetryGroup.get()
            params += ' --max_resolution %0.3f' % self.maxRes.get()
            params += ' --padding 2'
            params += ' --thr %d' % self.numberOfThreads.get()
            params += ' --sampling %f' % Ts
            self.runJob('xmipp_reconstruct_fourier', params)
            
            

        self.runJob('xmipp_resolution_fsc', "--ref %s -i %s -o %s --sampling_rate %f"%\
                    (fnRoot+"_00_%02d.vol"%iteration,fnRoot+"_01_%02d.vol"%iteration,fnRoot+"_fsc_%02d.xmd"%iteration,Ts), numberOfMpi=1)
        
        
        mdFSC = xmipp.MetaData(fnRoot+"_fsc_%02d.xmd"%iteration)
        for id in mdFSC:
            fscValue = mdFSC.getValue(xmipp.MDL_RESOLUTION_FRC,id)
            maxFreq = mdFSC.getValue(xmipp.MDL_RESOLUTION_FREQREAL,id)
            if fscValue<0.5:
                break
        fh = open(fnRoot+"_freq.txt","a")
        fh.write("%f\n"%maxFreq)
        fh.close()
        
        
        #for noise
        self.runJob('xmipp_resolution_fsc', "--ref %s -i %s -o %s --sampling_rate %f"%\
                    (fnRootN+"_00_%02d.vol"%iteration,fnRootN+"_01_%02d.vol"%iteration,fnRootN+"_fsc_%02d.xmd"%iteration,Ts), numberOfMpi=1)
        
        cleanPattern(fnRoot+"_noises_0?_0?.xmd")
        cleanPattern(fnRoot+"_noisesOld_0?.xmd")
        cleanPattern(fnRoot+"_noisesL_0?.xmd")
        cleanPattern(fnRoot+"_noises2_0?.stk")
        
        
        mdFSCN = xmipp.MetaData(fnRootN+"_fsc_%02d.xmd"%iteration)
        for id in mdFSCN:
            fscValueN = mdFSCN.getValue(xmipp.MDL_RESOLUTION_FRC,id)
            maxFreqN = mdFSCN.getValue(xmipp.MDL_RESOLUTION_FREQREAL,id)
            if fscValueN<0.5:
                break
        fhN = open(fnRootN+"_freq.txt","a")
        fhN.write("%f\n"%maxFreqN)
        fhN.close()
       
        if not debugging:
            cleanPattern(fnRoot+"_0?_0?.vol")
            cleanPattern(fnRoot+"_images_0?_0?.xmd") 
            cleanPattern(fnRoot+"_fsc_0?.xmd")
            cleanPattern(fnRootN+"_0?_0?.vol")
            cleanPattern(fnRoot+"_noises_0?.stk")
            cleanPattern(fnRootN+"_fsc_0?.xmd")
            cleanPattern(fnImgsAlign + "_0?_0?.xmd")
            
    def gatherResultsStep(self, debugging):
        fnFreqs = sorted(glob.glob(self._getExtraPath("fraction*_freq.txt")))
        subset = 0
        
        numberOfParticles=getFloatListFromValues(self.numberOfParticles.get())
        validationMd = xmipp.MetaData()

        for fnFreq in fnFreqs:
            print fnFreq
            data = []
            fnFreqOpen = open(fnFreq, "r")
            for line in fnFreqOpen:
                fields = line.split()
                rowdata = map(float, fields)
                data.extend(rowdata)
            meanRes = (sum(data)/len(data))
            data[:] = [(x-meanRes)**2 for x in data]
            varRes = (sum(data)/(len(data)-1))
            stdRes = sqrt(varRes)
            
            objId = validationMd.addObject()
            validationMd.setValue(xmipp.MDL_COUNT,long(numberOfParticles[subset]),objId)
            validationMd.setValue(xmipp.MDL_AVG,meanRes, objId)  
            validationMd.setValue(xmipp.MDL_STDDEV,stdRes,objId)

            subset += 1

        validationMd.write(self._defineResultsName())
        
        #for noise
        fnFreqsN = sorted(glob.glob(self._getExtraPath("Nfraction*_freq.txt")))
        subset = 0
        
        numberOfParticles=getFloatListFromValues(self.numberOfParticles.get())
        validationMdN = xmipp.MetaData()

        for fnFreq in fnFreqsN:
            data = []
            fnFreqOpen = open(fnFreq, "r")
            for line in fnFreqOpen:
                fields = line.split()
                rowdata = map(float, fields)
                data.extend(rowdata)
            meanRes = (sum(data)/len(data))
            data[:] = [(x-meanRes)**2 for x in data]
            varRes = (sum(data)/(len(data)-1))
            stdRes = sqrt(varRes)
            
            objId = validationMdN.addObject()
            validationMdN.setValue(xmipp.MDL_COUNT,long(numberOfParticles[subset]),objId)
            validationMdN.setValue(xmipp.MDL_AVG,meanRes, objId)  
            validationMdN.setValue(xmipp.MDL_STDDEV,stdRes,objId)

            subset += 1

        validationMdN.write(self._defineResultsNoiseName())
        
        if not debugging:
            cleanPattern(self._getExtraPath("fraction*_freq.txt"))
            cleanPattern(self._getExtraPath("Nfraction*_freq.txt"))
            cleanPattern(self._getExtraPath('Ref_Projections*'))
        
    #--------------------------- INFO functions -------------------------------------------- 
    def _summary(self):
        """ Should be overriden in subclasses to 
        return summary message for NORMAL EXECUTION. 
        """
        msg=[]
        msg.append("Number of particles: "+self.numberOfParticles.get())
        msg.append("Number of times that reconstruction is performed per each particles subset: %2d" % self.numberOfIterations.get())
        return msg
    
    def _methods(self):
        messages = []
        messages.append('B. Heymann "Validation of 3D EM Reconstructions"')
        return messages
    
    def _citations(self):
        return ['B.Heymann2015']
    
    def _validate(self):
        errors=[]
        maxNumberOfParticles=max(getFloatListFromValues(self.numberOfParticles.get()))
        if maxNumberOfParticles>0.5*self.inputParticles.get().getSize():
            errors.append("The number of tested particles should not be larger than 1/2 of the input set of particles")
         
        return errors
          
    #--------------------------- UTILS functions --------------------------------------------
    def _defineResultsName(self):
        return self._getExtraPath('results.xmd')
    
    def _defineResultsNoiseName(self):
        return self._getExtraPath('resultsNoise.xmd')
