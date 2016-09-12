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

from pyworkflow.protocol.params import (PointerParam, FloatParam, 
                                        NumericListParam, IntParam,
                                        StringParam, BooleanParam,
                                        LEVEL_ADVANCED)
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
    Check how the resolution changes with the number of projections used for 
    3D reconstruction. 
    
    NOTE:
    Using the output plot, with the reconstruction of aligned gaussian noise, 
    you can assess the validity of the reconstruction from your micrograph 
    images. Practically, if the resolution of reconstruction based on your 
    images is not considerably different from aligned gaussian noise one 
    (for less number of particles),your images may not produce a valid 
    reconstruction.
    
    This method has been proposed by:
    B. Heymann "Validation of 3D EM Reconstructions", 2015. 
    (see References)
    """
    _label = 'validate overfitting'    
    #--------------------------- DEFINE param functions --------------------------------------------   
   
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputParticles', PointerParam, 
                      pointerClass='SetOfParticles',
                      pointerCondition='hasAlignmentProj',
                      label="Input particles",  
                      help='Select the input images from the project.')     
        form.addParam('doResize', BooleanParam, default=False,
                      label='Resize input particles and volume?',
                      help="If obtaining the best possible reconstruction "
                           "is not your goal, you can resize your input "
                           "particales and volume to reduce running time "
                           "of the protocl")
        form.addParam('newSize', FloatParam, default=64, condition="doResize",
                      label='New size (px)',
                      help="Resizing input particles and volume"
                           "using fourier method")        
        form.addParam('input3DReference', PointerParam,
                 pointerClass='Volume,SetOfVolumes',
                 label='Initial 3D reference volume', 
                 help="Input 3D reference reconstruction to "
                      "align gaussian noise.") 
        form.addParam('symmetryGroup', StringParam, default='c1',
                      label="Symmetry group", 
                      help="See [[Xmipp Symmetry][http://www2.mrc-lmb.cam.ac.uk/Xmipp/index.php/Conventions_%26_File_formats#Symmetry]] page"
                           "for a description of the symmetry format"
                           "accepted by Xmipp")
        form.addParam('numberOfParticles', NumericListParam,
                      default="10 20 50 100 200 500 1000 1500 2000 3000 5000",
                      expertLevel=LEVEL_ADVANCED,
                      label="Number of particles",
                      help="Number of particles in each subset and consequently "
                           "number of subsets (for instance, in default values,"
                           "a number of 6 subsets with given values are chosen)\n"
                           "Note:\n"
                           "The number of particles in each subset should not "
                           "be larger than 1/2 of the input set of particles. "
                           "The protocol consider this issue automatically. It "
                           "means that if the input set of particles are lower "
                           "than 10,000, you could leave default values unchanged.") 
        form.addParam('numberOfIterations', IntParam, default=10,
                      expertLevel=LEVEL_ADVANCED,
                      label="Number of times the randomization is performed") 
        form.addParam('maxRes', FloatParam, default = 0.5,
                      expertLevel=LEVEL_ADVANCED,
                      label="Maximum resolution (dig.freq)",  
                      help='Nyquist is 0.5') 
        form.addParam('angSampRate', FloatParam, default = 5,
                      expertLevel=LEVEL_ADVANCED,
                      label="Angular sampling rate")  
                      
        form.addParallelSection(threads=0, mpi=4)
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
        
        inputNameRefVol = self.input3DReference.get().getFileName()
        fnNewVol = self._getExtraPath('newVolume.vol')
        fnNewImgStk = self._getExtraPath('newImages.stk')
        fnNewImgMd = self._getExtraPath('newImages.xmd')
        #do resizing
        if self.doResize.get():            
            args = "-i %s ""-o %s --fourier %f " % (inputNameRefVol, 
                                                    fnNewVol,
                                                    self.newSize)
            self.runJob("xmipp_image_resize", args)                                   
            
            args = "-i %s -o %s --fourier %f" % (particlesMd,
                                               fnNewImgStk,
                                               self.newSize)
            args += " --save_metadata_stack %s" % fnNewImgMd
            args += " --keep_input_columns"
            
            self.runJob("xmipp_image_resize", args)
                                    
            oldSize = self.inputParticles.get().getDim()[0]
            scaleFactor = oldSize/self.newSize.get()
            
            args = "-i %s" % fnNewImgMd
            args += " --operate modify_values 'shiftX=shiftX*%f'" % scaleFactor 
            self.runJob('xmipp_metadata_utilities', args)
            
            args =  "-i %s" % fnNewImgMd
            args += " --operate modify_values 'shiftY=shiftY*%f'" % scaleFactor           
            self.runJob('xmipp_metadata_utilities', args)
                        
        #projections from reference volume        
        if self.doResize.get():
            args = "-i %s -o %s" % (fnNewVol,
                                    self._getExtraPath('Ref_Projections.stk'))
            args += " --experimental_images %s" % fnNewImgMd
        else:
            args = "-i %s -o %s" % (self.input3DReference.get().getFileName(),
                                    self._getExtraPath('Ref_Projections.stk'))
            args += " --experimental_images %s" % particlesMd    
        args += " --sampling_rate %f --sym %s" % (self.angSampRate,
                                                  self.symmetryGroup.get())
        args += " --min_tilt_angle 0 --max_tilt_angle 90"                
        args += " --compute_neighbors --angular_distance -1"            
        self.runJob("xmipp_angular_project_library",
                    args,
                    numberOfMpi = self.numberOfMpi.get() * self.numberOfThreads.get())
        
        numberOfParticles = getFloatListFromValues(self.numberOfParticles.get())
        fractionCounter = 0
        maxNumberOfParticles = 0.5 * self.inputParticles.get().getSize()        
        for number in numberOfParticles:
            if number <= maxNumberOfParticles:                
                for iteration in range(0,self.numberOfIterations.get()):
                    self._insertFunctionStep('reconstructionStep', number, 
                                             fractionCounter, iteration, 
                                             debugging, fnNewImgMd,
                                             particlesMd)
                fractionCounter+=1     
        self._insertFunctionStep('gatherResultsStep', debugging)        
    #--------------------------- STEPS functions --------------------------------------------
    
    def reconstructionStep(self, numberOfImages, fractionCounter,
                           iteration, debugging,
                           fnNewImgMd, particlesMd):
        fnRoot = self._getExtraPath("fraction%02d"%fractionCounter)
        Ts = self.inputParticles.get().getSamplingRate()
        
        #for noise
        fnRootN = self._getExtraPath("Nfraction%02d"%fractionCounter)            
        
        for i in range(0,2):
            fnImgs = fnRoot+"_images_%02d_%02d.xmd"%(i, iteration)
            
            if self.doResize.get():
                self.runJob("xmipp_metadata_utilities",
                            "-i %s -o %s --operate random_subset %d" % (
                            fnNewImgMd, fnImgs,numberOfImages), 
                            numberOfMpi = 1)
            else:
                self.runJob("xmipp_metadata_utilities",
                            "-i %s -o %s --operate random_subset %d" % (
                            particlesMd,fnImgs,numberOfImages),
                            numberOfMpi = 1)
        
            params =  '  -i %s' % fnImgs
            params += '  -o %s' % fnRoot+"_%02d_%02d.vol" % (i, iteration)
            params += ' --sym %s' % self.symmetryGroup.get()
            params += ' --max_resolution %0.3f' % self.maxRes
            params += ' --padding 2'
            params += ' --thr 1'
            #params += ' --thr %d' % self.numberOfThreads.get()
            params += ' --sampling %f' % Ts
            self.runJob('xmipp_reconstruct_fourier', params)
            
            #for noise
            noiseStk = fnRoot+"_noises_%02d.stk"%i
            self.runJob ("xmipp_image_convert",
                         "-i %s -o %s" % (fnImgs, noiseStk), numberOfMpi = 1)
            self.runJob("xmipp_image_operate", 
                        "-i %s --mult 0" % noiseStk)
            self.runJob("xmipp_transform_add_noise", 
                        "-i %s --type gaussian 3" % noiseStk, numberOfMpi = 1)
            fnImgsNL = fnRoot+"_noisesL_%02d.xmd" % i
            noiseStk2 = fnRoot+"_noises2_%02d.stk" % i
            self.runJob ("xmipp_image_convert", 
                         "-i %s -o %s --save_metadata_stack %s" % (
                         noiseStk, noiseStk2, fnImgsNL), numberOfMpi = 1)
            fnImgsNoiseOld = fnRoot+"_noisesOld_%02d.xmd" % i
            fnImgsN = fnRoot+"_noises_%02d_%02d.xmd" % (i, iteration)
            self.runJob("xmipp_metadata_utilities",
                        '-i %s -o %s --operate drop_column "image"' % ( 
                        fnImgs,fnImgsNoiseOld), numberOfMpi = 1)
            self.runJob("xmipp_metadata_utilities",
                        "-i %s  --set merge %s -o %s" % (
                        fnImgsNL, fnImgsNoiseOld, fnImgsN), numberOfMpi = 1)               
            
            #alignment gaussian noise
            fnImgsAlign = self._getExtraPath("Nfraction_alignment%02d" % fractionCounter)
            fnImgsAlignN = fnImgsAlign + "_%02d_%02d.xmd" % (i,
                                                             iteration)
                             
            args="-i %s -o %s -r %s --Ri 0 --Ro -1 --mem 2 --append " % (
                 fnImgsN,fnImgsAlignN,
                 self._getExtraPath("Ref_Projections.stk"))
            
            self.runJob('xmipp_angular_projection_matching',
                        args, 
                        numberOfMpi = self.numberOfMpi.get())
                        #numberOfMpi = self.numberOfMpi.get() * self.numberOfThreads.get()) 
           
            params =  '  -i %s' % fnImgsAlignN
            params += '  -o %s' % fnRootN+"_%02d_%02d.vol"%(i, iteration)
            params += ' --sym %s' % self.symmetryGroup.get()
            params += ' --max_resolution %0.3f' % self.maxRes
            params += ' --padding 2'
            params += ' --thr 1'
            #params += ' --thr %d' % self.numberOfThreads.get()
            params += ' --sampling %f' % Ts
            self.runJob('xmipp_reconstruct_fourier', params)     
            
        self.runJob('xmipp_resolution_fsc', 
                    "--ref %s -i %s -o %s --sampling_rate %f" % \
                    (fnRoot + "_00_%02d.vol" % iteration,
                     fnRoot + "_01_%02d.vol" % iteration,
                     fnRoot + "_fsc_%02d.xmd" % iteration,Ts),
                     numberOfMpi = 1)
                
        mdFSC = xmipp.MetaData(fnRoot + "_fsc_%02d.xmd" % iteration)
        for id in mdFSC:
            fscValue = mdFSC.getValue(xmipp.MDL_RESOLUTION_FRC,id)
            maxFreq = mdFSC.getValue(xmipp.MDL_RESOLUTION_FREQREAL,id)
            if fscValue < 0.5:
                break
        fh = open(fnRoot + "_freq.txt","a")
        fh.write("%f\n" % maxFreq)
        fh.close()
                
        #for noise
        self.runJob('xmipp_resolution_fsc',
                    "--ref %s -i %s -o %s --sampling_rate %f" % (
                     fnRootN + "_00_%02d.vol" % iteration,
                     fnRootN + "_01_%02d.vol" % iteration,
                     fnRootN + "_fsc_%02d.xmd" % iteration,Ts),
                     numberOfMpi = 1)
        
        cleanPattern(fnRoot + "_noises_0?_0?.xmd")
        cleanPattern(fnRoot + "_noisesOld_0?.xmd")
        cleanPattern(fnRoot + "_noisesL_0?.xmd")
        cleanPattern(fnRoot + "_noises2_0?.stk")        
        
        mdFSCN = xmipp.MetaData(fnRootN + "_fsc_%02d.xmd" % iteration)
        for id in mdFSCN:
            fscValueN = mdFSCN.getValue(xmipp.MDL_RESOLUTION_FRC, id)
            maxFreqN = mdFSCN.getValue(xmipp.MDL_RESOLUTION_FREQREAL, id)
            if fscValueN < 0.5: 
                break
        fhN = open(fnRootN + "_freq.txt", "a")
        fhN.write("%f\n" % maxFreqN)
        fhN.close()
       
        if not debugging:
            cleanPattern(fnRoot + "_0?_0?.vol")
            cleanPattern(fnRoot + "_images_0?_0?.xmd") 
            cleanPattern(fnRoot + "_fsc_0?.xmd")
            cleanPattern(fnRootN + "_0?_0?.vol")
            cleanPattern(fnRoot + "_noises_0?.stk")
            cleanPattern(fnRootN + "_fsc_0?.xmd")
            cleanPattern(fnImgsAlign + "_0?_0?.xmd")
            
    
    def gatherResultsStep(self, debugging):
        self._writeFreqsMetaData("fraction*_freq.txt", 
                                 self._defineResultsName())
        self._writeFreqsMetaData("Nfraction*_freq.txt", 
                                 self._defineResultsNoiseName())        
        
        if not debugging:
            cleanPattern(self._getExtraPath("fraction*_freq.txt"))
            cleanPattern(self._getExtraPath("Nfraction*_freq.txt"))
            cleanPattern(self._getExtraPath('Ref_Projections*'))
            cleanPattern(self._getExtraPath('newImages.stk'))
            cleanPattern(self._getExtraPath('newImages.xmd'))
            cleanPattern(self._getExtraPath('newVolume.vol'))        
    #--------------------------- INFO functions -------------------------------------------- 
    
    def _summary(self):
        """ Should be overriden in subclasses to 
        return summary message for NORMAL EXECUTION. 
        """
        numberOfParticles = getFloatListFromValues(self.numberOfParticles.get())
        maxNumberOfParticles = 0.5 * self.inputParticles.get().getSize() 
        intNumberOfParticles = [int(float(x)) for x in numberOfParticles]        
        particlesNumber = ''
        for number in intNumberOfParticles:
            if number <= maxNumberOfParticles: 
                particlesNumber += str(number)+'  '         
        msg = []
        msg.append("Number of particles: " + particlesNumber)        
        msg.append("Number of times that reconstruction is performed per each "
                   "particles subset: %d" % self.numberOfIterations)
        return msg
    
    def _methods(self):
        messages = []
        messages.append('B. Heymann "Validation of 3D EM Reconstructions"')
        return messages
    
    def _citations(self):
        return ['B.Heymann2015']
    
    def _validate(self):
        errors=[]
        if (self.doResize.get() and 
            self.newSize.get() > self.inputParticles.get().getDim()[0]):
            errors.append("Fourier resize method cannot be used "
                          "to increase the dimensions") 
        if (self.doResize.get() and 
            self.newSize.get() == self.inputParticles.get().getDim()[0]):
            errors.append("The new chosen size is equal to the "
                          "recent particles size")
        return errors 
             
    #--------------------------- UTILS functions --------------------------------------------
    def _writeFreqsMetaData(self, pattern, outputFn):
        """ Glob and read frequencies either from reconstruction
        or noise and write the proper metadata file. 
        """
        fnFreqs = sorted(glob.glob(self._getExtraPath(pattern)))
        subset = 0
        
        numberOfParticles = getFloatListFromValues(self.numberOfParticles.get())
        validationMd = xmipp.MetaData()

        for fnFreq in fnFreqs:
            print fnFreq
            data = []
            fnFreqOpen = open(fnFreq, "r")
            for line in fnFreqOpen:
                fields = line.split()
                rowdata = map(float, fields)
                data.extend(rowdata)
            meanRes = (sum(data) / len(data))
            data[:] = [(x-meanRes) ** 2 for x in data]
            varRes = (sum(data) / (len(data) - 1))
            stdRes = sqrt(varRes)
            
            objId = validationMd.addObject()
            validationMd.setValue(xmipp.MDL_COUNT, 
                                  long(numberOfParticles[subset]),
                                  objId)
            validationMd.setValue(xmipp.MDL_AVG,meanRes, objId)  
            validationMd.setValue(xmipp.MDL_STDDEV,stdRes,objId)
            subset += 1

        validationMd.write(outputFn)        
    
    def _defineResultsName(self):
        return self._getExtraPath('results.xmd')
    
    def _defineResultsNoiseName(self):
        return self._getExtraPath('resultsNoise.xmd')
