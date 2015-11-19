# **************************************************************************
# *
# * Authors:     C.O.S. Sorzano (coss@cnb.csic.es)
# *              Mohsen Kazemi  (mkazemi@cnb.csic.es)
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
from pyworkflow.utils.path import cleanPattern, cleanPath
import os
import xmipp
import glob
from pyworkflow.object import Float, String
from math import sqrt
from plotter import XmippPlotter


class XmippProtValidateOverfitting(ProtReconstruct3D):
    """    
    Check how the FSC changes with the number of projections used for 3D reconstruction. This method
    has been proposed by B. Heymann at ***
    """
    _label = 'validate overfitting'
    
    #--------------------------- DEFINE param functions --------------------------------------------   
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputParticles', PointerParam, pointerClass='SetOfParticles',pointerCondition='hasAlignmentProj',
                      label="Input particles",  
                      help='Select the input images from the project.')     
        form.addParam('symmetryGroup', StringParam, default='c1',
                      label="Symmetry group", 
                      help='See [[Xmipp Symmetry][http://www2.mrc-lmb.cam.ac.uk/Xmipp/index.php/Conventions_%26_File_formats#Symmetry]] page '
                           'for a description of the symmetry format accepted by Xmipp') 
        form.addParam('numberOfParticles', NumericListParam, default="100 200 500 1000 2000 5000", expertLevel=LEVEL_ADVANCED,
                      label="Number of particles") 
        form.addParam('numberOfIterations', IntParam, default=10, expertLevel=LEVEL_ADVANCED,
                      label="Number of times the randomization is performed.") 
        form.addParam('maxRes', FloatParam, default=-1, expertLevel=LEVEL_ADVANCED,
                      label="Maximum resolution (A)",  
                      help='Maximum resolution (in Angstrom) to consider \n'
                           'in Fourier space (default Nyquist).\n'
                           'Param *--maxres* in Xmipp.') 
        form.addParam('pad', FloatParam, default=2, expertLevel=LEVEL_ADVANCED,
                      label="Padding factor")

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
        self._insertFunctionStep('convertInputStep')
        numberOfParticles=getFloatListFromValues(self.numberOfParticles.get())
        fractionCounter=0
        for number in numberOfParticles:
            if number<self.inputParticles.get().getSize():
                for iteration in range(0,self.numberOfIterations.get()):
                    self._insertFunctionStep('reconstructionStep',number,fractionCounter,iteration)
                fractionCounter+=1     
        self._insertFunctionStep('gatherResultsStep')
        #self._insertFunctionStep('createOutputStep')
        
    
    #--------------------------- STEPS functions --------------------------------------------
    def convertInputStep(self):
        particlesMd = self._getFileName('input_xmd')
        imgSet = self.inputParticles.get()
        writeSetOfParticles(imgSet, particlesMd)

    def reconstructionStep(self, numberOfImages, fractionCounter, iteration):
        fnRoot = self._getExtraPath("fraction%02d"%fractionCounter)
        Ts = self.inputParticles.get().getSamplingRate()
        for i in range(0,2):
            fnImgs=fnRoot+"_images_%02d.xmd"%i
            self.runJob("xmipp_metadata_utilities","-i %s -o %s --operate random_subset %d"%\
                        (self._getFileName('input_xmd'),fnImgs,numberOfImages),numberOfMpi=1)
        
            params =  '  -i %s' % fnImgs
            params += '  -o %s' % fnRoot+"_%02d.vol"%i
            params += ' --sym %s' % self.symmetryGroup.get()
            params += ' --max_resolution %0.3f' % self.maxRes.get()
            params += ' --padding %0.3f' % self.pad.get()
            params += ' --thr %d' % self.numberOfThreads.get()
            params += ' --sampling %f' % Ts
            self.runJob('xmipp_reconstruct_fourier', params)

        self.runJob('xmipp_resolution_fsc', "--ref %s -i %s -o %s --sampling_rate %f"%\
                    (fnRoot+"_00.vol",fnRoot+"_01.vol",fnRoot+"_fsc.xmd",Ts), numberOfMpi=1)
        cleanPattern(fnRoot+"_0?.vol")
        cleanPattern(fnRoot+"_images_0?.xmd")
        
        mdFSC = xmipp.MetaData(fnRoot+"_fsc.xmd")
        for id in mdFSC:
            fscValue = mdFSC.getValue(xmipp.MDL_RESOLUTION_FRC,id)
            maxFreq = mdFSC.getValue(xmipp.MDL_RESOLUTION_FREQREAL,id)
            if fscValue<0.5:
                break
        fh = open(fnRoot+"_freq.txt","a")
        fh.write("%f\n"%maxFreq)
        fh.close()
        cleanPath(fnRoot+"_fsc.xmd")
        
    def gatherResultsStep(self):
        fnFreqs = sorted(glob.glob(self._getExtraPath("fraction*_freq.txt")))
        subset = 0
        
        numberOfParticles=getFloatListFromValues(self.numberOfParticles.get())
        validationMd = xmipp.MetaData()

        for fnFreq in fnFreqs:
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

        validationMd.write(self._getExtraPath('results.xmd'))
        #cleanPattern(self._getExtraPath("fraction*_freq.txt"))
        
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
