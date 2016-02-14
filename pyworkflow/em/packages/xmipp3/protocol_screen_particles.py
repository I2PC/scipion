# **************************************************************************
# *
# * Authors:     Laura del Cano (laura.cano@cnb.csic.es)
# *              Jose Gutierrez (jose.gutierrez@cnb.csic.es)
# *              I. Foche (ifoche@cnb.csic.es)
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

import pyworkflow.em.metadata as md
from pyworkflow.object import String
from pyworkflow.protocol.params import (EnumParam, IntParam, Positive, Range,
                                        LEVEL_ADVANCED, FloatParam)
from pyworkflow.em.protocol import ProtProcessParticles
from convert import writeSetOfParticles, setXmippAttributes

class XmippProtScreenParticles(ProtProcessParticles):
    """ Classify particles according their similarity to the others in order to detect outliers. """
    _label = 'screen particles'

    # Automatic Particle rejection enum
    REJ_NONE = 0
    REJ_MAXZSCORE = 1
    REJ_PERCENTAGE =2
    REJ_PERCENTAGE_SSNR =1
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineProcessParams(self, form):
        
        form.addParam('autoParRejection', EnumParam, choices=['None', 'MaxZscore', 'Percentage'],
                      label="Automatic particle rejection based on Zscore", default=self.REJ_NONE,
                      display=EnumParam.DISPLAY_COMBO, expertLevel=LEVEL_ADVANCED,
                      help='How to automatically reject particles. It can be none (no rejection), '
                      'maxZscore (reject a particle if its Zscore [a similarity index] is larger than this value), '
                      'Percentage (reject a given percentage in each one of the screening criteria). ')
        form.addParam('maxZscore', FloatParam, default=3, condition='autoParRejection==1',
                      label='Maximum Zscore', expertLevel=LEVEL_ADVANCED,
                      help='Maximum Zscore.', validators=[Positive])      
        form.addParam('percentage', IntParam, default=5, condition='autoParRejection==2',
                      label='Percentage (%)', expertLevel=LEVEL_ADVANCED,
                      help='The worse percentage of particles according to metadata labels: ZScoreShape1, ZScoreShape2, ZScoreSNR1, ZScoreSNR2, ZScoreHistogram are automatically disabled. Therefore, the total number of disabled particles belongs to [percetage, 5*percentage]', validators=[Range(0, 100, error="Percentage must be between 0 and 100.")])

        form.addParam('autoParRejectionSSNR', EnumParam, choices=['None', 'Percentage'],
                      label="Automatic particle rejection based on SSNR", default=self.REJ_NONE,
                      display=EnumParam.DISPLAY_COMBO, expertLevel=LEVEL_ADVANCED,
                      help='How to automatically reject particles. It can be none (no rejection), '
                      'Percentage (reject a given percentage of the lowest SSNRs). ')
        form.addParam('percentageSSNR', IntParam, default=5, condition='autoParRejectionSSNR==1',
                      label='Percentage (%)', expertLevel=LEVEL_ADVANCED,
                      help='The worse percentage of particles according to SSNR are automatically disabled.', 
                      validators=[Range(0, 100, error="Percentage must be between 0 and 100.")])
        form.addParallelSection(threads=0, mpi=0)
        
    def _getDefaultParallel(self):
        """This protocol doesn't have mpi version"""
        return (0, 0)
     
    #--------------------------- INSERT steps functions --------------------------------------------            
    def _insertAllSteps(self):
        """ Mainly prepare the command line for call cl2d program"""
        # Convert input images if necessary
        self._insertFunctionStep('sortImages', self.inputParticles.getObjId()) 
        self._insertFunctionStep('sortImagesSSNR', self.inputParticles.getObjId()) 
        
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions --------------------------------------------
    def sortImages(self, inputId):
        imagesMd = self._getPath('images.xmd')
        writeSetOfParticles(self.inputParticles.get(), imagesMd)
        args = "-i Particles@%s --addToInput " % imagesMd
        
        if self.autoParRejection == self.REJ_MAXZSCORE:
            args += "--zcut " + str(self.maxZscore.get())
        
        elif self.autoParRejection == self.REJ_PERCENTAGE:
            args += "--percent " + str(self.percentage.get())

        self.runJob("xmipp_image_sort_by_statistics", args)
        
        self.outputMd = String(imagesMd)

    def sortImagesSSNR(self, inputId):
        imagesMd = self._getPath('images.xmd')
        args = "-i Particles@%s " % imagesMd
        
        if self.autoParRejectionSSNR == self.REJ_PERCENTAGE_SSNR:
            args += "--ssnrpercent " + str(self.percentageSSNR.get())

        self.runJob("xmipp_image_ssnr", args)
        
    def createOutputStep(self):
        imgSet = self.inputParticles.get()
        partSet = self._createSetOfParticles()
        print "self.outputMd.get(): ", self.outputMd.get()
        partSet.copyInfo(imgSet)
        partSet.copyItems(imgSet,
                            updateItemCallback=self._updateParticle,
                            itemDataIterator=md.iterRows(self.outputMd.get(), sortByLabel=md.MDL_ITEM_ID))
        
        self._defineOutputs(outputParticles=partSet)
        self._defineSourceRelation(imgSet, partSet)

    #--------------------------- INFO functions --------------------------------------------                
    def _summary(self):
        import os
        summary = []
        if not hasattr(self, 'outputParticles'):
            summary.append("Output particles not ready yet.")
        else:
            zscores = [p._xmipp_zScore.get() for p in self.outputParticles]
            summary.append("The minimum ZScore is %.2f" % min(zscores))
            summary.append("The maximum ZScore is %.2f" % max(zscores))
            summary.append("The mean ZScore is %.2f" % (sum(zscores)*1.0/len(self.outputParticles)))
        return summary
    
    def _validate(self):
        pass
        
    def _citations(self):
        return ['Vargas2013b']
    
    def _methods(self):
        methods = []
        if hasattr(self, 'outputParticles'):
            outParticles = len(self.outputParticles) if self.outputParticles is not None else None
            particlesRejected = len(self.inputParticles.get())-outParticles if outParticles is not None else None
            particlesRejectedText = ' ('+str(particlesRejected)+')' if particlesRejected is not None else ''
            rejectionText = [
                             '',# REJ_NONE
                             ' and removing those not reaching %s%s' % (str(self.maxZscore.get()), particlesRejectedText),# REJ_MAXZSCORE
                             ' and removing worst %s percent%s' % (str(self.percentage.get()), particlesRejectedText)# REJ_PERCENTAGE
                             ]
            methods.append('Input dataset %s of %s particles was sorted by'
                           ' its ZScore using xmipp_image_sort_by_statistics'
                           ' program%s. ' % (self.getObjectTag('inputParticles'), len(self.inputParticles.get()), rejectionText[self.autoParRejection.get()]))
            methods.append('Output set is %s.'%self.getObjectTag('outputParticles'))
        return methods
    
    #--------------------------- UTILS functions -------------------------------------------- 
    def _updateParticle(self, item, row):
        setXmippAttributes(item, row, md.MDL_ZSCORE, md.MDL_ZSCORE_SHAPE1, md.MDL_ZSCORE_SHAPE2, md.MDL_ZSCORE_SNR1, md.MDL_ZSCORE_SNR2, md.MDL_CUMULATIVE_SSNR)
        if row.getValue(md.MDL_ENABLED) <= 0:
            item._appendItem = False
        else:
            item._appendItem = True