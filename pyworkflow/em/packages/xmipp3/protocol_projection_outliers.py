# **************************************************************************
# *
# * Authors:     Josue Gomez Blanco (jgomez@cnb.csic.es)
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
This sub-package contains wrapper around Projection Outliers Xmipp program
"""

from math import floor
from os.path import exists

from pyworkflow.protocol.params import IntParam, PointerParam
from pyworkflow.utils.path import cleanPath, createLink
from pyworkflow.em.protocol import ProtAnalysis2D
from pyworkflow.em.data import SetOfClasses2D, Image
from pyworkflow.em.packages.xmipp3.convert import setXmippAttributes, xmippToLocation
import pyworkflow.em.metadata as md
from pyworkflow.protocol.constants import LEVEL_ADVANCED

import xmipp
from xmipp3 import ProjMatcher

        
class XmippProtProjectionOutliers(ProtAnalysis2D, ProjMatcher):
    """Compares a set of classes or averages with the corresponding projections of a reference volume.
    The set of images must have a 3D angular assignment and the protocol computes the residues
    (the difference between the experimental images and the reprojections). The zscore of the mean
    and variance of the residues are computed. Large values of these scores may indicate outliers.
    The protocol also analyze the covariance matrix of the residual and computes the logarithm of
    its determinant [Cherian2013]. The extremes of this score (called zScoreResCov), that is
    values particularly low or high, may indicate outliers."""
    _label = 'projection outliers'
    
    def __init__(self, **args):
        ProtAnalysis2D.__init__(self, **args)
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSet', PointerParam, label="Input images", important=True, 
                      pointerClass='SetOfParticles', pointerCondition='hasAlignmentProj')
        form.addParam('inputVolume', PointerParam, label="Volume to compare images to", important=True,
                      pointerClass='Volume',
                      help='Volume to be used for class comparison')
        form.addParam('residualSize', IntParam, label="Size of the residuals to analyze", expertLevel=LEVEL_ADVANCED, default=64,
                      help='The input volume and the input images will be resized to this size before analyzing the residuals.')
        form.addParallelSection(threads=0, mpi=1)
    
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep("convertStep")
        self._insertFunctionStep("produceResiduals")
        self._insertFunctionStep("evaluateResiduals")
        self._insertFunctionStep("createOutputStep")

    #--------------------------- STEPS functions ---------------------------------------------------
    def convertStep(self):
        from convert import writeSetOfParticles
        fnImgs = self._getExtraPath('images.xmd')
        writeSetOfParticles(self.inputSet.get(), fnImgs)
        
        currentDim=self._getDimensions()
        if currentDim!=self.residualSize.get():
            self.runJob("xmipp_image_resize","-i %s -o %s --fourier %d --save_metadata_stack %s"%(fnImgs,self._getExtraPath("images.stk"),self.residualSize.get(),fnImgs))
            scaleFactor=float(self.residualSize.get())/currentDim
            self.runJob("xmipp_metadata_utilities",'-i %s --operate modify_values "shiftX=shiftX*%f"'%(fnImgs,scaleFactor),numberOfMpi=1)
            self.runJob("xmipp_metadata_utilities",'-i %s --operate modify_values "shiftY=shiftY*%f"'%(fnImgs,scaleFactor),numberOfMpi=1)
        
        fnVol=self.inputVolume.get().getFileName()
        if self.inputVolume.get().getDim()[0]!=self.residualSize.get() or fnVol.endswith('.mrc'):
            if fnVol.endswith(".mrc"):
                fnVol+=":mrc"
            self.runJob("xmipp_image_resize","-i %s -o %s --fourier %d"%(fnVol,self._getExtraPath("volume.vol"),self.residualSize.get()))
        else:
            createLink(fnVol,self._getExtraPath("volume.vol"))
    
    def produceResiduals(self):
        anglesOutFn=self._getExtraPath("anglesCont.stk")
        residualsOutFn=self._getExtraPath("residuals.stk")
        projectionsOutFn=self._getExtraPath("projections.stk")
        xdim=self.inputVolume.get().getDim()[0]
        originalTs=self.inputVolume.get().getSamplingRate()
        Ts=originalTs*xdim/float(self.residualSize.get())
        self.runJob("xmipp_angular_continuous_assign2", "-i %s -o %s --ref %s --optimizeAngles --optimizeGray --optimizeShift --max_shift %d --oresiduals %s --oprojections %s --sampling %f" %\
                    (self._getExtraPath('images.xmd'),anglesOutFn,self._getExtraPath("volume.vol"),floor(xdim*0.05),residualsOutFn,projectionsOutFn,Ts))
        fnNewParticles=self._getExtraPath("images.stk")
        if exists(fnNewParticles):
            cleanPath(fnNewParticles)
    
    def evaluateResiduals(self):
        # Evaluate each image
        fnAutoCorrelations = self._getExtraPath("autocorrelations.xmd")
        stkAutoCorrelations = self._getExtraPath("autocorrelations.stk")
        stkResiduals = self._getExtraPath("residuals.stk")
        anglesOutFn=self._getExtraPath("anglesCont.xmd")
        self.runJob("xmipp_image_residuals", " -i %s -o %s --save_metadata_stack %s" % (stkResiduals, stkAutoCorrelations, fnAutoCorrelations), numberOfMpi=1)
        self.runJob("xmipp_metadata_utilities", '-i %s --operate rename_column "image imageResidual"' % fnAutoCorrelations, numberOfMpi=1)
        self.runJob("xmipp_metadata_utilities", '-i %s --set join %s imageResidual' % (anglesOutFn, fnAutoCorrelations), numberOfMpi=1)
        cleanPath(fnAutoCorrelations)
    
    def createOutputStep(self):
        fnImgs = self._getExtraPath('images.stk')
        if exists(fnImgs):
            cleanPath(fnImgs)

        outputSet = self._createSetOfParticles()
        imgSet = self.inputSet.get()
        imgFn = self._getExtraPath("anglesCont.xmd")
        outputSet.copyInfo(imgSet)
        outputSet.setAlignmentProj()
        outputSet.copyItems(imgSet,
                            updateItemCallback=self._processRow,
                            itemDataIterator=md.iterRows(imgFn, sortByLabel=md.MDL_ITEM_ID))
        self._defineOutputs(outputParticles=outputSet)
        self._defineSourceRelation(self.inputSet, outputSet)

    def _processRow(self, particle, row):
        setXmippAttributes(particle, row,
                           xmipp.MDL_ZSCORE_RESVAR, xmipp.MDL_ZSCORE_RESMEAN, xmipp.MDL_ZSCORE_RESCOV, xmipp.MDL_IMAGE_ORIGINAL,
                           xmipp.MDL_COST, xmipp.MDL_CONTINUOUS_GRAY_A, xmipp.MDL_CONTINUOUS_GRAY_B, xmipp.MDL_CONTINUOUS_X, xmipp.MDL_CONTINUOUS_Y)
        def __setXmippImage(label):
            attr = '_xmipp_' + xmipp.label2Str(label)
            if not hasattr(particle, attr):
                img = Image()
                setattr(particle, attr, img)
                img.setSamplingRate(particle.getSamplingRate())
            else:
                img = getattr(particle, attr)
            img.setLocation(xmippToLocation(row.getValue(label)))
        
        __setXmippImage(xmipp.MDL_IMAGE)
        __setXmippImage(xmipp.MDL_IMAGE_REF)
        __setXmippImage(xmipp.MDL_IMAGE_RESIDUAL)
        __setXmippImage(xmipp.MDL_IMAGE_COVARIANCE)

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        summary = []
        summary.append("Images evaluated: %i" % self.inputSet.get().getSize())
        summary.append("Volume: %s" % self.inputVolume.getNameId())
        return summary
    
    def _methods(self):
        methods = []
        if hasattr(self, 'outputParticles'):
            methods.append("We evaluated %i input images %s regarding to volume %s."\
                           %(self.inputSet.get().getSize(), self.getObjectTag('inputSet'), self.getObjectTag('inputVolume')) )
            methods.append("The residuals were evaluated according to their mean, variance and covariance structure [Cherian2013].")
        return methods
    
    #--------------------------- UTILS functions --------------------------------------------
    def _getDimensions(self):
        imgSet = self.inputSet.get()
        if isinstance(imgSet, SetOfClasses2D):
            xDim = imgSet.getImages().getDim()[0]
        else:
            xDim = imgSet.getDim()[0]
        return xDim
