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

from pyworkflow.object import Float
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.protocol.params import PointerParam, StringParam, FloatParam, BooleanParam
from pyworkflow.utils.path import cleanPath
from pyworkflow.em.constants import ALIGN_PROJ
from pyworkflow.em.protocol import ProtAnalysis2D
from pyworkflow.em.data import Class2D, SetOfClasses2D
from pyworkflow.em.packages.xmipp3.utils import iterMdRows
from pyworkflow.em.packages.xmipp3.convert import rowToAlignment
from math import floor

import xmipp
from xmipp3 import ProjMatcher

        
class XmippProtProjectionOutliers(ProtAnalysis2D, ProjMatcher):
    """Compares a set of classes or averages with the corresponding projections of a reference volume """
    _label = 'projection outliers'
    
    def __init__(self, **args):
        ProtAnalysis2D.__init__(self, **args)
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSet', PointerParam, label="Input averages", important=True, 
                      pointerClass='SetOfParticles', pointerCondition='hasAlignmentProj')
        form.addParam('inputVolume', PointerParam, label="Volume to compare classes to", important=True,
                      pointerClass='Volume',
                      help='Volume to be used for class comparison')
        form.addParallelSection(threads=0, mpi=1)
    
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        # Convert input images if necessary
        self.imgsFn = self._getExtraPath('input_imgs.xmd')
        vol = self.inputVolume.get()
        
        self._insertFunctionStep("convertStep", self.imgsFn)
        self._insertFunctionStep("produceResiduals", vol.getFileName(), self.imgsFn, vol.getSamplingRate())
        self._insertFunctionStep("evaluateResiduals", self.imgsFn)
        #self._insertFunctionStep("createOutputStep", self.imgsFn)

    #--------------------------- STEPS functions ---------------------------------------------------
    def convertStep(self, imgsFn):
        from convert import writeSetOfParticles
        writeSetOfParticles(self.inputSet.get(), self.imgsFn)
    
    def produceResiduals(self,volumeFn,anglesFn,Ts):
        if volumeFn.endswith(".mrc"):
            volumeFn+=":mrc"
        anglesOutFn=self._getExtraPath("anglesCont.stk")
        residualsOutFn=self._getExtraPath("residuals.stk")
        projectionsOutFn=self._getExtraPath("projections.stk")
        xdim=self.inputVolume.get().getDim()[0]
        self.runJob("xmipp_angular_continuous_assign2", "-i %s -o %s --ref %s --optimizeAngles --optimizeGray --optimizeShift --max_shift %d --oresiduals %s --oprojections %s --sampling %f" %\
                    (anglesFn,anglesOutFn,volumeFn,floor(xdim*0.05),residualsOutFn,projectionsOutFn,Ts))
    
    def evaluateResiduals(self, outImgsFn):
        # Evaluate each image
        fnAutoCorrelations = self._getExtraPath("autocorrelations.xmd")
        stkAutoCorrelations = self._getExtraPath("autocorrelations.stk")
        stkResiduals = self._getExtraPath("residuals.stk")
        anglesOutFn=self._getExtraPath("anglesCont.xmd")
        self.runJob("xmipp_image_residuals", " -i %s -o %s --save_metadata_stack %s" % (stkResiduals, stkAutoCorrelations, fnAutoCorrelations), numberOfMpi=1)
        self.runJob("xmipp_metadata_utilities", '-i %s --operate rename_column "image imageResidual"' % fnAutoCorrelations, numberOfMpi=1)
        self.runJob("xmipp_metadata_utilities", '-i %s --set join %s imageResidual' % (anglesOutFn, fnAutoCorrelations), numberOfMpi=1)
        cleanPath(fnAutoCorrelations)
    
    def createOutputStep(self, outImgsFn):
        inputSet = self.inputSet.get()
        if isinstance(inputSet, SetOfClasses2D):
            outputSet = self._createSetOfClasses2D(inputSet.getImages())
            outputName = 'outputClasses'
        else: # SetOfAverages
            outputSet = self._createSetOfAverages()
            outputSet.setAlignment3D()
            outputName = 'outputAverages'
            
        md = xmipp.MetaData(outImgsFn)
        outputSet.copyInfo(inputSet)
        outputSet.copyItems(inputSet, 
                            updateItemCallback=self.updateItem,
                            itemDataIterator=iterMdRows(md))

        self._defineOutputs(**{outputName: outputSet})
        self._defineTransformRelation(inputSet, outputSet)
    
    #--------------------------- INFO functions --------------------------------------------
    def _validate(self):
        errors = []
        vol = self.inputVolume.get()
        xDim = self._getDimensions()
        volDim = vol.getDim()[0]
        
        if volDim != xDim:
            errors.append("Make sure that the volume and the images have the same size")
        return errors    
    
    def _summary(self):
        summary = []
        summary.append("Images evaluated: %i" % self.inputSet.get().getSize())
        summary.append("Volume: %s" % self.inputVolume.getNameId())
        summary.append("symmetry: %s" % self.symmetryGroup.get())
        return summary
    
    def _methods(self):
        methods = []
        if hasattr(self, 'outputClasses') or hasattr(self, 'outputAverages'):
            methods.append("We evaluated %i input images %s regarding to volume %s"
                           " using %s symmetry" %(self.inputSet.get().getSize(), self.getObjectTag('inputSet'), \
                                                  self.getObjectTag('inputVolume'), self.symmetryGroup.get()) )
        return methods
    
    #--------------------------- UTILS functions --------------------------------------------
    def _getDimensions(self):
        imgSet = self.inputSet.get()
        if isinstance(imgSet, SetOfClasses2D):
            xDim = imgSet.getImages().getDim()[0]
        else:
            xDim = imgSet.getDim()[0]
        return xDim
    
    def updateItem(self, item, row):
        from convert import locationToXmipp
        # ToDo: uncomment this lines when the output metadata has ITEM_ID
#         if item.getObjId() != row.getValue(xmipp.MDL_ITEM_ID):
#             raise Exception("The objId is not equal to ITEM_ID. Please, sort the metadata.")
        if isinstance(item, Class2D):
            img = item.getRepresentative()
            index, fn = img.getLocation()
        else:
            index, fn = item.getLocation()
            
        objLoc = locationToXmipp(index, fn)
        mdLoc = row.getValue(xmipp.MDL_IMAGE)
        if objLoc != mdLoc:
            print objLoc+" "+mdLoc
            raise Exception("The image isn't the same. Please, sort the metadata.")

        item._xmipp_maxCC = Float(row.getValue(xmipp.MDL_MAXCC))
        item._xmipp_zScoreResCov = Float(row.getValue(xmipp.MDL_ZSCORE_RESCOV))
        item._xmipp_zScoreResMean = Float(row.getValue(xmipp.MDL_ZSCORE_RESMEAN))
        item._xmipp_zScoreResVar = Float(row.getValue(xmipp.MDL_ZSCORE_RESVAR))
        if isinstance(item, Class2D):
            particle = item.getRepresentative()
        else:
            particle = item
        particle.setTransform(rowToAlignment(row, alignType=ALIGN_PROJ))
