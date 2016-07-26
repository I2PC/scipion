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

from math import floor
import os

from pyworkflow.protocol.params import PointerParam, StringParam, FloatParam, BooleanParam
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.em.constants import ALIGN_PROJ
from pyworkflow.utils.path import cleanPath
from pyworkflow.em.protocol import ProtAnalysis3D
from pyworkflow.em.data import SetOfClasses2D, Image, SetOfAverages, SetOfParticles, Class2D
from pyworkflow.em.packages.xmipp3.convert import setXmippAttributes, xmippToLocation
import pyworkflow.em.metadata as md

import xmipp
from xmipp3 import ProjMatcher
from pyworkflow.em.packages.xmipp3.convert import rowToAlignment

        
class XmippProtCompareReprojections(ProtAnalysis3D, ProjMatcher):
    """Compares a set of classes or averages with the corresponding projections
     of a reference volume.
    """
    _label = 'compare reprojections'
    
    def __init__(self, **args):
        ProtAnalysis3D.__init__(self, **args)
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSet', PointerParam, label="Input images", important=True, 
                      pointerClass='SetOfClasses2D, SetOfAverages, SetOfParticles')
        form.addParam('inputVolume', PointerParam, label="Volume to compare images to", important=True,
                      pointerClass='Volume',
                      help='Volume to be used for class comparison')
        
        form.addParam('useAssignment', BooleanParam, default=True,
                      label='Use input angular assignment (if available)')
        form.addParam('symmetryGroup', StringParam, default="c1",
                      label='Symmetry group', 
                      help='See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Symmetry for a description of the symmetry groups format'
                        'If no symmetry is present, give c1')
        form.addParam('angularSampling', FloatParam, default=5, expertLevel=LEVEL_ADVANCED,
                      label='Angular sampling rate',
                      help='In degrees.'
                      ' This sampling defines how fine the projection gallery from the volume is explored.')
        form.addParallelSection(threads=0, mpi=8)
    
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        # Convert input images if necessary
        self.imgsFn = self._getExtraPath('input_imgs.xmd')
        vol = self.inputVolume.get()
        
        self._insertFunctionStep("convertStep", self.imgsFn)
        imgSet = self.inputSet.get()
        if not self.useAssignment or isinstance(imgSet, SetOfClasses2D) or isinstance(imgSet, SetOfAverages) or (isinstance(imgSet, SetOfParticles) and not imgSet.hasAlignmentProj()):
            anglesFn = self._getExtraPath('angles.xmd')
            self._insertFunctionStep("projMatchStep", self.inputVolume.get().getFileName(), self.angularSampling.get(), self.symmetryGroup.get(), self.imgsFn,
                                     anglesFn, self.inputVolume.get().getDim()[0])
        else:
            anglesFn=self.imgsFn
        self._insertFunctionStep("produceResiduals", vol.getFileName(), anglesFn, vol.getSamplingRate())
        self._insertFunctionStep("evaluateResiduals")
        self._insertFunctionStep("createOutputStep")

    #--------------------------- STEPS functions ---------------------------------------------------
    def convertStep(self, imgsFn):
        from convert import writeSetOfClasses2D, writeSetOfParticles
        imgSet = self.inputSet.get()
        if isinstance(imgSet, SetOfClasses2D):
            writeSetOfClasses2D(imgSet, self.imgsFn, writeParticles=True)
        else:
            writeSetOfParticles(imgSet, self.imgsFn)
    
    def produceResiduals(self,volumeFn,anglesFn,Ts):
        if volumeFn.endswith(".mrc"):
            volumeFn+=":mrc"
        anglesOutFn=self._getExtraPath("anglesCont.stk")
        residualsOutFn=self._getExtraPath("residuals.stk")
        projectionsOutFn=self._getExtraPath("projections.stk")
        xdim=self.inputVolume.get().getDim()[0]
        self.runJob("xmipp_angular_continuous_assign2", "-i %s -o %s --ref %s --optimizeAngles --optimizeGray --optimizeShift --max_shift %d --oresiduals %s --oprojections %s --sampling %f" %\
                    (anglesFn,anglesOutFn,volumeFn,floor(xdim*0.05),residualsOutFn,projectionsOutFn,Ts))
    
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
        outputSet = self._createSetOfParticles()
        imgSet = self.inputSet.get()
        imgFn = self._getExtraPath("anglesCont.xmd")
        self.newAssignmentPerformed = os.path.exists(self._getExtraPath("angles.xmd"))
        self.samplingRate = self.inputSet.get().getSamplingRate()
        if isinstance(imgSet, SetOfClasses2D):
            outputSet = self._createSetOfClasses2D(imgSet)
            outputSet.copyInfo(imgSet.getImages())
        elif isinstance(imgSet, SetOfAverages):
            outputSet = self._createSetOfAverages()
            outputSet.copyInfo(imgSet)
        else:
            outputSet = self._createSetOfParticles()
            outputSet.copyInfo(imgSet)
            if not self.newAssignmentPerformed:
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
