# **************************************************************************
# *
# * Authors:     Amaya Jimenez (ajimenez@cnb.csic.es)
# *              Javier Mota (jmota@cnb.csic.es)
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

from math import floor
import os

from pyworkflow import VERSION_1_1
from pyworkflow.protocol.params import PointerParam, StringParam, FloatParam, BooleanParam
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.em.constants import ALIGN_PROJ
from pyworkflow.utils.path import cleanPath
from pyworkflow.em.protocol import ProtAnalysis3D
from pyworkflow.em.data import SetOfClasses2D, Image, SetOfAverages, SetOfParticles, Class2D
from pyworkflow.em.packages.xmipp3.convert import setXmippAttributes, xmippToLocation
import pyworkflow.em.metadata as md
from pyworkflow.protocol.constants import LEVEL_ADVANCED

from shutil import copy
import xmipp
from xmipp3 import ProjMatcher
from pyworkflow.em.packages.xmipp3.convert import rowToAlignment

        
class XmippProtGenerateReprojections(ProtAnalysis3D):
    """Compares a set of classes or averages with the corresponding projections of a reference volume.
    The set of images must have a 3D angular assignment and the protocol computes the residues
    (the difference between the experimental images and the reprojections). The zscore of the mean
    and variance of the residues are computed. Large values of these scores may indicate outliers.
    The protocol also analyze the covariance matrix of the residual and computes the logarithm of
    its determinant [Cherian2013]. The extremes of this score (called zScoreResCov), that is
    values particularly low or high, may indicate outliers."""
    _label = 'generate reprojections'
    _lastUpdateVersion = VERSION_1_1
    
    def __init__(self, **args):
        ProtAnalysis3D.__init__(self, **args)
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSet', PointerParam, label="Input images", important=True, 
                      pointerClass='SetOfParticles')
        form.addParam('inputVolume', PointerParam, label="Volume to compare images to", important=True,
                      pointerClass='Volume',
                      help='Volume to be used for class comparison')
        form.addParallelSection(threads=0, mpi=8)
    
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        # Convert input images if necessary
        self.imgsFn = self._getExtraPath('input_imgs.xmd')
        vol = self.inputVolume.get()
        
        self._insertFunctionStep("convertStep", self.imgsFn)
        imgSet = self.inputSet.get()
        anglesFn=self.imgsFn
        self._insertFunctionStep("produceProjections", vol.getFileName(),
                                 anglesFn,
                                 vol.getSamplingRate())
        self._insertFunctionStep("createOutputStep")

    #--------------------------- STEPS functions ---------------------------------------------------
    def convertStep(self, imgsFn):
        from convert import writeSetOfParticles
        imgSet = self.inputSet.get()
        writeSetOfParticles(imgSet, self.imgsFn)

        from pyworkflow.em.convert import ImageHandler
        img = ImageHandler()
        fnVol = self._getTmpPath("volume.vol")
        img.convert(self.inputVolume.get(), fnVol)
        xdim=self.inputVolume.get().getDim()[0]

        imgXdim=imgSet.getDim()[0]
        if xdim!=imgXdim:
            self.runJob("xmipp_image_resize","-i %s --dim %d"%(fnVol,imgXdim),numberOfMpi=1)
    
    def produceProjections(self, fnVol, fnAngles, Ts):
        fnVol = self._getTmpPath("volume.vol")
        anglesOutFn=self._getExtraPath("anglesCont.stk")
        projectionsOutFn=self._getExtraPath("projections.stk")
        args="-i %s -o %s --ref %s --oprojections %s --sampling %f " \
             "--max_angular_change 90"%(fnAngles,anglesOutFn,fnVol,
                                  projectionsOutFn,Ts)
        self.runJob("xmipp_angular_continuous_assign2", args)
        fnNewParticles=self._getExtraPath("images.stk")
        if os.path.exists(fnNewParticles):
            cleanPath(fnNewParticles)
    
    def createOutputStep(self):

        fnImgs = self._getExtraPath('images.stk')
        if os.path.exists(fnImgs):
            cleanPath(fnImgs)

        imgSet = self.inputSet.get()
        imgFn = self._getExtraPath("anglesCont.xmd")
        self.newAssignmentPerformed = os.path.exists(self._getExtraPath("angles.xmd"))
        self.samplingRate = self.inputSet.get().getSamplingRate()
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

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        summary = []
        summary.append("Images evaluated: %i" % self.inputSet.get().getSize())
        summary.append("Volume: %s" % self.inputVolume.getNameId())
        return summary

    def _methods(self):
        methods = []
        if hasattr(self, 'outputParticles'):
            methods.append(
                "We projected the volume %s following the directions in %i input images %s." \
                % (
                self.getObjectTag('inputVolume'), self.inputSet.get().getSize(),
                self.getObjectTag('inputSet')))
        return methods

