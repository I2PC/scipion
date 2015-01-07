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
from pyworkflow.protocol.constants import LEVEL_EXPERT
from pyworkflow.protocol.params import PointerParam, StringParam, FloatParam, BooleanParam
from pyworkflow.utils.path import cleanPath
from pyworkflow.em.constants import ALIGN_PROJ
from pyworkflow.em.protocol import ProtAnalysis2D
from pyworkflow.em.data import Class2D, SetOfClasses2D
from pyworkflow.em.packages.xmipp3.utils import iterMdRows
from pyworkflow.em.packages.xmipp3.convert import rowToAlignment

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
                      pointerClass='SetOfClasses2D, SetOfAverages',
                      help='Select the input set to compare.'
                           'It should be a SetOfClasses2D or a SetOfAverages')
        form.addParam('inputVolume', PointerParam, label="Volume to compare classes to", important=True,
                      pointerClass='Volume',
                      help='Volume to be used for class comparison')
        form.addParam('hasVolCTFCorrected', BooleanParam, default=False,
                      label='Volume has been corrected by CTF?')        

        form.addParam('symmetryGroup', StringParam, default="c1",
                      label='Symmetry group', 
                      help='See http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Symmetry for a description of the symmetry groups format'
                        'If no symmetry is present, give c1')
        form.addParam('angularSampling', FloatParam, default=5, expertLevel=LEVEL_EXPERT,
                      label='Angular sampling rate',
                      help='In degrees.'
                      ' This sampling defines how fine the projection gallery from the volume is explored.')
        
        form.addParallelSection(threads=0, mpi=1)
    
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        """ Mainly prepare the command line for call cl2d program"""
        # Convert input images if necessary
        imgsFn = self._getPath('input_imgs.xmd')
        outImgsFn = self._getExtraPath('output_imgs.xmd')
        anglesFn = self._getExtraPath('angles.xmd')
        vol = self.inputVolume.get()
        
        angSampling = self.angularSampling.get()
        sym = self.symmetryGroup.get()
        
        self._insertFunctionStep("convertStep", imgsFn)
        
        self._insertFunctionStep("projMatchStep", vol.getFileName(), angSampling, sym, imgsFn, anglesFn, self._getDimensions())
        self._insertFunctionStep("joinStep", imgsFn, anglesFn)
        self._insertFunctionStep("produceAlignedImagesStep", self.hasVolCTFCorrected.get(), outImgsFn, imgsFn)
        self._insertFunctionStep("evaluateStep", outImgsFn)
#         self._insertFunctionStep("sortStep", outImgsFn)
        self._insertFunctionStep("createOutputStep", outImgsFn)

    #--------------------------- STEPS functions ---------------------------------------------------
    def convertStep(self, imgsFn):
        from convert import writeSetOfClasses2D, writeSetOfParticles
        imgSet = self.inputSet.get()
        if isinstance(imgSet, SetOfClasses2D):
            writeSetOfClasses2D(imgSet, imgsFn, writeParticles=False)
        else:
            writeSetOfParticles(imgSet, imgsFn)
    
    def joinStep(self, imgsFn, anglesFn):
        imgSet = self.inputSet.get()
        if isinstance(imgSet, SetOfClasses2D):
            self.runJob("xmipp_metadata_utilities", "-i classes@%s --set join %s --mode append" % (imgsFn, anglesFn), numberOfMpi=1)
        else:
            self.runJob("xmipp_metadata_utilities", "-i Particles@%s --set join %s --mode append" % (imgsFn, anglesFn), numberOfMpi=1)
    
    def evaluateStep(self, outImgsFn):
        # Evaluate each image
        fnAutoCorrelations = self._getExtraPath("autocorrelations.xmd")
        stkAutoCorrelations = self._getExtraPath("autocorrelations.stk")
        stkDiff = self._getExtraPath("diff.stk")
        args1 = " -i %s -o %s --save_metadata_stack %s"
        args2 = " -i %s --set merge %s"
        outClasses = 'classes_aligned@' + outImgsFn
        self.runJob("xmipp_image_residuals", args1 % (stkDiff, stkAutoCorrelations, fnAutoCorrelations), numberOfMpi=1)
        self.runJob("xmipp_metadata_utilities", args2 % (outClasses, fnAutoCorrelations), numberOfMpi=1)
        cleanPath(fnAutoCorrelations)
    
    def sortStep(self, outImgsFn):
        self.runJob("xmipp_metadata_utilities", "-i classes_aligned@%s --operate sort zScoreResCov desc" % (outImgsFn), numberOfMpi=1)
    
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
            methods.append("We evaluated %i images regarding to volume %s"
                           " using %s symmetry" %(self.inputSet.get().getSize(),\
                                                  self.inputVolume.getNameId(), self.symmetryGroup.get()) )
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
            raise Exception("The the image isn't the same. Please, sort the metadata.")

        item._xmipp_maxCC = Float(row.getValue(xmipp.MDL_MAXCC))
        item._xmipp_zScoreResCov = Float(row.getValue(xmipp.MDL_ZSCORE_RESCOV))
        item._xmipp_zScoreResMean = Float(row.getValue(xmipp.MDL_ZSCORE_RESMEAN))
        item._xmipp_zScoreResVar = Float(row.getValue(xmipp.MDL_ZSCORE_RESVAR))
        if isinstance(item, Class2D):
            particle = item.getRepresentative()
        else:
            particle = item
        particle.setTransform(rowToAlignment(row, alignType=ALIGN_PROJ))
