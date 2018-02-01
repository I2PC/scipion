# **************************************************************************
# *
# * Authors:     C.O.S. Sorzano (coss@cnb.csic.es)
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

from os.path import join, exists
import math

import pyworkflow.protocol.params as params
from pyworkflow import VERSION_1_2
from pyworkflow.utils.path import makePath
from pyworkflow.em.convert import ImageHandler, ALIGN_PROJ
from pyworkflow.em.data import Image
from pyworkflow.em.protocol import ProtAnalysis3D
from pyworkflow.em.packages.xmipp3.convert import writeSetOfParticles
import pyworkflow.em.metadata as md

import xmipp
from convert import rowToAlignment, setXmippAttributes, xmippToLocation
from xmipp3 import findRow
from constants import SYM_URL


class XmippProtCompareAngles(ProtAnalysis3D):
    """    
    Compare two sets of angles. The output is a list of all common particles with
    the angular difference between both assignments. The output is constructed by 
    keeping the information from the Set 1 and adding the shiftDiff and angularDiff.
    """

    _label = 'compare angles'
    _lastUpdateVersion = VERSION_1_2
    
    def __init__(self, *args, **kwargs):
        ProtAnalysis3D.__init__(self, *args, **kwargs)
        
    #--------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        
        form.addParam('inputParticles1', params.PointerParam,
                      pointerClass='SetOfParticles',
                      pointerCondition='hasAlignmentProj',
                      label="Input particles 1",  
                      help='Select the input experimental images with an '
                           'angular assignment.')

        form.addParam('inputParticles2', params.PointerParam,
                      pointerClass='SetOfParticles',
                      pointerCondition='hasAlignmentProj',
                      label="Input particles 2",  
                      help='Select the input experimental images with an '
                           'angular assignment.')

        form.addParam('symmetryGroup', params.StringParam, default='c1',
                      label="Symmetry group", 
                      help='See %s page for a description of the symmetries '
                           'accepted by Xmipp' % SYM_URL)

    
    #--------------------------- INSERT steps functions ------------------------

    def _insertAllSteps(self):        
        self._insertFunctionStep('convertInputStep',
                                 self.inputParticles1.get().getObjId(),
                                 self.inputParticles2.get().getObjId())

        self._insertFunctionStep('analyzeDistanceStep',
                                 self.inputParticles1.get().getObjId(),
                                 self.inputParticles2.get().getObjId(),
                                 self.symmetryGroup.get())

        self._insertFunctionStep('createOutputStep')
    
    #--------------------------- STEPS functions -------------------------------

    def convertInputStep(self, particlesId1, particlesId2):
        """ Write the input images as a Xmipp metadata file. 
        particlesId: is only need to detect changes in
        input particles and cause restart from here.
        """
        writeSetOfParticles(self.inputParticles1.get(), self._getExtraPath("angles1.xmd"))
        writeSetOfParticles(self.inputParticles2.get(), self._getExtraPath("angles2.xmd"))

    def analyzeDistanceStep(self, particlesId1, particlesId2, symmetryGroup):
        self.runJob("xmipp_metadata_utilities","-i %s -o %s --operate keep_column itemId"%\
                    (self._getExtraPath("angles1.xmd"),self._getTmpPath("ids1.xmd")))
        self.runJob("xmipp_metadata_utilities","-i %s -o %s --operate keep_column itemId"%\
                    (self._getExtraPath("angles2.xmd"),self._getTmpPath("ids2.xmd")))
        self.runJob("xmipp_metadata_utilities","-i %s --set intersection %s itemId -o %s"%\
                    (self._getTmpPath("ids1.xmd"),self._getTmpPath("ids2.xmd"),self._getTmpPath("ids.xmd")))
        self.runJob("xmipp_metadata_utilities","-i %s --set intersection %s itemId -o %s"%\
                    (self._getExtraPath("angles1.xmd"),self._getTmpPath("ids.xmd"),self._getExtraPath("angles1_common.xmd")))
        self.runJob("xmipp_metadata_utilities","-i %s --set intersection %s itemId -o %s"%\
                    (self._getExtraPath("angles2.xmd"),self._getTmpPath("ids.xmd"),self._getExtraPath("angles2_common.xmd")))
        self.runJob("xmipp_metadata_utilities","-i %s --operate sort itemId"%self._getExtraPath("angles1_common.xmd"))
        self.runJob("xmipp_metadata_utilities","-i %s --operate sort itemId"%self._getExtraPath("angles2_common.xmd"))

        self.runJob("xmipp_angular_distance","--ang1 %s --ang2 %s --sym %s --check_mirrors --oroot %s"%\
                    (self._getExtraPath("angles1_common.xmd"),self._getExtraPath("angles2_common.xmd"),
                     self.symmetryGroup,self._getTmpPath("angular_distance")))
        self.runJob("xmipp_metadata_utilities",'-i %s -o %s --operate keep_column "angleDiff shiftDiff"'%\
                    (self._getTmpPath("angular_distance.xmd"),self._getTmpPath("diffs.xmd")))
        self.runJob("xmipp_metadata_utilities","-i %s --set merge %s"%\
                    (self._getExtraPath("angles1_common.xmd"),self._getTmpPath("diffs.xmd")))


    def createOutputStep(self):
        imgSet1 = self.inputParticles1.get()
        imgSetOut = self._createSetOfParticles()
        imgSetOut.copyInfo(imgSet1)
        imgSetOut.setAlignmentProj()
        self.iterMd = md.iterRows(self._getExtraPath("angles1_common.xmd"), md.MDL_ITEM_ID)
        self.lastRow = next(self.iterMd) 
        imgSetOut.copyItems(imgSet1,
                            updateItemCallback=self._updateItem)
        self._defineOutputs(outputParticles=imgSetOut)
        self._defineSourceRelation(self.inputParticles1, imgSetOut)
        self._defineSourceRelation(self.inputParticles2, imgSetOut)

    def _updateItem(self, particle, row):
        count = 0
        
        while self.lastRow and particle.getObjId() == self.lastRow.getValue(md.MDL_ITEM_ID):
            count += 1
            if count:
                self._createItemMatrix(particle, self.lastRow)
            try:
                self.lastRow = next(self.iterMd)
            except StopIteration:
                self.lastRow = None
                    
        particle._appendItem = count > 0
        
    def _createItemMatrix(self, particle, row):
        setXmippAttributes(particle, row, xmipp.MDL_SHIFT_DIFF, xmipp.MDL_ANGLE_DIFF)
        #createItemMatrix(particle, row, align=em.ALIGN_PROJ)

    # --------------------------- INFO functions -------------------------------

    def _validate(self):
        validateMsgs = []
        if self.inputParticles1.get() and not self.inputParticles1.hasValue():
            validateMsgs.append('Please provide input particles 1.')            
        if self.inputParticles2.get() and not self.inputParticles2.hasValue():
            validateMsgs.append('Please provide input particles 2.')            
        return validateMsgs
    
    def _summary(self):
        summary = []
        return summary

