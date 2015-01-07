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
This sub-package contains wrapper around resolution3D Xmipp protocol
"""


from pyworkflow.em import *  
from pyworkflow.utils import *  
from math import floor
from xmipp3 import XmippProtocol
from convert import createXmippInputVolumes, readSetOfVolumes, locationToXmipp
from xmipp import MetaData, MDL_RESOLUTION_FREQ, MDL_RESOLUTION_FRC, MDL_RESOLUTION_DPR

class XmippProtResolution3D(ProtAnalysis3D):
    """ Computes resolution by several methods """
    _label = 'resolution 3D'
      
    #--------------------------- DEFINE param functions --------------------------------------------   
    def _defineParams(self, form):
        form.addSection(label='Input')
        # Volumes to process
        form.addParam('inputVolume', PointerParam, label="Volume to compare", important=True, 
                      pointerClass='Volume',
                      help='This volume will be compared to the reference volume.')  
        form.addParam('doFSC', BooleanParam, default=True,
                      label="Calculate FSC and DPR?", 
                      help='If set True calculate FSC and DPR.')
        form.addParam('referenceVolume', PointerParam, label="Reference volume", condition='doFSC', 
                      pointerClass='Volume',
                      help='Input volume will be compared to this volume.')  
        form.addParam('doStructureFactor', BooleanParam, default=True,
                      label="Calculate Structure Factor?", 
                      help='If set True calculate Structure Factor')
    
    #--------------------------- INSERT steps functions --------------------------------------------  
    def _insertAllSteps(self):
        """Insert all steps to calculate the resolution of a 3D reconstruction. """
        
        self.inputVol = locationToXmipp(*self.inputVolume.get().getLocation())
        self.refVol = locationToXmipp(*self.referenceVolume.get().getLocation())
        
        if self.doFSC:
            self._insertFunctionStep('calculateFscStep')
        if self.doStructureFactor:
            self._insertFunctionStep('structureFactorcStep')
        
        
        self._insertFunctionStep('createSummaryStep')

    #--------------------------- STEPS steps functions --------------------------------------------
    def calculateFscStep(self):
        """ Calculate the FSC between two volumes"""
        samplingRate = self.inputVolume.get().getSamplingRate()
        fscFn = self._defineFscName()
        args="--ref %s -i %s -o %s --sampling_rate %f --do_dpr" % (self.refVol, self.inputVol, fscFn, samplingRate)
        self.runJob("xmipp_resolution_fsc", args)
    
    def structureFactorcStep(self):
        """ Calculate the structure factors of the volume"""
        samplingRate = self.inputVolume.get().getSamplingRate()
        structureFn = self._defineStructFactorName()
        args = "-i %s -o %s --sampling %f" % (self.inputVol, structureFn, samplingRate)
        self.runJob("xmipp_volume_structure_factor", args)
    
    def calculateFSCResolution(self, md, threshold):
        xl=-1
        xr=-1
        yl=-1
        yr=-1
        leftSet=False
        rightSet=False
        for objId in md:
            freq=md.getValue(MDL_RESOLUTION_FREQ,objId)
            FSC=md.getValue(MDL_RESOLUTION_FRC,objId)
            if FSC>threshold:
                xl=freq
                yl=FSC
                leftSet=True
            else:
                xr=freq
                yr=FSC
                rightSet=True
                break
        if leftSet and rightSet:
            x=xl+(threshold-yl)/(yr-yl)*(xr-xl)
            return 1.0/x
        else:
            return -1
    
    def calculateDPRResolution(self, md, threshold):
        xl=-1
        xr=-1
        yl=-1
        yr=-1
        leftSet=False
        rightSet=False
        for objId in md:
            freq=md.getValue(MDL_RESOLUTION_FREQ,objId)
            DPR=md.getValue(MDL_RESOLUTION_DPR,objId)
            if DPR<threshold:
                xl=freq
                yl=DPR
                leftSet=True
            else:
                xr=freq
                yr=DPR
                rightSet=True
                break
        if leftSet and rightSet:
            x=xl+(threshold-yl)/(yr-yl)*(xr-xl)
            return 1.0/x
        else:
            return -1

    def createSummaryStep(self):
        summary=""
        if self.doFSC.get():
            summary+="FSC file: [%s]\n" % self._defineFscName()
            md=MetaData(self._defineFscName())
            summary+="   Resolution FSC=0.5: %3.2f Angstroms\n"%self.calculateFSCResolution(md,0.5)
            summary+="   Resolution FSC=0.143: %3.2f Angstroms\n"%self.calculateFSCResolution(md,0.143)
            summary+="   Resolution DPR=45: %3.2f Angstroms\n"%self.calculateDPRResolution(md,45)
        if self.doStructureFactor:
            summary+="Structure factor file: [%s]\n" % self._defineStructFactorName()
        self.summaryVar.set(summary)
    
    #--------------------------- INFO steps functions --------------------------------------------
    def _validate(self):
        validateMsgs = []
        
        if not self.inputVolume.hasValue():
            validateMsgs.append('Please provide an initial volume.')
        if self.doFSC.get() and not self.referenceVolume.hasValue():
            validateMsgs.append('Please provide a reference volume.')
            
        return validateMsgs
    
    def _summary(self):
        return [self.summaryVar.get()]
    
    def _methods(self):
        messages = []
        if self.doFSC.get():
            messages.append('We obtained the FSC of the volume %s' % self.inputVolume.get().getNameId())
            messages.append('taking the volume %s' % self.referenceVolume.get().getNameId() + ' as reference')
        return messages
              
    def _defineStructFactorName(self):
        return self._getPath('structureFactor.xmd')
    
    def _defineFscName(self):
        return self._getPath('fsc.xmd')
