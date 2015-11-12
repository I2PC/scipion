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

import pyworkflow.em as em
import pyworkflow.protocol.params as params
from pyworkflow.em.packages.xmipp3.convert import getImageLocation
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.utils.path import cleanPattern, createLink, moveFile, copyFile, makePath, cleanPath
from os.path import basename, exists
from pyworkflow.object import String
import os
from pyworkflow.em.data import SetOfNormalModes
from pyworkflow.em.packages.xmipp3 import XmippMdRow
from pyworkflow.em.packages.xmipp3.pdb.protocol_pseudoatoms_base import XmippProtConvertToPseudoAtomsBase
import xmipp
from pyworkflow.em.packages.xmipp3.nma.protocol_nma_base import XmippProtNMABase, NMA_CUTOFF_REL
#from convert import rowToMode, getNMAEnviron
#from pyworkflow.utils import redStr

ALIGN_MASK_CIRCULAR = 0
ALIGN_MASK_BINARY_FILE = 1

ALIGN_ALGORITHM_EXHAUSTIVE = 0
ALIGN_ALGORITHM_LOCAL = 1
ALIGN_ALGORITHM_EXHAUSTIVE_LOCAL = 2
ALIGN_ALGORITHM_FAST_FOURIER = 3

NMA_CUTOFF_REL = 1
NMA_CUTOFF_ABS = 0



class XmippProtStructureMapping(XmippProtConvertToPseudoAtomsBase,XmippProtNMABase):
    """ 
    Multivariate distance analysis of elastically aligned electron microscopy density maps
    for exploring pathways of conformational changes    
     """
    _label = 'structure mapping'
    
    #def __init__(self, **args):
    #    em.ProtAlignVolume.__init__(self, **args)
    #    self.stepsExecutionMode = em.STEPS_PARALLEL
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputVolumes', params.MultiPointerParam, pointerClass='SetOfVolumes,Volume',  
                      label="Input volume(s)", important=True, 
                      help='Select one or more volumes (Volume or SetOfVolumes)\n'
                           'to be aligned againt the reference volume.')
        
        form.addSection(label='Pseudoatom')
        XmippProtConvertToPseudoAtomsBase._defineParams(self,form)
        self.getParam("pseudoAtomRadius").setDefault(2)
        self.getParam("pseudoAtomTarget").setDefault(2)
        
        form.addSection(label='Normal Mode Analysis')
        XmippProtNMABase._defineParamsCommon(self,form)
        form.addParallelSection(threads=8, mpi=1)
        
    #--------------------------- INSERT steps functions --------------------------------------------    
    def _insertAllSteps(self):
        cutoffStr=''
        if self.cutoffMode == NMA_CUTOFF_REL:
            cutoffStr = 'Relative %f'%self.rcPercentage.get()
        else:
            cutoffStr = 'Absolute %f'%self.rc.get()

        volList = [vol.clone() for vol in self._iterInputVolumes()]
                    
        for voli in volList:
            fnIn = getImageLocation(voli)
            fnMask = self._insertMaskStep(fnIn)
            prefix="_%d"%voli.getObjId()        
            self._insertFunctionStep('convertToPseudoAtomsStep', fnIn, fnMask, voli.getSamplingRate(), prefix)
            fnPseudoAtoms = self._getPath("pseudoatoms_%d.pdb"%voli.getObjId())
            self._insertFunctionStep('computeModesStep', fnPseudoAtoms, self.numberOfModes, cutoffStr)
            self._insertFunctionStep('reformatOutputStep', os.path.basename(fnPseudoAtoms), prefix)
#         
#             #self._insertFunctionStep('qualifyModesStep', n, self.collectivityThreshold.get(), self.structureEM)
#             for volj in volList:
#                 if volj is not voli:
#                     refFn = getImageLocation(voli)
#                     volFn = getImageLocation(volj)
#                     stepId = self._insertFunctionStep('alignVolumeStep', refFn, volFn, vol.outputName, 
#                                               maskArgs, alignArgs, prerequisites=[])
#                     alignSteps.append(stepId)
            
        
                   
    #    self._insertFunctionStep('createOutputStep', prerequisites=alignSteps)
        
    #--------------------------- STEPS functions --------------------------------------------
    
   
        
        
    
    
    
    
    
    
    
    
    
    
    
    def alignVolumeStep(self, refFn, inVolFn, outVolFn, maskArgs, alignArgs):
        
        args = "--i1 %s --i2 %s --apply %s" % (refFn, inVolFn, outVolFn)
        args += maskArgs
        args += alignArgs
        
        self.runJob("xmipp_volume_align", args)
        
        if self.alignmentAlgorithm == ALIGN_ALGORITHM_EXHAUSTIVE_LOCAL:
            args = "--i1 %s --i2 %s --apply --local" % (refFn, outVolFn)
            self.runJob("xmipp_volume_align", args)
            
            
            
            
            
            
            
            
         
            
      
    #def createOutputStep(self):
    
    """ for pseudoatom"""
    
    #inputVol = self.inputStructure.get()
    #pdb = PdbFile(self._getPath('pseudoatoms.pdb'), pseudoatoms=True)
    #self.createChimeraScript(inputVol, pdb)
    #self._defineOutputs(outputPdb=pdb)
    #self._defineSourceRelation(self.inputStructure, pdb)
    
    
    
    """ for mode analysis""" 
    
    #fnSqlite = self._getPath('modes.sqlite')
    #nmSet = SetOfNormalModes(filename=fnSqlite)

    #md = xmipp.MetaData(self._getPath('modes.xmd'))
    #row = XmippMdRow()
        
    #for objId in md:
    #    row.readFromMd(md, objId)
    #    nmSet.append(rowToMode(row))
    #inputPdb = self.inputStructure.get()
    #nmSet.setPdb(inputPdb)
    #self._defineOutputs(outputModes=nmSet)
    #self._defineSourceRelation(self.inputStructure, nmSet)
    
    
    
    
    
    
    """ for align volume"""
    
    #    vols = []
    #    for vol in self._iterInputVolumes():
    #        outVol = em.Volume()
    #        outVol.setLocation(vol.outputName)
    #        outVol.setObjLabel(vol.getObjLabel())
    #        outVol.setObjComment(vol.getObjComment())
    #        vols.append(outVol) 
    #        
    #    if len(vols) > 1:
    #        volSet = self._createSetOfVolumes()
    #        volSet.setSamplingRate(self.inputReference.get().getSamplingRate())
    #        for vol in vols:
    #            volSet.append(vol)
    #        outputArgs = {'outputVolumes': volSet}
    #    else:
    #        vols[0].setSamplingRate(self.inputReference.get().getSamplingRate())
    #        outputArgs = {'outputVolume': vols[0]}
    #        
    #    self._defineOutputs(**outputArgs)
    #    for pointer in self.inputVolumes:
    #        self._defineSourceRelation(pointer, outputArgs.values()[0])

    
    #--------------------------- INFO functions --------------------------------------------
    
    def _validate(self):
        errors = []
        for pointer in self.inputVolumes:
            if pointer.pointsNone():
                errors.append('Invalid input, pointer: %s' % pointer.getObjValue())
                errors.append('              extended: %s' % pointer.getExtended())
        return errors    
    
    #def _summary(self):
    #    summary = []
    #    nVols = self._getNumberOfInputs()
    #        
    #    if nVols > 0:
    #        summary.append("Volumes to align: *%d* " % nVols)
    #    else:
    #        summary.append("No volumes selected.")
    #    summary.append("Alignment method: %s" % self.getEnumText('alignmentAlgorithm'))
    #            
    #    return summary
    
    #def _methods(self):
    #    nVols = self._getNumberOfInputs()
    #    
    #    if nVols > 0:
    #        methods = 'We aligned %d volumes against a reference volume using ' % nVols
            #TODO: Check a more descriptive way to add the reference and 
            # all aligned volumes to the methods (such as obj.getNameId())
            # also to show the number of volumes from each set in the input.
            # This approach implies to consistently include also the outputs
            # ids to be tracked in all the workflow's methods.
    #       if self.alignmentAlgorithm == ALIGN_ALGORITHM_FAST_FOURIER:
    #            methods += ' the Fast Fourier alignment described in [Chen2013].' 
    #            
    #        elif self.alignmentAlgorithm == ALIGN_ALGORITHM_LOCAL:
    #            methods += ' a local search of the alignment parameters.'
    #        elif self.alignmentAlgorithm == ALIGN_ALGORITHM_EXHAUSTIVE:
    #            methods += ' an exhaustive search.'
    #        elif self.alignmentAlgorithm == ALIGN_ALGORITHM_EXHAUSTIVE_LOCAL:
    #            methods += ' an exhaustive search followed by a local search.'
    #    else:
    #        methods = 'No methods available yet.'
    #        
    #    return [methods]
        
    #def _citations(self):
    #    if self.alignmentAlgorithm == ALIGN_ALGORITHM_FAST_FOURIER:
    #        return ['Chen2013']
        
    #--------------------------- UTILS functions --------------------------------------------
    def _iterInputVolumes(self):
        """ Iterate over all the input volumes. """
        for pointer in self.inputVolumes:
            item = pointer.get()
            if item is None:
                break
            itemId = item.getObjId()
            if isinstance(item, em.Volume):
                item.outputName = self._getExtraPath('output_vol%06d.vol' % itemId)
                yield item
            elif isinstance(item, em.SetOfVolumes):
                for vol in item:
                    vol.outputName = self._getExtraPath('output_vol%06d_%03d.vol' % (itemId, vol.getObjId()))
                    yield vol
                    
    def _getAlignArgs(self):
        alignArgs = ''
        
        if self.alignmentAlgorithm == ALIGN_ALGORITHM_FAST_FOURIER:
            alignArgs += " --frm"
            
        elif self.alignmentAlgorithm == ALIGN_ALGORITHM_LOCAL:
            alignArgs += " --local --rot %f %f 1 --tilt %f %f 1 --psi %f %f 1 -x %f %f 1 -y %f %f 1 -z %f %f 1 --scale %f %f 0.005" %\
               (self.initialRotAngle, self.initialRotAngle,
                self.initialTiltAngle, self.initialTiltAngle,
                self.initialInplaneAngle, self.initialInplaneAngle,
                self.initialShiftX, self.initialShiftX,
                self.initialShiftY, self.initialShiftY,
                self.initialShiftZ,self.initialShiftZ,
                self.initialScale, self.initialScale)
        else: # Exhaustive or Exhaustive+Local
            alignArgs += " --rot %f %f %f --tilt %f %f %f --psi %f %f %f -x %f %f %f -y %f %f %f -z %f %f %f --scale %f %f %f" %\
               (self.minRotationalAngle, self.maxRotationalAngle, self.stepRotationalAngle,
                self.minTiltAngle, self.maxTiltAngle, self.stepTiltAngle,
                self.minInplaneAngle, self.maxInplaneAngle, self.stepInplaneAngle,
                self.minimumShiftX, self.maximumShiftX, self.stepShiftX,
                self.minimumShiftY, self.maximumShiftY, self.stepShiftY,
                self.minimumShiftZ, self.maximumShiftZ, self.stepShiftZ,
                self.minimumScale, self.maximumScale, self.stepScale)
        return alignArgs
     
  