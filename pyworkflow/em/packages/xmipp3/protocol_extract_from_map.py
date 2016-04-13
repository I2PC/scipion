# **************************************************************************
# *
# * Authors:         Josue Gomez Blanco (jgomez@cnb.csic.es)
# *
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

from os.path import join, isfile
from shutil import copyfile
from pyworkflow.object import Float, String
from pyworkflow.protocol.params import (PointerParam, FloatParam, STEPS_PARALLEL,
                                        StringParam, BooleanParam, LEVEL_ADVANCED)
from pyworkflow.em.data import Volume
from pyworkflow.em import Viewer
import pyworkflow.em.metadata as md
from pyworkflow.em.protocol import ProtAnalysis3D
from pyworkflow.utils.path import moveFile, makePath
from pyworkflow.em.packages.xmipp3.convert import (writeSetOfParticles,
                                                   writeSetOfVolumes,
                                                   getImageLocation)


class XmippProtExtractFromVolume(ProtAnalysis3D):
    """    
    Extract the information contained in a volume to the experimental
    particles. The particles must be projection alignment in order to
    properly generate volume projections to extract the information.
    A typical case of use, is in the deletion of the capsid to the
    experimental image to only refine the genetic material.
    """
    _label = 'extract from volume'
    
    def __init__(self, *args, **kwargs):
        ProtAnalysis3D.__init__(self, *args, **kwargs)
        
    #--------------------------- DEFINE param functions --------------------------------------------   
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', PointerParam,
                      pointerClass='SetOfParticles', pointerCondition='hasAlignment',
                      label="Input particles", important=True,
                      help='Select the experimental particles.')
        form.addParam('inputVolume', PointerParam, pointerClass='Volume',
                      label="Input volume",
                      help='Select the input volume. Is desirable that the '
                           'volume was generated with the input particles.')
        form.addParallelSection(threads=0, mpi=0)
    
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        partSetId = self.inputParticles.get().getObjId()
        volume = self.inputVolume.get()
        volFn = getImageLocation(volume)
        
        self._insertFunctionStep('convertInputStep', partSetId)
        self._insertFunctionStep('phantomProject', volFn)
        self._insertFunctionStep('createOutputStep')
    
    #--------------------------- STEPS functions --------------------------------------------
    def convertInputStep(self, particlesId):
        """ Write the input images as a Xmipp metadata file. 
        particlesId: is only need to detect changes in
        input particles and cause restart from here.
        """
        writeSetOfParticles(self.inputParticles.get(), self._getInputParticlesFn())
    
    def phantomProject(self,volName):
        from convert import createParamPhantomFile
        
        imgSet = self.inputParticles.get()
        phantomFn = self._getExtraPath('params')
        pathParticles = self._getInputParticlesFn()
        dimX, _, _ = self.inputParticles.get().getDim()
        
        createParamPhantomFile(imgSet, phantomFn, dimX, pathParticles)
        
        params =  ' -i %s' % volName
        params += ' --params %s' % phantomFn
        params += ' -o %s' % self._getOutputRefsFn()
        params += ' --sampling_rate % 0.3f' % imgSet.getSamplingRate()
        
        self.runJob('xmipp_phantom_project', params)

    def createOutputStep(self):
        pass
#         outputVols = self._createSetOfVolumes()
#         imgSet = self.inputParticles.get()
#         for i, vol in enumerate(self._iterInputVols()):
#             volume = vol.clone()               
#             volDir = self._getVolDir(i+1)
#             volPrefix = 'vol%03d_' % (i+1)
#             validationMd = self._getExtraPath(volPrefix + 'validation.xmd')
#             moveFile(join(volDir, 'validation.xmd'), 
#                      validationMd)
#             clusterMd = self._getExtraPath(volPrefix + 'clusteringTendency.xmd')
#             moveFile(join(volDir, 'clusteringTendency.xmd'), clusterMd)
#             
#             outImgSet = self._createSetOfParticles(volPrefix)
#             
#             outImgSet.copyInfo(imgSet)
# 
#             outImgSet.copyItems(imgSet,
#                                 updateItemCallback=self._setWeight,
#                                 itemDataIterator=md.iterRows(clusterMd,
#                                                              sortByLabel=md.MDL_ITEM_ID))
#                         
#             mdValidatoin = md.MetaData(validationMd)
#             weight = mdValidatoin.getValue(md.MDL_WEIGHT, mdValidatoin.firstObject())
#             volume.weight = Float(weight)
#             volume.clusterMd = String(clusterMd)
#             volume.cleanObjId() # clean objects id to assign new ones inside the set            
#             outputVols.append(volume)
#             self._defineOutputs(outputParticles=outImgSet)
#         
#         outputVols.setSamplingRate(volume.getSamplingRate())
#         self._defineOutputs(outputVolumes=outputVols)
#         #self._defineTransformRelation(self.inputVolumes.get(), volume)
        
    #--------------------------- INFO functions -------------------------------------------- 
    def _validate(self):
        validateMsgs = []
        return validateMsgs

    def _summary(self):
        summary = ["Input particles:  %s" % self.inputParticles.get().getNameId()]
        summary.append("-----------------")
        return summary
    
    def _methods(self):
        messages = []
        return messages
    
    #--------------------------- UTILS functions --------------------------------------------
    def _getInputParticlesFn(self):
        return self._getPath('input_particles.xmd')
    
    def _getOutputRefsFn(self):
        return self._getPath('reference_particles.xmd')
    
    def _getVolDir(self, volIndex):
        return self._getExtraPath('vol%03d' % volIndex)
    
    def _defineMetadataRootName(self, mdrootname,volId):
        
        if mdrootname=='P':
            VolPrefix = 'vol%03d_' % (volId)
            return self._getExtraPath(VolPrefix+'clusteringTendency.xmd')
        if mdrootname=='Volume':

            VolPrefix = 'vol%03d_' % (volId)
            return self._getExtraPath(VolPrefix+'validation.xmd')
            
    def _definePName(self):
        fscFn = self._defineMetadataRootName('P')
        return fscFn
    
    def _defineVolumeName(self,volId):
        fscFn = self._defineMetadataRootName('Volume',volId)
        return fscFn
    
    def _setWeight(self, item, row):  
        item._xmipp_weightClusterability = Float(row.getValue(md.MDL_VOLUME_SCORE1))
