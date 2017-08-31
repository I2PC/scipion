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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os
from glob import glob
import numpy as np
import pyworkflow.em.metadata as md
import pyworkflow.em as em
from pyworkflow.em.protocol import EMProtocol
import pyworkflow.protocol.params as params
from pyworkflow import VERSION_1_1
from pyworkflow.em.packages.xmipp3.convert import getImageLocation
from pyworkflow.protocol.constants import LEVEL_ADVANCED, STEPS_PARALLEL
from pyworkflow.utils.path import cleanPattern, createLink, moveFile, copyFile, makePath, cleanPath
from pyworkflow.object import String
from pyworkflow.em.data import SetOfNormalModes
from pyworkflow.em.packages.xmipp3 import XmippMdRow
from pyworkflow.em.packages.xmipp3.pdb.protocol_pseudoatoms_base import XmippProtConvertToPseudoAtomsBase
import xmipp
from pyworkflow.em.packages.xmipp3.nma.protocol_nma_base import XmippProtNMABase, NMA_CUTOFF_REL

def mds(d, dimensions = 2):
    """
    Multidimensional Scaling - Given a matrix of interpoint distances,
    find a set of low dimensional points that have similar interpoint
    distances.
    """

    (n,n) = d.shape
    E = (-0.5 * d**2)

    # Use mat to get column and row means to act as column and row means.
    Er = np.mat(np.mean(E,1))
    Es = np.mat(np.mean(E,0))

    # From Principles of Multivariate Analysis: A User's Perspective (page 107).
    F = np.array(E - np.transpose(Er) - Es + np.mean(E))

    [U, S, V] = np.linalg.svd(F)

    Y = U * np.sqrt(S)

    return (Y[:,0:dimensions], S)

class XmippProtStructureMapping(XmippProtConvertToPseudoAtomsBase,
                                XmippProtNMABase):
    """ 
    A quantitive analysis of dissimilarities (distances) among the EM maps
    that placing the entire set of density maps in to a common space of
    comparison.The approach is based on statistical analysis of distance
    among elastically aligned EM maps, and results in visualizing those maps
    as points in a lower dimensional distance space.    
    """
    _label = 'structure mapping'
    _lastUpdateVersion = VERSION_1_1
           
    #--------------------------- DEFINE param functions --------------------------------------------
    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL

    def _defineParams(self, form):
        
        form.addSection(label='Input')
        form.addParam('inputVolumes', params.MultiPointerParam, 
                      pointerClass='SetOfVolumes,Volume',  
                      label="Input volume(s)", important=True, 
                      help='Select one or more volumes (Volume or SetOfVolumes)\n'
                           'to be aligned againt the reference volume.')
        form.addParam('keepingOutputFiles', params.BooleanParam, 
                      default=False, expertLevel=LEVEL_ADVANCED,
                      label="Keeping intermediate output files",
                      help="Set to true if you want to keep all intermediate "
                           "produced files during rigid and elastic alignment.")        
        
        form.addSection(label='Pseudoatom')
        XmippProtConvertToPseudoAtomsBase._defineParams(self,form)
        self.getParam("pseudoAtomRadius").setDefault(1.5)
        self.getParam("pseudoAtomTarget").setDefault(2)
        
        form.addSection(label='Normal Mode Analysis')
        form.addParam('rigidAlignment', params.BooleanParam, default=True, expertLevel=LEVEL_ADVANCED,
                      label='Rigid alignment',
                      help='Perform a rigid alignment before elastic alignment')
        form.addParam('optimizeScale', params.BooleanParam, default=False, expertLevel=LEVEL_ADVANCED,
                      label="Optimize scale", condition="rigidAlignment")
        XmippProtNMABase._defineParamsCommon(self,form)
        
               
        form.addParallelSection(threads=4, mpi=1)
        
    #--------------------------- INSERT steps functions --------------------------------------------    
    def _insertAllSteps(self):
        
        cutoffStr=''
        if self.cutoffMode == NMA_CUTOFF_REL:
            cutoffStr = 'Relative %f'%self.rcPercentage.get()
        else:
            cutoffStr = 'Absolute %f'%self.rc.get()

        alignArgs = self._getAlignArgs()
        self.alignmentAlgorithm = 1 # Local alignment
                                                
        volList = [vol.clone() for vol in self._iterInputVolumes()]
        nVoli = 1
                                   
        for voli in volList:
            fnIn = getImageLocation(voli)
            fnMask = self._insertMaskStep(fnIn)
            suffix = "_%d"%nVoli        
            
            self._insertFunctionStep('convertToPseudoAtomsStep', 
                                     fnIn, fnMask, voli.getSamplingRate(), suffix)
            fnPseudoAtoms = self._getPath("pseudoatoms_%d.pdb"%nVoli)
            
            self._insertFunctionStep('computeModesStep', fnPseudoAtoms,
                                     self.numberOfModes, cutoffStr)
            self._insertFunctionStep('reformatOutputStep', os.path.basename(fnPseudoAtoms))
                                            
            stepQualify = self._insertFunctionStep('qualifyModesStep', self.numberOfModes, 
                                     self.collectivityThreshold.get(), 
                                     self._getPath("pseudoatoms_%d.pdb"%nVoli), suffix)
            
            #rigid alignment            
            nVolj = 1
            deps = []
            for volj in volList:
                if nVolj != nVoli:
                    inVolFn = getImageLocation(volj)
                    if self.rigidAlignment:
                        refFn = getImageLocation(voli)
                        volId = volj.getObjId()
                        outVolFn = self._getPath('outputRigidAlignment_vol_%d_to_%d.vol' % (nVolj, nVoli))
                        stepId=self._insertFunctionStep('alignVolumeStep', refFn, inVolFn, outVolFn,
                                                        volId, prerequisites=[stepQualify])
                    else:
                        outVolFn = inVolFn
                        stepId = stepQualify
                    deps.append(self._insertFunctionStep('elasticAlignmentStep', nVoli, voli.getSamplingRate(), nVolj, inVolFn, prerequisites=[stepId]))
                nVolj += 1
            self._insertFunctionStep('gatherSingleVolumeStep',prerequisites=deps) # This is a synchronization step, does not do any real work
            nVoli += 1
               
        self._insertFunctionStep('gatherResultsStep')
        self._insertFunctionStep('managingOutputFilesStep')
                                        
    #--------------------------- STEPS functions --------------------------------------------
    def alignVolumeStep(self, refFn, inVolFn, outVolFn, volId):
        args = "--i1 %s --i2 %s --apply %s" % (refFn, inVolFn, outVolFn)
        args += " --local --rot 0 0 1 --tilt 0 0 1 --psi 0 0 1 -x 0 0 1 -y 0 0 1 -z 0 0 1 --dontScale" 
        args += " --copyGeo %s" % (
                self._getExtraPath('transformation-matrix_vol%06d.txt'%volId))        
        self.runJob("xmipp_volume_align", args)

    def elasticAlignmentStep(self, nVoli, Ts, nVolj, fnAlignedVolj):
        fnVolOut = self._getExtraPath('DeformedVolume_Vol_%d_To_Vol_%d' % (nVolj, nVoli))
        if os.path.exists(fnVolOut+".pdb"):
            return
        
        makePath(self._getExtraPath("modes%d"%nVoli))
        
        for i in range(self.numberOfModes.get() + 1):
            if i==0:
                i += 1 
            copyFile (self._getPath("modes/vec.%d"%i),
                      self._getExtraPath("modes%d/vec.%d"%(nVoli, i)))
            
        mdVol = xmipp.MetaData()
        fnOutMeta = self._getExtraPath('RigidAlignVol_%d_To_Vol_%d.xmd' % (nVolj, nVoli))
        mdVol.setValue(xmipp.MDL_IMAGE, fnAlignedVolj, mdVol.addObject())      
        mdVol.write(fnOutMeta)
                                              
        fnPseudo = self._getPath("pseudoatoms_%d.pdb"%nVoli)
        fnModes = self._getPath("modes_%d.xmd"%nVoli)
        fnDeform = self._getExtraPath('compDeformVol_%d_To_Vol_%d.xmd' % (nVolj, nVoli))
        sigma = Ts * self.pseudoAtomRadius.get()
        fnPseudoOut = self._getExtraPath('PseudoatomsDeformedPDB_Vol_%d_To_Vol_%d.pdb' % 
                                         (nVolj, nVoli))
        self.runJob('xmipp_nma_alignment_vol', 
                    "-i %s --pdb %s --modes %s --sampling_rate %s -o %s --fixed_Gaussian %s --opdb %s"%\
                   (fnOutMeta, fnPseudo, fnModes, Ts, fnDeform, sigma, fnPseudoOut))
        
        self.runJob('xmipp_volume_from_pdb', "-i %s -o %s --sampling %s --fixed_Gaussian %s" % 
                    (fnPseudoOut, fnVolOut, Ts, sigma))
    
    def gatherSingleVolumeStep(self):
        pass
    
    def gatherResultsStep(self):
                         
        volList = [vol.clone() for vol in self._iterInputVolumes()]            
        #score and distance matrix calculation         
        score = [[0 for i in volList] for i in volList]
        nVoli = 1
        
        for voli in volList:
            nVolj = 1
            for volj in volList:
                if nVolj == nVoli:
                    score[(nVoli-1)][(nVolj-1)] = 0
                else:
                    elasticRow = xmipp.MetaData(self._getExtraPath('compDeformVol_%d_To_Vol_%d.xmd' %
                                                                   (nVolj, nVoli)))
                    maxCc = elasticRow.getValue(md.MDL_MAXCC,1)
                    score[(nVoli-1)][(nVolj-1)] = (1 - maxCc)
                nVolj += 1
            nVoli += 1     
                              
        fnRoot = self._getExtraPath ("DistanceMatrix.txt")   
        distance = [[0 for i in volList] for i in volList]
        nVoli = 1
        
        for i in volList:
            nVolj = 1
            for j in volList:
                distance[(nVoli-1)][(nVolj-1)] = (score[(nVoli-1)][(nVolj-1)]+score[(nVolj-1)][(nVoli-1)])/2
                fh = open(fnRoot,"a")
                fh.write("%f\t"%distance[(nVoli-1)][(nVolj-1)])
                fh.close()
                nVolj += 1  
            fh = open(fnRoot,"a")
            fh.write("\n")
            fh.close()
            nVoli += 1                     
        
        distance = np.asarray(distance)
        for i in range(1, 4):
            embed,_ = mds(distance,i)
            embedExtended = np.pad(embed,((0,0),(0,i-embed.shape[1])),"constant",constant_values=0)
            print(embedExtended)
            np.savetxt(self._defineResultsName(i),embedExtended)        
       
    def managingOutputFilesStep(self): 
        cleanPattern(self._getPath('pseudoatoms*'))
        cleanPattern(self._getPath('modes'))
        cleanPattern(self._getExtraPath('vec_ani.pkl'))  
        
        if not self.keepingOutputFiles.get():
            cleanPattern(self._getPath('warnings*'))
            cleanPattern(self._getPath('outputRigid*'))
            cleanPattern(self._getPath('modes*.xmd'))
            cleanPattern(self._getExtraPath('RigidAlign*'))
            cleanPattern(self._getExtraPath('pseudoatoms*'))
            cleanPattern(self._getExtraPath('Pseudoatoms*'))
            cleanPattern(self._getExtraPath('comp*'))
            cleanPattern(self._getExtraPath('Deform*'))
            cleanPattern(self._getExtraPath('transformation-matrix*'))
            cleanPattern(self._getExtraPath('modes*'))
                                      
    #--------------------------- INFO functions --------------------------------------------
    def _validate(self):
        errors = []
        for pointer in self.inputVolumes:
            if pointer.pointsNone():
                errors.append('Invalid input, pointer: %s' % pointer.getObjValue())
                errors.append('              extended: %s' % pointer.getExtended())                
        return errors 
           
    def _summary(self):
        summary = []
        nVols = self._getNumberOfInputs()
            
        if nVols > 0:
            summary.append("Volumes to calculate StructMap: *%d* " % nVols)
        else:
            summary.append("No volumes selected.")
                        
        return summary 
       
    def _citations(self):
        return ['Sorzano2016']
    
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
                    vol.outputName = self._getExtraPath('output_vol%06d_%03d.vol' %
                                                         (itemId, vol.getObjId()))
                    yield vol
    
    def _getNumberOfInputs(self):
        """ Return the total number of input volumes. """
        nVols = 0
        
        for _ in self._iterInputVolumes():
            nVols += 1    
            
        return nVols
                    
    def _getAlignArgs(self):
        alignArgs = " --local --rot 0 0 1 --tilt 0 0 1 --psi 0 0 1 -x 0 0 1 -y 0 0 1 -z 0 0 1"
        if self.optimizeScale:
            alignArgs += " --scale 1 1 0.005"
        else:
            alignArgs += " --dontScale"                    
        return alignArgs
        
    def _defineResultsName(self,i):
        return self._getExtraPath('CoordinateMatrix%d.txt'%i)
        
    
