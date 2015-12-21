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

import os
from glob import glob
import pyworkflow.em.metadata as md
import pyworkflow.em as em
import pyworkflow.protocol.params as params
from pyworkflow.em.packages.xmipp3.convert import getImageLocation
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.utils.path import cleanPattern, createLink, moveFile, copyFile, makePath, cleanPath
from pyworkflow.object import String
from pyworkflow.em.data import SetOfNormalModes
from pyworkflow.em.packages.xmipp3 import XmippMdRow
from pyworkflow.em.packages.xmipp3.pdb.protocol_pseudoatoms_base import XmippProtConvertToPseudoAtomsBase
import xmipp
from pyworkflow.em.packages.xmipp3.nma.protocol_nma_base import XmippProtNMABase, NMA_CUTOFF_REL
from pyworkflow.em.packages.xmipp3.protocol_align_volume import XmippProtAlignVolume
from sklearn import manifold

class XmippProtStructureMapping(XmippProtConvertToPseudoAtomsBase,XmippProtNMABase, XmippProtAlignVolume):
    """ 
    Multivariate distance analysis of elastically aligned electron microscopy density maps
    for exploring pathways of conformational changes    
     """
    _label = 'structure mapping'
           
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
        
               
        form.addParallelSection(threads=4, mpi=1)
        
    #--------------------------- INSERT steps functions --------------------------------------------    
    def _insertAllSteps(self):
        
        cutoffStr=''
        if self.cutoffMode == NMA_CUTOFF_REL:
            cutoffStr = 'Relative %f'%self.rcPercentage.get()
        else:
            cutoffStr = 'Absolute %f'%self.rc.get()

        maskArgs = ''
        alignArgs = self._getAlignArgs()
        ALIGN_ALGORITHM_EXHAUSTIVE_LOCAL = 1
        self.alignmentAlgorithm = ALIGN_ALGORITHM_EXHAUSTIVE_LOCAL       
                                
        volList = [vol.clone() for vol in self._iterInputVolumes()]
        nVoli = 1
                                   
        for voli in volList:
            fnIn = getImageLocation(voli)
            fnMask = self._insertMaskStep(fnIn)
            suffix = "_%d"%nVoli        
            
            self._insertFunctionStep('convertToPseudoAtomsStep', fnIn, fnMask, voli.getSamplingRate(), suffix)
            fnPseudoAtoms = self._getPath("pseudoatoms_%d.pdb"%nVoli)
            
            self._insertFunctionStep('computeModesStep', fnPseudoAtoms, self.numberOfModes, cutoffStr)
            self._insertFunctionStep('reformatOutputStep', os.path.basename(fnPseudoAtoms))
                                            
            self._insertFunctionStep('qualifyModesStep', self.numberOfModes, self.collectivityThreshold.get(), 
                                        self._getPath("pseudoatoms_%d.pdb"%nVoli), suffix)
            #rigid alignment            
            nVolj = 1
            for volj in volList:
                if nVolj != nVoli:
                    refFn = getImageLocation(voli)
                    inVolFn = getImageLocation(volj)
                    outVolFn = self._getPath('outputRigidAlignment_vol_%d_to_%d.vol' % (nVolj, nVoli))
                    self._insertFunctionStep('alignVolumeStep', refFn, inVolFn, outVolFn, maskArgs, alignArgs)
                nVolj += 1   
            
                 
            #elastic alignment
            self._insertFunctionStep('elasticAlignmentStep',nVoli, voli )
            nVoli += 1
               
        self._insertFunctionStep('gatherResultsStep')
                                        
    #--------------------------- STEPS functions --------------------------------------------
    def elasticAlignmentStep(self, nVoli, voli):
            
        makePath(self._getExtraPath("modes%d"%nVoli))
        
        for i in range(self.numberOfModes.get() + 1):
            if i == 0 :
                i += 1 
            copyFile (self._getPath("modes/vec.%d"%i), self._getExtraPath("modes%d/vec.%d"%(nVoli, i)))
            
        mdVols = xmipp.MetaData()
        files = glob(self._getPath('outputRigidAlignment_vol_*_to_%d.vol')%nVoli)
        fnOutMeta = self._getExtraPath('RigidAlignToVol_%d.xmd')%nVoli
        for f in files:
            mdVols.setValue(xmipp.MDL_IMAGE, f, mdVols.addObject())      
        mdVols.write(fnOutMeta)
                                              
        fnPseudo = self._getPath("pseudoatoms_%d.pdb"%nVoli)
        fnModes = self._getPath("modes_%d.xmd"%nVoli)
        Ts = voli.getSamplingRate()
        fnDeform = self._getExtraPath("compDeformVol_%d.xmd"%nVoli)
        sigma = Ts * self.pseudoAtomRadius.get() 
        self.runJob('xmipp_nma_alignment_vol', "-i %s --pdb %s --modes %s --sampling_rate %s -o %s --fixed_Gaussian %s"%\
                (fnOutMeta, fnPseudo, fnModes, Ts, fnDeform, sigma))
        
    
    
    def gatherResultsStep(self):
                
        volList = [vol.clone() for vol in self._iterInputVolumes()]
            
        #score and distance matrix calculation         
        score = [[0 for i in volList] for i in volList]
        nVoli = 1
        
        for voli in volList:
            elastAlign = xmipp.MetaData(self._getExtraPath("compDeformVol_%d.xmd"%nVoli))
            mdIter = md.iterRows(elastAlign)
            nVolj = 1
            for volj in volList:
                if nVolj == nVoli:
                    score[(nVoli-1)][(nVolj-1)] = 0
                    
                else:
                    elasticRow = next(mdIter)
                    maxCc = elasticRow.getValue(md.MDL_MAXCC)
                    score[(nVoli-1)][(nVolj-1)] = (1 - maxCc)
                nVolj += 1
            nVoli += 1 
                                           
        #print "score matrix is: "
        #print score
                
        fnRoot = self._getExtraPath ("DistanceMatrix.txt")   
        distance = [[0 for i in volList] for i in volList]
        nVoli = 1
        for i in volList:
            nVolj = 1
            for j in volList:
                distance[(nVoli-1)][(nVolj-1)] = (score[(nVoli-1)][(nVolj-1)] + score[(nVolj-1)][(nVoli-1)])/2
                fh = open(fnRoot,"a")
                fh.write("%f\t"%distance[(nVoli-1)][(nVolj-1)])
                fh.close()
                nVolj += 1  
            fh = open(fnRoot,"a")
            fh.write("\n")
            fh.close()
            nVoli += 1                     
               
        for i in range(1, 4):
               
            mds = manifold.MDS(n_components=i, metric=True, max_iter=3000, eps=1e-9, random_state=0, dissimilarity="precomputed", n_jobs=1)
            embed3d = mds.fit(distance).embedding_ 
            
            print i
            print mds
            print embed3d  
               
            nVoli = 1
            for x in volList:
                for y in range(i):
                    fh = open(self._getExtraPath ("CoordinateMatrix%d.txt"%i),"a")
                    fh.write("%f\t"%embed3d[(nVoli - 1)][(y)])
                    fh.close()
                    fh = open(self._getExtraPath ("CoordinateMatrixColumnF%d.txt"%i),"a")
                    fh.write("%f\n"%embed3d[(nVoli - 1)][(y)])
                    fh.close()
                fh = open(self._getExtraPath ("CoordinateMatrix%d.txt"%i),"a")
                fh.write("\n")
                fh.close()
                nVoli += 1 
        
        copyFile (self._getExtraPath ("CoordinateMatrixColumnF1.txt"), self._defineResultsName1())
        copyFile (self._getExtraPath ("CoordinateMatrixColumnF2.txt"), self._defineResultsName2())
        copyFile (self._getExtraPath ("CoordinateMatrixColumnF3.txt"), self._defineResultsName3())
                       
        cleanPattern(self._getExtraPath('pseudoatoms*'))
        cleanPattern(self._getExtraPath('vec_ani.pkl'))  
        cleanPattern(self._getExtraPath('CoordinateMatrixColumnF*'))                      
        
    #--------------------------- INFO functions --------------------------------------------
    
    def _validate(self):
        errors = []
        for pointer in self.inputVolumes:
            if pointer.pointsNone():
                errors.append('Invalid input, pointer: %s' % pointer.getObjValue())
                errors.append('              extended: %s' % pointer.getExtended())
        numberOfDimensions=self.numberOfDimensions.get()
        if numberOfDimensions > 3 or numberOfDimensions < 1:
            errors.append("The number of dimensions should be 1, 2 or, at most, 3.")
        
        return errors 
       
    def _summary(self):
        summary = []
        nVols = self._getNumberOfInputs()
            
        if nVols > 0:
            summary.append("Volumes to calculate StructMap: *%d* " % nVols)
        else:
            summary.append("No volumes selected.")
                        
        return summary
    
    def _methods(self):
        messages = []
        messages.append('C.O.S. Sorzano et. al. "StructMap: Multivariate distance analysis of elastically aligned electron microscopy density maps\n'
                        '                         for exploring pathways of conformational changes"')
        return messages
        
    def _citations(self):
        return ['C.O.S.Sorzano2015']
    
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
    
    def _getNumberOfInputs(self):
        """ Return the total number of input volumes. """
        nVols = 0
        for _ in self._iterInputVolumes():
            nVols += 1    
            
        return nVols
                    
    def _getAlignArgs(self):
        alignArgs = ''
        alignArgs += " --local --rot %f %f 1 --tilt %f %f 1 --psi %f %f 1 -x %f %f 1 -y %f %f 1 -z %f %f 1 --scale %f %f 0.005" %\
                (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1)
                       
        return alignArgs
        
    def _defineResultsName1(self):
        return self._getExtraPath('CoordinateMatrixColumn1.txt')
    
    def _defineResultsName2(self):
        return self._getExtraPath('CoordinateMatrixColumn2.txt')
    
    def _defineResultsName3(self):
        return self._getExtraPath('CoordinateMatrixColumn3.txt')
     
    