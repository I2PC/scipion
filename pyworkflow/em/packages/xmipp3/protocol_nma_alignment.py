# **************************************************************************
# *
# * Authors:  Carlos Oscar Sanchez Sorzano (coss@cnb.csic.es), May 2013
# *           Qiyu Jin
# *           Slavica Jonic                (jonic@impmc.upmc.fr)
# * Ported to Scipion:
# *           J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es), Jan 2014
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

import math
from glob import glob

from pyworkflow.em import *  
from pyworkflow.utils import * 
from pyworkflow.protocol.constants import LEVEL_EXPERT, LEVEL_ADVANCED
from convert import createXmippInputImages
import xmipp

NMA_ALIGNMENT_WAV = 0
NMA_ALIGNMENT_PROJ = 1    

        
class XmippProtAlignmentNMA(EMProtocol):
    """ Protocol for flexible angular alignment. """
    _label = 'nma analysis'
    _references = ["[[http://www.ncbi.nlm.nih.gov/pubmed/15885434][Jonic, et.al, Ultramic (2005)]]",
                   "[[http://www.ncbi.nlm.nih.gov/pubmed/15099579][Sorzano, et.al, JSB (2004)]]"
                   ]
    
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', PointerParam, label="Input particles", 
                      pointerClass='SetOfParticles',
                      help='Select the set of particles that you want to use for flexible analysis.')  
        form.addParam('inputPdb', PointerParam, label="Input PDB",  
                      pointerClass='PdbFile',
                      help='Atomic or pseudo-atomic structure to apply the normal modes.')
        form.addParam('inputModes', PointerParam, label="Normal modes",  
                      pointerClass='NormalModes',
                      help='Set of normal modes to explore.')
        
        form.addSection(label='Angular assignment and mode detection')
        form.addParam('trustRegionScale', IntParam, default=1,
                      expertLevel=LEVEL_ADVANCED,
                      label='Trust region scale',
                      help='For elastic alignment, this parameter scales the initial \n'
                           'value of the trust region radius for optimization purposes.\n'
                           'Use larger values for larger expected deformation amplitudes.')    
        form.addParam('alignmentMethod', EnumParam, default=NMA_ALIGNMENT_WAV,
                      choices=['wavelets & splines', 'projection matching'],
                      label='Alignment method',
                      help='For rigid-body alignment, use Projection Matching (faster) instead\n'
                           'of Wavelets and Splines (more accurate). In the case of Wavelets \n'
                           'and Splines, the size of images should be a power of 2.')
        form.addParam('discreteAngularSampling', FloatParam, default=10,
                      label="Discrete angular sampling (deg)", 
                      help='This parameter is used in Projection Matching and Wavelets methods\n'
                           'for a rough rigid-body alignment. It is the angular step (in degrees)\n'
                           'with which the library of reference projections is computed. This \n'
                           'alignment is refined with Splines method if Wavelets and Splines \n'
                           'alignment is chosen.')
                      
        form.addParallelSection(threads=1, mpi=8)    
             
    def _printWarnings(self, *lines):
        """ Print some warning lines to 'warnings.xmd', 
        the function should be called inside the working dir."""
        fWarn = open("warnings.xmd",'w')
        for l in lines:
            print >> fWarn, l
        fWarn.close()
        
    def _defineSteps(self):
        atomsFn = self.inputPdb.get().getFileName()
        modesFn = self.inputModes.get().getFileName()
        # Convert input images if necessary
        imgFn = createXmippInputImages(self, self.inputParticles.get())
        self._insertFunctionStep('copyFilesStep', atomsFn, modesFn)
        self._insertFunctionStep("performNmaStep", imgFn,
                        self._getBasePath(atomsFn), self._getBasePath(modesFn))
        self._insertFunctionStep("extractDeformationsStep")
        
    def copyFilesStep(self, atomsFn, modesFn):
        """ Copy the input files to the local working dir. """
        for fn in [modesFn, atomsFn]:
            copyFile(fn, self._getBasePath(fn))
        
    def performNmaStep(self, imgFn, atomsFn, modesFn):
        sampling = self.inputParticles.get().getSamplingRate()
        discreteAngularSampling = self.discreteAngularSampling.get()
        trustRegionScale = self.trustRegionScale.get()
        odir = self._getTmpPath()
        
        args = "-i %(imgFn)s --pdb %(atomsFn)s --modes %(modesFn)s --sampling_rate %(sampling)f "
        args += "--discrAngStep %(discreteAngularSampling)f --odir %(odir)s --centerPDB "
        args += "--trustradius_scale %(trustRegionScale)d --resume "
        
        if self.inputPdb.get().getPseudoAtoms():
            args += "--fixed_Gaussian "
        
        if self.alignmentMethod == NMA_ALIGNMENT_PROJ:
            args += "--projMatch "
    
        self.runJob(None, "xmipp_nma_alignment", args % locals())
        cleanPath(self._getPath('nmaTodo.xmd'))
    
    def extractDeformationsStep(self, imgFn):
        md = xmipp.MetaData(imgFn)
        defFn = self._getExtraPath("deformations.txt")
        fhDef = open(defFn, 'w')
        for objId in md:
            deformation = md.getValue(xmipp.MDL_NMA, objId)
            for coef in deformation:
                fhDef.write("%f " % coef)
            fhDef.write("\n")
        fhDef.close()
        return defFn
    
    def createOutputStep(self):
        if self.structureEM:
            pdb = PdbFile(self._getPath('pseudoatoms.pdb'))
            self._defineOutputs(outputPdb=pdb)
        modes = NormalModes()
        self._defineOutputs(outputModes=modes)

    def _summary(self):
        summary = []
        return summary
    
    def _validate(self):
        validateMsgs = []
        return validateMsgs
    
