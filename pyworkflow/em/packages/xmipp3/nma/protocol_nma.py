# **************************************************************************
# *
# * Authors:  Carlos Oscar Sanchez Sorzano (coss@cnb.csic.es), May 2013
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os
import math
from os.path import basename, exists, join

from pyworkflow.utils import redStr
from pyworkflow.utils.path import copyFile, createLink, makePath, cleanPath, moveFile
from pyworkflow.protocol.params import (PointerParam, IntParam, FloatParam, 
                                        LEVEL_ADVANCED)
from pyworkflow.em.data import SetOfNormalModes
from pyworkflow.em.packages.xmipp3 import XmippMdRow
import xmipp
from protocol_nma_base import XmippProtNMABase, NMA_CUTOFF_REL
from convert import rowToMode, getNMAEnviron
        
class XmippProtNMA(XmippProtNMABase):
    """ Flexible angular alignment using normal modes """
    _label = 'nma analysis'
    
    def _defineParams(self, form):
        form.addSection(label='Normal Mode Analysis')
        form.addParam('inputStructure', PointerParam, label="Input structure", important=True, 
                      pointerClass='PdbFile',
                      help='The input structure can be an atomic model (true PDB) or a pseudoatomic model\n'
                           '(an EM volume converted into pseudoatoms)')
        XmippProtNMABase._defineParamsCommon(self,form)
        form.addParam('rtbBlockSize', IntParam, default=10,
                      expertLevel=LEVEL_ADVANCED,
                      label='Number of residues per RTB block',
                      help='This is the RTB block size for the RTB NMA method. \n'
                           'When calculating the normal modes, aminoacids are grouped\n'
                           'into blocks of this size that are moved translationally\n'
                           'and rotationally together.') 
        form.addParam('rtbForceConstant', FloatParam, default=10,
                      expertLevel=LEVEL_ADVANCED,
                      label='Interaction force constant',
                      help='This is the RTB block size for the RTB NMA method. \n'
                           'When calculating the normal modes, aminoacids are grouped\n'
                           'into blocks of this size that are moved translationally\n'
                           'and rotationally together.')
              
        form.addSection(label='Animation')        
        form.addParam('amplitude', FloatParam, default=50,
                      label="Amplitude") 
        form.addParam('nframes', IntParam, default=10,
                      expertLevel=LEVEL_ADVANCED,
                      label='Number of frames')
        form.addParam('downsample', FloatParam, default=1,
                      expertLevel=LEVEL_ADVANCED,
                      # condition=isEm
                      label='Downsample pseudoatomic structure',
                      help='Downsample factor 2 means removing one half of the atoms or pseudoatoms.')
        form.addParam('pseudoAtomThreshold', FloatParam, default=0,
                      expertLevel=LEVEL_ADVANCED,
                      # cond
                      label='Pseudoatom mass threshold',
                      help='Remove pseudoatoms whose mass is below this threshold. This value should be between 0 and 1.\n'
                           'A threshold of 0 implies no atom removal.')

                                   
    def _insertAllSteps(self):
        # Some steps will differ if the input is a volume or a pdb file
        self.structureEM = self.inputStructure.get().getPseudoAtoms()
        n = self.numberOfModes.get()
        
        # Link the input
        inputFn = self.inputStructure.get().getFileName()
        localFn = self._getPath(basename(inputFn))
        self._insertFunctionStep('copyPdbStep', inputFn, localFn, self.structureEM)
        
        # Construct string for relative-absolute cutoff
        # This is used to detect when to reexecute a step or not
        cutoffStr=''
        if self.cutoffMode == NMA_CUTOFF_REL:
            cutoffStr = 'Relative %f'%self.rcPercentage.get()
        else:
            cutoffStr = 'Absolute %f'%self.rc.get()

        # Compute modes
        self.pseudoAtomRadius=1
        if self.structureEM:
            with open(inputFn, 'r') as fh:
                first_line = fh.readline()
                second_line = fh.readline()
                self.pseudoAtomRadius = float(second_line.split()[2])
            self._insertFunctionStep('computeModesStep', self.inputStructure.get().getFileName(), n, cutoffStr)
            self._insertFunctionStep('reformatOutputStep',"pseudoatoms.pdb")
        else:
            if self.cutoffMode == NMA_CUTOFF_REL:
                params = '-i %s --operation distance_histogram %s' % (localFn, self._getExtraPath('atoms_distance.hist'))
                self._insertRunJobStep("xmipp_pdb_analysis", params)
            self._insertFunctionStep('computePdbModesStep', n, self.rtbBlockSize.get(), self.rtbForceConstant.get(), cutoffStr)
            self._insertFunctionStep('reformatPdbOutputStep', n)
            self.PseudoAtomThreshold=0.0
        
        self._insertFunctionStep('qualifyModesStep', n, self.collectivityThreshold.get(), self.structureEM)
        self._insertFunctionStep('animateModesStep', n,
                                 self.amplitude.get(), self.nframes.get(), self.downsample.get(), 
                                 self.pseudoAtomThreshold.get(), self.pseudoAtomRadius)
        self._insertFunctionStep('computeAtomShiftsStep', n)
        self._insertFunctionStep('createOutputStep')
        
    def copyPdbStep(self, inputFn, localFn, isEM):
        """ Copy the input pdb file and also create a link 'atoms.pdb' """
        copyFile(inputFn, localFn)
        if isEM:
            fnOut=self._getPath('pseudoatoms.pdb')
        else:
            fnOut=self._getPath('atoms.pdb')
        if not os.path.exists(fnOut):
            createLink(localFn, fnOut)
        
    def computePdbModesStep(self, numberOfModes, RTBblockSize, RTBForceConstant, cutoffStr):
        rc = self._getRc(self._getExtraPath('atoms_distance.hist'))
                
        self._enterWorkingDir()
        self.runJob('nma_record_info_PDB.py', "%d %d atoms.pdb %f %f" % (numberOfModes, RTBblockSize, rc, RTBForceConstant),
                    env=getNMAEnviron())
        self.runJob("nma_elnemo_pdbmat","",env=getNMAEnviron())
        self.runJob("nma_diagrtb","",env=getNMAEnviron())

        if not exists("diagrtb.eigenfacs"):
            msg = "Modes cannot be computed. Check the number of modes you asked to compute and/or consider "
            msg += "increasing cut-off distance. The maximum number of modes allowed by the method for atomic "
            msg += "normal mode analysis is 6 times the number of RTB blocks but the protocol allows only up "
            msg += "to 200 modes as 20-100 modes are usually enough. If the number of modes is below the minimum "
            msg += "between 200 and 6 times the number of RTB blocks, consider increasing cut-off distance."
            self._printWarnings(redStr(msg) + '\n')
        self.runJob("rm","-f *.dat_run diagrtb.dat pdbmat.xyzm pdbmat.sdijf pdbmat.dat")
        
        self._leaveWorkingDir()
        
    def reformatPdbOutputStep(self, numberOfModes):
        self._enterWorkingDir()
        
        makePath('modes')
        Natoms = self._countAtoms("atoms.pdb")
        fhIn = open('diagrtb.eigenfacs')
        fhAni = open('vec_ani.txt','w')
        
        for n in range(numberOfModes):
            # Skip two lines
            fhIn.readline()
            fhIn.readline()
            fhOut=open('modes/vec.%d'%(n+1),'w')
            for i in range(Natoms):
                line=fhIn.readline()
                fhOut.write(line)
                fhAni.write(line.rstrip().lstrip()+" ")
            fhOut.close()
            if n!=(numberOfModes-1):
                fhAni.write("\n")
        fhIn.close()
        fhAni.close()
        self.runJob("nma_prepare_for_animate.py","",env=getNMAEnviron())
        cleanPath("vec_ani.txt")
        moveFile('vec_ani.pkl', 'extra/vec_ani.pkl')

        self._leaveWorkingDir()
        
    def animateModesStep(self, numberOfModes,amplitude,nFrames,downsample,pseudoAtomThreshold,pseudoAtomRadius):
        makePath(self._getExtraPath('animations'))
        self._enterWorkingDir()
        
        if self.structureEM:
            fn = "pseudoatoms.pdb"
            self.runJob("nma_animate_pseudoatoms.py","%s extra/vec_ani.pkl 7 %d %f extra/animations/animated_mode %d %d %f"%\
                      (fn,numberOfModes,amplitude,nFrames,downsample,pseudoAtomThreshold),env=getNMAEnviron())
        else:
            fn="atoms.pdb"
            self.runJob("nma_animate_atoms.py","%s extra/vec_ani.pkl 7 %d %f extra/animations/animated_mode %d"%\
                      (fn,numberOfModes,amplitude,nFrames),env=getNMAEnviron())
        
        for mode in range(7,numberOfModes+1):
            fnAnimation = join("extra", "animations", "animated_mode_%03d" % mode)
            fhCmd=open(fnAnimation+".vmd",'w')
            fhCmd.write("mol new %s.pdb\n" % self._getPath(fnAnimation))
            fhCmd.write("animate style Loop\n")
            fhCmd.write("display projection Orthographic\n")
            if self.structureEM:
                fhCmd.write("mol modcolor 0 0 Beta\n")
                fhCmd.write("mol modstyle 0 0 Beads %f 8.000000\n"%(pseudoAtomRadius))
            else:
                fhCmd.write("mol modcolor 0 0 Index\n")
                #fhCmd.write("mol modstyle 0 0 Beads 1.000000 8.000000\n")
                fhCmd.write("mol modstyle 0 0 NewRibbons 1.800000 6.000000 2.600000 0\n")
            fhCmd.write("animate speed 0.5\n")
            fhCmd.write("animate forward\n")
            fhCmd.close();
        
        self._leaveWorkingDir()
  
    def computeAtomShiftsStep(self, numberOfModes):
        fnOutDir = self._getExtraPath("distanceProfiles")
        makePath(fnOutDir)
        maxShift=[]
        maxShiftMode=[]
        
        for n in range(7, numberOfModes+1):
            fnVec = self._getPath("modes", "vec.%d" % n)
            if exists(fnVec):
                fhIn = open(fnVec)
                md = xmipp.MetaData()
                atomCounter = 0
                for line in fhIn:
                    x, y, z = map(float, line.split())
                    d = math.sqrt(x*x+y*y+z*z)
                    if n==7:
                        maxShift.append(d)
                        maxShiftMode.append(7)
                    else:
                        if d>maxShift[atomCounter]:
                            maxShift[atomCounter]=d
                            maxShiftMode[atomCounter]=n
                    atomCounter+=1
                    md.setValue(xmipp.MDL_NMA_ATOMSHIFT,d,md.addObject())
                md.write(join(fnOutDir,"vec%d.xmd" % n))
                fhIn.close()
        md = xmipp.MetaData()
        for i, _ in enumerate(maxShift):
            fnVec = self._getPath("modes", "vec.%d" % (maxShiftMode[i]+1))
            if exists(fnVec):
                objId = md.addObject()
                md.setValue(xmipp.MDL_NMA_ATOMSHIFT, maxShift[i],objId)
                md.setValue(xmipp.MDL_NMA_MODEFILE, fnVec, objId)
        md.write(self._getExtraPath('maxAtomShifts.xmd'))
                                                      
    def createOutputStep(self):
        fnSqlite = self._getPath('modes.sqlite')
        nmSet = SetOfNormalModes(filename=fnSqlite)

        md = xmipp.MetaData(self._getPath('modes.xmd'))
        row = XmippMdRow()
        
        for objId in md:
            row.readFromMd(md, objId)
            nmSet.append(rowToMode(row))
        inputPdb = self.inputStructure.get()
        nmSet.setPdb(inputPdb)
        self._defineOutputs(outputModes=nmSet)
        self._defineSourceRelation(self.inputStructure, nmSet)
