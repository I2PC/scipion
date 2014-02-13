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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************

import math
from glob import glob

from pyworkflow.em import *  
from pyworkflow.utils import * 
from pyworkflow.protocol.constants import LEVEL_EXPERT, LEVEL_ADVANCED
import xmipp

#from xmipp3 import XmippProtocol
NMA_MASK_NONE = 0
NMA_MASK_THRE = 1
NMA_MASK_FILE = 2

NMA_CUTOFF_ABS = 0
NMA_CUTOFF_REL = 1    

        
class XmippProtNMA(EMProtocol):
    """ Protocol for flexible analysis using NMA. """
    _label = 'nma analysis'
    
    def _defineParams(self, form):
        form.addSection(label='Normal Mode Analysis')
        form.addParam('inputStructure', PointerParam, label="Input structure", important=True, 
                      pointerClass='PdbFile',
                      help='The input structure can be an atomic model (true PDB) or a pseudoatomic model\n'
                           '(an EM volume converted into pseudoatoms)')
        form.addParam('numberOfModes', IntParam, default=20,
                      label='Number of modes',
                      help='The maximum number of modes allowed by the method for \n'
                           'atomic normal mode analysis is 6 times the number of  \n'
                           'RTB blocks and for pseudoatomic normal mode analysis 3\n'
                           'times the number of pseudoatoms. However, the protocol\n'
                           'allows only up to 200 modes as 20-100 modes are usually\n'
                           'enough. The number of modes given here should be below \n'
                           'the minimum between these two numbers.')    
        form.addParam('cutoffMode', EnumParam, choices=['absolute', 'relative'],
                      default=NMA_CUTOFF_REL,
                      label='Cut-off mode',
                      help='Absolute distance allows specifying the maximum distance (in Angstroms) for which it\n'
                           'is considered that two atoms are connected. Relative distance allows to specify this distance\n'
                           'as a percentile of all the distances between an atom and its nearest neighbors.')
        form.addParam('rc', FloatParam, default=8,
                      label="Cut-off distance (A)", condition='cutoffMode==%d' % NMA_CUTOFF_ABS,
                      help='Atoms or pseudoatoms beyond this distance will not interact.')
        form.addParam('rcPercentage', FloatParam, default=95,
                      label="Cut-off percentage", condition='cutoffMode==%d' % NMA_CUTOFF_REL,
                      help='The interaction cutoff distance is calculated as the distance\n'
                           'below which is this percentage of interatomic or interpseudoatomic\n'
                           'distances. \n'
                           'Atoms or pseudoatoms beyond this distance will not interact.')      
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
                           'into blocks of this size that are moved translationally  \n'
                           'and rotationally together.')
        form.addParam('collectivityThreshold', FloatParam, default=0.15,
                      label='Threshold on collectivity',
                      help='Collectivity degree is related to the number of atoms or \n'
                           'pseudoatoms that are affected by the mode, and it is normalized\n'
                           'between 0 and 1. Modes below this threshold are deselected in  \n'
                           'the modes metadata file. Set to 0 for no deselection. You can  \n'
                           'always modify the selection manually after the modes metadata  \n'
                           'file is created. The modes metadata file can be used with      \n'
                           'Flexible fitting protocol. Modes 1-6 are always deselected as  \n'
                           'they are related to rigid-body movements.')
              
        form.addSection(label='Animation')        
        form.addParam('amplitude', FloatParam, default=50,
                      label="Amplitud") 
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
                      label='Pseudoatom mass thresold',
                      help='Remove pseudoatoms whose mass is below this threshold. This value should be between 0 and 1.\n'
                           'A threshold of 0 implies no atom removal.')
                      
        form.addParallelSection(threads=1, mpi=8)    
             
    def _printWarnings(self, *lines):
        """ Print some warning lines to 'warnings.xmd', 
        the function should be called inside the working dir."""
        fWarn = open("warnings.xmd",'wa')
        for l in lines:
            print >> fWarn, l
        fWarn.close()
        
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
            rccutoffStr = 'Absolute %f'%self.rc.get()

        # Compute modes
        self.pseudoAtomRadius=1
        if self.structureEM:
            with open(inputFn, 'r') as fh:
                first_line = fh.readline()
                second_line = fh.readline()
                self.pseudoAtomRadius = float(second_line.split()[2])
            self._insertFunctionStep('computeModesStep', n, cutoffStr)
            self._insertFunctionStep('reformatOutputStep')
        else:
            if self.cutoffMode == NMA_CUTOFF_REL:
                params = '-i %s --operation distance_histogram %s' % (localFn, self._getExtraPath('atoms_distance.hist'))
                self._insertRunJobStep("xmipp_pdb_analysis", params)
            self._insertFunctionStep('computePdbModesStep', n, self.rtbBlockSize.get(), self.rtbForceConstant.get(), cutoffStr)
            self._insertFunctionStep('reformatPdbOutputStep', n)
            self.PseudoAtomThreshold=0.0
        
        self._insertFunctionStep('qualifyModesStep', n, self.collectivityThreshold.get())
        self._insertFunctionStep('animateModesStep', n,
                                 self.amplitude.get(), self.nframes.get(), self.downsample.get(), 
                                 self.pseudoAtomThreshold.get(), self.pseudoAtomRadius)
        self._insertFunctionStep('computeAtomShiftsStep', n)
        self._insertFunctionStep('createOutputStep')
        
    def computeModesStep(self, numberOfModes, cutoffStr):
        fnDistanceHist=os.path.join(os.path.split(self.inputStructure.get().getFileName())[0],'extra','pseudoatoms_distance.hist')
        rc = self._getRc(fnDistanceHist)
        self._enterWorkingDir()
        self.runJob("nma_record_info.py","%d pseudoatoms.pdb %d" % (numberOfModes, rc))
        self.runJob("nma_pdbmat.pl","pdbmat.dat")
        self.runJob("nma_diag_arpack","")
        if not exists("fort.11"):
            self._printWarnings(redStr("Modes cannot be computed. Check the number of modes you asked to compute and/or consider increasing cut-off distance. The maximum number of modes allowed by the method for pseudoatomic normal mode analysis is 3 times the number of pseudoatoms but the protocol allows only up to 200 modes as 20-100 modes are usually enough.  If the number of modes is below the minimum between 200 and 3 times the number of pseudoatoms, consider increasing cut-off distance."))
        cleanPath("diag_arpack.in", "pdbmat.dat")
        self._leaveWorkingDir()
        
    def _getRc(self, fnDistanceHist):
        if self.cutoffMode == NMA_CUTOFF_REL:
            rc = self._computeCutoff(fnDistanceHist, self.rcPercentage.get())
        else:
            rc = self.rc.get()
        return rc
        
    def _computeCutoff(self, fnHist, rcPercentage):
        mdHist = xmipp.MetaData(fnHist)
        distances = mdHist.getColumnValues(xmipp.MDL_X)
        distanceCount = mdHist.getColumnValues(xmipp.MDL_COUNT)
        # compute total number of distances
        nCounts = 0
        for count in distanceCount:
            nCounts+=count
        # Compute threshold
        NcountThreshold = nCounts*rcPercentage/100.0
        nCounts = 0
        for i in range(len(distanceCount)):
            nCounts+=distanceCount[i]
            if nCounts>NcountThreshold:
                rc = distances[i]
                break
        msg = "Cut-off distance = %s A" % rc
        print msg
        self._enterWorkingDir()
        self._printWarnings(msg)
        self._leaveWorkingDir()

        return rc    
    
    def reformatOutputStep(self):
        self._enterWorkingDir()
        n = self._countAtoms("pseudoatoms.pdb")
        self.runJob("nma_reformat_vector_foranimate.pl","%d fort.11" % n)
        self.runJob("cat","vec.1* > vec_ani.txt")
        self.runJob("rm","-f vec.1*")
        self.runJob("nma_reformat_vector.pl","%d fort.11" % n)
        makePath("modes")
        self.runJob("mv","-f vec.* modes")
        self.runJob("nma_prepare_for_animate.py","")
        self.runJob("rm","-f vec_ani.txt fort.11 matrice.sdijf")
        moveFile('vec_ani.pkl','extra/vec_ani.pkl')
        self._leaveWorkingDir()
        
    def _countAtoms(self, fnPDB):
        fh = open(fnPDB, 'r')
        n = 0
        for line in fh:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                n += 1
        fh.close()
        return n
    
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
        
        self.runJob("nma_record_info_PDB.py","%d %d atoms.pdb %f %f" % (numberOfModes,RTBblockSize,rc,RTBForceConstant))
        self.runJob("nma_elnemo_pdbmat","")
        self.runJob("nma_diagrtb","")

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
        runJob(log,"nma_prepare_for_animate.py","")
        cleanPath("vec_ani.txt")
        moveFile('vec_ani.pkl', 'extra/vec_ani.pkl')

        self._leaveWorkingDir()
        
    def qualifyModesStep(self, numberOfModes, collectivityThreshold):
        self._enterWorkingDir()
        
        fnVec = glob("modes/vec.*")
        
        if len(fnVec) < numberOfModes:
            msg = "There are only %d modes instead of %d. "
            msg += "Check the number of modes you asked to compute and/or consider increasing cut-off distance."
            msg += "The maximum number of modes allowed by the method for atomic normal mode analysis is 6 times"
            msg += "the number of RTB blocks and for pseudoatomic normal mode analysis 3 times the number of pseudoatoms. "
            msg += "However, the protocol allows only up to 200 modes as 20-100 modes are usually enough. If the number of"
            msg += "modes is below the minimum between these two numbers, consider increasing cut-off distance." 
            self._printWarnings(redStr(msg % (len(fnVec), numberOfModes)))
    
        fnDiag = "diagrtb.eigenfacs"
        
        if self.structureEM:
            self.runJob("nma_reformatForElNemo.sh","%d" % numberOfModes)
            fnDiag = "diag_arpack.eigenfacs"
            
        self.runJob("echo","%s | nma_check_modes" % fnDiag)
        cleanPath(fnDiag)
        
        fh = open("Chkmod.res")
        MDout = xmipp.MetaData()
        collectivityList = []
        
        for n in range(numberOfModes):
            line = fh.readline()
            collectivity = float(line.split()[1])
            collectivityList.append(collectivity)
    
            id = MDout.addObject()
            modefile = self._getPath("modes", "vec.%d" % (n+1))
            MDout.setValue(xmipp.MDL_NMA_MODEFILE, modefile, id)
            MDout.setValue(xmipp.MDL_ORDER, long(n+1), id)
            
            if n >= 6:
                MDout.setValue(xmipp.MDL_ENABLED, 1, id)
            else:
                MDout.setValue(xmipp.MDL_ENABLED, -1, id)
            MDout.setValue(xmipp.MDL_NMA_COLLECTIVITY, collectivity, id)
            
            if collectivity < collectivityThreshold:
                MDout.setValue(xmipp.MDL_ENABLED,-1,id)
        fh.close()
        idxSorted = [i[0] for i in sorted(enumerate(collectivityList), key=lambda x:x[1])]
        score = [0]*numberOfModes
        for i in range(numberOfModes):
            score[i] += i+1
            score[idxSorted[i]] += numberOfModes - i
        i = 0
        for id in MDout:
            score_i = float(score[i])/(2.0*numberOfModes)
            MDout.setValue(xmipp.MDL_NMA_SCORE,score_i,id)
            i+=1
        MDout.write("modes.xmd")
        cleanPath("Chkmod.res")
        self._leaveWorkingDir()

    def animateModesStep(self, numberOfModes,amplitude,nFrames,downsample,pseudoAtomThreshold,pseudoAtomRadius):
        makePath(self._getExtraPath('animations'))
        self._enterWorkingDir()
        
        if self.structureEM:
            fn = "pseudoatoms.pdb"
            self.runJob("nma_animate_pseudoatoms.py","%s extra/vec_ani.pkl 7 %d %f extra/animations/animated_mode %d %d %f"%\
                      (fn,numberOfModes,amplitude,nFrames,downsample,pseudoAtomThreshold))
        else:
            fn="atoms.pdb"
            runJob(log,"nma_animate_atoms.py","%s extra/vec_ani.pkl 7 %d %f extra/animations/animated_mode %d"%\
                      (fn,numberOfModes,amplitude,nFrames))
        
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
                fhCmd.write("mol modstyle 0 0 Beads 1.000000 8.000000\n")
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
            fhIn = open(self._getPath("modes", "vec.%d" % n))
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
            id = md.addObject()
            md.setValue(xmipp.MDL_NMA_ATOMSHIFT, maxShift[i],id)
            md.setValue(xmipp.MDL_NMA_MODEFILE, self._getPath("modes", "vec.%d" % (maxShiftMode[i]+1)), id)
        md.write(self._getExtraPath('maxAtomShifts.xmd'))
                                                      
    def createOutputStep(self):
        modes = NormalModes(filename=self._getPath('modes.xmd'))
        self._defineOutputs(outputModes=modes)
        self._defineSourceRelation(self.inputStructure.get(), self.outputModes)

    def _summary(self):
        summary = []
        return summary
    
    def _validate(self):
        validateMsgs = []
        if which('nma_diag_arpack') is '':
            validateMsgs.append('Check that the programs nma_* are in $XMIPP_HOME/bin')
        return validateMsgs
    
    def _citations(self):
        return ['Nogales2013','Jin2014']
