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
    _references = ['[[http://www.ncbi.nlm.nih.gov/pubmed/23671335][Nogales-Cadenas, et.al, NAR (2013)]]']
    
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputStructure', PointerParam, label="Input structure", important=True, 
                      pointerClass='PdbFile, SetOfVolumes',
                      help='You can choose either a PDB atomic structure or EM volume')  
        form.addParam('maskMode', EnumParam, choices=['none', 'threshold', 'file'], 
                      default=NMA_MASK_NONE, 
                      label='Mask mode', display=EnumParam.DISPLAY_COMBO,
                      help='')        
        form.addParam('maskThreshold', FloatParam, default=0.01, 
                      condition='maskMode==%d' % NMA_MASK_THRE,
                      label='Threshold value',
                      help='TODO: More help?')
        form.addParam('volumeMask', PointerParam, pointerClass='VolumeMask',
                      label='Mask volume', condition='maskMode==%d' % NMA_MASK_FILE,
                      )          
        form.addParam('pseudoAtomRadius', IntParam, default=1, 
                      label='Pseudoatom radius (vox)',
                      help='Pseudoatoms are defined as Gaussians whose \n'
                           'standard deviation is this value in voxels') 
        form.addParam('pseudoAtomTarget', FloatParam, default=5, 
                      condition='maskMode==%d' % NMA_MASK_THRE,
                      label='Threshold value',
                      help='This value is a percentage (between 0.001 and 100) \n'
                           'specifying how fine you want to approximate the EM \n'
                           'volume by the pseudoatomic structure. Lower values \n'
                           'imply lower approximation error, and consequently, \n'
                           'more pseudoatoms.')        
              
        form.addSection(label='Normal Mode Analysis')
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
                      help='TODO: More help.')
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
                           'into blocks of this size that are moved translationally  \n'
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
              
        form.addSection(label='Normal Mode Analysis')        
        form.addParam('amplitud', FloatParam, default=50,
                      label="Amplitud") 
        form.addParam('nframes', IntParam, default=10,
                      expertLevel=LEVEL_ADVANCED,
                      label='Number of frames')
        form.addParam('downsample', FloatParam, default=1,
                      expertLevel=LEVEL_ADVANCED,
                      label='Interaction force constant',
                      help='This is the RTB block size for the RTB NMA method. \n'
                           'When calculating the normal modes, aminoacids are grouped\n'
                           'into blocks of this size that are moved translationally  \n'
                           'and rotationally together.')
        form.addParam('pseudoAtomThreshold', FloatParam, default=0,
                      expertLevel=LEVEL_ADVANCED,
                      label='Threshold on collectivity',
                      help='Remove pseudoatoms whose mass is below this threshold. \n'
                           'This value should be between 0 and 1.')
        
                      
        form.addParallelSection(threads=1, mpi=8)    
             
    def _printWarnings(self, *lines):
        """ Print some warning lines to 'warnings.xmd', 
        the function should be called inside the working dir."""
        fWarn = open("warnings.xmd",'w')
        for l in lines:
            print >> fWarn, l
        fWarn.close()
        
    def _defineSteps(self):
        # Some steps will differ if the input is a volume or a pdb file
        self.structureEM = not isinstance(self.inputStructure.get(), PdbFile)
        n = self.numberOfModes.get()
        
        if self.structureEM:
            inputStructure = self.inputStructure.get().getFirstItem()
            fnMask = self._insertMaskStep()
            self.sampling = inputStructure.getSamplingRate()
            self._insertFunctionStep('convertToPseudoAtomsStep', 
                                     inputStructure.getFileName(), fnMask)
            self._insertFunctionStep('computeModesStep', n)
            self._insertFunctionStep('reformatOutputStep')
            self._insertFunctionStep('createChimeraScriptStep')
        else:
            # Link the input
            inputFn = self.inputStructure.get().getFileName()
            localFn = self._getPath(basename(inputFn))
            self._insertFunctionStep('copyPdbStep', inputFn, localFn)
            
            # Compute modes
            if self.cutoffMode == NMA_CUTOFF_REL:
                params = '-i %s --operation distance_histogram %s' % (localFn, self._getExtraPath('atoms_distance.hist'))
                self._insertRunJobStep("xmipp_pdb_analysis", params)
            self._insertFunctionStep('computePdbModesStep', n, self.rtbBlockSize.get(), self.rtbForceConstant.get())
            self._insertFunctionStep('reformatPdbOutputStep', n)
            self.PseudoAtomThreshold=0.0
        
        self._insertFunctionStep('qualifyModesStep', n, self.collectivityThreshold.get())
        self._insertFunctionStep('animateModesStep', n,
                                 self.amplitud.get(), self.nframes.get(), self.downsample.get(), 
                                 self.pseudoAtomThreshold.get(), self.pseudoAtomRadius.get())
        self._insertFunctionStep('computeAtomShiftsStep', n)
        self._insertFunctionStep('createOutputStep')
        
    def _insertMaskStep(self):
        """ Check the mask selected and insert the necessary steps.
        Return the mask filename if needed.
        """
        fnMask = ''
        if self.maskMode == NMA_MASK_THRE:
            fnMask = self._getExtraPath('mask.vol')
            maskParams = '-i %s -o %s --select below %f --substitute binarize' % (self.inputStructure.get().getFirstItem().getFileName(), fnMask, self.maskThreshold.get())
            self._insertRunJobStep('xmipp_transform_threshold', maskParams)
        elif self.maskMode == NMA_MASK_FILE:
            fnMask = self.volumeMask.get().getFileName()
        return fnMask
        
    def convertToPseudoAtomsStep(self, inputFn, fnMask):
        prefix = 'pseudoatoms'
        outputFn = self._getPath(prefix)
        sampling = self.sampling
        sigma = sampling * self.pseudoAtomRadius.get() 
        targetErr = self.pseudoAtomTarget.get()
        params = "-i %(inputFn)s -o %(outputFn)s --sigma %(sigma)f "
        params += "--targetError %(targetErr)f --sampling_rate %(sampling)f -v 2 --intensityColumn Bfactor"
        if fnMask:
            params += " --mask binary_file %(fnMask)s"
        self.runJob(None, "xmipp_volume_to_pseudoatoms", params=params % locals())
        for suffix in ["_approximation.vol", "_distance.hist"]:
            moveFile(self._getPath(prefix+suffix), self._getExtraPath(prefix+suffix))
        cleanPattern(self._getPath(prefix+'_*'))
     
    def computeModesStep(self, numberOfModes):
        rc = self._getRc('pseudoatoms')
        self._enterWorkingDir()
        self.runJob(None, "nma_record_info.py","%d pseudoatoms.pdb %d" % (numberOfModes, rc))
        self.runJob(None, "nma_pdbmat.pl","pdbmat.dat")
        self.runJob(None, "nma_diag_arpack","")
        if not exists("fort.11"):
            self._printWarnings(redStr("Modes cannot be computed. Check the number of modes you asked to compute and/or consider increasing cut-off distance. The maximum number of modes allowed by the method for pseudoatomic normal mode analysis is 3 times the number of pseudoatoms but the protocol allows only up to 200 modes as 20-100 modes are usually enough.  If the number of modes is below the minimum between 200 and 3 times the number of pseudoatoms, consider increasing cut-off distance."))
        cleanPath("diag_arpack.in", "pdbmat.dat")
        self._leaveWorkingDir()
        
    def _getRc(self, prefix):
        if self.cutoffMode == NMA_CUTOFF_REL:
            rc = self._computeCutoff(self._getExtraPath('%s_distance.hist' % prefix), self.rcPercentage.get())
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
        self.runJob(None,"nma_reformat_vector_foranimate.pl","%d fort.11" % n)
        self.runJob(None,"cat","vec.1* > vec_ani.txt")
        self.runJob(None,"rm","-f vec.1*")
        self.runJob(None,"nma_reformat_vector.pl","%d fort.11" % n)
        makePath("modes")
        self.runJob(None,"mv","-f vec.* modes")
        self.runJob(None,"nma_prepare_for_animate.py","")
        self.runJob(None,"rm","-f vec_ani.txt fort.11 matrice.sdijf")
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
    
    def createChimeraScriptStep(self):
        radius = self.sampling * self.pseudoAtomRadius.get() 
        input = self.inputStructure.get().getFirstItem()
        inputFn = input.getFileName()
        fhCmd = open(self._getPath("chimera.cmd"),'w')
        fhCmd.write("open pseudoatoms.pdb\n")
        fhCmd.write("rangecol bfactor,a 0 white 1 red\n")
        fhCmd.write("setattr a radius %f\n" % radius)
        fhCmd.write("represent sphere\n")
        fhCmd.write("open %s\n" % inputFn)
        
        threshold = 0.01
        if self.maskMode == NMA_MASK_THRE:
            self.maskThreshold.get()
        xdim, _, _, _ = input.getDim()
        origin = xdim / 2
        fhCmd.write("volume #1 level %f transparency 0.5 voxelSize %f originIndex %d\n" % (threshold, self.sampling, origin))
        fhCmd.close()
        
    def copyPdbStep(self, inputFn, localFn):
        """ Copy the input pdb file and also create a link 'atoms.pdb' """
        copyFile(inputFn, localFn)
        createLink(localFn, self._getPath('atoms.pdb'))
        
    def computePdbModesStep(self, numberOfModes, RTBblockSize, RTBForceConstant):
        rc = self._getRc('atoms')
                
        self._enterWorkingDir()
        
        self.runJob(None,"nma_record_info_PDB.py","%d %d atoms.pdb %f %f" % (numberOfModes,RTBblockSize,rc,RTBForceConstant))
        self.runJob(None,"nma_elnemo_pdbmat","")
        self.runJob(None,"nma_diagrtb","")

        if not exists("diagrtb.eigenfacs"):
            msg = "Modes cannot be computed. Check the number of modes you asked to compute and/or consider "
            msg += "increasing cut-off distance. The maximum number of modes allowed by the method for atomic "
            msg += "normal mode analysis is 6 times the number of RTB blocks but the protocol allows only up "
            msg += "to 200 modes as 20-100 modes are usually enough. If the number of modes is below the minimum "
            msg += "between 200 and 6 times the number of RTB blocks, consider increasing cut-off distance."
            self._printWarnings(redStr(msg) + '\n')
        self.runJob(None,"rm","-f *.dat_run diagrtb.dat pdbmat.xyzm pdbmat.sdijf pdbmat.dat")
        
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
            self.runJob(None,"nma_reformatForElNemo.sh","%d" % numberOfModes)
            fnDiag = "diag_arpack.eigenfacs"
            
        self.runJob(None,"echo","%s | nma_check_modes" % fnDiag)
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
            self.runJob(None,"nma_animate_pseudoatoms.py","%s extra/vec_ani.pkl 7 %d %f extra/animations/animated_mode %d %d %f"%\
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
                fhCmd.write("mol modstyle 0 0 Beads %f 8.000000\n"%(pseudoAtomRadius*self.sampling))
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
        if self.structureEM:
            pdb = PdbFile(self._getPath('pseudoatoms.pdb'), pseudoatoms=True)
            self._defineOutputs(outputPdb=pdb)
        modes = NormalModes(filename=self._getPath('modes.xmd'))
        self._defineOutputs(outputModes=modes)

    def _summary(self):
        summary = []
        return summary
    
    def _validate(self):
        validateMsgs = []
        return validateMsgs
    
