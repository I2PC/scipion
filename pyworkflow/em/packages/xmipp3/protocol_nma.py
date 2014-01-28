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


from pyworkflow.em import *  
from pyworkflow.utils import *  
import xmipp
from convert import createXmippInputImages, readSetOfVolumes, createXmippInputVolumes
from pyworkflow.protocol.constants import LEVEL_EXPERT, LEVEL_ADVANCED

#from xmipp3 import XmippProtocol
NMA_MASK_NONE = 0
NMA_MASK_THRE = 1
NMA_MASK_FILE = 2

NMA_CUTOFF_ABS = 0
NMA_CUTOFF_REL = 1    
        
class XmippProtNMA(EMProtocol):
    """ Protocol for flexible analysis using NMA. """
    _label = 'ml3d'
    _reference = ['[[http://www.ncbi.nlm.nih.gov/pubmed/23671335][Nogales-Cadenas, et.al, NAR (2013)]]'
                  ]
    
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
                      label="Cut-off distance", condition='cutoffMode==%d' % NMA_CUTOFF_ABS,
                      help='Atoms or pseudoatoms beyond this distance will not interact.')
        form.addParam('rcPercentage', FloatParam, default=95,
                      label="Cut-off distance", condition='cutoffMode==%d' % NMA_CUTOFF_ABS,
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
        form.addParam('PseudoAtomThreshold', FloatParam, default=0,
                      expertLevel=LEVEL_ADVANCED,
                      label='Threshold on collectivity',
                      help='Remove pseudoatoms whose mass is below this threshold. \n'
                           'This value should be between 0 and 1.')
        
                      
        form.addParallelSection(threads=1, mpi=8)    
             
    def _printWarnings(self, *lines):
        """ Print some warning lines to 'warnings.xmd' """
        fWarn = open(self._getPath("warnings.xmd"),'w')
        for l in lines:
            print >> fWarn, l
        fWarn.close()
        
    def _defineSteps(self):
        # Some steps will differ if the input is a volume or a pdb file
        inputStructure = self.inputStructure.get()
        if isinstance(inputStructure, Volume):
            fnMask = self._insertMaskStep()
            self._insertFunctionStep('convertToPseudoAtomsStep', 
                                     inputStructure.getFileName(), fnMask)
            self._insertFunctionStep('computeModesStep')
            self._insertFunctionStep('reformatOutputStep')
            self._insertFunctionStep('createChimeraScriptStep')
        else:
            pass
        
    def _insertMaskStep(self):
        """ Check the mask selected and insert the necessary steps.
        Return the mask filename if needed.
        """
        fnMask = ''
        if self.maskMode == NMA_MASK_THRE:
            fnMask = self._getExtraPath('mask.vol')
            maskParams = '-i %s -o %s --select below %f --substitute binarize' % (inputStructure.getFileName(), fnMask, self.maskThreshold.get())
            self._insertRunJobStep('xmipp_transform_threshold', maskParams)
        elif self.maskMode == NMA_MASK_FILE:
            fnMask = self.volumeMask.get().getFileName()
        return fnMask
        
    def convertToPseudoAtomsStep(self, inputFn, fnMask):
        prefix = 'pseudoatoms'
        outputFn = self._getPath(prefix)
        sampling = self.inputStructure.getSamplingRate()
        sigma = sampling * self.pseudoAtomRadius.get() 
        targetErr = self.pseudoAtomTarget.get()
        params = "-i %(inputFn)s -o %(outputFn)s --sigma %(sigma)f "
        params += "--targetError %(targetErr)f --sampling_rate %(sampling)f -v 2 --intensityColumn Bfactor"
        if fnMask:
            params += " --mask binary_file %(fnMask)s"
        self.insertRunJobStep("xmipp_volume_to_pseudoatoms", params=params % locals())
        for suffix in ["_approximation.vol", "_distance.hist"]:
            moveFile(self._getPath(prefix+suffix), self._getExtraPath(prefix+suffix))
        cleanPattern(self._getPath(prefix+'_*'))
     
    def computeModesStep(self):
        if self.cutoffMode == NMA_CUTOFF_REL:
            rc = self._computeCutoff(self._getExtraPath('pseudoatoms_distance.hist'), self.rcPercentage.get())
            
        self._enterWorkingDir()
        runJob(None, "nma_record_info.py","%d pseudoatoms.pdb %d" % (self.numberOfModes.get(), rc))
        runJob(None, "nma_pdbmat.pl","pdbmat.dat")
        runJob(None, "nma_diag_arpack","")
        if not os.path.exists("fort.11"):
            self._printWarnings(redStr("Modes cannot be computed. Check the number of modes you asked to compute and/or consider increasing cut-off distance. The maximum number of modes allowed by the method for pseudoatomic normal mode analysis is 3 times the number of pseudoatoms but the protocol allows only up to 200 modes as 20-100 modes are usually enough.  If the number of modes is below the minimum between 200 and 3 times the number of pseudoatoms, consider increasing cut-off distance."))
        cleanPath("diag_arpack.in", "pdbmat.dat")
        self._leaveWorkingDir()
        
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
        self._printWarnings(msg)

        return rc    
    
    def reformatOutputStep(self):
        self._enterWorkingDir()
        n = self._countAtoms("pseudoatoms.pdb")
        runJob(log,"nma_reformat_vector_foranimate.pl","%d fort.11" % n)
        runJob(log,"cat","vec.1* > vec_ani.txt")
        runJob(log,"rm","-f vec.1*")
        runJob(log,"nma_reformat_vector.pl","%d fort.11" % n)
        makePath("modes")
        runJob(log,"mv","-f vec.* modes")
        runJob(log,"nma_prepare_for_animate.py","")
        runJob(log,"rm","vec_ani.txt fort.11 matrice.sdijf")
        moveFile(log,'vec_ani.pkl','extra/vec_ani.pkl')
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
        sampling = self.inputStructure.getSamplingRate()
        radius = sampling * self.pseudoAtomRadius.get() 
        inputFn = self.inputStructure.getFileName()
        fhCmd = open(self._getPath("chimera.cmd"),'w')
        fhCmd.write("open pseudoatoms.pdb\n")
        fhCmd.write("rangecol bfactor,a 0 white 1 red\n")
        fhCmd.write("setattr a radius %f\n" % radius)
        fhCmd.write("represent sphere\n")
        fhCmd.write("open %s\n" % inputFn)
        
        threshold = 0.01
        if self.maskMode == NMA_MASK_THRE:
            self.maskThreshold.get()
        xdim, _, _, _ = self.inputStructure.getDim()
        origin = xdim / 2
        fhCmd.write("volume #1 level %f transparency 0.5 voxelSize %f originIndex %d\n" % (threshold, sampling, origin))
        fhCmd.close()
                                                    
    def createOutput(self):
        #lastIter = self._lastIteration()
        lastIter = 'iter%03d' % self._lastIteration()
        md = xmipp.MetaData(self._getExtraPath(lastIter, 'iter_volumes.xmd'))
        md.addItemId()
        fn = self._getPath('output_volumes.xmd')
        md.write('Volumes@%s' % fn)
        volumes = self._createSetOfVolumes()
        volumes.setSamplingRate(self.inputImages.get().getSamplingRate())
        readSetOfVolumes(fn, volumes)
        volumes.write()
        self._defineOutputs(outputVolumes=volumes)

    def _summary(self):
        summary = []
        return summary
    
    def _validate(self):
        validateMsgs = []
        return validateMsgs
    
