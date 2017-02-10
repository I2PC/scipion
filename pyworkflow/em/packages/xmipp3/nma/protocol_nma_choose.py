# **************************************************************************
# *
# * Authors:  Carlos Oscar Sanchez Sorzano (coss@cnb.csic.es), March 2014
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

from protocol_convert_to_pseudoatoms_base import *
from protocol_nma_base import *
from pyworkflow.utils.path import createLink, cleanPath
from pyworkflow.protocol.params import BooleanParam
from xmipp import MetaData, MDL_NMA, MDL_ENABLED, MDL_NMA_MINRANGE, MDL_NMA_MAXRANGE

class XmippProtNMAChoose(XmippProtConvertToPseudoAtomsBase, XmippProtNMABase):
    """ Protocol for choosing a volume to construct an NMA analysis """
    _label = 'choose NMA'
    def __init__(self, **args):
        XmippProtConvertToPseudoAtomsBase.__init__(self, **args)
        XmippProtNMABase.__init__(self, **args)
        self.stepsExecutionMode = STEPS_PARALLEL
    
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputStructures', PointerParam, label="Input volumes", important=True, 
                      pointerClass='SetOfVolumes')
        form.addParam('alignVolumes', BooleanParam, label="Align volumes", default=False,
                      help="Align deformed PDBs to volume to maximize match")
        XmippProtConvertToPseudoAtomsBase._defineParams(self,form)
        form.addParallelSection(threads=4, mpi=1)    

        form.addSection(label='Normal Mode Analysis')
        XmippProtNMABase._defineParamsCommon(self,form)
             
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        inputStructures = self.inputStructures.get()
        self.sampling = inputStructures.getSamplingRate()
        filenames=[]
        for inputStructure in inputStructures:
            filenames.append(getImageLocation(inputStructure))
          
        deps = []  
        for volCounter in range(1,len(filenames)+1):
            fnIn=filenames[volCounter-1]
            prefix="_%02d"%volCounter
            fnMask = self._insertMaskStep(fnIn, prefix)
            
            self._insertFunctionStep('convertToPseudoAtomsStep', inputStructure, fnIn, fnMask, prefix, prerequisites=deps)
            parentId=self._insertFunctionStep('computeNMAStep',self._getPath("pseudoatoms%s.pdb"%prefix), prefix)
            deps=[]
            for volCounter2 in range(1,len(filenames)+1):
                if volCounter2!=volCounter:
                    args="-i %s --pdb %s --modes %s --sampling_rate %f -o %s --fixed_Gaussian %f --opdb %s"%\
                         (filenames[volCounter2-1],self._getPath("pseudoatoms%s.pdb"%prefix), \
                          self._getPath("modes%s.xmd"%prefix),self.sampling,\
                          self._getExtraPath('alignment_%02d_%02d.xmd'%(volCounter,volCounter2)),\
                          self.sampling*self.pseudoAtomRadius.get(),
                          self._getExtraPath('alignment_%02d_%02d.pdb'%(volCounter,volCounter2)))
                    if self.alignVolumes.get():
                        args+=" --alignVolumes"
                    stepId=self._insertRunJobStep("xmipp_nma_alignment_vol",args,prerequisites=[parentId])
                    deps.append(stepId)
            
        self._insertFunctionStep('evaluateDeformationsStep',prerequisites=deps)
        

    #--------------------------- Step functions --------------------------------------------
    def convertToPseudoAtomsStep(self, inputStructure, fnIn, fnMask, prefix):
        XmippProtConvertToPseudoAtomsBase.convertToPseudoAtomsStep(self, fnIn, fnMask, prefix)
        self.createChimeraScriptStep(inputStructure, fnIn, prefix)
        createLink(self._getPath("pseudoatoms%s.pdb"%prefix),self._getPath("pseudoatoms.pdb"))

    def computeNMAStep(self, fnIn, prefix):
        cutoffStr=''
        if self.cutoffMode == NMA_CUTOFF_REL:
            cutoffStr = 'Relative %f'%self.rcPercentage.get()
        else:
            cutoffStr = 'Absolute %f'%self.rc.get()
        self.computeModesStep(fnIn, self.numberOfModes.get(), cutoffStr)
        self.reformatOutputStep("pseudoatoms.pdb")
        self.qualifyModesStep(self.numberOfModes.get(), self.collectivityThreshold.get(), True)
        fnModes=self._getPath("modes.xmd")
        fnModesPrefix=self._getPath("modes%s.xmd"%prefix)
        self.runJob("xmipp_metadata_utilities",
                    "-i %s --operate modify_values \"nmaModeFile=replace(nmaModeFile,'/modes/','/modes%s/')\" -o %s"%
                    (fnModes,prefix,fnModesPrefix))
        self.runJob("mv","%s %s"%(self._getPath('modes'),self._getPath('modes%s'%prefix)))

        # Remove intermediate files
        cleanPath(self._getPath("pseudoatoms.pdb"), fnModes, self._getExtraPath('vec_ani.pkl'))
    
    def evaluateDeformationsStep(self):
        N = self.inputStructures.get().getSize()
        import numpy
        distances=numpy.zeros([N,N])
        for volCounter in range(1,N+1):
            pdb1=open(self._getPath('pseudoatoms_%02d.pdb'%volCounter)).readlines()
            for volCounter2 in range(1,N+1):
                if volCounter!=volCounter2:
                    davg=0.
                    Navg=0.
                    pdb2=open(self._getExtraPath('alignment_%02d_%02d.pdb'%(volCounter,volCounter2))).readlines()
                    for i in range(len(pdb1)):
                        line1=pdb1[i]
                        if line1.startswith("ATOM"):
                            line2=pdb2[i]
                            x1=float(line1[30:37])
                            y1=float(line1[38:45])
                            z1=float(line1[46:53])
                            x2=float(line2[30:37])
                            y2=float(line2[38:45])
                            z2=float(line2[46:53])
                            dx=x1-x2
                            dy=y1-y2
                            dz=z1-z2
                            d=math.sqrt(dx*dx+dy*dy+dz*dz)
                            davg+=d
                            Navg+=1
                    if Navg>0:
                        davg/=Navg
                    distances[volCounter-1,volCounter2-1]=davg
        distances=0.5*(distances+numpy.transpose(distances))
        numpy.savetxt(self._getPath('distances.txt'),distances)
        distances1D=numpy.mean(distances,axis=0)
        print("Average distance to rest of volumes=",distances1D)
        imin=numpy.argmin(distances1D)
        print("The volume in the middle is pseudoatoms_%02d.pdb"%(imin+1))
        createLink(self._getPath("pseudoatoms_%02d.pdb"%(imin+1)),self._getPath("pseudoatoms.pdb"))
        createLink(self._getPath("modes_%02d.xmd"%(imin+1)),self._getPath("modes.xmd"))
        createLink(self._getExtraPath("pseudoatoms_%02d_distance.hist"%(imin+1)),self._getExtraPath("pseudoatoms_distance.hist"))

        # Measure range
        minDisplacement= 1e38*numpy.ones([self.numberOfModes.get(),1])
        maxDisplacement=-1e38*numpy.ones([self.numberOfModes.get(),1])
        mdNMA=MetaData(self._getPath("modes.xmd"))
        for volCounter in range(1,N+1):
            if volCounter!=imin+1:
                md=MetaData(self._getExtraPath("alignment_%02d_%02d.xmd"%(imin+1,volCounter)))
                displacements=md.getValue(MDL_NMA, md.firstObject())
                idx1=0
                idx2=0
                for idRow in mdNMA:
                    if mdNMA.getValue(MDL_ENABLED,idRow)==1:
                        minDisplacement[idx2]=min(minDisplacement[idx2],displacements[idx1])
                        maxDisplacement[idx2]=max(maxDisplacement[idx2],displacements[idx1])
                        idx1+=1
                    else:
                        minDisplacement[idx2]=0
                        maxDisplacement[idx2]=0
                    idx2+=1
        idx2=0
        for idRow in mdNMA:
            mdNMA.setValue(MDL_NMA_MINRANGE,float(minDisplacement[idx2]),idRow)
            mdNMA.setValue(MDL_NMA_MAXRANGE,float(maxDisplacement[idx2]),idRow)
            idx2+=1
        mdNMA.write(self._getPath("modes.xmd"))

        # Create output
        volCounter=0
        for inputStructure in self.inputStructures.get():
            if volCounter==imin:
                print("The corresponding volume is %s"%(getImageLocation(inputStructure)))
                finalStructure=inputStructure
                break
            volCounter+=1

        pdb = PdbFile(self._getPath('pseudoatoms.pdb'), pseudoatoms=True)
        self._defineOutputs(outputPdb=pdb)
        modes = NormalModes(filename=self._getPath('modes.xmd'))
        self._defineOutputs(outputModes=modes)
        
        self._defineSourceRelation(self.inputStructures, self.outputPdb)
        # ToDo: the self.outputPdb should be a Pointer, not an object
#         self._defineSourceRelation(self.outputPdb, self.outputModes)

    #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        summary = []
        summary.append('Pseudoatom radius (voxels): %f'%self.pseudoAtomRadius.get())
        summary.append('Approximation target error (%%): %f'%self.pseudoAtomTarget.get())
        return summary

    def _methods(self):
        summary = []
#        summary.append('We converted the volume %s into a pseudoatomic representation with Gaussian atoms (sigma=%f A and a target error'\
#                       ' of %f%%) [Nogales2013].'%(self.inputStructure.get().getNameId(),
#                                     self.pseudoAtomRadius.get()*self.inputStructure.get().getSamplingRate(),
#                                     self.pseudoAtomTarget.get()));
#        if self.hasAttribute('outputPdb'):
#            summary.append('We refer to the pseudoatomic model as %s.'%self.outputPdb.getNameId())
        return summary

    def _citations(self):
        return ['Nogales2013','Jin2014']
        