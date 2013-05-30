#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Protocol for Normal Mode analysis of atomic and EM structures
# Author: Carlos Oscar Sanchez Sorzano, May 2013
#         Slavica Jonic
#

import glob,os,re,sys,shutil,time
from protlib_base import *
from config_protocols import protDict
from protlib_utils import runJob, which, runShowJ, runChimera, runVMD
from protlib_filesystem import changeDir, createLink, getExt, moveFile, createDir, deleteFile
from xmipp import MetaData, MDL_X, MDL_COUNT, MDL_NMA_MODEFILE, MDL_ORDER, MDL_ENABLED, MDL_NMA_COLLECTIVITY, MDL_NMA_SCORE, SingleImgSize

class ProtNMA(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.nma.name, scriptname, project)
        self.Import = 'from protocol_nma import *'    

    def createFilenameTemplates(self):
        return {
            'normal_modes': "%(WorkingDir)s/normal_modes.xmd",
            'pseudoatoms':  '%(WorkingDir)s/pseudoatoms',
            'extra_pseudoatoms':  '%(ExtraDir)s/pseudoatoms'
            }

    def defineSteps(self):
        self.insertStep('createDir',path=self.ExtraDir)
        if self.StructureType=="EM":
            # Link the input
            self.insertStep("createLink",source=self.InputStructure,dest=self.workingDirPath(os.path.basename(self.InputStructure)))
            
            # Mask or link input structure
            fnMask=""
            if self.MaskMode=="Threshold":
                fnMask=self.extraPath('mask.vol')
                self.insertRunJobStep("xmipp_transform_threshold", params="-i %s -o %s --select below %f --substitute binarize"%\
                                      (self.InputStructure,fnMask,self.Threshold),verifyFiles=[fnMask])
            else:
                fnMask=self.MaskFile
            
            # Convert to pseudoatoms
            fnOut=self.getFilename("pseudoatoms")
            params="-i %s -o %s --sigma %f --targetError %f --sampling_rate %f -v 2 --intensityColumn Bfactor"%\
                (self.InputStructure,fnOut,self.PseudoAtomRadius*self.Sampling,self.PseudoAtomTarget/100.0,self.Sampling)
            if fnMask!="":
                params+=" --mask binary_file %s"%fnMask
            self.insertRunJobStep("xmipp_volume_to_pseudoatoms", params=params,verifyFiles=[fnOut+".pdb"])
            self.insertStep('moveFile',source=fnOut+"_approximation.vol",dest=self.getFilename("extra_pseudoatoms")+"_approximation.vol")
            self.insertStep('moveFile',source=fnOut+"_distance.hist",dest=self.getFilename("extra_pseudoatoms")+"_distance.hist")
            self.insertRunJobStep("rm", params=fnOut+"_*")
            self.insertStep('computeModes',WorkingDir=self.WorkingDir, NumberOfModes=self.NumberOfModes,
                            CutoffMode=self.CutoffMode, Rc=self.Rc, RcPercentage=self.RcPercentage)
            self.insertStep('reformatOutput',WorkingDir=self.WorkingDir)
            self.insertStep('qualifyModes',WorkingDir=self.WorkingDir,NumberOfModes=self.NumberOfModes)
            self.insertStep('animateModes',WorkingDir=self.WorkingDir,LastMode=self.NumberOfModes,
                            Amplitude=self.Amplitude,NFrames=self.NFrames,Downsample=self.Downsample,
                            PseudoAtomThreshold=self.PseudoAtomThreshold, PseudoAtomRadius=self.PseudoAtomRadius,
                            Sampling=self.Sampling)
            if self.StructureType=="EM":
                self.insertStep('generateChimeraScript',WorkingDir=self.WorkingDir,MaskMode=self.MaskMode,Threshold=self.Threshold, \
                                InputStructure=self.InputStructure, PseudoAtomRadius=self.PseudoAtomRadius, Sampling=self.Sampling)
    
    def summary(self):
        message=[]
        message.append('NMA of [%s]'%self.InputStructure)
        message.append('%d flexible modes: [%s]'%(self.NumberOfModes,self.workingDirPath('modes.xmd')))
        return message
    
    def validate(self):
        errors = []
        if self.MaskMode=="Binary mask" and not os.path.exists(self.MaskFile):
            errors.append(self.MaskFile+" does not exist")
        if self.CutoffMode=="Relative" and self.StructureType=="PDB":
            errors.append("Cut-off in percentage is not meant for PDB structures. It can only be used with EM.")
            # TODO: relative cut-off for PDBs
        if which("diag_arpack")=="":
            errors.append("Cannot find diag_arpack in the PATH")
        return errors
    
    def visualize(self):
        plots = [k for k in ['DisplayPseudoAtom'
                            ,'DisplayPseudoApproximation'
                            ,'DisplayModes'
                           ] if self.ParamsDict[k]]
        if len(plots):
            self.launchVisualize(plots)
        if self.DisplaySingleMode!="":
            os.system("vmd -e %s"%self.extraPath("animations/animated_mode_%03d.vmd"%self.DisplaySingleMode))

    def visualizeVar(self, varName):
        self.launchVisualize([varName])
    
    def launchVisualize(self,selectedPlots):
        def doPlot(plotName):
            return plotName in selectedPlots
        if doPlot('DisplayPseudoAtom'):
            runChimera(self.workingDirPath("chimera.cmd"))
        if doPlot('DisplayPseudoApproximation'):
            runShowJ(self.InputStructure+" "+self.getFilename("extra_pseudoatoms")+"_approximation.vol")
        if doPlot('DisplayModes'):
            runShowJ(self.workingDirPath('modes.xmd'))

def countAtoms(fnPDB):
    fh = open(fnPDB, 'r')
    Natoms=0
    for line in fh:
        if line.find('ATOM')!=-1:
            Natoms+=1
    fh.close()
    return Natoms

def computeModes(log, WorkingDir, NumberOfModes, CutoffMode, Rc, RcPercentage):
    if CutoffMode=="Relative":
        MDhist=MetaData(os.path.join(WorkingDir,"extra/pseudoatoms_distance.hist"))
        distances=MDhist.getColumnValues(MDL_X)
        distanceCount=MDhist.getColumnValues(MDL_COUNT)
        
        # compute total number of distances
        Ncounts=0
        for count in distanceCount:
            Ncounts+=count
        
        # Compute threshold
        NcountThreshold=Ncounts*RcPercentage/100.0
        Ncounts=0
        for i in range(len(distanceCount)):
            Ncounts+=distanceCount[i]
            if Ncounts>NcountThreshold:
                Rc=distances[i]
                break
        print("Cut-off distance="+str(Rc))
    currentDir=os.getcwd()
    changeDir(log,WorkingDir)
    runJob(log,"nma_record_info.py","%d pseudoatoms.pdb %d"%(NumberOfModes,int(Rc)))
    runJob(log,"nma_pdbmat.pl","pdbmat.dat")
    runJob(log,"diag_arpack","")
    runJob(log,"rm","diag_arpack.in pdbmat.dat")
    changeDir(log,currentDir)

def reformatOutput(log,WorkingDir):
    currentDir=os.getcwd()
    changeDir(log,WorkingDir)
    Natoms=countAtoms("pseudoatoms.pdb")
    runJob(log,"nma_reformat_vector_foranimate.pl","%d fort.11"%Natoms)
    runJob(log,"cat","vec.1* > vec_ani.txt")
    runJob(log,"rm","-f vec.1*")
    runJob(log,"nma_reformat_vector.pl","%d fort.11"%Natoms)
    createDir(log,"modes")
    runJob(log,"mv","-f vec.* modes")
    runJob(log,"nma_prepare_for_animate.py","")
    runJob(log,"rm","vec_ani.txt fort.11 matrice.sdijf")
    moveFile(log,'vec_ani.pkl','extra/vec_ani.pkl')
    changeDir(log,currentDir)

def qualifyModes(log,WorkingDir,NumberOfModes):
    currentDir=os.getcwd()
    changeDir(log,WorkingDir)
    
    runJob(log,"nma_reformatForElNemo.sh","%d"%NumberOfModes)
    runJob(log,"echo","diag_arpack.eigenfacs | check_modes")
    deleteFile(log,"diag_arpack.eigenfacs")
    
    fh=open("Chkmod.res")
    MDout=MetaData()
    collectivityList=[]
    for n in range(NumberOfModes):
        line=fh.readline()
        collectivity=float(line.split()[1])
        collectivityList.append(collectivity)

        id=MDout.addObject()
        MDout.setValue(MDL_NMA_MODEFILE,"modes/vec.%d"%(n+1),id)
        MDout.setValue(MDL_ORDER,long(n+1),id)
        MDout.setValue(MDL_ENABLED,1,id)
        MDout.setValue(MDL_NMA_COLLECTIVITY,collectivity,id)
    fh.close()
    idxSorted=[i[0] for i in sorted(enumerate(collectivityList), key=lambda x:x[1])]
    i=0;
    for id in MDout:
        MDout.setValue(MDL_NMA_SCORE,float(idxSorted[i]+1+(NumberOfModes-i))/(2.0*NumberOfModes),id)
        i+=1
    MDout.write("modes.xmd")
    deleteFile(log,"Chkmod.res")
    changeDir(log,currentDir)

def animateModes(log,WorkingDir,LastMode,Amplitude,NFrames,Downsample,PseudoAtomThreshold,PseudoAtomRadius,Sampling):
    createDir(log,os.path.join(WorkingDir,"extra/animations"))
    currentDir=os.getcwd()
    changeDir(log,WorkingDir)
    runJob(log,"nma_animate.py","pseudoatoms.pdb extra/vec_ani.pkl 7 %d %f extra/animations/animated_mode %d %d %f"%\
                  (LastMode,Amplitude,NFrames,Downsample,PseudoAtomThreshold))
    
    for mode in range(7,LastMode+1):
        fnAnimation="extra/animations/animated_mode_%03d"%mode
        fhCmd=open(fnAnimation+".vmd",'w')
        fhCmd.write("mol new "+os.path.join(WorkingDir,fnAnimation)+".pdb\n")
        fhCmd.write("animate style Loop\n")
        fhCmd.write("display projection Orthographic\n")
        fhCmd.write("mol modcolor 0 0 Beta\n")
        fhCmd.write("mol modstyle 0 0 Beads %f 8.000000\n"%(PseudoAtomRadius*Sampling))
        fhCmd.write("animate speed 0.5\n")
        fhCmd.write("animate forward\n")
        fhCmd.close();
    
    changeDir(log,currentDir)

def generateChimeraScript(log,WorkingDir,MaskMode,Threshold,InputStructure,PseudoAtomRadius,Sampling):
    fhCmd=open(os.path.join(WorkingDir,"chimera.cmd"),'w')
    fhCmd.write("open pseudoatoms.pdb\n")
    fhCmd.write("rangecol bfactor,a 0 white 1 red\n")
    fhCmd.write("setattr a radius "+str(PseudoAtomRadius*Sampling)+"\n")
    fhCmd.write("represent sphere\n")
    fhCmd.write("open "+InputStructure+"\n")
    if MaskMode!="Threshold":
        Threshold=0.01
    [xdim,ydim,zdim,ndim]=SingleImgSize(InputStructure)
    fhCmd.write("volume #1 level "+str(Threshold)+" transparency 0.5 voxelSize "+str(Sampling)+" originIndex "+str(xdim/2)+"\n")
    fhCmd.close()
