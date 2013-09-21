#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Protocol for Normal Mode analysis of atomic and EM structures
# Author: Carlos Oscar Sanchez Sorzano, May 2013
#         Slavica Jonic
#

import glob,os,re,sys,shutil,time,math
from protlib_base import *
from config_protocols import protDict
from protlib_utils import runJob, which, runShowJ, runChimera, runVMD
from protlib_filesystem import changeDir, createLink, getExt, moveFile, createDir, deleteFile
from xmipp import MetaData, MDL_X, MDL_COUNT, MDL_NMA_MODEFILE, MDL_ORDER, MDL_ENABLED, MDL_NMA_COLLECTIVITY, MDL_NMA_SCORE, \
    MDL_NMA_ATOMSHIFT, getImageSize

class ProtNMA(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.nma.name, scriptname, project)
        self.Import = 'from protocol_nma import *'    

    def createFilenameTemplates(self):
        return {
            'normal_modes': "%(WorkingDir)s/normal_modes.xmd",
            'pseudoatoms':  '%(WorkingDir)s/pseudoatoms',
            'extra_pseudoatoms':  '%(ExtraDir)s/pseudoatoms',
            'atoms':  '%(WorkingDir)s/atoms',
            'extra_atoms':  '%(ExtraDir)s/atoms'
            }

    def defineSteps(self):
        self.insertStep('createDir',path=self.ExtraDir)
        if self.StructureType=="EM":
            # Link the input
            fnLocalInputStructure=self.workingDirPath(os.path.basename(self.InputStructure))
            self.insertStep("createLink",source=self.InputStructure,dest=fnLocalInputStructure)
            
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
                (self.InputStructure,fnOut,self.PseudoAtomRadius*self.Sampling,self.PseudoAtomTarget,self.Sampling)
            if fnMask!="":
                params+=" --mask binary_file %s"%fnMask
            self.insertRunJobStep("xmipp_volume_to_pseudoatoms", params=params,verifyFiles=[fnOut+".pdb"])
            self.insertStep('moveFile',source=fnOut+"_approximation.vol",dest=self.getFilename("extra_pseudoatoms")+"_approximation.vol")
            self.insertStep('moveFile',source=fnOut+"_distance.hist",dest=self.getFilename("extra_pseudoatoms")+"_distance.hist")
            self.insertRunJobStep("rm", params=fnOut+"_*")
            
            # Compute modes
            self.insertStep('computeModes',WorkingDir=self.WorkingDir, NumberOfModes=self.NumberOfModes,
                            CutoffMode=self.CutoffMode, Rc=self.Rc, RcPercentage=self.RcPercentage)
            self.insertStep('reformatOutput',WorkingDir=self.WorkingDir)
            
            if self.StructureType=="EM":
                self.insertStep('generateChimeraScript',WorkingDir=self.WorkingDir,MaskMode=self.MaskMode,Threshold=self.Threshold, \
                                InputStructure=fnLocalInputStructure, PseudoAtomRadius=self.PseudoAtomRadius, Sampling=self.Sampling)
        else: # PDB
            # Link the input
            self.insertStep("copyFile",source=self.InputStructure,dest=self.workingDirPath(os.path.basename(self.InputStructure)))
            self.insertStep("createLink",source=self.workingDirPath(os.path.basename(self.InputStructure)),dest=self.workingDirPath("atoms.pdb"))
            
            # Compute modes
            if self.CutoffMode=='Relative':
                params="-i "+self.InputStructure+" --operation distance_histogram "+self.getFilename("extra_atoms")+"_distance.hist"
                self.insertRunJobStep("xmipp_pdb_analysis",params=params)
            self.insertStep('computeModesPDB',WorkingDir=self.WorkingDir, NumberOfModes=self.NumberOfModes,
                            CutoffMode=self.CutoffMode, Rc=self.Rc, RcPercentage=self.RcPercentage, RTBblockSize=self.RTBblockSize,
                            RTBForceConstant=self.RTBForceConstant)
            self.insertStep('reformatOutputPDB',WorkingDir=self.WorkingDir, NumberOfModes=self.NumberOfModes)
            self.PseudoAtomThreshold=0.0
        self.insertStep('qualifyModes',WorkingDir=self.WorkingDir,NumberOfModes=self.NumberOfModes,StructureType=self.StructureType,
                        CollectivityThreshold=self.CollectivityThreshold)
        self.insertStep('animateModes',WorkingDir=self.WorkingDir,LastMode=self.NumberOfModes,
                        Amplitude=self.Amplitude,NFrames=self.NFrames,Downsample=self.Downsample,
                        PseudoAtomThreshold=self.PseudoAtomThreshold, PseudoAtomRadius=self.PseudoAtomRadius,
                        Sampling=self.Sampling,StructureType=self.StructureType)
        self.insertStep('computeAtomShifts',WorkingDir=self.WorkingDir,NumberOfModes=self.NumberOfModes)

    def summary(self):
        message=[]
        message.append('NMA of [%s]'%self.InputStructure)
        message.append('%d flexible modes: [%s]'%(self.NumberOfModes,self.workingDirPath('modes.xmd')))
        return message
    
    def validate(self):
        errors = []
        if self.MaskMode=="Binary mask" and not os.path.exists(self.MaskFile):
            errors.append(self.MaskFile+" does not exist")
        if which("nma_diag_arpack")=="":
            errors.append("Cannot find nma_diag_arpack in the PATH")
        if which("nma_diagrtb")=="":
            errors.append("Cannot find nma_diagrtb in the PATH")
        return errors
    
    def visualize(self):
        plots = [k for k in ['DisplayPseudoAtom'
                            ,'DisplayPseudoApproximation'
                            ,'DisplayModes'
                            ,'DisplayMaxDistanceProfile'
                            ,'DisplayDistanceProfile'
                           ] if self.ParamsDict[k]]
        if len(plots):
            self.launchVisualize(plots)
        if self.DisplaySingleMode!="":
            os.system("vmd -e %s"%self.extraPath("animations/animated_mode_%03d.vmd"%int(self.DisplaySingleMode)))

    def visualizeVar(self, varName):
        self.launchVisualize([varName])
    
    def launchVisualize(self,selectedPlots):
        def doPlot(plotName):
            return plotName in selectedPlots
        if self.StructureType=="EM":
            if doPlot('DisplayPseudoAtom'):
                runChimera(self.workingDirPath("chimera.cmd"))
            if doPlot('DisplayPseudoApproximation'):
                runShowJ(self.InputStructure+" "+self.getFilename("extra_pseudoatoms")+"_approximation.vol")
        if doPlot('DisplayModes'):
            runShowJ(self.workingDirPath('modes.xmd'))
        if doPlot('DisplayMaxDistanceProfile'):
            fnProfile=os.path.join("extra","maxAtomShifts.xmd")
            os.system("xmipp_metadata_plot -i "+self.workingDirPath(fnProfile)+' -y nmaAtomShift --title "Maximum atom shifts" &')
        if doPlot('DisplayDistanceProfile'):
            fnProfile=os.path.join("extra","distanceProfiles","vec%d"%int(self.DisplaySingleMode)+".xmd")
            os.system("xmipp_metadata_plot -i "+self.workingDirPath(fnProfile)+' -y nmaAtomShift --title "Atom shifts for mode '+\
                      self.DisplaySingleMode+'" &')

def countAtoms(fnPDB):
    fh = open(fnPDB, 'r')
    Natoms=0
    for line in fh:
        if line.find('ATOM')!=-1:
            Natoms+=1
    fh.close()
    return Natoms

def computeCutoff(fnHist, RcPercentage):
    MDhist=MetaData(fnHist)
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
    return Rc

def computeModes(log, WorkingDir, NumberOfModes, CutoffMode, Rc, RcPercentage):
    if CutoffMode=="Relative":
        Rc=computeCutoff(os.path.join(WorkingDir,"extra/pseudoatoms_distance.hist"),RcPercentage)
    currentDir=os.getcwd()
    changeDir(log,WorkingDir)
    runJob(log,"nma_record_info.py","%d pseudoatoms.pdb %d"%(NumberOfModes,int(Rc)))
    runJob(log,"nma_pdbmat.pl","pdbmat.dat")
    runJob(log,"nma_diag_arpack","")
    runJob(log,"rm","diag_arpack.in pdbmat.dat")
    changeDir(log,currentDir)

def computeModesPDB(log, WorkingDir, NumberOfModes, CutoffMode, Rc, RcPercentage, RTBblockSize, RTBForceConstant):
    if CutoffMode=="Relative":
        Rc=computeCutoff(os.path.join(WorkingDir,"extra/atoms_distance.hist"),RcPercentage)
    currentDir=os.getcwd()
    changeDir(log,WorkingDir)
    runJob(log,"nma_record_info_PDB.py","%d %d atoms.pdb %f %f"%(NumberOfModes,RTBblockSize,Rc,RTBForceConstant))
    runJob(log,"nma_elnemo_pdbmat","")
    runJob(log,"nma_diagrtb","")
    runJob(log,"rm","*.dat_run diagrtb.dat pdbmat.xyzm pdbmat.sdijf pdbmat.dat")
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

def reformatOutputPDB(log,WorkingDir,NumberOfModes):
    currentDir=os.getcwd()
    changeDir(log,WorkingDir)
    createDir(log,"modes")
    Natoms=countAtoms("atoms.pdb")
    fhIn=open('diagrtb.eigenfacs')
    fhAni=open('vec_ani.txt','w')
    for n in range(NumberOfModes):
        # Skip two lines
        fhIn.readline()
        fhIn.readline()
        fhOut=open('modes/vec.%d'%(n+1),'w')
        for i in range(Natoms):
            line=fhIn.readline()
            fhOut.write(line)
            fhAni.write(line.rstrip().lstrip()+" ")
        fhOut.close()
        if n!=(NumberOfModes-1):
            fhAni.write("\n")
    fhIn.close()
    fhAni.close()
    runJob(log,"nma_prepare_for_animate.py","")
    runJob(log,"rm","vec_ani.txt")
    moveFile(log,'vec_ani.pkl','extra/vec_ani.pkl')
    changeDir(log,currentDir)

def qualifyModes(log,WorkingDir,NumberOfModes,StructureType,CollectivityThreshold):
    currentDir=os.getcwd()
    changeDir(log,WorkingDir)
    
    fnDiag="diagrtb.eigenfacs"
    if StructureType=="EM":
        runJob(log,"nma_reformatForElNemo.sh","%d"%NumberOfModes)
        fnDiag="diag_arpack.eigenfacs"
        
    runJob(log,"echo","%s | nma_check_modes"%fnDiag)
    deleteFile(log,fnDiag)
    
    fh=open("Chkmod.res")
    MDout=MetaData()
    collectivityList=[]
    for n in range(NumberOfModes):
        line=fh.readline()
        collectivity=float(line.split()[1])
        collectivityList.append(collectivity)

        id=MDout.addObject()
        MDout.setValue(MDL_NMA_MODEFILE,os.path.join(WorkingDir,"modes/vec.%d"%(n+1)),id)
        MDout.setValue(MDL_ORDER,long(n+1),id)
        if n>=6:
            MDout.setValue(MDL_ENABLED,1,id)
        else:
            MDout.setValue(MDL_ENABLED,-1,id)
        MDout.setValue(MDL_NMA_COLLECTIVITY,collectivity,id)
        if collectivity<CollectivityThreshold:
            MDout.setValue(MDL_ENABLED,-1,id)
    fh.close()
    idxSorted=[i[0] for i in sorted(enumerate(collectivityList), key=lambda x:x[1])]
    score=[0]*NumberOfModes
    for i in range(NumberOfModes):
       score[i]+=i+1
       score[idxSorted[i]]+=NumberOfModes-i
    i=0
    for id in MDout:
        score_i=float(score[i])/(2.0*NumberOfModes)
        MDout.setValue(MDL_NMA_SCORE,score_i,id)
        i+=1
    MDout.write("modes.xmd")
    deleteFile(log,"Chkmod.res")
    changeDir(log,currentDir)

def animateModes(log,WorkingDir,LastMode,Amplitude,NFrames,Downsample,PseudoAtomThreshold,PseudoAtomRadius,Sampling,StructureType):
    createDir(log,os.path.join(WorkingDir,"extra/animations"))
    currentDir=os.getcwd()
    changeDir(log,WorkingDir)
    if StructureType=="EM":
        fn="pseudoatoms.pdb"
    else:
        fn="atoms.pdb"
    runJob(log,"nma_animate.py","%s extra/vec_ani.pkl 7 %d %f extra/animations/animated_mode %d %d %f"%\
                  (fn,LastMode,Amplitude,NFrames,Downsample,PseudoAtomThreshold))
    
    for mode in range(7,LastMode+1):
        fnAnimation="extra/animations/animated_mode_%03d"%mode
        fhCmd=open(fnAnimation+".vmd",'w')
        fhCmd.write("mol new "+os.path.join(WorkingDir,fnAnimation)+".pdb\n")
        fhCmd.write("animate style Loop\n")
        fhCmd.write("display projection Orthographic\n")
        if StructureType=="EM":
            fhCmd.write("mol modcolor 0 0 Beta\n")
            fhCmd.write("mol modstyle 0 0 Beads %f 8.000000\n"%(PseudoAtomRadius*Sampling))
        else:
            fhCmd.write("mol modcolor 0 0 Index\n")
            fhCmd.write("mol modstyle 0 0 Beads 1.000000 8.000000\n")
        fhCmd.write("animate speed 0.5\n")
        fhCmd.write("animate forward\n")
        fhCmd.close();
    
    changeDir(log,currentDir)

def computeAtomShifts(log,WorkingDir,NumberOfModes):
    fnOutDir=os.path.join(WorkingDir,"extra/distanceProfiles")
    createDir(log,fnOutDir)
    maxShift=[]
    maxShiftMode=[]
    for n in range(7,NumberOfModes+1):
        fhIn=open(os.path.join(WorkingDir,"modes/vec."+str(n)))
        md=MetaData()
        atomCounter=0;
        for line in fhIn:
            coord=line.split()
            x=float(coord[0])
            y=float(coord[1])
            z=float(coord[2])
            d=math.sqrt(x*x+y*y+z*z)
            if n==7:
                maxShift.append(d)
                maxShiftMode.append(7)
            else:
                if d>maxShift[atomCounter]:
                    maxShift[atomCounter]=d
                    maxShiftMode[atomCounter]=n
            atomCounter+=1
            md.setValue(MDL_NMA_ATOMSHIFT,d,md.addObject())
        md.write(os.path.join(fnOutDir,"vec"+str(n)+".xmd"))
        fhIn.close()
    md=MetaData()
    for i in range(len(maxShift)):
        id=md.addObject()
        md.setValue(MDL_NMA_ATOMSHIFT,maxShift[i],id)
        md.setValue(MDL_NMA_MODEFILE,os.path.join(WorkingDir,"modes/vec."+str(maxShiftMode[i]+1)),id)
    md.write(os.path.join(WorkingDir,'extra','maxAtomShifts.xmd'))

def generateChimeraScript(log,WorkingDir,MaskMode,Threshold,InputStructure,PseudoAtomRadius,Sampling):
    fhCmd=open(os.path.join(WorkingDir,"chimera.cmd"),'w')
    fhCmd.write("open pseudoatoms.pdb\n")
    fhCmd.write("rangecol bfactor,a 0 white 1 red\n")
    fhCmd.write("setattr a radius "+str(PseudoAtomRadius*Sampling)+"\n")
    fhCmd.write("represent sphere\n")
    fhCmd.write("open "+os.path.basename(InputStructure)+"\n")
    if MaskMode!="Threshold":
        Threshold=0.01
    [xdim,ydim,zdim,ndim]=getImageSize(InputStructure)
    fhCmd.write("volume #1 level "+str(Threshold)+" transparency 0.5 voxelSize "+str(Sampling)+" originIndex "+str(xdim/2)+"\n")
    fhCmd.close()
