#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
#
# General script for Xmipp-based pre-processing of single-particles: 
#  - phase flipping
#  - extraction of particles
#  - normalization
#  - sort_junk

# Author: Carlos Oscar, August 2011
#
from protlib_base import *
import xmipp
import os
from protlib_utils import runJob
from protlib_filesystem import deleteFile

class ProtExtractParticles(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.extract_particles.name, scriptname, project)
        self.Import = 'from protocol_extract_particles import *'
        # COSS: Falta incluir una ejecucion de automaticos
        # COSS: Falta tiltpairs
        self.pickingRunname=self.PickingRun.replace(protDict.particle_pick.name,"")
        if self.pickingRunname[0]=="_":
            self.pickingRunname=self.pickingRunname[1:]
        self.pickingDir= os.path.join(protDict.particle_pick.dir,self.pickingRunname)
        self.familyFile = os.path.join(self.pickingDir, "families.xmd")

    def defineSteps(self):
        families=xmipp.MetaData(self.familyFile)
        self.familyList=[]
        for id in families:
            familyName=families.getValue(xmipp.MDL_PICKING_FAMILY,id)
            particleSize=families.getValue(xmipp.MDL_PICKING_PARTICLE_SIZE,id)
            self.familyList.append((familyName,particleSize))
            familyDir=os.path.join(self.WorkingDir,familyName)
            self.Db.insertStep('createDir',verifyfiles=[familyDir],path=familyDir)
        
        idMPI=self.Db.insertStep('runStepGapsMpi',passDb=True, script=self.scriptName, NumberOfMpi=self.NumberOfMpi)
        verifyFiles=[]        
        mD=xmipp.MetaData(os.path.join(self.pickingDir,"micrographs.sel"))
        for id in mD:
            micrograph=mD.getValue(xmipp.MDL_IMAGE,id)
            micrographName=os.path.split(os.path.split(micrograph)[0])[1]
            posFile=os.path.join(self.pickingDir,micrographName+".pos")
            if os.path.exists(posFile):
                parent_id=None
                micrographToExtract=micrograph
                if self.DoFlip:
                    ctf=mD.getValue(xmipp.MDL_CTFMODEL,id)
                    micrographToExtract=os.path.join(self.TmpDir,micrographName+"_flipped.xmp")
                    parent_id=self.Db.insertStep('phaseFlip',verifyfiles=[micrographToExtract],execution_mode=SqliteDb.EXEC_GAP,
                                      micrograph=micrograph,ctf=ctf,fnOut=micrographToExtract)
                (tasks,outputFiles)=self.whichTasks(posFile,micrographName)
                verifyFiles+=outputFiles
                parent_id=self.Db.insertStep('extractTasks',verifyfiles=outputFiles,execution_mode=SqliteDb.EXEC_GAP,parent_step_id=parent_id,
                                      micrographToExtract=micrographToExtract,tasks=tasks,
                                      doNorm=self.doNorm,doLog=self.DoLog,doInvert=self.DoInvert,
                                      bgRadius=self.BackGroundRadius, doRemoveDust=self.DoRemoveDust,
                                      dustRemovalThreshold=self.DustRemovalThreshold)
                if self.DoFlip:
                    self.Db.insertStep('deleteFile',execution_mode=SqliteDb.EXEC_GAP,parent_step_id=parent_id,filename=micrographToExtract,verbose=True)
        self.Db.updateVerifyFiles(idMPI,verifyFiles)

        # COSS Falta unir en un unico selfile (por familia) y sort_junk
                
    def validate(self):
        errors = []
        if not os.path.exists(self.familyFile):
            errors.append("Cannot find "+self.familyFile)
        fnMicrographs=os.path.join(self.pickingDir,"micrographs.sel")
        if not os.path.exists(fnMicrographs):
            errors.append("Cannot find "+fnMicrographs)
        else:
            mD=xmipp.MetaData(fnMicrographs)
            if self.DoFlip and not mD.containsLabel(xmipp.MDL_CTFMODEL):
                errors.append(fnMicrographs+" does not contain CTF information for phase flipping")
        return errors

    def whichTasks(self,posFile,micrographName):
        fileList=[]
        tasks=[]
        for family in self.familyList:
            familyName=family[0]
            particleSize=family[1]
            blockName="%s@%s"%(familyName,posFile)
            mD=xmipp.MetaData(blockName)
            if mD.size()>0:
                fnOut=os.path.join(self.WorkingDir,familyName,micrographName+".stk")
                tasks.append((blockName,fnOut,particleSize))
                fileList.append(fnOut)
        return (tasks,fileList)

def phaseFlip(log,micrograph,ctf,fnOut):
    runJob(log,"xmipp_ctf_phase_flip"," -i "+micrograph+" --ctf "+ctf+" -o "+fnOut)

def extractTasks(log,micrographToExtract,tasks,
                    doNorm,doLog,doInvert,bgRadius,doRemoveDust,dustRemovalThreshold):
    for task in tasks:
        (blockName,fnOut,particleSize)=task
        extract(log,micrographToExtract,blockName,particleSize,fnOut,
                doNorm,doLog,doInvert,bgRadius,doRemoveDust,dustRemovalThreshold)

def extract(log,micrographToExtract,blockName,particleSize,fnOut,
            doNorm,doLog,doInvert,bgRadius,doRemoveDust,dustRemovalThreshold):
        # Extract 
        rootname=os.path.splitext(fnOut)[0]
        arguments="-i "+micrographToExtract+" --pos "+blockName+" --oroot "+rootname+" --Xdim "+str(particleSize)
        if doInvert:
            arguments+=" --invert"
        if doLog:
            arguments+=" --log"
        runJob(log,"xmipp_micrograph_scissor",arguments)
        
        # Normalize 
        if doNorm:
            if bgRadius==0:
                bgRadius=int(particleSize/2)
            arguments="-i "+fnOut+' --method Ramp --background circle '+str(bgRadius)
            if doRemoveDust:
                arguments+=' --thr_black_dust -' + str(dustRemovalThreshold)+' --thr_white_dust ' + str(dustRemovalThreshold)
            runJob(log,"xmipp_transform_normalize",arguments)
