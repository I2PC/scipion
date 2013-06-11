#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# General script for Xmipp-based particle picking
#
# Author: Carlos Oscar Sorzano, August, 2011
#
#------------------------------------------------------------------------------------------------

from config_protocols import protDict
from os.path import exists
from protlib_base import *
from protlib_filesystem import createLink, deleteFiles, replaceFilenameExt
from protlib_utils import runJob
from protocol_particle_pick import launchParticlePickingGUI, getTemplateFiles, PM_READONLY, \
    PM_REVIEW
from xmipp import MetaData, MD_APPEND, MDL_IMAGE, \
    MDL_PICKING_PARTICLE_SIZE, MDL_PICKING_MICROGRAPH_STATE, MDL_ENABLED, \
    MDL_COST, MDL_MICROGRAPH, MDValueRange, getBlocksInMetaDataFile
import glob

# Create a GUI automatically from a selfile of micrographs
class ProtParticlePickingAuto(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.particle_pick_auto.name, scriptname, project)
        self.Import = "from protocol_particle_pick_auto import *"
        self.setPreviousRun(self.SupervisedRun)
        self.pickingDir = self.PrevRun.WorkingDir
        self.keysToImport = ['micrographs', 'config', 'acquisition']
       
        self.inputFilename(*self.keysToImport)
        self.inputProperty('TiltPairs', 'MicrographsMd', 'Fast')
        self.inputProperty('TiltPairs', 'MicrographsMd')
        self.micrographs = self.getFilename('micrographs')
        self.model="model"
        
        
    def createFilenameTemplates(self):
        return {
            'training': join('%(ExtraDir)s', '%(model)s_training.txt'),
                 'pos': join('%(ExtraDir)s', '%(micrograph)s.pos'),
#                    'mask': join('%(ExtraDir)s', '%(model)s_mask.xmp'),
                'pca': join('%(ExtraDir)s', '%(model)s_pca_model.stk'),
                'rotpca': join('%(ExtraDir)s', '%(model)s_rotpca_model.stk'),
                'svm': join('%(ExtraDir)s', '%(model)s_svm.txt'),
                'svm2': join('%(ExtraDir)s', '%(model)s_svm2.txt'),
                'average': join('%(ExtraDir)s', '%(model)s_particle_avg.xmp'),
                'config': join('%(ExtraDir)s','config.xmd')
                }
        
   
       
    def defineSteps(self):
        self.insertStep("createDir",verifyfiles=[self.ExtraDir],path=self.ExtraDir)
        md = MetaData(self.Input['config'])
        for objId in md:
            self.particleSizeForAuto = md.getValue(MDL_PICKING_PARTICLE_SIZE, objId)
        filesToImport = [self.Input[k] for k in self.keysToImport]
        filesToImport += getTemplateFiles(self.PrevRun)
        
        filesToImport.append(self.PrevRun.getFilename('training', model=self.model))
        filesToImport.append(self.PrevRun.getFilename('pca', model=self.model))
        filesToImport.append(self.PrevRun.getFilename('rotpca', model=self.model))
        filesToImport.append(self.PrevRun.getFilename('svm', model=self.model))
        #filesToImport.append(self.PrevRun.getFilename('svm2', model=model))
        filesToImport.append(self.PrevRun.getFilename('average', model=self.model))
        filesToImport.append(self.PrevRun.getFilename('config', model=self.model))
        self.insertImportOfFiles(filesToImport)

        md = MetaData(self.Input['micrographs'])
        particleSize = self.particleSizeForAuto
        modelRoot = self.extraPath(self.model)
        for objId in md:
            # Get micrograph path and name
            path = md.getValue(MDL_MICROGRAPH, objId)
            micrographName = removeBasenameExt(path)
            proceed = True
            fnPos = self.PrevRun.getFilename('pos', micrograph=micrographName)
            if xmippExists(fnPos):
                mdheader = MetaData("header@" + fnPos)
                for idheader in mdheader:
                    state = mdheader.getValue(MDL_PICKING_MICROGRAPH_STATE, idheader)
                if state != "Available":
                    proceed = False
                    break
            if proceed:
                oroot = self.extraPath(micrographName)
                cmd = "-i %(path)s --particleSize %(particleSize)d --model %(modelRoot)s --outputRoot %(oroot)s --mode autoselect" % locals()
                if self.Fast:
                    cmd += " --fast "
                self.insertParallelRunJobStep("xmipp_micrograph_automatic_picking", cmd)
                                             
        self.insertStep('gatherResults', WorkingDir=self.WorkingDir, PickingDir=self.pickingDir)
        self.insertStep('deleteTempFiles', ExtraDir=self.ExtraDir)

    def summary(self):
        
        summary = ["Supervised picking RUN: <%s> " % self.SupervisedRun,
                   "Input directory: [%s] " % self.pickingDir]
        autoFiles = glob.glob(self.workingDirPath("extra/*auto.pos"))
        if len(autoFiles) > 0:
            Nparticles = 0
            Nmicrographs = len(autoFiles)
            first = True
            maxTime = 0
            minTime = 0
            for f in autoFiles:
                md = MetaData(f)
                Nparticles += md.size()
                if first:
                    first = False
                    minTime = os.stat(f).st_mtime
                    maxTime = minTime
                else:
                    t = os.stat(f).st_mtime
                    minTime = min(minTime, t)
                    maxTime = max(maxTime, t)
            summary.append("<%d> particles automatically picked from <%d> micrographs in <%d> minutes" % (Nparticles, Nmicrographs, int((maxTime - minTime) / 60.0)))
        else:
            
                fnExtractList = self.getFilename('extract_list', model=self.model)
                msg = "model: "
                if os.path.exists(fnExtractList):
                    MD = MetaData("mic.*@" + fnExtractList)
                    MDauto = MetaData()
                    MDauto.importObjects(MD, MDValueRange(MDL_COST, 0., 1.))
                    Nparticles = MDauto.size()
                    #minTime = os.stat(self.workingDirPath(family+"_mask.xmp")).st_mtime
                    #maxTime = os.stat(fnExtractList).st_mtime
                    #Nmicrographs = len(getBlocksInMetaDataFile(fnExtractList))
                    #msg += "<%d> particles automatically picked from <%d> micrographs in <%d> minutes"%(Nparticles,Nmicrographs,int((maxTime-minTime)/60.0))
                summary.append(msg)
        return summary
    
    def validate(self):
        errors = []
        return errors
    
    def visualize(self):
        mode = PM_REVIEW + " %s"
        for f in glob.glob(self.extraPath("*extract_list.xmd")):
            launchParticlePickingGUI(None, self.Input['micrographs'], self.ExtraDir, mode % f, self.TiltPairs, self.Memory)
            
            

def gatherResults(log, WorkingDir, PickingDir):
    md = MetaData(getProtocolFilename("micrographs",WorkingDir=WorkingDir))
    mdpos = MetaData()
    mdposAux = MetaData()
    fnExtractList = getProtocolFilename("extract_list", ExtraDir=join(WorkingDir,"extra"))
    particlesblock = 'particles@'
    for objId in md:
        fullname = md.getValue(MDL_MICROGRAPH, objId)
        name = os.path.split(replaceFilenameExt(fullname, ''))[1]        
        
        fn = join(PickingDir, "extra", name + ".pos")
        mdpos.clear()
        print "Looking for "+fn
        if exists(fn):
            try:
                print "Reading: "+ fn          
                mdpos.read(particlesblock + fn)
                mdpos.setValueCol(MDL_COST, 2.0)
                print "Found : "+str(mdpos.size())
            except:
                # Skip this metadata if it is corrupted
                pass
        else:
            print "    It does not exist"
        for path in [PickingDir, WorkingDir]:
            fn = join(path, "extra", name + "_auto.pos")
            print "Looking for "+fn
            if exists(fn):
                try:
                    print "Reading: %s" % fn          
                    mdposAux.read(particlesblock + fn)
                    mdposAux.removeDisabled();
                    mdposAux.removeLabel(MDL_ENABLED)
                    print "Found : "+str(mdposAux.size())
                    mdpos.unionAll(mdposAux) 
                    print "Altogether : "+str(mdpos.size())
                except:
                    # Skip this metadata if it is corrupted
                    pass
            else:
                print "    It does not exist"
        # Append alphanumeric prefix to help identifying the block s
        mdpos.write("mic_%s@%s" % (name, fnExtractList), MD_APPEND)

def deleteTempFiles(log, ExtraDir):
    deleteFiles(log, glob.glob(join(ExtraDir, "*_auto.pos")), True)

#		
# Main
#     
if __name__ == '__main__':
    protocolMain(ProtParticlePickingAuto)
