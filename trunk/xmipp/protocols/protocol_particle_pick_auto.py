#!/usr/bin/env python
#------------------------------------------------------------------------------------------------
# General script for Xmipp-based particle picking
#
# Author: Carlos Oscar Sorzano, August, 2011
#
#------------------------------------------------------------------------------------------------

from config_protocols import protDict
from protlib_base import *
from protlib_utils import runJob
from protlib_filesystem import createLink, deleteFiles, replaceFilenameExt
from xmipp import MetaData, MD_APPEND,MDL_IMAGE, MDL_PICKING_FAMILY, \
                  MDL_PICKING_MICROGRAPH_FAMILY_STATE, MDL_PICKING_PARTICLE_SIZE,\
                  MDL_ENABLED, MDL_COST, MDL_MICROGRAPH
import glob
from os.path import exists

# Create a GUI automatically from a selfile of micrographs
class ProtParticlePickingAuto(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.particle_pick_auto.name, scriptname, project)
        self.Import="from protocol_particle_pick_auto import *"
        pickingProt = self.getProtocolFromRunName(self.PickingRun)
        self.pickingDir = pickingProt.WorkingDir
        self.pickingMicrographs = pickingProt.getFilename('micrographs')
        self.pickingFamilies = pickingProt.getFilename('families')
        self.micrographs = self.getFilename('micrographs')
        self.families = self.getEquivalentFilename(pickingProt, self.pickingFamilies)

        # Families to process
        self.familiesForAuto = []
        self.particleSizeForAuto = []
        md = MetaData(self.pickingFamilies)
        for objId in md:
            family = md.getValue(MDL_PICKING_FAMILY, objId)
            if exists(join(self.pickingDir,family+"_training.txt")):
                self.familiesForAuto.append(family)
                particleSize = md.getValue(MDL_PICKING_PARTICLE_SIZE, objId)
                self.particleSizeForAuto.append(particleSize)


    def createFilenameTemplates(self):
        return {
                'training': join('%(dir)s', '%(family)s_training.txt'),
                    'mask': join('%(dir)s', '%(family)s_mask.xmp'),
                     'pos': join('%(dir)s', '%(micrograph)s.pos')
                }
        
    def defineSteps(self):
        self.insertStep('createLink', source=self.pickingMicrographs, dest=self.micrographs)
        self.insertStep('createLink', source=self.pickingFamilies, dest=self.families)
        for family in self.familiesForAuto:
            getFilename = lambda k, d: self.getFilename(k, dir=d, family=family)
            self.insertStep('createLink',source=getFilename('training', self.pickingDir),
                            dest=getFilename('training', self.WorkingDir))
            self.insertStep('createLink', source=getFilename('mask', self.pickingDir),
                            dest=getFilename('mask', self.WorkingDir))

        md = MetaData(self.pickingMicrographs)
        for i, family in enumerate(self.familiesForAuto):
            particleSize = self.particleSizeForAuto[i]
            modelRoot = self.workingDirPath(family)
            for objId in md:
                # Get micrograph path and name
                path = md.getValue(MDL_MICROGRAPH, objId)
                micrographName = replaceFilenameExt(os.path.basename(path), '')
                proceed = True
                fnPos = self.getFilename('pos', dir=self.pickingDir, micrograph=micrographName)
                if exists(fnPos):
                    mdaux = MetaData("families@" + fnPos)
                    for idaux in mdaux:
                        familyAux = mdaux.getValue(MDL_PICKING_FAMILY,idaux)
                        if family == familyAux:
                            state = mdaux.getValue(MDL_PICKING_MICROGRAPH_FAMILY_STATE,idaux)
                            if state != "Available":
                                proceed = False
                                break
                if proceed:
                    oroot = join(self.WorkingDir, micrographName)
                    cmd = "-i %(path)s --particleSize %(particleSize)d --model %(modelRoot)s --outputRoot %(oroot)s --mode autoselect" % locals()
                    self.insertParallelRunJobStep("xmipp_micrograph_automatic_picking", cmd)
                                             
        for family in self.familiesForAuto:
            self.insertStep('gatherResults', Family=family, WorkingDir=self.WorkingDir, PickingDir=self.pickingDir)
        self.insertStep('deleteTempFiles',WorkingDir=self.WorkingDir)

    def summary(self):
        summary = []
        summary.append("Manual picking run: <%s> " % self.PickingRun)
        summary.append("Automatic picking of the following models: "+",".join(self.familiesForAuto))
        autoFiles = glob.glob(self.workingDirPath("*auto.pos"))
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
                minTime = min(minTime,t)
                maxTime = max(maxTime,t)
        summary.append("<%d> particles automatically picked from <%d> micrographs in <%d> minutes"%(Nparticles,Nmicrographs,int((maxTime-minTime)/60.0)))
        return summary
    
    def validate(self):
        errors = []
        return errors
    
    def visualize(self):
        for f in glob.glob(self.workingDirPath("*extract_list.xmd")):
            params="-i %s -o %s --mode review %s &" % (self.pickingMicrographs, self.WorkingDir, f)
            os.system("xmipp_micrograph_particle_picking " + params)

def gatherResults(log,Family,WorkingDir,PickingDir):
    familyFn = lambda fn: "%s@%s" % (Family, fn)
    md = MetaData(join(WorkingDir,"micrographs.xmd"))
    mdpos = MetaData()
    mdposAux = MetaData()
    fnExtractList = join(WorkingDir, Family + "_extract_list.xmd")
    for objId in md:
        fullname = md.getValue(MDL_MICROGRAPH, objId)
        name = os.path.split(replaceFilenameExt(fullname, ''))[1]        
        
        fn = join(PickingDir, name + ".pos")
        mdpos.clear()
        if exists(fn):
            try:            
                mdpos.read(familyFn(fn))
                mdpos.setValueCol(MDL_COST,2.0)
            except:
                pass
        for path in [PickingDir, WorkingDir]:
            fn = join(path, name + "_auto.pos")
            if exists(fn):
                try:
                    mdposAux.read( familyFn(fn))
                    mdposAux.removeDisabled();
                    mdposAux.removeLabel(MDL_ENABLED)
                    mdpos.unionAll(mdposAux) 
                except:
                    pass
        # Append alphanumeric prefix to help identifying the block 
        mdpos.write("mic_%s@%s" % (name, fnExtractList), MD_APPEND)

def deleteTempFiles(log,WorkingDir):
    deleteFiles(log,glob.glob(join(WorkingDir,"*_auto.pos")),True)

#		
# Main
#     
if __name__ == '__main__':
    protocolMain(ProtParticlePicking)
