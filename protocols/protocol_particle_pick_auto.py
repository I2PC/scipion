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
from protlib_particles import getMetadataWithPickedParticles
from protlib_filesystem import createLink, deleteFiles, replaceFilenameExt
from protlib_utils import runJob
from protocol_particle_pick import launchParticlePickingGUI, getTemplateFiles, PM_READONLY, PM_REVIEW
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
        self.keysToImport = ['micrographs', 'config', 'acquisition', 'templates']
       
        self.inputFilename(*self.keysToImport)
        self.inputProperty('TiltPairs', 'MicrographsMd', 'Fast')
        self.inputProperty('TiltPairs', 'MicrographsMd')
        self.micrographs = self.getFilename('micrographs')
        self.model="model"
        
    def createFilenameTemplates(self):
        return {
            'training': join('%(ExtraDir)s', '%(model)s_training.txt'),
                 'pos': join('%(ExtraDir)s', '%(micrograph)s.pos'),
                 'templates': join('%(ExtraDir)s', 'templates.stk'),
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
        filesToImport += [self.PrevRun.getFilename(k, model=self.model) for k in ['training', 'pca', 'rotpca', 'svm', 'average', 'config']]
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
                blocks = getBlocksInMetaDataFile(fnPos)
                copy = True
                if 'header' in blocks:
                    mdheader = MetaData("header@" + fnPos)
                    state = mdheader.getValue(MDL_PICKING_MICROGRAPH_STATE, mdheader.firstObject())
                    if state == "Available":
                        copy = False
                if copy:
                    # Copy manual .pos file of this micrograph
                    self.insertCopyFile(fnPos, self.getFilename('pos', micrograph=micrographName))
                    proceed = False
            if proceed:
                oroot = self.extraPath(micrographName)
                cmd = "-i %(path)s --particleSize %(particleSize)d --model %(modelRoot)s --outputRoot %(oroot)s --mode autoselect" % locals()
                if self.Fast:
                    cmd += " --fast "
                self.insertParallelRunJobStep("xmipp_micrograph_automatic_picking", cmd)
                                             
    def summary(self):
        summary = ["Input directory: [%s] " % self.pickingDir]
        totalCount=0
        for fnPos in glob.glob(self.extraPath('*.pos')):
            md=getMetadataWithPickedParticles(fnPos)
            totalCount+=md.size()
        summary.append("Number of particles: %d"%totalCount)
        return summary
    
    def validate(self):
        errors = []
        return errors
    
    def visualize(self):
        mode = PM_REVIEW
        launchParticlePickingGUI(None, self.Input['micrographs'], self.ExtraDir, mode, self.TiltPairs, self.Memory)
            
#		
# Main
#     
if __name__ == '__main__':
    protocolMain(ProtParticlePickingAuto)
