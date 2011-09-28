#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Protocol for Xmipp-based 2D alignment and classification,
# using maximum-likelihood principles
#
#   Author:  Sjors Scheres, January 2008
#  Updated:  J. M. de la Rosa Trevin July 2011
#

from os.path import join
from protlib_base import XmippProtocol, protocolMain
from config_protocols import protDict
from xmipp import MetaData, Image, MDL_IMAGE

class ProtML3D(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.ml3d.name, scriptname, project)
        self.Import = 'from protocol_ml3d import *'
        self.ParamsStr = ''
        
    def defineSteps(self):        
        restart = False
        if restart:            #Not yet implemented
            pass
#        # Restarting a previous run...
#        else:
#            # Execute protocol in the working directory
#            os.chdir(self.WorkingDir)
#            self.restart_MLrefine3D(RestartIter)
        else:
            
            if self.DoMlf:
                self.ParamsDict['ORoot'] = 'mlf3d'
            else:
                self.ParamsDict['ORoot'] = 'ml3d'
                
            if self.DoMlf and self.DoCorrectAmplitudes:
                ctfFile = self.ParamsDict['ctfFile'] = self.workingDirPath('my.ctfdat')
                # Copy CTFdat to the workingdir as well
                self.insertStep('copyFile', [ctfFile], 
                                   source=self.InCtfDatFile, dest=ctfFile )
            
            self.ParamsDict['InitialVols'] = initVols = self.workingDirPath('initial_volumes.stk')
            self.insertStep('copyVolumes', [initVols], 
                               inputMd=self.RefMd, outputStack=initVols)
        
            if self.DoCorrectGreyScale:
                self.insertCorrectGreyScaleSteps()
                
            if self.DoLowPassFilterReference:
                self.insertFilterStep()
                
            if self.DoML3DClassification:
                self.insertML3DStep()
            
    # Insert the step of launch some program
    def insertRunJob(self, prog, vf=[]): 
        self.insertStep('runJob', verifyfiles=[self.ParamsDict[k] for k in vf],
                        programname=prog, params=self.ParamsStr % self.ParamsDict,
                        NumberOfMpi=self.NumberOfMpi, NumberOfThreads=self.NumberOfThreads)
    # Crude correction of grey-scale, by performing a single iteration of 
    # projection matching and fourier reconstruction
    def insertCorrectGreyScaleSteps(self):
    #    print '*********************************************************************'
    #    print '*  Correcting absolute grey scale of initial reference:'
        cgsDir = self.workingDirPath('CorrectGreyscale')
        self.insertStep('createDir', path=cgsDir)
        volStack = self.ParamsDict['InitialVols'] = join(cgsDir, 'corrected_volumes.stk')
        #FIXME: better to create an empty file
        #self.insertStep('copyVolumes', [volStack], 
        #                       inputMd=self.RefMd, outputStack=volStack)
        md  = MetaData(self.RefMd)
        index = 1
        for idx in md:
            volDir = join(cgsDir, 'vol%03d' % index)
            projs = join(volDir, 'projections')
            self.insertStep('createDir', path=volDir)
            self.ParamsDict.update({
                'inputVol': md.getValue(MDL_IMAGE, idx),
                'outputVol': "%(index)d@%(volStack)s" % locals(),
                'projRefs': projs + ".stk",
                'docRefs': projs + ".doc",
                'corrRefs': join(volDir, 'corrected_refs.stk'),
                'projMatch': join(volDir, "proj_match.doc")
                })
            
            self.ParamsStr = ' -i %(inputVol)s --experimental_images %(ImgMd)s -o %(projRefs)s' + \
                    ' --sampling_rate %(ProjMatchSampling)d --sym %(Symmetry)s' + \
                    'h --compute_neighbors --angular_distance -1' 
                       
            self.insertRunJob('xmipp_angular_project_library', ['projRefs', 'docRefs'])

            self.ParamsStr = '-i %(ImgMd)s -o %(projMatch)s --ref %(projRefs)s' 
            self.insertRunJob('xmipp_angular_projection_matching', ['projMatch'])
 
#FIXME: COMMENTED THIS STEP UNTIL COMPLETION BY ROBERTO    
#            params = '-i %(projMatch)s --lib %(docRefs)s -o %(corrRefs)s'
#            insertRunJob('xmipp_angular_class_average', [corrRefs])

            self.ParamsStr = '-i %(projMatch)s -o %(outputVol)s --sym %(Symmetry)s --weight --thr %(NumberOfThreads)d'
            self.insertRunJob('xmipp_reconstruct_fourier')
        
    def insertFilterStep(self):
        self.ParamsDict['FilteredVols'] = 'filtered_volumes.stk'
        self.ParamsStr = '-i %(InitialVols)s -o %(FilteredVols)s --fourier low_pass %(LowPassFilter) --sampling %(PixelSize)s'
        self.insertRunJob('xmipp_transform_filter', ['FilteredVols'])
        self.ParamsDict['InitialVols'] = self.ParamsDict['FilteredVols']
        
    def insertML3DStep(self):
        self.ParamsStr = "-i %(InitialVols)s --oroot %(ORoot)s --ref %(InitialVols)s --iter %(NumberOfIterations) " + \
                         "--sym %(Symmetry) --ang %(AngularSampling)s %(ExtraParams)s"
        if self.NumberOfThreads > 1:
            self.ParamsStr += " --thr %(NumberOfThreads)d"
        if self.DoNorm:
            self.ParamsStr += " --norm"
        
        if self.DoMlf:
            if self.DoCorrectAmplitudes:
                self.ParamsStr += ' --ctfdat %(ctfFile)d'
            else:
                self.ParamsStr += ' --no_ctf --pixel_size %(PixelSize)f'
            self.ParamsStr += ""
            if not self.ImagesArePhaseFlipped:
                self.ParamsStr += " --not_phase_flipped"
                
''' This function will copy input references into a stack in working directory'''
def copyVolumes(log, inputMd, outputStack):
    from protlib_filesystem import deleteFile
    deleteFile(log, outputStack)
    md = MetaData(inputMd)
    img = Image()
    i = 1
    for idx in md:
        img.read(md.getValue(MDL_IMAGE, idx))
        img.write('%d@%s' % (i, outputStack))
        i += 1
