#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Protocol for Xmipp-based 2D alignment and classification,
# using maximum-likelihood principles
#
#   Author:  Sjors Scheres, January 2008
#  Updated:  J. M. de la Rosa Trevin July 2011
#

import os
from os.path import join
from protlib_base import XmippProtocol, protocolMain
from config_protocols import protDict
from protlib_filesystem import copyFile, createDir, deleteDir, deleteFile
from xmipp import MetaData, Image, MDL_IMAGE
from protocol__mltomo import WorkingDir
from protlib_utils import runJob

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
                
            if (self.DoMlf and self.DoCorrectAmplitudes):
                ctfFile = self.workingDirPath('my.ctfdat')
                # Copy CTFdat to the workingdir as well
                self.insertStep('copyFile', [ctfFile], 
                                   source=self.InCtfDatFile, dest=ctfFile )
            
            self.ParamsDict['InitialVols'] = initVols = self.workingDirPath('initial_volumes.stk')
            self.insertStep('copyInputReferences', [initVols],  workingDir=self.WorkingDir,
                               referenceMd=self.RefMd, outRefs=initVols)
        
            if self.DoCorrectGreyScale:
                self.insertCorrectGreyScaleSteps()
#                self.insertStep('correctGreyScale', [],
#                                    workingDir=self.WorkingDir, imgMd=self.ImgMd, initRefs=self.RefMd, 
#                                    samplingRate=self.ProjMatchSampling, sym=self.Symmetry, 
#                                    numberOfMpi=self.NumberOfMpi, numberOfThreads=self.NumberOfThreads )

    # Insert the step of launch some program
    def insertRunJob(self, prog, vf=[]): 
        self.insertStep('runJob', verifyfiles=[self.ParamsDict[k] for k in vf],
                        programname=prog, params=self.ParamsStr % self.ParamsDict(),
                        NumberOfMpi=self.NumberOfMpi, NumberOfThreads=self.NumberOfThreads)
    # Crude correction of grey-scale, by performing a single iteration of 
    # projection matching and fourier reconstruction
    def insertCorrectGreyScaleSteps(self):
    #    print '*********************************************************************'
    #    print '*  Correcting absolute grey scale of initial reference:'
        cgsDir = self.workingDirPath('CorrectGreyscale')
        self.insertStep('createDir', path=cgsDir)
        volStack = self.ParamsDict['InitialVols'] = self.workingDirPath(cgsDir, 'corrected_volumes.stk')
        md  = MetaData(self.RefMd)
        index = 1
        for idx in md:
            volDir = join(cgsDir, 'vol%03d' % index)
            projs = join(volDir, 'projections')
            self.insertStep('createDir', path=volDir)
            self.ParamsDict.update({
                'InputVol': md.getValue(MDL_IMAGE, idx),
                'OutputVol': "%(index)d@%(volStack)s" % locals(),
                'projRefs': projs + ".stk",
                'docRefs': projs + ".doc",
                'corrRefs': join(volDir, 'corrected_refs.stk'),
                'projMatch': join(volDir, "proj_match.doc")
                })
            
            self.ParamsStr = ' -i %(InputVol)s --experimental_images %(ImgMd)s -o %(projRefs)s' + \
                    ' --sampling_rate %(ProjMatchSampling)d --sym %(Symmetry)sh' + \
                    ' --compute_neighbors --angular_distance -1' 
                       
            self.insertRunJob('xmipp_angular_project_library', ['projRefs', 'docRefs'])

            self.params = '-i %(imgMd)s -o %(projMatch)s --ref %(projRefs)s' 
            self.insertRunJob('xmipp_angular_projection_matching', ['projMatch'])
 
#FIXME: COMMENTED THIS STEP UNTIL COMPLETION BY ROBERTO    
#            params = '-i %(projMatch)s --lib %(docRefs)s -o %(corrRefs)s'
#            insertRunJob('xmipp_angular_class_average', [corrRefs])

            self.params = '-i %(projMatch)s -o %(outputVol)s --sym %(Symmetry)s --weight --thr %(NumberOfThreads)'
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
            self.ParamsStr += ""
    
def initRefsFileName(workingDir):
    return join(workingDir, 'initial_refs.stk')

''' This function will copy input references into a stack in working directory'''
def copyInputReferences(log, workingDir, referenceMd, outRefs):
    deleteFile(log, outRefs)
    md = MetaData(referenceMd)
    img = Image()
    i = 1
    for idx in md:
        img.read(md.getValue(MDL_IMAGE, idx))
        img.write('%d@%s' % (i, outRefs))
        i += 1

#
## Perform the actual classification
#def execute_ML3D_classification(self):
#    import os
#    import launch_job
#
#    dirname='RunML3D/'
#    if not os.path.exists(dirname):
#        os.makedirs(dirname)
#
#    if (self.DoJustRefine):
#        if os.path.exists('filtered_reference.vol'):
#            reference='filtered_reference.vol'
#        elif os.path.exists('corrected_reference.vol'):
#            reference='corrected_reference.vol'
#        else:
#            reference='initial_reference.vol'
#    else:
#        reference='ml3d_seeds.sel'
#
#    self.execute_MLrefine3D(self.InSelFile,
#                            dirname+'ml3d',
#                            reference,
#                            self.AngularSampling,
#                            self.NumberOfIterations,
#                            self.Symmetry,
#                            self.ImagesArePhaseFlipped,
#                            self.SeedsAreAmplitudeCorrected,
#                            self.ExtraParamsMLrefine3D)
#
#
## Either for seeds generation or for ML3D-classification
#def execute_MLrefine3D(self,inselfile,outname,
#                       volname,sampling,iter,symmetry,phase_flipped,amplitude_corrected,extraparam):
#    import os
#    import launch_job
#
#    print '*********************************************************************'
#    print '*  Executing ml_refine3d program :' 
#    params= ' -i '    + str(inselfile) + \
#            ' -o '    + str(outname) + \
#            ' -vol '  + str(volname) + \
#            ' -iter ' + str(iter) + \
#            ' -sym '  + symmetry +\
#            ' -ang '  + str(sampling)
#    params+=' '+extraparam
#    if (self.NumberOfThreads > 1):
#        params+=' -thr ' + str(self.NumberOfThreads)
#    if (self.DoNorm):
#        params+=' -norm '
#    if (self.DoFourier):
#        params+=' -fourier '
#    if (self.DoMlf):
#        if (self.DoCorrectAmplitudes):
#            params+= ' -ctfdat my.ctfdat'
#        else:
#            params+= ' -no_ctf -pixel_size ' + str(self.PixelSize)
#        if (not phase_flipped):
#            params+= ' -not_phase_flipped'
#        if (not amplitude_corrected):
#            params+= ' -ctf_affected_refs'
#        if (self.HighResLimit > 0):
#            params += ' -high ' + str(self.HighResLimit)
#
#    if (self.DoMlf):
#        program="xmipp_mlf_refine3d"
#    else:
#        program="xmipp_ml_refine3d"
#       
#    launchJob(program,
#                          params,
#                          self.log,
#                          self.DoParallel,
#                          self.NumberOfMpi,
#                          self.NumberOfThreads,
#                          self.SystemFlavour)
#
#
## Either for seeds generation or for ML3D-classification
#def restart_MLrefine3D(self,iter):
#    import os
#    import launch_job, utils_xmipp
#
#    print '*********************************************************************'
#    print '*  Restarting ml(f)_refine3d program :' 
#    restartname = utils_xmipp.composeFileName('RunML3D/ml3d_it',iter,'log')
#    params= ' -restart ' + restartname
#
#    if (self.DoMlf):
#        program="xmipp_mlf_refine3d"
#    else:
#        program="xmipp_ml_refine3d"
#
#    launchJob(program,
#                          params,
#                          self.log,
#                          self.DoParallel,
#                          self.NumberOfMpi,
#                          self.NumberOfThreads,
#                          self.SystemFlavour)
#
#def close(self):
#    message='Done!'
#    print '*',message
#    print '*********************************************************************'
#
