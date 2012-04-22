#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Protocol for Xmipp-based 2D alignment and classification,
# using maximum-likelihood principles
#
#   Author:  Sjors Scheres, January 2008
#  Updated:  J. M. de la Rosa Trevin July 2011
#

from os.path import join, exists
from protlib_base import XmippProtocol, protocolMain
from config_protocols import protDict
from xmipp import MetaData, Image, MDL_IMAGE, MDL_ITER, MDL_LL

class ProtML3D(XmippProtocol):
    def __init__(self, scriptname, project):
        XmippProtocol.__init__(self, protDict.ml3d.name, scriptname, project)
        self.Import = 'from protocol_ml3d import *'
        self.ParamsStr = ''
        if self.DoMlf:
            self.progId = 'mlf'
        else:
            self.progId = 'ml'
        self.ParamsDict['ProgId'] = self.progId
        self.ORoot = self.ParamsDict['ORoot'] = self.workingDirPath('%s3d' % self.progId)
                        
    def createFilenameTemplates(self):
        return {
                'iter_logs': '%(ORoot)s_%(ProgId)s2d_iter_logs.xmd',
                'iter_refs': '%(ORoot)s_%(ProgId)s2d_iter_refs.xmd',
                'vols': 'iter%(iter)06d@%(ORoot)s_vols.xmd'
                }
        
    def summary(self):
        md = MetaData(self.ImgMd)
        lines = ["Input images:  [%s] (<%u>)" % (self.ImgMd, md.size())]
        if self.DoMlf:
            if self.DoCorrectAmplitudes:
                suffix = "with CTF correction "
            else:
                suffix = "ignoring CTF effects "
            lines.append("Using a ML in <Fourier-space> " + suffix)

        lines.append("Reference volumes(s): [%s]" % self.RefMd)
        
        if self.NumberOfReferences > 1:
            lines.append("Number of references per volume: <%d>" % self.NumberOfReferences)
        
        logs = self.getFilename('iter_logs')    
        
        if exists(logs):
            md = MetaData(logs)
            id = md.lastObject()
            iteration = md.getValue(MDL_ITER, id)
            lines.append("Last iteration:  <%d>" % iteration)
            LL = md.getValue(MDL_LL, id)
            lines.append("LogLikelihood:  %f" % LL)
            mdRefs = self.getFilename('iter_refs')
            lines.append("Last 2D classes: [iter%06d@%s]" % (iteration, mdRefs))
            fnVols = self.getFilename('vols', iter=iteration)
            lines.append("Last 3D classes: [%s]" % fnVols)
        
        return lines
    
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
            initVols = self.ParamsDict['InitialVols'] = self.workingDirPath('initial_volumes.stk')
            self.mdVols = MetaData(self.RefMd)
            
            self.insertStep('copyVolumes', [initVols], 
                               inputMd=self.RefMd, outputStack=initVols)
            
            if self.DoCorrectGreyScale:
                self.insertCorrectGreyScaleSteps()
                
            if self.DoLowPassFilterReference:
                self.insertFilterStep()
                
            if self.NumberOfReferences > 1:
                self.insertGenerateRefSteps()
                
            if self.DoML3DClassification:
                self.insertML3DStep(self.ImgMd, self.ORoot, self.ParamsDict['InitialVols'], self.NumberOfIterations)
            
    # Insert the step of launch some program
    def insertRunJob(self, prog, vf=[], files=None, useProcs=True):
        if files is None:
            files = [self.ParamsDict[k] for k in vf]
        NumberOfMpi = NumberOfThreads = 1
        if useProcs:
            NumberOfMpi = self.NumberOfMpi
            NumberOfThreads = self.NumberOfThreads
        self.insertStep('runJob', programname=prog, 
                        params = self.ParamsStr % self.ParamsDict, 
                        verifyfiles = files, 
                        NumberOfMpi = NumberOfMpi, 
                        NumberOfThreads = NumberOfThreads)
    
    # Crude correction of grey-scale, by performing a single iteration of 
    # projection matching and fourier reconstruction
    def insertCorrectGreyScaleSteps(self):
        ''' Correct the initial reference greyscale '''
        cgsDir = self.workingDirPath('CorrectGreyscale')
        self.insertStep('createDir', path=cgsDir)
        volStack = self.ParamsDict['InitialVols'] = join(cgsDir, 'corrected_volumes.stk')
        index = 1
        outputVol = ''
        for idx in self.mdVols:
            volDir = join(cgsDir, 'vol%03d' % index)
            projs = join(volDir, 'projections')
            self.insertStep('createDir', path=volDir)
            outputVol = "%(index)d@%(volStack)s" % locals()
            self.ParamsDict.update({
                'inputVol': self.mdVols.getValue(MDL_IMAGE, idx),
                'outputVol': outputVol,
                'projRefs': projs + ".stk",
                'docRefs': projs + ".doc",
                'corrRefs': join(volDir, 'corrected_refs.stk'),
                'projMatch': join(volDir, "proj_match.doc")
                })
            self.mdVols.setValue(MDL_IMAGE, outputVol, idx)
            self.ParamsStr = ' -i %(inputVol)s --experimental_images %(ImgMd)s -o %(projRefs)s' + \
                    ' --sampling_rate %(ProjMatchSampling)d --sym %(Symmetry)s' + \
                    'h --compute_neighbors --angular_distance -1' 
                       
            self.insertRunJob('xmipp_angular_project_library', ['projRefs', 'docRefs'])

            self.ParamsStr = '-i %(ImgMd)s -o %(projMatch)s --ref %(projRefs)s' 
            self.insertRunJob('xmipp_angular_projection_matching', ['projMatch'])
 
#FIXME: COMMENTED THIS STEP UNTIL COMPLETION BY ROBERTO    
#            self.ParamsStr = '-i %(projMatch)s --lib %(docRefs)s -o %(corrRefs)s'
#            self.insertRunJob('xmipp_angular_class_average', ['corrRefs'])

            self.ParamsStr = '-i %(projMatch)s -o %(outputVol)s --sym %(Symmetry)s --weight --thr %(NumberOfThreads)d'
            self.insertRunJob('xmipp_reconstruct_fourier', ['outputVol'])
            index += 1
        
    def insertFilterStep(self):
        volStack = self.ParamsDict['FilteredVols'] = self.workingDirPath('filtered_volumes.stk')
        index = 1
        outputVol = ''
        for idx in self.mdVols:
            outputVol = "%(index)d@%(volStack)s" % locals()
            self.mdVols.setValue(MDL_IMAGE, outputVol, idx)
            index += 1
        self.ParamsStr = '-i %(InitialVols)s -o %(FilteredVols)s --fourier low_pass %(LowPassFilter) --sampling %(PixelSize)s'
        self.insertRunJob('xmipp_transform_filter', ['FilteredVols'], useProcs=False)
        self.ParamsDict['InitialVols'] = self.ParamsDict['FilteredVols']
            
    def insertGenerateRefSteps(self):
        ''' Generete more reference volumes than provided in input reference '''
        grDir = self.workingDirPath('GeneratedReferences')
        # Create dir for seeds generation
        self.insertStep('createDir', path=grDir)
        # Split images metadata
        nvols = self.ParamsDict['NumberOfVols'] = self.mdVols.size() * self.NumberOfReferences
        sroot = self.ParamsDict['SplitRoot'] = join(grDir, 'images')
        self.ParamsStr = '-i %(ImgMd)s -n %(NumberOfVols)d --oroot %(SplitRoot)s'
        files = ['%s%06d.xmd' % (sroot, i) for i in range(1, nvols+1)]        
        self.insertRunJob('xmipp_metadata_split', files=files, useProcs=False)
        
        
        volStack = self.ParamsDict['InitialVols'] = self.workingDirPath('generated_volumes.stk')
        index = 1
        copyVols = []
        for idx in self.mdVols:
            for i in range(self.NumberOfReferences):
                outputVol = "%d@%s" % (index, volStack)
                generatedVol = join(grDir, "vol%03d_iter%06d_vol%06d.vol" % (index, 1, 1))
                copyVols.append((outputVol, generatedVol))
                self.insertML3DStep(files[index-1], join(grDir, 'vol%03d' % index), self.mdVols.getValue(MDL_IMAGE, idx), 1)
                #self.mdVols.setValue(MDL_IMAGE, outputVol, idx)
                index += 1
                
        for outVol, genVol in copyVols:
            self.ParamsDict.update({'outVol': outVol, 'genVol':genVol})
            self.ParamsStr = '-i %(genVol)s -o %(outVol)s'
            self.insertRunJob('xmipp_image_convert', files=[volStack], useProcs=False)
        
        
    def insertML3DStep(self, inputImg, oRoot, initialVols, numberOfIters):
        self.ParamsDict.update({
                         '_ImgMd': inputImg,
                         '_ORoot': oRoot,
                         '_InitialVols': initialVols,
                         '_NumberOfIterations': numberOfIters       
                                })
        self.ParamsStr = "-i %(_ImgMd)s --oroot %(_ORoot)s --ref %(_InitialVols)s --iter %(_NumberOfIterations)d " + \
                         "--sym %(Symmetry)s --ang %(AngularSampling)s %(ExtraParams)s"
#        if self.NumberOfReferences > 1:
#            self.ParamsStr += " --nref %(NumberOfReferences)s"
        if self.NumberOfThreads > 1:
            self.ParamsStr += " --thr %(NumberOfThreads)d"
        if self.DoNorm:
            self.ParamsStr += " --norm"
        
        if self.DoMlf:
            if not self.DoCorrectAmplitudes:
                self.ParamsStr += ' --no_ctf %(PixelSize)f'
            if not self.ImagesArePhaseFlipped:
                self.ParamsStr += " --not_phase_flipped"

        self.ParamsStr += " --recons %(ReconstructionMethod)s "
        if self.ReconstructionMethod == 'wslART':
            self.ParamsStr += " %(ARTExtraParams)s"
        else:
            self.ParamsStr += " %(FourierExtraParams)s" 
        self.insertRunJob('xmipp_%s_refine3d' % self.progId, [])
                
''' This function will copy input references into a stack in working directory'''
def copyVolumes(log, inputMd, outputStack):
    from protlib_filesystem import deleteFile
    deleteFile(log, outputStack)
    md = MetaData(inputMd)
    img = Image()
    for i, idx in enumerate(md):
        img.read(md.getValue(MDL_IMAGE, idx))
        img.write('%d@%s' % (i + 1, outputStack))
