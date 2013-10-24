#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
#Wrapper for relion 1.2

#   Author: Roberto Marabini
#
from protlib_base import XmippProtocol, protocolMain
from xmipp import MetaData, MDL_SAMPLINGRATE, MDL_IMAGE, MDL_CTF_MODEL, FileName,\
                 MDValueEQ, MDL_REF3D, MD_APPEND
import sys
from os.path import join, exists
from protlib_filesystem import xmippExists, findAcquisitionInfo, moveFile, \
                               replaceBasenameExt, replaceFilenameExt
from protlib_utils import runShowJ, getListFromVector, getListFromRangeString, \
                          runJob, runChimera, which, runChimeraClient
from protlib_import import exportReliontoMetadataFile
from glob import glob

class ProtRelionBase(XmippProtocol):
    relionFiles=['data','optimiser','sampling']
    def __init__(self, protDictName, scriptname, project):
        XmippProtocol.__init__(self, protDictName, scriptname, project)
        self.ParamsStr = ''
        if self.NumberOfMpi > 1:
            self.program = 'relion_refine_mpi'
        else:
            self.program = 'relion_refine'
        self.ParamsDict['program'] = self.program
        
        self.addParam('ORoot', self.WorkingDir + '/')        
        acquisionInfo = self.findAcquisitionInfo(self.ImgMd)

        if not acquisionInfo is None: 
            md = MetaData(acquisionInfo)
            self.addParam('SamplingRate', md.getValue(MDL_SAMPLINGRATE, md.firstObject()))

    def summary(self,totalNumberOfIter=0):
        lines = ["Wrapper implemented for RELION **1.2**"]
        lastIteration=self.lastIter()
        if totalNumberOfIter < 1:
            lines += ['Performed <%d> iterations ' % lastIteration]
        else:
            lines += ['Performed <%d/%d> iterations ' % (lastIteration,totalNumberOfIter)]
            
        return lines 

    def lastIter(self):
        NumberOfIterations=0
        relionDataTemplate = self.getFilename('dataRe', iter=0 )
        for i in range (0,1000):
            fileName = relionDataTemplate.replace('000',("%03d"%i))
            if exists(fileName):
                NumberOfIterations = i
            else:
                break
        return NumberOfIterations

    def validate(self):
        errors = []
        md = MetaData(self.ImgMd)
        if md.containsLabel(MDL_IMAGE):
            # Check that have same size as images:
            from protlib_xmipp import validateInputSize
            if not xmippExists(self.Ref3D):
               errors.append("Reference: <%s> doesn't exists" % ref)
            if len(errors) == 0:
                validateInputSize([self.Ref3D], self.ImgMd, md, errors)
        else:
            errors.append("Input metadata <%s> doesn't contains image label" % self.ImgMd)
            
        if self.DoCTFCorrection and not md.containsLabel(MDL_CTF_MODEL):
            errors.append("CTF correction selected and input metadata <%s> doesn't contains CTF information" % self.ImgMd)
            
        # Check relion is installed
        if len(which('relion_refine')) == 0:
            errors.append('<relion_refine> was not found.') 
        if len(which('relion_movie_handler')) == 0:
            errors.append('''program "relion_movie_handler" is missing. 
                             Are you sure you have relion version 1.2?.''') 
        return errors 
    
    def defineSteps(self): 
        self.doContinue = len(self.ContinueFrom) > 0
        #in the future we need to deal with continue
        restart = False
        #temporary normalized files
        tmpFileNameSTK = join(self.ExtraDir, 'norRelion.stk')
        tmpFileNameXMD = join(self.ExtraDir, 'norRelion.xmd')
        #create extra output directory
        self.insertStep('createDir', verifyfiles=[self.ExtraDir], path=self.ExtraDir)
        #normalize data

        if self.DoNormalizeInputImage:
            self.insertStep('runNormalizeRelion', verifyfiles=[tmpFileNameSTK,tmpFileNameXMD],
                            inputMd  = self.ImgMd,
                            outputMd = tmpFileNameSTK,
                            normType = 'NewXmipp',
                            bgRadius = int(self.MaskDiameterA/2.0/self.SamplingRate),
                            Nproc    = self.NumberOfMpi*self.NumberOfThreads
                            )
        else:
            self.Db.insertStep('createLink',
                               verifyfiles=[tmpFileNameXMD],
                               source=self.ImgMd,
                               dest=tmpFileNameXMD)
        # convert input metadata to relion model
        self.ImgStar = self.extraPath(replaceBasenameExt(tmpFileNameXMD, '.star'))
        self.insertStep('convertImagesMd', verifyfiles=[self.ImgStar],
                        inputMd=tmpFileNameXMD, 
                        outputRelion=self.ImgStar                            
                        )
    def defineSteps2(self,lastIteration,NumberOfClasses,ExtraInputs=None,ExtraOutputs=None):
        verifyFiles=[]
        inputs=[]
        for i in range (0,lastIteration+1):
            for v in self.relionFiles:
                 verifyFiles += [self.getFilename(v+'Xm', iter=i )]
#            for v in self.relionFiles:
#                 inputs += [self.getFilename(v+'Re', iter=i )]
        if not ExtraOutputs is None:
#            inputs      += ExtraInputs
            verifyFiles += ExtraOutputs
        #data to store in Working dir so it can be easily accessed by users
        lastIterationVolumeFns = []
        standardOutputClassFns = []
        for ref3d in range (1,NumberOfClasses+1):
            try:
                lastIterationVolumeFns += [self.getFilename('volume', iter=lastIteration, ref3d=ref3d )]
            except:
                lastIterationVolumeFns += [self.getFilename('volumeFinal', iter=lastIteration, ref3d=ref3d )]
            
            standardOutputClassFns += ["images_ref3d%06d@"%ref3d + self.workingDirPath("classes_ref3D.xmd")]
        lastIterationMetadata = "images@"+self.getFilename('data'+'Xm', iter=lastIteration )
        
        extra = self.workingDirPath('extra')
        self.insertStep('convertRelionMetadata', verifyfiles=verifyFiles,
                        inputs  = inputs,
                        outputs = verifyFiles,
                        lastIterationVolumeFns = lastIterationVolumeFns,
                        lastIterationMetadata  = lastIterationMetadata,
                        NumberOfClasses        = NumberOfClasses,
                        standardOutputClassFns = standardOutputClassFns,
                        standardOutputImageFn  = "images@" + self.workingDirPath("images.xmd"),
                        standardOutputVolumeFn = "volumes@" + self.workingDirPath("volumes.xmd"),
                        extraDir=extra
                        )
#Do NOT change relion binary files, let us use MRC
#                 for ref3d in range(1,NumberOfClasses+1):
#                       inputs  += [self.getFilename('volumeMRC', iter=it, ref3d=ref3d )]
#                       outputs += [self.getFilename('volume', iter=it, ref3d=ref3d )]
#           
#        self.insertStep('convertRelionBinaryData',
#                        inputs  = inputs,
#                        outputs = outputs
#                        )

    def createFilenameTemplates(self):
        extra = self.workingDirPath('extra')
        self.extraIter = join(extra, 'relion_it%(iter)03d_')
        myDict={}
        return myDict

    def visualize(self):
        
        plots = [k for k in ['TableImagesPerClass'
                           , 'DisplayResolutionPlotsFSC'
                           , 'DisplayResolutionPlotsSSNR'
                           , 'TableChange'
                           , 'Likelihood'
                           , 'DisplayReconstruction'
                           , 'DisplayAngularDistribution'
                           ] if self.ParamsDict[k]]

        if len(plots):
            self.launchRelionPlots(plots)

    def visualizeVar(self, varName):
        self.launchRelionPlots([varName])
        
    def launchRelionPlots(self, selectedPlots):
        ############TABLES STOP ANALIZE
        ''' Launch some plots for a Projection Matching protocol run '''
        import numpy as np
        _log = self.Log

        xplotter=None
        self._plot_count = 0
        
        def doPlot(plotName):
            return plotName in selectedPlots
        DisplayRef3DNo = self.parser.getTkValue('DisplayRef3DNo')
        VisualizeIter = self.parser.getTkValue('VisualizeIter')
        
        if DisplayRef3DNo == 'all':
            ref3Ds = range(1,self.NumberOfClasses+1)
        else:
            ref3Ds = map(int, getListFromVector(self.parser.getTkValue('SelectedRef3DNo')))

        #Get last iteration
        self.NumberOfIterations = self.lastIter()

        if VisualizeIter == 'last':
            iterations = [self.NumberOfIterations]
        elif VisualizeIter == 'all':
            iterations = range(1,self.NumberOfIterations+1)
        else:
            iterations = map(int, getListFromVector(self.parser.getTkValue('SelectedIters')))
        #check if last iteration exists

        runShowJExtraParameters = ' --dont_wrap --view '+ self.parser.getTkValue('DisplayVolumeSlicesAlong') 
        self.DisplayVolumeSlicesAlong=self.parser.getTkValue('DisplayVolumeSlicesAlong')
        lastIteration = iterations[-1]
        lastRef3D = ref3Ds[-1]
        lastVolume = self.getFilename('volumeMRC', iter=lastIteration, ref3d=lastRef3D )
        if not xmippExists(lastVolume):
            print "does not exists"
            message = "No data available for <iteration %d> and <class %d>"%\
                       (int(lastIteration),int(lastRef3D))
            showError("File not Found", message, self.master)
            return
        
        #if relion has not finished metadata has not been converted to xmipp
        #do it now if needed
        #=============
        outputs=[]#output files
        inputs=[]
        for i in iterations:
            for v in self.relionFiles:
                fileName =self.getFilename(v+'Xm', iter=i )
                if not xmippExists(fileName):
                    outputs += [fileName]
                    inputs += [self.getFilename(v+'Re', iter=i )]
        lastIterationVolumeFns = []
        standardOutputClassFns = []
        for ref3d in range (1,self.NumberOfClasses+1):
            lastIterationVolumeFns += [self.getFilename('volume', iter=lastIteration, ref3d=ref3d )]
            standardOutputClassFns += ["images_ref3d%06d@"%ref3d + self.workingDirPath("classes_ref3D.xmd")]
        lastIterationMetadata = "images@"+self.getFilename('data'+'Xm', iter=lastIteration )

import os
def convertRelionMetadata(log, inputs,
                          outputs,
                          lastIterationVolumeFns,
                          lastIterationMetadata,
                          NumberOfClasses,
                          standardOutputClassFns,
                          standardOutputImageFn,
                          standardOutputVolumeFn,
                          extraDir
                          ):
    """ Convert the relion style MetaData to one ready for xmipp.
    Main differences are: STAR labels are named different and
    optimiser.star -> changes in orientation, offset. number images assigned to each class
    data.star -> loglikelihood (per image) may be used to delete worst images (10%)
                 orientation.shift per particle
    model.star:average_P_max (plot per iteration)
               block: model_classes
                  class distribution, number of particles per class
                  estimated error in orientation and translation
               block: model_class_N: 
                    resolution-dependent SNR: report resol where it drops below 1? (not that important?)
                                    (only in auto-refine) Gold-std FSC: make plot!
    """
    for i in glob(join(extraDir,"relion*star")):
        o = replaceFilenameExt(i,'.xmd')
        exportReliontoMetadataFile(i,o)
    #create images. xmd and class metadata
    #lastIteration = self.NumberOfIterations
    #this images cames from relion
    md = MetaData(lastIterationMetadata)
    #total number Image
    numberImages = md.size()


    #total number volumes 
    comment  = " numberImages=%d..................................................... "%numberImages
    comment += " numberRef3D=%d........................................................."%NumberOfClasses
    md.setComment(comment)
    md.write(standardOutputImageFn)
    #data_images_ref3d000001
    mdOut = MetaData()
    mdOut.setComment(comment)
    f = FileName(standardOutputClassFns[0])
    f=f.removeBlockName()
    if exists(f):
        os.remove(f)
    for i in range (0,NumberOfClasses):
        mdOut.clear()
        mdOut.importObjects(md, MDValueEQ(MDL_REF3D, i+1))
        mdOut.write(standardOutputClassFns[i],MD_APPEND)
        
    #volume.xmd, metada with volumes
    mdOut.clear()
    for lastIterationVolumeFn in lastIterationVolumeFns:
        objId = mdOut.addObject()
        mdOut.setValue(MDL_IMAGE, lastIterationVolumeFn, objId)
    mdOut.write(standardOutputVolumeFn)
        
def convertRelionBinaryData(log, inputs,outputs):
    """Make sure mrc files are volumes properlly defined"""
    program = "xmipp_image_convert"
    for i,o in zip(inputs,outputs):
        args = "-i %s -o %s  --type vol"%(i,o)
        runJob(log, program, args )
        
def renameOutput(log, WorkingDir, ProgId):
    ''' Remove ml2d prefix from:
        ml2dclasses.stk, ml2dclasses.xmd and ml2dimages.xmd'''
    prefix = '%s2d' % ProgId
    for f in ['%sclasses.stk', '%sclasses.xmd', '%simages.xmd']:
        f = join(WorkingDir, f % prefix)
        nf = f.replace(prefix, '')
        moveFile(log, f, nf)
##        inputs=[]
##        outputs=[]
##        for it in iterations:
##             for ref3d in range(1,self.NumberOfClasses+1):
##                 fileName = self.getFilename('volume', iter=it, ref3d=ref3d )
##                 if not xmippExists(fileName):
##                     inputs  += [self.getFilename('volumeMRC', iter=it, ref3d=ref3d )]
##                     outputs += [fileName]
##        convertRelionBinaryData(None, inputs, outputs)
#                       
#        #============
#        if doPlot('TableImagesPerClass'):
#            _r = []
#            mdOut = MetaData()
#            maxRef3D = 1
#            _c = tuple(['Iter/Ref'] + ['Ref3D_%d' % ref3d for ref3d in ref3Ds])
#            for it in iterations:
#                md = MetaData("images@"+self.getFilename('data'+'Xm', iter=it ))
#                #agregar por ref3D
#                mdOut.aggregate(md, AGGR_COUNT, MDL_REF3D, MDL_REF3D, MDL_COUNT)
#                tmp = ['0']*(self.NumberOfClasses+1)
#                tmp[0]=('Iter_%d'%it)
#                for objId in mdOut:
#                    _ref3D = mdOut.getValue(MDL_REF3D,objId)
#                    _count = mdOut.getValue(MDL_COUNT,objId)
#                    if _ref3D not in ref3Ds:
#                        tmp [_ref3D] = '---'
#                    else:
#                        tmp [_ref3D] = str(_count)
#                _r.append(tuple(tmp))
#            showTable(_c,_r, title='Images per Ref3D', width=100)
#
#        if doPlot('DisplayResolutionPlotsFSC'):
#            self.ResolutionThreshold=float(self.parser.getTkValue('ResolutionThreshold'))
#
#            if(len(ref3Ds) == 1):
#                gridsize1 = [1, 1]
#            elif (len(ref3Ds) == 2):
#                gridsize1 = [2, 1]
#            else:
#                gridsize1 = [(len(ref3Ds)+1)/2, 2]
#            xplotter = XmippPlotter(*gridsize1,windowTitle="ResolutionFSC")
#            #print 'gridsize1: ', gridsize1
#            #print 'iterations: ', iterations
#            #print 'ref3Ds: ', ref3Ds
#            
#            for ref3d in ref3Ds:
#                plot_title = 'Ref3D_%s' % ref3d
#                a = xplotter.createSubPlot(plot_title, 'Armstrongs^-1', 'Fourier Shell Correlation', yformat=False)
#                legendName = []
#                blockName = 'model_class_%d@'%ref3d
#                for it in iterations:
#                    file_name = blockName + self.getFilename('model'+'Xm', iter=it )
#                    #file_name = self.getFilename('ResolutionXmdFile', iter=it, ref=ref3d)
#                    if xmippExists(file_name):
#                        #print 'it: ',it, ' | file_name:',file_name
#                        md = MetaData(file_name)
#                        resolution_inv = [md.getValue(MDL_RESOLUTION_FREQ, id) for id in md]
#                        frc = [md.getValue(MDL_RESOLUTION_FRC, id) for id in md]
#                        a.plot(resolution_inv, frc)
#                        legendName.append('Iter_'+str(it))
#                    xplotter.showLegend(legendName)
#                if (self.ResolutionThreshold < max(frc)):
#                    a.plot([min(resolution_inv), max(resolution_inv)],[self.ResolutionThreshold, self.ResolutionThreshold], color='black', linestyle='--')
#                    a.grid(True)
#            xplotter.draw()
#    
#        if doPlot('DisplayResolutionPlotsSSNR'):
#            if(len(ref3Ds) == 1):
#                gridsize1 = [1, 1]
#            elif (len(ref3Ds) == 2):
#                gridsize1 = [2, 1]
#            else:
#                gridsize1 = [(len(ref3Ds)+1)/2, 2]
#            xplotter = XmippPlotter(*gridsize1,windowTitle="ResolutionSSNR")
#            #print 'gridsize1: ', gridsize1
#            #print 'iterations: ', iterations
#            #print 'ref3Ds: ', ref3Ds
#            md = MetaData()
#            activateMathExtensions()
#            for ref3d in ref3Ds:
#                plot_title = 'Ref3D_%s' % ref3d
#                a = xplotter.createSubPlot(plot_title, 'Armstrongs^-1', 'log(SSNR)', yformat=False)
#                legendName=[]
#                blockName = 'model_class_%d@'%ref3d
#                for it in iterations:
#                    file_name = blockName + self.getFilename('model'+'Xm', iter=it )
#                    #file_name = self.getFilename('ResolutionXmdFile', iter=it, ref=ref3d)
#                    #print file_name
#                    if xmippExists(file_name):
#                        #print 'it: ',it, ' | file_name:',file_name
#                        mdOut = MetaData(file_name)
#                        md.clear()
#                        # only cross by 1 is important
#                        md.importObjects(mdOut, MDValueGT(MDL_RESOLUTION_SSNR, 0.9))
#                        md.operate("resolutionSSNR=log(resolutionSSNR)")
#                        resolution_inv = [md.getValue(MDL_RESOLUTION_FREQ, id) for id in md]
#                        frc = [md.getValue(MDL_RESOLUTION_SSNR, id) for id in md]
#                        a.plot(resolution_inv, frc)
#                        legendName.append('Iter_'+str(it))
#                    xplotter.showLegend(legendName)
#            xplotter.draw()
#
#        if doPlot('TableChange'):
#            _r = []
#            mdOut = MetaData()
#            _c = ('Iter','offset','orientation','classAssigment')
#            
#            for it in iterations:
#                md = MetaData(self.getFilename('optimiser'+'Xm', iter=it ))
#                #agregar por ref3D
#                firstObject = md.firstObject()
#                cOrientations = md.getValue(MDL_AVG_CHANGES_ORIENTATIONS,firstObject)
#                cOffsets = md.getValue(MDL_AVG_CHANGES_OFFSETS,firstObject)
#                cClasses = md.getValue(MDL_AVG_CHANGES_CLASSES,firstObject)
#                _r.append(("Iter_%d" % it, cOrientations,cOffsets,cClasses))
#            showTable(_c,_r, title='Changes per Iteration', width=100)
#
#        if doPlot('Likelihood'):
#            for it in iterations:
#                fileName = 'images@'+ self.getFilename('data'+'Xm', iter=it )
#                md = MetaData(fileName)
#                md.sort(MDL_LL, False)
#                md.write(fileName)
#                runShowJ(fileName)
#                xplotter = XmippPlotter(windowTitle="max Likelihood particles sorting Iter_%d"%it)
#                xplotter.createSubPlot("Particle sorting: Iter_%d"%it, "Particle number", "maxLL")
#                xplotter.plotMd(md, False, mdLabelY=MDL_LL)
#
#
##        if doPlot('AvgPMAX'):
##            plotMd=MetaData()
##            for it in iterations:
##                fileName = 'model_general@'+ self.getFilename('model'+'Xm', iter=it )
##                md = MetaData(fileName)
##                pmax = md.getValue(MDL_AVGPMAX,md.firstObject())
##                objId=plotMd.addObject()
##                plotMd.setValue(MDL_COUNT,long(it),objId)
##                plotMd.setValue(MDL_AVGPMAX,pmax,objId)
##            print plotMd
##            xplotter = XmippPlotter(windowTitle="Avg PMax")
##            xplotter.createSubPlot("Avg PMax /per iteration", 'Iteration', 'AvgPmax')
##            xplotter.plotMd(plotMd, mdLabelX=MDL_COUNT, mdLabelY=MDL_AVGPMAX)
#
#        if doPlot('AvgPMAX'):
#            _r = []
#            mdOut = MetaData()
#            _c = ('Iteration','Avg PMax')
#            
#            for it in iterations:
#                fileName = 'model_general@'+ self.getFilename('model'+'Xm', iter=it )
#                md = MetaData(fileName)
#                pmax = md.getValue(MDL_AVGPMAX,md.firstObject())
#                _r.append(("Iter_%d" % it, pmax))
#            showTable(_c,_r, title='Avg PMax /per iteration', width=100)
#
#
#        if doPlot('DisplayReconstruction'):
#
#            for ref3d in ref3Ds:
#                for it in iterations:
#                    #relion save volume either in spider or mrc therefore try first one of them and then the other
#                    file_name = self.getFilenameAlternative('volume','volumeMRC', iter=it, ref3d=ref3d )
#                    if xmippExists(file_name):
#                        #Chimera
#                        if(self.DisplayVolumeSlicesAlong == 'surface'):
#                            runChimeraClient(file_name)
#                        else:
#                        #Xmipp_showj (x,y and z shows the same)
#                            try:
#                                runShowJ(file_name, extraParams = runShowJExtraParameters)
#                            except Exception, e:
#                                showError("Error launching java app", str(e), self.master)
#                    else:
#                        print "file ", file_name, "those not exits"
#        from tempfile import NamedTemporaryFile    
#        if doPlot('DisplayAngularDistribution'):
#            self.DisplayAngularDistributionWith = self.parser.getTkValue('DisplayAngularDistributionWith')
#            if(self.DisplayAngularDistributionWith == '3D'):
#                    #relion save volume either in spider or mrc therefore try first one of them and then the other
#                fileNameVol = self.getFilenameAlternative('volume', 'volumeMRC', iter=iterations[-1], ref3d=1 )
#                (Xdim, Ydim, Zdim, Ndim) = getImageSize(fileNameVol)
#                for ref3d in ref3Ds:
#                    #relion save volume either in spider or mrc therefore try first one of them and then the other
#                    fileNameVol = self.getFilenameAlternative('volume', 'volumeMRC',iter=iterations[-1], ref3d=ref3d )
#                    md=MetaData('images@'+ self.getFilename('data'+'Xm', iter=iterations[-1] ))
#                    mdOut = MetaData()
#                    mdOut.importObjects(md, MDValueEQ(MDL_REF3D, ref3d))
#                    md.clear()
#                    md.aggregateMdGroupBy(mdOut, AGGR_COUNT, [MDL_ANGLE_ROT, MDL_ANGLE_TILT], MDL_ANGLE_ROT, MDL_WEIGHT)
#                    md.setValue(MDL_ANGLE_PSI,0.0, md.firstObject())
#                    _OuterRadius = int(float(self.parser.getTkValue('MaskDiameterA'))/2.0/self.SamplingRate)
#                    #do not delete file since chimera needs it
#                    ntf = NamedTemporaryFile(dir="/tmp/", suffix='.xmd',delete=False)
#                    md.write("angularDist@"+ntf.name)
#                    parameters = ' --mode projector 256 -a ' + "angularDist@"+ ntf.name + " red "+\
#                         str(float(_OuterRadius) * 1.1)
#                    runChimeraClient(fileNameVol,parameters)
#            else: #DisplayAngularDistributionWith == '2D'
#                if(len(ref3Ds) == 1):
#                    gridsize1 = [1, 1]
#                elif (len(ref3Ds) == 2):
#                    gridsize1 = [2, 1]
#                else:
#                    gridsize1 = [(len(ref3Ds)+1)/2, 2]
#                
#                xplotter = XmippPlotter(*gridsize1, mainTitle='Iteration_%d' % iterations[-1], windowTitle="AngularDistribution")
#                
#                for ref3d in ref3Ds:
#                    #relion save volume either in spider or mrc therefore try first one of them and then the other
#                    fileNameVol = self.getFilenameAlternative('volume', 'volumeMRC', iter=iterations[-1], ref3d=ref3d )
#                    md=MetaData('images@'+ self.getFilename('data'+'Xm', iter=iterations[-1] ))
#                    mdOut = MetaData()
#                    mdOut.importObjects(md, MDValueEQ(MDL_REF3D, ref3d))
#                    md.clear()
#                    md.aggregateMdGroupBy(mdOut, AGGR_COUNT, [MDL_ANGLE_ROT, MDL_ANGLE_TILT], MDL_ANGLE_ROT, MDL_WEIGHT)
#                    #keep only direction projections from the right ref3D
#                    plot_title = 'Ref3D_%d' % ref3d
#                    xplotter.plotAngularDistribution(plot_title, md)
#                xplotter.draw()
#
#        if xplotter:
#            xplotter.show()

def runNormalizeRelion(log,inputMd,outputMd,normType,bgRadius,Nproc,):
    stack = inputMd
    if exists(outputMd):
        remove(outputMd)
    program = "xmipp_transform_normalize"
    args = "-i %(stack)s -o %(outputMd)s "

    if bgRadius <= 0:
        particleSize = MetaDataInfo(stack)[0]
        bgRadius = int(particleSize/2)

    if normType=="OldXmipp":
        args += "--method OldXmipp"
    elif normType=="NewXmipp":
        args += "--method NewXmipp --background circle %(bgRadius)d"
    else:
        args += "--method Ramp --background circle %(bgRadius)d"
    runJob(log, program, args % locals(), Nproc)

def convertImagesMd(log, inputMd, outputRelion):
    """ Convert the Xmipp style MetaData to one ready for Relion.
    Main differences are: STAR labels are named different and
    in Xmipp the images rows contains a path to the CTFModel, 
    while Relion expect the explict values in the row.
    Params:
     input: input filename with the Xmipp images metadata
     output: output filename for the Relion metadata
    """
    from protlib_import import exportMdToRelion
    
    md = MetaData(inputMd)
    # Get the values (defocus, magnification, etc) from each 
    # ctfparam files and put values in the row
    md.fillExpand(MDL_CTF_MODEL)
    # Create the mapping between relion labels and xmipp labels
    exportMdToRelion(md, outputRelion)