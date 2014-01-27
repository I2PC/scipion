#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Module with classes and functions needed for Relion protocols
#
#  Author: Roberto Marabini            roberto@cnb.csic.es    
# Updated: J. M. de la Rosa Trevin     jmdelarosa@cnb.csic.es (Jan 2014)
#
#------------------------------------------------------------------------------------------------
#
import sys
import re
from os.path import join, exists
from glob import glob

from xmipp import *
from protlib_base import XmippProtocol, protocolMain
from protlib_filesystem import xmippExists, findAcquisitionInfo, moveFile, \
                               replaceBasenameExt, replaceFilenameExt
from protlib_utils import runShowJ, getListFromRangeString, \
                          runJob, runChimera, which, runChimeraClient
from protlib_import import exportReliontoMetadataFile, addRelionLabels
from protlib_gui_ext import showWarning, showTable, showError
from protlib_gui_figure import XmippPlotter
from protlib_xmipp import getSampling, validateInputSize


class ProtRelionBase(XmippProtocol):
    def __init__(self, protDictName, scriptname, project):
        self.relionFiles=['data','optimiser','sampling'] 
        XmippProtocol.__init__(self, protDictName, scriptname, project)
        self.ParamsStr = ''
        if self.NumberOfMpi > 1:
            self.program = 'relion_refine_mpi'
        else:
            self.program = 'relion_refine'
        self.ParamsDict['program'] = self.program
        if self.DoContinue:
            self.setPreviousRunFromFile(self.optimiserFileName)
            #if optimizer has not been properly selected this will 
            #fail, let us go ahead and handle the situation in verify
            try:
                self.inputProperty('ImgMd')
            except:
                print "Cannot find the input file with experimental images in protocol %s" % self.PrevRunName

        self.addParam('ORoot', self.WorkingDir + '/')
        self.addParam('SamplingRate', getSampling(self.ImgMd))

    def summary(self):
        lines = ["Wrapper implemented for RELION **1.2**"]
        return lines 

    def papers(self):
        papers = ['Bayesian view: Scheres, JMB (2012) [http://www.ncbi.nlm.nih.gov/pubmed/22100448]',
                  'RELION implementation: Scheres, JSB (2012) [http://www.ncbi.nlm.nih.gov/pubmed/23000701]'
                  ]
        return papers

    def lastIter(self):
        fileNameTemplate = self.getFilename('dataRe', iter=0)
        fileNameTemplate = fileNameTemplate.replace('000','???')
        a = sorted(glob(fileNameTemplate))
        if not a:
            return 0 
        return int (re.findall(r'_it\d+_',a[-1])[0][3:6])

    def firstIter(self):
        fileNameTemplate = self.getFilename('dataRe', iter=0)
        fileNameTemplate = fileNameTemplate.replace('000','???')
        a = sorted(glob(fileNameTemplate))[0]
        if not a:
            return 1 
        i = int (re.findall(r'_it\d+_',a)[0][3:6])
        if i == 0:
             i = 1
        return i

    def _containsCTF(self, md):
        """ Validate where the metadata contains information of the CTF. """
        def containsCTFLabels():
            for l in [MDL_CTF_CS, MDL_CTF_DEFOCUS_ANGLE, MDL_CTF_DEFOCUSU, MDL_CTF_DEFOCUSV, MDL_CTF_Q0, MDL_CTF_VOLTAGE]:
                if not md.containsLabel(l):
                    return False
            return True                
        return (md.containsLabel(MDL_CTF_MODEL) or containsCTFLabels())
        
    def _validateInputSize(self, inputParam, inputLabel, mandatory, md, errors):
        """ Validate, if the mask value is non-empty,
        that the mask file exists and have the same 
        dimensions than the input images. 
        """
        inputValue = getattr(self, inputParam, '')
        errMsg = "Param <%s> file should exists" % inputLabel
        if not mandatory:
            errMsg += " if not empty"
            
        if mandatory or inputValue.strip():
            if not xmippExists(inputValue):
                errors.append(errMsg)
            else:
                validateInputSize([inputValue], self.ImgMd, md, errors, inputLabel)
            
    def validate(self):
        errors = []
        md = MetaData()
        print "reading md"
        md.read(self.ImgMd, 1) # Read only the first object in md
        print "done"
        
        if md.containsLabel(MDL_IMAGE):
            self._validateInputSize('Ref3D', 'Reference volume', True, md, errors)
            self._validateInputSize('ReferenceMask', 'Reference mask',False, md, errors)
            self._validateInputSize('SolventMask','Solvent mask', False, md, errors)
        else:
            errors.append("Input metadata <%s> doesn't contains image label" % self.ImgMd)
            
        if self.DoCTFCorrection and not self._containsCTF(md):
            errors.append("CTF correction selected and input metadata <%s> doesn't contains CTF information" % self.ImgMd)
            
        # Check relion is installed
        if len(which('relion_refine')) == 0:
            errors.append('<relion_refine> was not found.') 
        if len(which('relion_movie_handler')) == 0:
            errors.append('''program "relion_movie_handler" is missing. 
                             Are you sure you have relion version 1.2?.''') 
        return errors 
    
    def defineSteps(self): 
        #remove ContinueFrom
        #self.doContinue = len(self.ContinueFrom) > 0
        #in the future we need to deal with continue
        restart = False
        #temporary normalized files
        tmpFileNameSTK = join(self.ExtraDir, 'norRelion.stk')
        tmpFileNameXMD = join(self.ExtraDir, 'norRelion.xmd')
        #create extra output directory
        self.insertStep('createDir', verifyfiles=[self.ExtraDir], path=self.ExtraDir)
        self.insertStep("linkAcquisitionInfo", InputFile=self.ImgMd, dirDest=self.WorkingDir) 
        #normalize data

        if self.DoNormalizeInputImage:
            self.insertStep('runNormalizeRelion', verifyfiles=[tmpFileNameSTK,tmpFileNameXMD],
                            inputMd  = self.ImgMd,
                            outputMd = tmpFileNameSTK,
                            normType = 'NewXmipp',
                            bgRadius = int(self.MaskRadiusA/self.SamplingRate),
                            Nproc    = self.NumberOfMpi*self.NumberOfThreads
                            )
        else:
            self.insertStep('createLink',
                               verifyfiles=[tmpFileNameXMD],
                               source=self.ImgMd,
                               dest=tmpFileNameXMD)
        # convert input metadata to relion model
        self.ImgStar = self.extraPath(replaceBasenameExt(tmpFileNameXMD, '.star'))
        self.insertStep('convertImagesMd', verifyfiles=[self.ImgStar],
                        inputMd=tmpFileNameXMD, 
                        outputRelion=self.ImgStar                            
                        )
        
    def defineSteps2(self, firstIteration
                         , lastIteration
                         , NumberOfClasses
                         , ExtraInputs=None
                         , ExtraOutputs=None):
        verifyFiles=[]
        inputs=[]
        for i in range (firstIteration, lastIteration+1):
            for v in self.relionFiles:
                 verifyFiles += [self.getFilename(v+'Xm', iter=i, workingDir=self.WorkingDir )]
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
                lastIterationVolumeFns += [self.getFilename('volume', iter=lastIteration, ref3d=ref3d, workingDir=self.WorkingDir )]
            except:
                lastIterationVolumeFns += [self.getFilename('volumeFinal', iter=lastIteration, ref3d=ref3d, workingDir=self.WorkingDir )]
            
            standardOutputClassFns += ["images_ref3d%06d@"%ref3d + self.workingDirPath("classes_ref3D.xmd")]
        lastIterationMetadata = "images@"+self.getFilename('data'+'Xm', iter=lastIteration, workingDir=self.WorkingDir )
        
        extra = self.workingDirPath('extra')
        inputMetadatas = []
        for it in range (firstIteration,self.NumberOfIterations+1): #always list all iterations
            inputMetadatas += ["images@"+self.getFilename('data'+'Xm', iter=it, workingDir=self.WorkingDir )]

        self.insertStep('convertRelionMetadata', verifyfiles=verifyFiles,
                        inputs  = inputs,
                        outputs = verifyFiles,
                        lastIterationVolumeFns = lastIterationVolumeFns,
                        lastIterationMetadata  = lastIterationMetadata,
                        NumberOfClasses        = NumberOfClasses,
                        standardOutputClassFns = standardOutputClassFns,
                        standardOutputImageFn  = "images@" + self.workingDirPath("images.xmd"),
                        standardOutputVolumeFn = "volumes@" + self.workingDirPath("volumes.xmd"),
                        extraDir=extra,
                        inputMetadatas=inputMetadatas
                        #,imagesAssignedToClassFn=self.getFilename('imagesAssignedToClass', workingDir=self.WorkingDir)
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
        self.extraIter = self.extraPath('relion_it%(iter)03d_')
        myDict = {
                  'iter_data': self.extraIter + 'data.star'
                  }

        return myDict

    def visualize(self):
        
        plots = [k for k in ['TableImagesPerClass'
                           , 'DisplayResolutionPlotsFSC'
                           , 'DisplayResolutionPlotsSSNR'
                           , 'TableChange'
                           , 'Likelihood'
                           , 'DisplayReconstruction'
                           , 'DisplayAngularDistribution'
                           , 'DisplayImagesClassification'
                           , 'AvgPMAX'
                           ] if self.ParamsDict[k]]

        if len(plots):
            self.launchRelionPlots(plots)

    def visualizeVar(self, varName):
        print "visualizeVar, varName: ", varName
        self.launchRelionPlots([varName])
        
    def _loadVisualizeItersAndRefs(self):
        """ Load selected iterations and classes 3D for visualization mode. """
        DisplayRef3DNo = self.parser.getTkValue('DisplayRef3DNo')
        VisualizeIter = self.parser.getTkValue('VisualizeIter')
        
        if DisplayRef3DNo == 'all':
            self._visualizeRef3Ds = range(1, self.NumberOfClasses + 1)
        else:
            self._visualizeRef3Ds = getListFromRangeString(self.parser.getTkValue('SelectedRef3DNo'))
        lastIteration = self.lastIter()
        firstIter = self.firstIter()
        
        if VisualizeIter == 'last':
            self._visualizeIterations = [lastIteration]
        else:
            self._visualizeIterations = getListFromRangeString(self.parser.getTkValue('SelectedIters'))

        self._visualizeLastIteration = self._visualizeIterations[-1]
        self._visualizeLastRef3D = self._visualizeRef3Ds[-1]
        
        self._visualizeVolumesMode = self.parser.getTkValue('DisplayReconstruction')
        
        
    def launchRelionPlots(self, selectedPlots):
        ############TABLES STOP ANALIZE
        ''' Launch some plots for a Projection Matching protocol run '''
        from tempfile import NamedTemporaryFile    
        TmpDir=join(self.WorkingDir,"tmp")
        import numpy as np
        _log = self.Log

        xplotter=None
        self._plot_count = 0
        
        print "launchRelionPlots, before loading..."
        self._loadVisualizeItersAndRefs()
        
        def doPlot(plotName):
            return plotName in selectedPlots
        #check if last iteration exists
        
        if( self.relionType == 'classify' ):
            lastVolume = self.getFilename('volume', iter=self._visualizeLastIteration, ref3d=self._visualizeLastRef3D, workingDir=self.WorkingDir )
            _volume='volume'
        else:
            lastVolume = self.getFilename('volumeh1', iter=self._visualizeLastIteration, ref3d=self._visualizeLastRef3D, workingDir=self.WorkingDir )
            _volume='volumeFinal'
            fileNameRe =self.getFilename('modelXmFinalRe')
            fileNameXm =self.getFilename('modelXmFinalXm')
            if xmippExists(fileNameRe) and (not xmippExists(fileNameXm)):
                #since relion crasesh very often after finishing its job
                #I cannot assume that conversion is done
                exportReliontoMetadataFile(fileNameRe,fileNameXm)
                
        if not xmippExists(lastVolume):
            message = "No data available for <iteration %d> and <class %d>"%\
                       (int(self._visualizeLastIteration),int(lastRef3D))
            showError("File not Found", message, self.master)
            return
        
#        #if relion has not finished metadata has not been converted to xmipp
#        #do it now if needed
#        for i in iterations:
#            for v in self.relionFiles:
#                fileName =self.getFilename(v+'Xm', iter=i)
#                if not xmippExists(fileName):
#                    exportReliontoMetadataFile(self.getFilename(v+'Re', iter=i ),fileName)
        if doPlot('TableImagesPerClass'):
            mdOut = MetaData()
            #mdOut = self.imagesAssignedToClass2(firstIter,self._visualizeLastIteration)
            mdOut = self.imagesAssignedToClass2(iterations[0],iterations[-1])
            print mdOut
            _r = []
            maxRef3D = 1
            _c = tuple(['Iter/Ref'] + ['Ref3D_%d' % ref3d for ref3d in ref3Ds])

            oldIter = iterations[0]
            tmp = ['---']*(len(ref3Ds)+1)
            for objId in mdOut:
                _ref3D = mdOut.getValue(MDL_REF3D,objId)
                _count = mdOut.getValue(MDL_CLASS_PERCENTAGE,objId)
                _iter  = mdOut.getValue(MDL_ITER,objId)
                if oldIter != _iter:
                    tmp[0]=('Iter_%d'%(_iter-1))
                    _r.append(tuple(tmp))
                    oldIter = _iter
                    tmp = ['---']*(len(ref3Ds)+1)
                tmp [_ref3D] = str(_count)
            tmp[0]=('Iter_%d'%oldIter)
            _r.append(tuple(tmp))
            showTable(_c,_r, title='Images per Class', width=100)
            
        if doPlot('AvgPMAX'):
            _r = []
            mdOut = MetaData()
            
            if(self.relionType=='classify'):
                _c = ('Iteration','Avg PMax')
                
                for it in iterations: # range (firstIter,self._visualizeLastIteration+1): #alwaya list all iteration
                    fileName = 'model_general@'+ self.getFilename('model'+'Xm', iter=it , workingDir=self.WorkingDir)
                    md = MetaData(fileName)
                    pmax = md.getValue(MDL_AVGPMAX,md.firstObject())
                    _r.append(("Iter_%d" % it, pmax))
            else:
                _c = ('Iteration','Avg PMax h1', 'Avg PMax h2')
                
                for it in iterations:#range (firstIter,self._visualizeLastIteration+1): #alwaya list all iteration
                    fileName = 'model_general@'+ self.getFilename('half1_model'+'Xm', iter=it , workingDir=self.WorkingDir)
                    md = MetaData(fileName)
                    pmax1 = md.getValue(MDL_AVGPMAX,md.firstObject())
                    fileName = 'model_general@'+ self.getFilename('half2_model'+'Xm', iter=it , workingDir=self.WorkingDir)
                    md = MetaData(fileName)
                    pmax2 = md.getValue(MDL_AVGPMAX,md.firstObject())
                    _r.append(("Iter_%d" % it, pmax1,pmax2))
                fileName = 'model_general@'+ self.getFilename('modelXmFinalXm', workingDir=self.WorkingDir)
                if(self.relionType=='refine' and iterations[-1]==self.NumberOfIterations):
                    md = MetaData(fileName)
                    pmax1 = md.getValue(MDL_AVGPMAX,md.firstObject())
                    _r.append(("all_particles", pmax1,'---'))
                
            showTable(_c,_r, title='Avg PMax /per iteration', width=100)
            
        if doPlot('DisplayReconstruction'):
            volumeKeys = self._getVolumeKeys()
                
            for ref3d in self._visualizeRef3Ds:
                for it in self._visualizeIterations:
                  for volKey in volumeKeys:
                    self.display3D(self.getFilename(volKey, iter=it, ref3d=ref3d , workingDir=self.WorkingDir))
            
            if self.relionType=='refine' and iterations[-1]==self.NumberOfIterations:
                self.display3D(self.getFilename('volumeFinal',ref3d=1, workingDir=self.WorkingDir))
                
        if doPlot('DisplayImagesClassification'):
            self._visualizeDisplayImagesClassification()
                
#            for ref3d in ref3Ds:
#                for iter in iterations:
#                    md=MetaData('images@'+ self.getFilename('dataXm', iter=iter, workingDir=self.WorkingDir ))
#                    mdOut = MetaData()
#                    mdOut.importObjects(md, MDValueEQ(MDL_REF3D, ref3d))
#                    ntf = NamedTemporaryFile(dir=TmpDir, suffix='_it%03d_class%03d.xmd'%(iter,ref3d),delete=False)
#                    mdOut.write(ntf.name)
#                    self.display2D(ntf.name)
#                    md.clear()
                #save temporary file
        if doPlot('DisplayResolutionPlotsSSNR'):
            names = []
            windowTitle = {}
            if(self.relionType=='classify'):
                names = ['modelXm']
                windowTitle[names[0]]="ResolutionSSNR"
            else:
                names = ['half1_modelXm','half2_modelXm']
                windowTitle[names[0]]="ResolutionSSNR_half1"
                windowTitle[names[1]]="ResolutionSSNR_half2"
                
            if(len(ref3Ds) == 1):
                gridsize1 = [1, 1]
            elif (len(ref3Ds) == 2):
                gridsize1 = [2, 1]
            else:
                gridsize1 = [(len(ref3Ds)+1)/2, 2]
            md = MetaData()
            activateMathExtensions()
            for name in names:
              xplotter = XmippPlotter(*gridsize1,windowTitle=windowTitle[name])
              for ref3d in ref3Ds:
                plot_title = 'Class_%s' % ref3d
                a = xplotter.createSubPlot(plot_title, 'Armstrongs^-1', 'log(SSNR)', yformat=False)
                legendName=[]
                blockName = 'model_class_%d@'%ref3d
                for it in iterations:
                    file_name = blockName + self.getFilename(name, iter=it, workingDir=self.WorkingDir )
                    #file_name = self.getFilename('ResolutionXmdFile', iter=it, ref=ref3d)
                    if xmippExists(file_name):
                        mdOut = MetaData(file_name)
                        md.clear()
                        # only cross by 1 is important
                        md.importObjects(mdOut, MDValueGT(MDL_RESOLUTION_SSNR, 0.9))
                        md.operate("resolutionSSNR=log(resolutionSSNR)")
                        resolution_inv = [md.getValue(MDL_RESOLUTION_FREQ, id) for id in md]
                        frc = [md.getValue(MDL_RESOLUTION_SSNR, id) for id in md]
                        a.plot(resolution_inv, frc)
                        legendName.append('Iter_'+str(it))
                    xplotter.showLegend(legendName)
              a.grid(True)
              xplotter.draw()
              
            if(self.relionType=='refine'):
                file_name = blockName + self.getFilename('modelXmFinalXm')
                if xmippExists(file_name):
                    gridsize1 = [1, 1]
                    xplotter = XmippPlotter(*gridsize1,windowTitle="ResolutionSSNRAll")
                    plot_title = 'SSNR for all images, final iteration'
                    a = xplotter.createSubPlot(plot_title, 'Armstrongs^-1', 'log(SSNR)', yformat=False)
                    blockName = 'model_class_%d@'%1
                    mdOut = MetaData(file_name)
                    md.clear()
                    # only cross by 1 is important
                    md.importObjects(mdOut, MDValueGT(MDL_RESOLUTION_SSNR, 0.9))
                    md.operate("resolutionSSNR=log(resolutionSSNR)")
                    resolution_inv = [md.getValue(MDL_RESOLUTION_FREQ, id) for id in md]
                    frc = [md.getValue(MDL_RESOLUTION_SSNR, id) for id in md]
                    a.plot(resolution_inv, frc)
                    a.grid(True)
                    xplotter.draw()
                    
        if doPlot('DisplayResolutionPlotsFSC'):
            self.ResolutionThreshold=float(self.parser.getTkValue('ResolutionThreshold'))
            names = []
            windowTitle = {}
            #not used for classification but I may change my mind in the future
            if(self.relionType=='classify'):
                names = ['modelXm']
                windowTitle[names[0]]="ResolutionFSC"
            else:
                names = ['half1_modelXm','half2_modelXm']
                windowTitle[names[0]]="ResolutionFSC_half1"
                windowTitle[names[1]]="ResolutionFSC_half2"
            if(len(ref3Ds) == 1):
                gridsize1 = [1, 1]
            elif (len(ref3Ds) == 2):
                gridsize1 = [2, 1]
            else:
                gridsize1 = [(len(ref3Ds)+1)/2, 2]
            md = MetaData()
            activateMathExtensions()
            for name in names:
              xplotter = XmippPlotter(*gridsize1,windowTitle=windowTitle[name])
              for ref3d in ref3Ds:
                plot_title = 'Ref3D_%s' % ref3d
                a = xplotter.createSubPlot(plot_title, 'Armstrongs^-1', 'FSC', yformat=False)
                legendName=[]
                blockName = 'model_class_%d@'%ref3d
                for it in iterations:
                    file_name = blockName + self.getFilename(name, iter=it, workingDir=self.WorkingDir )
                    #file_name = self.getFilename('ResolutionXmdFile', iter=it, ref=ref3d)
                    if xmippExists(file_name):
                        mdOut = MetaData(file_name)
                        resolution_inv = [mdOut.getValue(MDL_RESOLUTION_FREQ, id) for id in mdOut]
                        frc = [mdOut.getValue(MDL_RESOLUTION_FRC, id) for id in mdOut]
                        a.plot(resolution_inv, frc)
                        legendName.append('Iter_'+str(it))
                    xplotter.showLegend(legendName)
                if (self.ResolutionThreshold < max(frc)):
                    a.plot([min(resolution_inv), max(resolution_inv)],[self.ResolutionThreshold, self.ResolutionThreshold], color='black', linestyle='--')
                    a.grid(True)
              xplotter.draw()
              
            if(self.relionType=='refine'):
                file_name = blockName + self.getFilename('modelXmFinalXm', workingDir=self.WorkingDir)
                print "final model", file_name
                if xmippExists(file_name):
                    print "final model exists"
                    gridsize1 = [1, 1]
                    xplotter = XmippPlotter(*gridsize1,windowTitle="ResolutionSSNRAll")
                    plot_title = 'FSC for all images, final iteration'
                    a = xplotter.createSubPlot(plot_title, 'Armstrongs^-1', 'log(SSNR)', yformat=False)
                    blockName = 'model_class_%d@'%1
                    mdOut = MetaData(file_name)
                    # only cross by 1 is important
                    resolution_inv = [mdOut.getValue(MDL_RESOLUTION_FREQ, id) for id in mdOut]
                    frc = [mdOut.getValue(MDL_RESOLUTION_FRC, id) for id in mdOut]
                    a.plot(resolution_inv, frc)
                    if (self.ResolutionThreshold < max(frc)):
                        a.plot([min(resolution_inv), max(resolution_inv)],[self.ResolutionThreshold, self.ResolutionThreshold], color='black', linestyle='--')
                        a.grid(True)
                    xplotter.draw()

        if doPlot('DisplayAngularDistribution'):
            self._visualizeDisplayAngularDistribution()

        if doPlot('TableChange'):
            _r = []
            mdOut = MetaData()
            if(self.relionType=='classify'):
                _c = ('Iter','offset','orientation','classAssigment')
            else:
                _c = ('Iter','offset','orientation')
            
            print " Computing average changes in: offset, abgles, and class belonging"
            for it in iterations:
                print "Computing data for iteration; %03d"% (it)
                fileName = self.getFilename('optimiser'+'Xm', iter=it )
                md = MetaData(fileName)
                #agregar por ref3D
                firstObject = md.firstObject()
                cOrientations = md.getValue(MDL_AVG_CHANGES_ORIENTATIONS,firstObject)
                cOffsets = md.getValue(MDL_AVG_CHANGES_OFFSETS,firstObject)
                cClasses = md.getValue(MDL_AVG_CHANGES_CLASSES,firstObject)
                if(self.relionType=='classify'):
                    _r.append(("Iter_%d" % it, cOrientations,cOffsets,cClasses))
                else:
                    _r.append(("Iter_%d" % it, cOrientations,cOffsets))
            showTable(_c,_r, title='Changes per Iteration', width=100)


        if doPlot('Likelihood'):
            for it in iterations:
                print "Computing data for iteration; %03d"% (it)
                fileName = 'images@'+ self.getFilename('data'+'Xm', iter=it, workingDir=self.WorkingDir )
                md = MetaData(fileName)
                md.sort(MDL_LL, False)
                md.write(fileName)
                runShowJ(fileName)
                xplotter = XmippPlotter(windowTitle="max Likelihood particles sorting Iter_%d"%it)
                xplotter.createSubPlot("Particle sorting: Iter_%d"%it, "Particle number", "maxLL")
                xplotter.plotMd(md, False, mdLabelY=MDL_LL)

        if xplotter:
            xplotter.show()
                
    def display3D(self, filename):
        if xmippExists(filename):
            #Chimera
            if self._visualizeVolumesMode == 'chimera':
                runChimeraClient(filename)
            else:
                try:
                    runShowJ(filename, extraParams=' --dont_wrap ')
                except Exception, e:
                    showError("Error launching xmipp_showj: ", str(e), self.master)
        else:
            showError("Can not visualize volume file <%s>\nit does not exists." % filename, self.master)

    def display2D(self,filename):
        if xmippExists(filename):
            try:
                os.system('xmipp_showj %s %s' % (filename, '--dont_wrap'))
            except Exception, e:
                showError("Error launching java app", str(e), self.master)
        else:
            showError("Can not visualize file <%s>\nit does not exists." % filename, self.master)

    def imagesAssignedToClass2(self,firstIter, lastIteration):
            inputMetadatas = []
            for it in range (firstIter,lastIteration+1): #always list all iterations
                inputMetadatas += ["model_classes@"+self.getFilename('model'+'Xm', iter=it, workingDir=self.WorkingDir )]
            return imagesAssignedToClass(inputMetadatas)
        
    def _visualizeDisplayAngularDistribution(self):        
        self.SpheresMaxradius = float(self.parser.getTkValue('SpheresMaxradius'))
        self.DisplayAngularDistributionWith = self.parser.getTkValue('DisplayAngularDistributionWith')
        
        volumeKeys = self._getVolumeKeys()
        
        if self.DisplayAngularDistribution == 'chimera':
            fileNameVol = self.getFilename(volumeKeys[0], iter=iterations[-1], ref3d=1, workingDir=self.WorkingDir )
            (Xdim, Ydim, Zdim, Ndim) = getImageSize(fileNameVol)
            for it in iterations:
                for ref3d in ref3Ds:
                    print "Computing data for iteration; %03d and class %03d"% (it, ref3d)
                    #relion save volume either in spider or mrc therefore try first one of them and then the other
                    #name[0] or name[1] as bakcground is an arbitrary selection
                    fileNameVol = self.getFilename(volumeKeys[0], iter=it, ref3d=ref3d, workingDir=self.WorkingDir )
                    fileNameMd  = 'images@'+ self.getFilename('data'+'Xm', iter=it, workingDir=self.WorkingDir)
                    print "fileNameVol", fileNameVol
                    print "fileNameMd", fileNameMd
                    md=MetaData(fileNameMd)
                    #md=MetaData('images@'+ self.getFilename('data'+'Xm', iter=iterations[-1] ))
                    mdOut = MetaData()
                    mdOut.importObjects(md, MDValueEQ(MDL_REF3D, ref3d))
                    md.clear()
                    md.aggregateMdGroupBy(mdOut, AGGR_COUNT, [MDL_ANGLE_ROT, MDL_ANGLE_TILT], MDL_ANGLE_ROT, MDL_WEIGHT)
                    md.setValue(MDL_ANGLE_PSI,0.0, md.firstObject())
                    _OuterRadius = int(float(self.MaskRadiusA)/self.SamplingRate)
                    #do not delete file since chimera needs it
                    ntf = NamedTemporaryFile(dir=TmpDir, suffix='_it%03d_class%03d.xmd'%(it,ref3d),delete=False)
                    md.write("angularDist@"+ntf.name)
                    parameters = ' --mode projector 256 -a ' + "angularDist@"+ ntf.name + " red "+\
                         str(float(_OuterRadius) * 1.1)
                    if self.SpheresMaxradius > 0:
                        parameters += ' %f ' % self.SpheresMaxradius
                    runChimeraClient(fileNameVol,parameters)
        else: #DisplayAngularDistributionWith == '2D'
            if(len(ref3Ds) == 1):
                gridsize1 = [1, 1]
            elif (len(ref3Ds) == 2):
                gridsize1 = [2, 1]
            else:
                gridsize1 = [(len(ref3Ds)+1)/2, 2]
            
            
            for it in iterations:
              xplotter = XmippPlotter(*gridsize1, mainTitle='Iteration_%d' % iterations[-1], windowTitle="AngularDistribution")
              for ref3d in ref3Ds:
                print "Computing data for iteration; %03d and class %03d"% (it, ref3d)
                #relion save volume either in spider or mrc therefore try first one of them and then the other
                fileNameMd  = 'images@'+ self.getFilename('data'+'Xm', iter=it, workingDir=self.WorkingDir)
                md=MetaData(fileNameMd)
                mdOut = MetaData()
                mdOut.importObjects(md, MDValueEQ(MDL_REF3D, ref3d))
                md.clear()
                md.aggregateMdGroupBy(mdOut, AGGR_COUNT, [MDL_ANGLE_ROT, MDL_ANGLE_TILT], MDL_ANGLE_ROT, MDL_WEIGHT)
                #keep only direction projections from the right ref3D
                plot_title = 'Class_%d' % ref3d
                xplotter.plotAngularDistribution(plot_title, md)
            xplotter.draw()
 
    def _visualizeDisplayImagesClassification(self):
        """ Read Relion _data.star images file and 
        generate a new metadata with the Xmipp classification standard:
        a 'classes' block and a 'class00000?_images' block per class.
        If the new metadata was already written, it is just shown.
        """
        for it in self._visualizeIterations:
            data_star = self.getFilename('iter_data', iter=it)
            data_classes = data_star.replace('.star', '_classes.star')
            if not xmippExists(data_classes):
                addRelionLabels(replace=True)
                md = MetaData(data_star)
                mdClasses = MetaData()
                mdClasses.aggregate(md, AGGR_COUNT, MDL_REF3D, MDL_IMAGE, MDL_CLASS_COUNT)
                #md2.aggregate(md, AGGR_COUNT, MDL_REF3D, MDL_IMAGE, MDL_COUNT)
                mdClasses.write('classes@' + data_classes)
                mdImages = MetaData()
                
                for objId in mdClasses:
                    ref3d = mdClasses.getValue(MDL_REF3D, objId)
                    mdImages.importObjects(md, MDValueEQ(MDL_REF3D, ref3d))
                    print "writing to '%s', %d images" % (('class%06d_images' % ref3d) + data_classes, mdImages.size())
                    mdImages.write(('class%06d_images@' % ref3d) + data_classes, MD_APPEND)
                
            self.display2D(data_classes)

def imagesAssignedToClass(inputMetadatas):
#        mdOut = MetaData()
#        mdAux = MetaData()
#        maxRef3D = 1
#        it=1
#        for mdName in inputMetadatas:
#            print "processing metadata", mdName
#            md = MetaData(mdName)
#            mdAux.aggregate(md, AGGR_COUNT, MDL_REF3D, MDL_REF3D, MDL_COUNT)
#            mdAux.fillConstant(MDL_ITER,it)
#            mdOut.unionAll(mdAux)
#            it += 1
#        return mdOut

        it=1
        mdOut = MetaData()
        # read percentages and store them in a single metadata table
        # with columns MDL_ITER MDL_REF3D and MDL_CLASS_PERCENTAGE
        for mdName in inputMetadatas:
            md = MetaData(mdName)
            mdOut.unionAll(md)
            md.clear()
        # change MDL_IMAGE_REF by 3DREF
        for id in mdOut:
            fileName=mdOut.getValue(MDL_IMAGE_REF,id)
            ref3d = getClass(fileName)
            mdOut.setValue(MDL_REF3D,ref3d,id)
            iter = getIteration(fileName)
            mdOut.setValue(MDL_ITER,iter,id)
        return mdOut

import os
def convertRelionMetadata(log, inputs,
                          outputs,
                          lastIterationVolumeFns,
                          lastIterationMetadata,
                          NumberOfClasses,
                          standardOutputClassFns,
                          standardOutputImageFn,
                          standardOutputVolumeFn,
                          extraDir,
                          inputMetadatas
                          #,imagesAssignedToClassFn
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
    comment = " numberImages=%d..................................................... " % numberImages
    comment += " numberRef3D=%d........................................................." % NumberOfClasses
    md.setComment(comment)
    md.write(standardOutputImageFn)
    #data_images_ref3d000001
    mdOut = MetaData()
    mdOut.setComment(comment)
    f = FileName(standardOutputClassFns[0])
    f=f.removeBlockName()
    if exists(f):
        os.remove(f)
    for i in range (1,NumberOfClasses):
        mdOut.clear()
        mdOut.importObjects(md, MDValueEQ(MDL_REF3D, i+1))
        mdOut.write(standardOutputClassFns[i],MD_APPEND)
        
    #volume.xmd, metada with volumes
    mdOut.clear()
    for lastIterationVolumeFn in lastIterationVolumeFns:
        objId = mdOut.addObject()
        mdOut.setValue(MDL_IMAGE, lastIterationVolumeFn, objId)
    mdOut.write(standardOutputVolumeFn)
#    if NumberOfClasses >1:
#        mdOut.clear()
#        mdOut = imagesAssignedToClass(inputMetadatas)
#        mdOut.write(imagesAssignedToClassFn)

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

def runNormalizeRelion(log,inputMd,outputMd,normType,bgRadius,Nproc,):
    stack = inputMd
    if exists(outputMd):
        os.remove(outputMd)
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
    if md.containsLabel(MDL_CTF_MODEL):
        md.fillExpand(MDL_CTF_MODEL)
    # Create the mapping between relion labels and xmipp labels
    exportMdToRelion(md, outputRelion)

def getIteration(fileName):
   return int (re.findall(r'_it\d+_',fileName)[0][3:6])

def getClass(fileName):
   return int (re.findall(r'_class\d+',fileName)[0][6:9])
