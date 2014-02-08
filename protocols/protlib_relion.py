#!/usr/bin/env xmipp_python
#------------------------------------------------------------------------------------------------
# Module with classes and functions needed for Relion protocols
#
#  Author: Roberto Marabini            roberto@cnb.csic.es    
# Updated: J. M. de la Rosa Trevin     jmdelarosa@cnb.csic.es (Jan 2014)
#
#------------------------------------------------------------------------------------------------
#
import os
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
                  'iter_data': self.extraIter + 'data.star',
                  'iter_classes': self.extraIter + 'classes.xmd',
                  'iter_angularDist': self.extraIter + 'angularDist.xmd',
                  'avgPmax': self.tmpPath('iterations_avgPmax.xmd'),
                  'changes': self.tmpPath('iterations_changes.xmd')
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
        
        if doPlot('AvgPMAX'):
            self._visualizeAvgPMAX()
                
        if doPlot('DisplayReconstruction'):
            self._visualizeDisplayReconstruction()
                
        if doPlot('DisplayImagesClassification'):
            self._visualizeDisplayImagesClassification()

        if doPlot('DisplayResolutionPlotsSSNR'):
            self._visualizeDisplayResolutionPlotsSSNR()
                    
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
            self._visualizeTableChange()

        if doPlot('Likelihood'):
            self._visualizeLikelihood()
#            
#        if xplotter:
#            xplotter.show()
                
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

    def display2D(self,fn, extraParams=''):
        if xmippExists(fn):
            try:
                runShowJ(fn, extraParams=extraParams + ' --dont_wrap')
            except Exception, e:
                showError("Error launching java app", str(e), self.master)
        else:
            showError("Can not visualize file <%s>\nit does not exists." % fn, self.master)

    def imagesAssignedToClass2(self,firstIter, lastIteration):
            inputMetadatas = []
            for it in range (firstIter,lastIteration+1): #always list all iterations
                inputMetadatas += ["model_classes@"+self.getFilename('model'+'Xm', iter=it, workingDir=self.WorkingDir )]
            return imagesAssignedToClass(inputMetadatas)
        
    def _getGridSize(self):
        """ Figure out the layout of the plots given the number of references. """
        nrefs = len(self._visualizeRef3Ds)
        
        if nrefs == 1:
            gridsize = [1, 1]
        elif nrefs == 2:
            gridsize = [2, 1]
        else:
            gridsize = [(nrefs+1)/2, 2]
            
        return gridsize
                               
    def _getIterClasses(self, it):
        """ Return the .star file with the classes for this iteration.
        If the file doesn't exists, it will be created. 
        """
        addRelionLabels()
        data_star = self.getFilename('iter_data', iter=it)
        data_classes = self.getFilename('iter_classes', iter=it)
        
        if not xmippExists(data_classes):
            md = MetaData(data_star)
            mdClasses = MetaData()
            mdClasses.aggregate(md, AGGR_COUNT, MDL_REF3D, MDL_IMAGE, MDL_CLASS_COUNT)
            #md2.aggregate(md, AGGR_COUNT, MDL_REF3D, MDL_IMAGE, MDL_COUNT)
            mdClasses.write('classes@' + data_classes)
            mdImages = MetaData()
            
            for objId in mdClasses:
                ref3d = mdClasses.getValue(MDL_REF3D, objId)
                volFn = self.getFilename('volume', iter=it, ref3d=ref3d)
                mdClasses.setValue(MDL_IMAGE, volFn, objId)
                mdImages.importObjects(md, MDValueEQ(MDL_REF3D, ref3d))
                print "Writing to '%s', %d images" % (('class%06d_images' % ref3d) + data_classes, mdImages.size())
                mdImages.write(('class%06d_images@' % ref3d) + data_classes, MD_APPEND)
            mdClasses.write('classes@' + data_classes, MD_APPEND)
        return data_classes
    
    def _getIterAngularDist(self, it):
        """ Return the .star file with the classes angular distribution
        for this iteration. If the file not exists, it will be written.
        """
        addRelionLabels()
        data_classes = self._getIterClasses(it)
        mdClasses = MetaData('classes@' + data_classes)
        print mdClasses
        data_angularDist = self.getFilename('iter_angularDist', iter=it)
        
        if not xmippExists(data_angularDist):
            print "Writing angular distribution to: ", data_angularDist
            for objId in mdClasses:
                ref3d = mdClasses.getValue(MDL_REF3D, objId)
                print "Computing data class %03d" % ref3d
                mdImages = MetaData(('class%06d_images@' % ref3d) + data_classes)
                mdDist = MetaData()
                mdDist.aggregateMdGroupBy(mdImages, AGGR_COUNT, [MDL_ANGLE_ROT, MDL_ANGLE_TILT], MDL_ANGLE_ROT, MDL_WEIGHT)
                mdDist.setValueCol(MDL_ANGLE_PSI, 0.0)
                mdDist.write(('class%06d_angularDist@' % ref3d) + data_angularDist, MD_APPEND)
                    
        return data_angularDist
 
    def _visualizeDisplayAngularDistribution(self):        
        self.SpheresMaxradius = float(self.parser.getTkValue('SpheresMaxradius'))
        self.DisplayAngularDistribution = self.parser.getTkValue('DisplayAngularDistribution')
        volumeKeys = self._getVolumeKeys()
        
        if self.DisplayAngularDistribution == 'chimera':
#            fileNameVol = self.getFilename(volumeKeys[0], iter=iterations[-1], ref3d=1, workingDir=self.WorkingDir )
#            (Xdim, Ydim, Zdim, Ndim) = getImageSize(fileNameVol)
            
            outerRadius = int(float(self.MaskRadiusA)/self.SamplingRate)
            radius = float(outerRadius) * 1.1

            for it in self._visualizeIterations:
                data_angularDist = self._getIterAngularDist(it)
                for ref3d in self._visualizeRef3Ds:
                    volFn = self.getFilename('volume', iter=it, ref3d=ref3d)
                    args = " --mode projector 256 -a class%06d_angularDist@%s red %f " % (ref3d, data_angularDist, radius)
                    if self.SpheresMaxradius > 0:
                        args += ' %f ' % self.SpheresMaxradius
                        print "runChimeraClient(%s, %s) " % (volFn, args)
                    runChimeraClient(volFn, args)
                    
        elif self.DisplayAngularDistribution == "2D plots": 
            gridsize = self._getGridSize()
            
            for it in self._visualizeIterations:
                data_angularDist = self._getIterAngularDist(it)
                xplotter = XmippPlotter(*gridsize, mainTitle='Iteration %d' % it, windowTitle="Angular Distribution")
                for ref3d in self._visualizeRef3Ds:
                    volFn = self.getFilename('volume', iter=it, ref3d=ref3d)
                    md = MetaData("class%06d_angularDist@%s" % (ref3d, data_angularDist))
                    plot_title = 'Class %d' % ref3d
                    xplotter.plotAngularDistribution(plot_title, md)
                xplotter.show()        
        
    def _visualizeDisplayImagesClassification(self):
        """ Read Relion _data.star images file and 
        generate a new metadata with the Xmipp classification standard:
        a 'classes' block and a 'class00000?_images' block per class.
        If the new metadata was already written, it is just shown.
        """
        for it in self._visualizeIterations:
            data_classes = self._getIterClasses(it)
            self.display2D(data_classes, extraParams='--mode metadata --render first')

    def _visualizeDisplayResolutionPlotsSSNR(self):
            names = []
            windowTitle = {}
            if(self.relionType=='classify'):
                names = ['modelXm']
                windowTitle[names[0]]="ResolutionSSNR"
            else:
                names = ['half1_modelXm','half2_modelXm']
                windowTitle[names[0]]="ResolutionSSNR_half1"
                windowTitle[names[1]]="ResolutionSSNR_half2"
                
            gridsize = self._getGridSize()
            md = MetaData()
            activateMathExtensions()
            for name in names:
              xplotter = XmippPlotter(*gridsize, windowTitle=windowTitle[name])
              for ref3d in self._visualizeRef3Ds:
                plot_title = 'SSNR Class %s' % ref3d
                a = xplotter.createSubPlot(plot_title, 'Armstrongs^-1', 'log(SSNR)', yformat=False)
                legendName=[]
                blockName = 'model_class_%d@'%ref3d
                for it in self._visualizeIterations:
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
              xplotter.show()
              
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
              
    def _getChangeLabels(self):
        """ This method should be redefined in each subclass. """
        pass
        
    def _visualizeTableChange(self): 
            _r = []
            mdIters = MetaData()
            labels = self._getChangeLabels()
            iterations = range(1, self._visualizeLastIteration+1)
#            if(self.relionType=='classify'):
#                labels = ('Iter','offset','orientation','classAssigment')
#            else:
#                _c = ('Iter','offset','orientation')
            
            print " Computing average changes in offset, angles, and class membership"
            for it in iterations:
                print "Computing data for iteration; %03d"% (it)
                objId = mdIters.addObject()
                mdIters.setValue(MDL_ITER, it, objId)
                #agregar por ref3D
                fn = self.getFilename('optimiser'+'Xm', iter=it )
                md = MetaData(fn)
                firstId = md.firstObject()
                for label in labels:
                    mdIters.setValue(label, md.getValue(label, firstId), objId)
            fn = self.getFilename('changes')
            mdIters.write(fn)
            self.display2D(fn)#                if(self.relionType=='classify'):
#                    _r.append(("Iter_%d" % it, cOrientations,cOffsets,cClasses))
#                else:
#                    _r.append(("Iter_%d" % it, cOrientations,cOffsets))
#            showTable(_c,_r, title='Changes per Iteration', width=100)
#                if(self.relionType=='classify'):
#                    _r.append(("Iter_%d" % it, cOrientations,cOffsets,cClasses))
#                else:
#                    _r.append(("Iter_%d" % it, cOrientations,cOffsets))
#            showTable(_c,_r, title='Changes per Iteration', width=100)
            
    def _visualizeLikelihood(self):  
        for it in self._visualizeIterations:
            print "Computing data for iteration %03d"% (it)
            fn = 'images@'+ self.getFilename('dataXm', iter=it)
            md = MetaData(fn)
            md.sort(MDL_LL, False)
            md.write(fn)
            runShowJ(fn)
            xplotter = XmippPlotter(windowTitle="max Likelihood particles sorting Iter_%d"%it)
            xplotter.createSubPlot("Particle sorting: Iter_%d"%it, "Particle number", "maxLL")
            xplotter.plotMd(md, False, mdLabelY=MDL_LL)
        xplotter.show()
        
    def _visualizeAvgPMAX(self):         
            _r = []
            mdOut = MetaData()
            
            if(self.relionType=='classify'):
                addRelionLabels()
                mdIters = MetaData()
                iterations = range(1, self._visualizeLastIteration+1)
                
                for it in iterations: # range (firstIter,self._visualizeLastIteration+1): #alwaya list all iteration
                    fn = 'model_general@'+ self.getFilename('modelRe', iter=it)
                    md = MetaData(fn)
                    pmax = md.getValue(MDL_AVGPMAX,md.firstObject())
                    objId = mdIters.addObject()
                    mdIters.setValue(MDL_ITER, it, objId)
                    mdIters.setValue(MDL_AVGPMAX, pmax, objId)
                fn = self.getFilename('avgPmax')
                mdIters.write(fn)
                self.display2D(fn)                    
                #self.display2D(fn, extraParams)
                xplotter = XmippPlotter()
                xplotter.createSubPlot("Avg PMax per Iterations", "Iterations", "Avg PMax")
                xplotter.plotMd(mdIters, MDL_ITER, MDL_AVGPMAX)
                xplotter.show()
            else:
                _c = ('Iteration','Avg PMax h1', 'Avg PMax h2')
                
                for it in self._visualizeIterations:#range (firstIter,self._visualizeLastIteration+1): #alwaya list all iteration
                    fn = 'model_general@'+ self.getFilename('half1_model'+'Xm', iter=it , workingDir=self.WorkingDir)
                    md = MetaData(fn)
                    pmax1 = md.getValue(MDL_AVGPMAX,md.firstObject())
                    fn = 'model_general@'+ self.getFilename('half2_model'+'Xm', iter=it , workingDir=self.WorkingDir)
                    md = MetaData(fn)
                    pmax2 = md.getValue(MDL_AVGPMAX,md.firstObject())
                    _r.append(("Iter_%d" % it, pmax1,pmax2))
                fn = 'model_general@'+ self.getFilename('modelXmFinalXm', workingDir=self.WorkingDir)
                if(self.relionType=='refine' and self._visualizeIterations[-1]==self.NumberOfIterations):
                    md = MetaData(fn)
                    pmax1 = md.getValue(MDL_AVGPMAX,md.firstObject())
                    _r.append(("all_particles", pmax1,'---'))

    def _visualizeDisplayReconstruction(self): 
            volumeKeys = self._getVolumeKeys()
                
            for ref3d in self._visualizeRef3Ds:
                for it in self._visualizeIterations:
                  for volKey in volumeKeys:
                    self.display3D(self.getFilename(volKey, iter=it, ref3d=ref3d , workingDir=self.WorkingDir))
            
            if self.relionType=='refine' and iterations[-1]==self.NumberOfIterations:
                self.display3D(self.getFilename('volumeFinal',ref3d=1, workingDir=self.WorkingDir))
                
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
