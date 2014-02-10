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
                               replaceBasenameExt, replaceFilenameExt, deleteFile
from protlib_utils import runShowJ, getListFromRangeString, \
                          runJob, runChimera, which, runChimeraClient
from protlib_import import exportReliontoMetadataFile, addRelionLabels
from protlib_gui_ext import showWarning, showTable, showError
from protlib_gui_figure import XmippPlotter
from protlib_xmipp import getSampling, validateInputSize


class ProtRelionBase(XmippProtocol):
    def __init__(self, protDictName, scriptname, project):
        self.FileKeys = ['data', 'optimiser', 'sampling'] 
        XmippProtocol.__init__(self, protDictName, scriptname, project)
        self.ParamsStr = ''
        
        program = 'relion_refine'
        if self.NumberOfMpi > 1:
            program += '_mpi'
        
        if self.DoContinue:
            if os.path.exists(self.optimiserFileName):
                self.setPreviousRunFromFile(self.optimiserFileName)
            try:
                self.inputProperty('ImgMd')
            except:
                print "Cannot find the input file with experimental images in protocol %s" % self.PrevRunName
        
        self.addParam('program', program)
        self.addParam('ORoot', self.WorkingDir + '/')
        self.addParam('SamplingRate', getSampling(self.ImgMd, 0.))
        
    def createFilenameTemplates(self):
        """ Define keywords for getting filenames for each iteration. """
        self.extraIter = self.extraPath('relion_it%(iter)03d_')
        myDict = {
                  #'data': self.extraIter + 'data.star',
                  #'volume': self.extraIter + 'class%(ref3d)03d.mrc:mrc',
                  
                  'data_sorted_xmipp': self.extraIter + 'data_sorted_xmipp.star',
                  'classes_xmipp': self.extraIter + 'classes_xmipp.xmd',
                  'angularDist_xmipp': self.extraIter + 'angularDist_xmipp.xmd',
                  'all_avgPmax_xmipp': self.tmpPath('iterations_avgPmax_xmipp.xmd'),
                  'all_changes_xmipp': self.tmpPath('iterations_changes_xmipp.xmd'),
                  'selected_volumes': self.tmpPath('selected_volumes_xmipp.xmd')
                  }
        
        # add to keys, data.star, optimiser.star and sampling.star
        for key in self.FileKeys:
            myDict[key] = self.extraIter + '%s.star' % key
            key_xmipp = key + '_xmipp'             
            myDict[key_xmipp] = self.extraIter + '%s.xmd' % key
        # add other keys that depends on prefixes
        for p in self._getPrefixes():            
            myDict['%smodel' % p] = self.extraIter + '%smodel.star' % p
            myDict['%svolume' % p] = self.extraIter + p + 'class%(ref3d)03d.mrc:mrc'

        return myDict

    def _summary(self):
        """ Should be overriden in subclasses to 
        return summary message for NORMAL EXECUTION. 
        """
        pass
    
    def _summaryContinue(self):
        """ Should be overriden in subclasses to
        return summary messages for CONTINUE EXECUTION.
        """
        pass
    
    def summary(self):
        if self.DoContinue:
            return self._summaryContinue()
        return self._summary()
        #return ["Wrapper implemented for RELION **1.2**"]

    def papers(self):
        papers = ['Bayesian view: Scheres, JMB (2012) [http://www.ncbi.nlm.nih.gov/pubmed/22100448]',
                  'RELION implementation: Scheres, JSB (2012) [http://www.ncbi.nlm.nih.gov/pubmed/23000701]'
                  ]
        return papers

    def lastIter(self):
        fileNameTemplate = self.getFilename('data', iter=0)
        fileNameTemplate = fileNameTemplate.replace('000','???')
        a = sorted(glob(fileNameTemplate))
        if not a:
            return 0 
        return int (re.findall(r'_it\d+_',a[-1])[0][3:6])

    def firstIter(self):
        fileNameTemplate = self.getFilename('data', iter=0)
        print "fileNameTemplate", fileNameTemplate
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
            
    def _validate(self):
        """ Should be overriden in subclasses to 
        return summary message for NORMAL EXECUTION. 
        """
        return []
    
    def _validateContinue(self):
        """ Should be overriden in subclasses to
        return summary messages for CONTINUE EXECUTION.
        """
        return []
               
    def validate(self):
        errors = []
        md = MetaData()
        md.read(self.ImgMd, 1) # Read only the first object in md
        
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
        
        if self.DoContinue:
            errors += self._validateContinue()
        else:
            errors += self._validate()
            
        return errors 
    
    
    
    def defineSteps(self): 
        #temporary normalized files
        images_xmd = join(self.ExtraDir, 'images.xmd')
        #create extra output directory
        self.insertStep('createDir', verifyfiles=[self.ExtraDir], path=self.ExtraDir)
        self.insertStep("linkAcquisitionInfo", InputFile=self.ImgMd, dirDest=self.WorkingDir) 
        #normalize data

        if self.DoNormalizeInputImage:
            images_stk = join(self.ExtraDir, 'images_normalized.stk')
            self.insertStep('runNormalizeRelion', verifyfiles=[images_stk, images_xmd],
                            inputMd  = self.ImgMd,
                            outputMd = images_stk,
                            normType = 'NewXmipp',
                            bgRadius = int(self.MaskRadiusA/self.SamplingRate),
                            Nproc    = self.NumberOfMpi*self.NumberOfThreads
                            )
        else:
            self.insertStep('createLink',
                               verifyfiles=[images_xmd],
                               source=self.ImgMd,
                               dest=images_xmd)
        # convert input metadata to relion model
        self.ImgStar = self.extraPath(replaceBasenameExt(images_xmd, '.star'))
        self.insertStep('convertImagesMd', verifyfiles=[self.ImgStar],
                        inputMd=images_xmd, 
                        outputRelion=self.ImgStar                            
                        )
        
        if self.DoContinue:
            self._insertStepsContinue()
            firstIteration = getIteration(self.optimiserFileName)
        else:
            self._insertSteps()
            firstIteration = 1
            
        self.defineSteps2(firstIteration, self.NumberOfIterations, self.NumberOfClasses)
        
    def defineSteps2(self, firstIteration, lastIteration, NumberOfClasses):
        verifyFiles = []
        inputs = []
        #TODO: check this
#        for i in range (firstIteration, lastIteration+1):
#            for v in self.FileKeys:
#                 verifyFiles += [self.getFilename(v+'Xm', iter=i, workingDir=self.WorkingDir )]
        verifyFiles += self._getExtraOutputs()
        #data to store in Working dir so it can be easily accessed by users
        lastIterVolFns = []
        standardOutputClassFns = []
        classesFn = "images_ref3d%06d@" + self.workingDirPath("classes_ref3D.xmd")
        finalSuffix = self._getFinalSuffix()
        
        for ref3d in range (1, NumberOfClasses+1):
            lastIterVolFns.append(self.getFilename('volume' + finalSuffix, iter=lastIteration, ref3d=ref3d))
            standardOutputClassFns.append(classesFn % ref3d)
        
        lastIterationMetadata = "images@" + self.getFilename('data_xmipp', iter=lastIteration)
        
        inputMetadatas = []
        for it in range (firstIteration, self.NumberOfIterations+1): #always list all iterations
            inputMetadatas.append("images@" + self.getFilename('data_xmipp', iter=it))

#        self.insertStep('convertRelionMetadata', verifyfiles=verifyFiles,
#                        inputs  = inputs,
#                        outputs = verifyFiles,
#                        lastIterationVolumeFns = lastIterVolFns,
#                        lastIterationMetadata  = lastIterationMetadata,
#                        NumberOfClasses        = NumberOfClasses,
#                        standardOutputClassFns = standardOutputClassFns,
#                        standardOutputImageFn  = "images@" + self.workingDirPath("images.xmd"),
#                        standardOutputVolumeFn = "volumes@" + self.workingDirPath("volumes.xmd"),
#                        extraDir=self.ExtraDir,
#                        inputMetadatas=inputMetadatas
#                        )
        
        self.insertStep('produceXmippOutputs', verifyfiles=[],
                        inputDataStar=self.getFilename('data' + finalSuffix, iter=lastIteration),
                        outputImages=self.workingDirPath('images.xmd'), 
                        outputClasses=self.workingDirPath('classes_ref3D.xmd'),
                        outputVolumes=self.workingDirPath('volumes.xmd')
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
        
        self._visualizeNrefs = len(self._visualizeRef3Ds)
        lastIteration = self.lastIter()
        firstIter = self.firstIter()
        
        if VisualizeIter == 'last':
            self._visualizeIterations = [lastIteration]
        else:
            self._visualizeIterations = getListFromRangeString(self.parser.getTkValue('SelectedIters'))

        self._visualizeLastIteration = self._visualizeIterations[-1]
        self._visualizeLastRef3D = self._visualizeRef3Ds[-1]
        
        self._visualizeVolumesMode = self.parser.getTkValue('DisplayReconstruction')
        
        
    def _checkIterData(self):
        """ check that for last visualized iteration, the data is produced. """
        volumes = []
        for prefix in self._getPrefixes():
            volumes.append(self.getFilename(prefix + 'volume', 
                                            iter=self._visualizeLastIteration, 
                                            ref3d=self._visualizeLastRef3D))
                
        for v in volumes:
            if not xmippExists(v):
                message = "No data available for <iteration %d> and <class %d>" % (self._visualizeLastIteration, self._visualizeLastRef3D)
                message += '\nVolume <%s> not found. ' % v
                showError("File not Found", message, self.master)
                return False
        return True
        
        
    def launchRelionPlots(self, selectedPlots):
        ############TABLES STOP ANALIZE
        ''' Launch some plots for a Projection Matching protocol run '''
        from tempfile import NamedTemporaryFile    
        TmpDir=join(self.WorkingDir,"tmp")
        import numpy as np
        _log = self.Log

        xplotter=None
        self._plot_count = 0
        
        self._loadVisualizeItersAndRefs()
        
        if not self._checkIterData():
            return
        
        def doPlot(plotName):
            return plotName in selectedPlots        
        
        if doPlot('AvgPMAX'):
            self._visualizeAvgPMAX()
                
        if doPlot('DisplayReconstruction'):
            self._visualizeDisplayReconstruction()
                
        if doPlot('DisplayImagesClassification'):
            self._visualizeDisplayImagesClassification()

        if doPlot('DisplayResolutionPlotsSSNR'):
            self._visualizeDisplayResolutionPlotsSSNR()
                    
        if doPlot('DisplayResolutionPlotsFSC'):
            self._visualizeDisplayResolutionPlotsFSC()

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
                # We use os.system to load the relion labels through showj
                os.system('xmipp_showj %(fn)s %(extraParams)s --dont_wrap &' % locals())
                #runShowJ(fn, extraParams=extraParams + ' --dont_wrap')
            except Exception, e:
                showError("Error launching java app", str(e), self.master)
        else:
            showError("Can not visualize file <%s>\nit does not exists." % fn, self.master)

    def imagesAssignedToClass2(self,firstIter, lastIteration):
            inputMetadatas = []
            for it in range (firstIter,lastIteration+1): #always list all iterations
                inputMetadatas += ["model_classes@"+self.getFilename('model'+'Xm', iter=it, workingDir=self.WorkingDir )]
            return imagesAssignedToClass(inputMetadatas)
        
    def _getGridSize(self, n=None):
        """ Figure out the layout of the plots given the number of references. """
        if n is None:
            n = self._visualizeNrefs
        
        if n == 1:
            gridsize = [1, 1]
        elif n == 2:
            gridsize = [2, 1]
        else:
            gridsize = [(n+1)/2, 2]
            
        return gridsize
                               
    def _getIterClasses(self, it):
        """ Return the .star file with the classes for this iteration.
        If the file doesn't exists, it will be created. 
        """
        addRelionLabels()
        data_star = self.getFilename('data', iter=it)
        data_classes = self.getFilename('classes_xmipp', iter=it)
        
        if not xmippExists(data_classes):
            createClassesFromImages(data_star, data_classes, it)
        return data_classes
    
    def _writeIterAngularDist(self, it):
        """ Write the angular distribution. Should be overriden in subclasses. """
        pass
    
    def _getIterAngularDist(self, it):
        """ Return the .star file with the classes angular distribution
        for this iteration. If the file not exists, it will be written.
        """
        data_angularDist = self.getFilename('angularDist_xmipp', iter=it)
        addRelionLabels()
        
        if not xmippExists(data_angularDist):
            self._writeIterAngularDist(it)
 
        return data_angularDist
    
    def _getIterSortedData(self, it):
        """ Sort the it??.data.star file by the maximum likelihood. """
        addRelionLabels(True)
        data_sorted = self.getFilename('data_sorted_xmipp', iter=it)
        if not xmippExists(data_sorted):
            print "Sorting particles by likelihood iteration %03d" % it
            fn = 'images@'+ self.getFilename('data', iter=it)
            md = MetaData(fn)
            md.sort(MDL_LL, False)
            md.write('images_sorted@' + data_sorted)
        
        return data_sorted
            
    def _visualizeDisplayReconstruction(self):
        """ Visualize selected volumes per iterations.
        If show in slices, create a temporal single metadata
        for avoid opening multiple windows. 
        """
        prefixes = self._getPrefixes()
        
        if self._visualizeVolumesMode == 'slices':
            fn = self.getFilename('selected_volumes')
            if xmippExists(fn):
                deleteFile(None, fn)
                
            for it in self._visualizeIterations:
                volMd = MetaData()
                for prefix in prefixes:
                    for ref3d in self._visualizeRef3Ds:
                        volFn = self.getFilename('%svolume' % prefix, iter=it, ref3d=ref3d)
                        volMd.setValue(MDL_IMAGE, volFn, volMd.addObject())
                block = 'volumes_iter%06d@' % it
                volMd.write(block + fn, MD_APPEND)
            self.display2D(fn)         
        else:
            for it in self._visualizeIterations:
              for prefix in prefixes:
                  for ref3d in self._visualizeRef3Ds:
                    self.display3D(self.getFilename('%svolume' % prefix, iter=it, ref3d=ref3d))
   

    def _visualizeDisplayResolutionPlotsFSC(self):
        self.ResolutionThreshold = float(self.parser.getTkValue('ResolutionThreshold'))
        n = self._visualizeNrefs * len(self._getPrefixes())
        gridsize = self._getGridSize(n)
        
        activateMathExtensions()
        addRelionLabels()
        
        xplotter = XmippPlotter(*gridsize, windowTitle='Resolution FSC')

        for prefix in self._getPrefixes():
          for ref3d in self._visualizeRef3Ds:
            plot_title = prefix + 'class %s' % ref3d
            a = xplotter.createSubPlot(plot_title, 'Armstrongs^-1', 'FSC', yformat=False)
            legendName=[]
            blockName = 'model_class_%d@' % ref3d
            for it in self._visualizeIterations:
                model_star = blockName + self.getFilename(prefix + 'model', iter=it)
                if xmippExists(model_star):
                    md = MetaData(model_star)
                    resolution_inv = [md.getValue(MDL_RESOLUTION_FREQ, id) for id in md]
                    frc = [md.getValue(MDL_RESOLUTION_FRC, id) for id in md]
                    a.plot(resolution_inv, frc)
                    legendName.append('iter %d' % it)
            xplotter.showLegend(legendName)
            if self.ResolutionThreshold < max(frc):
                a.plot([min(resolution_inv), max(resolution_inv)],[self.ResolutionThreshold, self.ResolutionThreshold], color='black', linestyle='--')
            a.grid(True)
            
        xplotter.show(True)
          
#        if(self.relionType=='refine'):
#            file_name = blockName + self.getFilename('modelXmFinalXm', workingDir=self.WorkingDir)
#            print "final model", file_name
#            if xmippExists(file_name):
#                print "final model exists"
#                gridsize1 = [1, 1]
#                xplotter = XmippPlotter(*gridsize1,windowTitle="ResolutionSSNRAll")
#                plot_title = 'FSC for all images, final iteration'
#                a = xplotter.createSubPlot(plot_title, 'Armstrongs^-1', 'log(SSNR)', yformat=False)
#                blockName = 'model_class_%d@'%1
#                md = MetaData(file_name)
#                # only cross by 1 is important
#                resolution_inv = [md.getValue(MDL_RESOLUTION_FREQ, id) for id in md]
#                frc = [md.getValue(MDL_RESOLUTION_FRC, id) for id in md]
#                a.plot(resolution_inv, frc)
#                if (self.ResolutionThreshold < max(frc)):
#                    a.plot([min(resolution_inv), max(resolution_inv)],[self.ResolutionThreshold, self.ResolutionThreshold], color='black', linestyle='--')
#                    a.grid(True)
#                xplotter.draw()
                                      
    def _visualizeDisplayAngularDistribution(self):        
        self.SpheresMaxradius = float(self.parser.getTkValue('SpheresMaxradius'))
        self.DisplayAngularDistribution = self.parser.getTkValue('DisplayAngularDistribution')
        
        if self.DisplayAngularDistribution == 'chimera':
            outerRadius = int(float(self.MaskRadiusA)/self.SamplingRate)
            radius = float(outerRadius) * 1.1

            for it in self._visualizeIterations:
                data_angularDist = self._getIterAngularDist(it)
                for ref3d in self._visualizeRef3Ds:
                    for prefix in self._getPrefixes():
                        volFn = self.getFilename(prefix + 'volume', iter=it, ref3d=ref3d)
                        args = " --mode projector 256 -a %sclass%06d_angularDist@%s red %f " % (prefix, ref3d, data_angularDist, radius)
                        if self.SpheresMaxradius > 0:
                            args += ' %f ' % self.SpheresMaxradius
                            print "runChimeraClient(%s, %s) " % (volFn, args)
                        runChimeraClient(volFn, args)
                    
        elif self.DisplayAngularDistribution == "2D plots": 
            n = self._visualizeNrefs * len(self._getPrefixes())
            gridsize = self._getGridSize(n)
            
            for it in self._visualizeIterations:
                data_angularDist = self._getIterAngularDist(it)
                xplotter = XmippPlotter(*gridsize, mainTitle='Iteration %d' % it, windowTitle="Angular Distribution")
                for ref3d in self._visualizeRef3Ds:
                    for prefix in self._getPrefixes():
                        md = MetaData("class%06d_angularDist@%s" % (ref3d, data_angularDist))
                        plot_title = '%sclass %d' % (prefix, ref3d)
                        xplotter.plotAngularDistribution(plot_title, md)
                xplotter.show(interactive=True)        
        
    def _getPrefixes(self):
        """ This should be redefined in sub-classes. """
        pass
    
    def _getIterFiles(self, iter):
        """ Return the files that should be produces for a give iteration. """
        return [self.getFilename(k, iter=iter) for k in self.FileKeys]
    
    def _visualizeDisplayResolutionPlotsSSNR(self):
        names = []
        windowTitle = {}
        prefixes = self._getPrefixes()
        
        n = len(prefixes) * self._visualizeNrefs
        gridsize = self._getGridSize(n)
        addRelionLabels()
        md = MetaData()
        activateMathExtensions()
        xplotter = XmippPlotter(*gridsize)
        
        for prefix in prefixes:
            for ref3d in self._visualizeRef3Ds:
              plot_title = 'Resolution SSNR %s, for Class %s' % (prefix, ref3d)
              a = xplotter.createSubPlot(plot_title, 'Armstrongs^-1', 'log(SSNR)', yformat=False)
              legendName=[]
              blockName = 'model_class_%d@' % ref3d
              for it in self._visualizeIterations:
                  file_name = blockName + self.getFilename(prefix + 'model', iter=it)
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
                      legendName.append('iter %d' % it)
                  xplotter.showLegend(legendName)
              a.grid(True)
        xplotter.show(True)
              
#            if(self.relionType=='refine'):
#                file_name = blockName + self.getFilename('modelXmFinalXm')
#                if xmippExists(file_name):
#                    gridsize1 = [1, 1]
#                    xplotter = XmippPlotter(*gridsize1,windowTitle="ResolutionSSNRAll")
#                    plot_title = 'SSNR for all images, final iteration'
#                    a = xplotter.createSubPlot(plot_title, 'Armstrongs^-1', 'log(SSNR)', yformat=False)
#                    blockName = 'model_class_%d@'%1
#                    mdOut = MetaData(file_name)
#                    md.clear()
#                    # only cross by 1 is important
#                    md.importObjects(mdOut, MDValueGT(MDL_RESOLUTION_SSNR, 0.9))
#                    md.operate("resolutionSSNR=log(resolutionSSNR)")
#                    resolution_inv = [md.getValue(MDL_RESOLUTION_FREQ, id) for id in md]
#                    frc = [md.getValue(MDL_RESOLUTION_SSNR, id) for id in md]
#                    a.plot(resolution_inv, frc)
#                    a.grid(True)
#                    xplotter.draw()
              
    def _getChangeLabels(self):
        """ This method should be redefined in each subclass. """
        pass
        
    def _visualizeTableChange(self): 
            mdIters = MetaData()
            labels = self._getChangeLabels()
            iterations = range(1, self._visualizeLastIteration+1)
            addRelionLabels()
            
            print " Computing average changes in offset, angles, and class membership"
            for it in iterations:
                print "Computing data for iteration; %03d"% (it)
                objId = mdIters.addObject()
                mdIters.setValue(MDL_ITER, it, objId)
                #agregar por ref3D
                fn = self.getFilename('optimiser', iter=it )
                md = MetaData(fn)
                firstId = md.firstObject()
                for label in labels:
                    mdIters.setValue(label, md.getValue(label, firstId), objId)
            fn = self.getFilename('all_changes_xmipp')
            mdIters.write(fn)
            self.display2D(fn, extraParams='--label_relion')              
            
    def _visualizeLikelihood(self):  
        for it in self._visualizeIterations:
            fn = self._getIterSortedData(it)
            md = MetaData(fn)
            self.display2D(fn)
            xplotter = XmippPlotter(windowTitle="max Likelihood particles sorting Iter_%d"%it)
            xplotter.createSubPlot("Particle sorting: Iter_%d" % it, "Particle number", "maxLL")
            xplotter.plotMd(md, False, mdLabelY=MDL_LL)
            xplotter.show(True)
        
    def _visualizeAvgPMAX(self):         
            addRelionLabels()            
            mdIters = MetaData()
            iterations = range(1, self._visualizeLastIteration+1)
            labels = [MDL_AVGPMAX, MDL_PMAX]
            colors = ['g', 'b']
            
            for it in iterations: # range (firstIter,self._visualizeLastIteration+1): #alwaya list all iteration
                objId = mdIters.addObject()
                mdIters.setValue(MDL_ITER, it, objId)
                for i, prefix in enumerate(self._getPrefixes()):
                    fn = 'model_general@'+ self.getFilename(prefix + 'model', iter=it)
                    md = MetaData(fn)
                    pmax = md.getValue(MDL_AVGPMAX, md.firstObject())
                    mdIters.setValue(labels[i], pmax, objId)
            fn = self.getFilename('all_avgPmax_xmipp')
            mdIters.write(fn)
            self.display2D(fn, extraParams='--label_relion')                    
            #self.display2D(fn, extraParams)
            xplotter = XmippPlotter()
            xplotter.createSubPlot("Avg PMax per Iterations", "Iterations", "Avg PMax")
            for label, color in zip(labels, colors):
                xplotter.plotMd(mdIters, MDL_ITER, label, color)
            xplotter.showLegend(self._getPrefixes())
            xplotter.show(True)

                
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
    f = f.removeBlockName()
    if exists(f):
        os.remove(f)
    #TODO: write classes block?
    for i in range (1, NumberOfClasses+1):
        mdOut.clear()
        mdOut.importObjects(md, MDValueEQ(MDL_REF3D, i+1))
        mdOut.write(standardOutputClassFns[i],MD_APPEND)
        
    #volume.xmd, metada with volumes
    mdOut.clear()
    for lastIterationVolumeFn in lastIterationVolumeFns:
        objId = mdOut.addObject()
        mdOut.setValue(MDL_IMAGE, lastIterationVolumeFn, objId)
    mdOut.write(standardOutputVolumeFn)
    
    
def produceXmippOutputs(log, inputDataStar, outputImages, outputClasses, outputVolumes, 
                        ):
    """ After Relion have finished, produce the files expected by Xmipp.
    Params:
        inputDataStar: the filename of the last data.star images file.
        outputImages: the filename to write the images as expected by xmipp.
        outputClasses: the filename to write the classes and images as expected by xmipp.
        outputVolumes: write a list of produced volumes.
    """
    addRelionLabels()
    # Write output images
    md = MetaData(inputDataStar)
    md.write(outputImages)
    # Create the classes metadata
    iteration = getIteration(inputDataStar)
    # If iteration is not None, means that is a classification
    # since for refinenment, the final images should not have iteration
    if iteration is not None:
        createClassesFromImages(inputDataStar, outputClasses, iteration)
        # Write volumes (the same as classes block)
        md.read('classes@' + outputClasses)
        md.write('volumes@' + outputVolumes)
    else:
        volFn = inputDataStar.replace('data.star', 'class001.mrc')
        md.clear()
        md.setValue(MDL_IMAGE, volFn, md.addObject())
    

def convertRelionBinaryData(log, inputs,outputs):
    """Make sure mrc files are volumes properlly defined"""
    program = "xmipp_image_convert"
    for i,o in zip(inputs,outputs):
        args = "-i %s -o %s  --type vol"%(i,o)
        runJob(log, program, args )

def runNormalizeRelion(log, inputMd, outputMd, normType, bgRadius, Nproc,):
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
    #md.dropColumn(MDL_CTF_MODEL)
    # Create the mapping between relion labels and xmipp labels
    #addRelionLabels(True, True)
    #md.write(outputRelion)
    exportMdToRelion(md, outputRelion)
    
def createClassesFromImages(imgsFn, classesFn, iter):
    """ Give a Relion data file, produces a classes metadata as expected by Xmipp.
    Params:
        imgsFn: the filename with relion star format.
        classesFn: filename where to write the classes.
    Note: it is assumed that addRelionLabels is done before calling this function.
    """
    md = MetaData(imgsFn)
    mdClasses = MetaData()
    mdClasses.aggregate(md, AGGR_COUNT, MDL_REF3D, MDL_IMAGE, MDL_CLASS_COUNT)
    #md2.aggregate(md, AGGR_COUNT, MDL_REF3D, MDL_IMAGE, MDL_COUNT)
    mdClasses.write('classes@' + classesFn)
    mdImages = MetaData()
    # We asume here that the volumes (classes3d) are in the same folder than imgsFn
    rootDir = os.path.dirname(imgsFn)
    volPattern = os.path.join(rootDir, 'relion_it%(iter)03d_class%(ref3d)03d.mrc:mrc')
    
    for objId in mdClasses:
        ref3d = mdClasses.getValue(MDL_REF3D, objId)
        volFn = volPattern % locals() # replace iter and class3d
        mdClasses.setValue(MDL_IMAGE, volFn, objId)
        mdImages.importObjects(md, MDValueEQ(MDL_REF3D, ref3d))
        print "Writing to '%s', %d images" % (('class%06d_images' % ref3d) + classesFn, mdImages.size())
        mdImages.write(('class%06d_images@' % ref3d) + classesFn, MD_APPEND)
    mdClasses.write('classes@' + classesFn, MD_APPEND)

def getIteration(fileName):
    matchList = re.findall(r'_it\d+_',fileName)
    if len(matchList):
        return int(matchList[0][3:6])
    return None

def getClass(fileName):
   return int(re.findall(r'_class\d+',fileName)[0][6:9])
