# **************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
# *              J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Josue Gomez Blanco (jgomez@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************
"""
This module implement the wrappers aroung Xmipp ML2D protocol
visualization program.
"""
import os
from os.path import exists

from pyworkflow.protocol.executor import StepExecutor
from pyworkflow.viewer import CommandView, Viewer, ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO
from pyworkflow.em.viewer import DataView, ClassesView, Classes3DView
from pyworkflow.gui.form import FormWindow
from pyworkflow.utils import createUniqueFileName, cleanPattern
from protocol_projmatch import XmippProtProjMatch
# from projmatch_initialize import createFilenameTemplates
from pyworkflow.em.packages.xmipp3.convert import * # change this
from pyworkflow.protocol.constants import LEVEL_EXPERT, LEVEL_ADVANCED
from pyworkflow.protocol.params import (PointerParam, BooleanParam, IntParam, 
                                        FloatParam, StringParam, Positive, GE,
                                        EnumParam, NumericRangeParam, TextParam,
                                        DigFreqParam)
import xmipp
from pyworkflow.em.packages.xmipp3.plotter import XmippPlotter
import numpy as np


ITER_LAST = 0
ITER_SELECTION = 1

ANGDIST_2DPLOT = 0
ANGDIST_CHIMERA = 1

VOLUME_SLICES = 0
VOLUME_CHIMERA = 1

REF_ALL = 0
REF_SEL = 1

CLASSES_ALL = 0
CLASSES_SEL = 1

class XmippProjMatchViewer(ProtocolViewer):
    """ Wrapper to visualize different type of data objects
    with the Xmipp program xmipp_showj
    """
    _targets = [XmippProtProjMatch]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    
    _label = 'viewer projection matching'
#     _plotVars = ['doShowLL', 'doShowPmax', 'doShowSignalChange', 'doShowMirror'] 
    
    def _defineParams(self, form):
        form.addSection(label='Visualization')
        group = form.addGroup('Overall results')
        form.addParam('viewIter', EnumParam, choices=['last', 'selection'], default=ITER_LAST, 
                      display=EnumParam.DISPLAY_LIST,
                      label="Iteration to visualize", 
                      help="""
*last*: only the last iteration will be visualized.
*selection*: you may specify a range of iterations.
Examples:
"1,5-8,10" -> [1,5,6,7,8,10]
"2,6,9-11" -> [2,6,9,10,11]
"2 5, 6-8" -> [2,5,6,7,8]                      
                           """)
        form.addParam('iterSelection', NumericRangeParam, 
                      condition='viewIter==%d' % ITER_SELECTION, 
                      label="Iterations list", 
                      help="Write the iteration list to visualize.")
        
        group = form.addGroup('Volumes display')
        group.addParam('showRef3DNo', EnumParam, choices=['all', 'selection'], default=REF_ALL,
                      display=EnumParam.DISPLAY_LIST,
                      label='Show results for reference 3D volumes?',
                      help='')
        group.addParam('ref3DSelection', NumericRangeParam, default='1',
                      condition='showRef3DNo == %d' % REF_SEL,
                      label='references list',
                      help='')
        group.addParam('matrixWidth', FloatParam, default=-1, 
                      expertLevel=LEVEL_EXPERT,
                      label='Width of projection galleries',
                      help='Usually a multiple of 2 is the right value. -1 => authomatic')
        group.addParam('showVolume', EnumParam, choices=['slices', 'chimera'], 
                      display=EnumParam.DISPLAY_LIST, default=VOLUME_SLICES,
                      label='Display volume with',
                      help='*slices*: display volumes as 2D slices along z axis.\n'
                           '*chimera*: display volumes as surface with Chimera.')
        group.addParam('showReference', BooleanParam, default=False,
                      label='Show reference volume?',
                      help='Volume after filtration and masking')
        group.addParam('showReconstruction', BooleanParam, default=False,
                      label='Show reconstructed volume?',
                      help='Volume as given by the reconstruction algorithm')
        group.addParam('showFilteredReconstruction', BooleanParam, default=False,
                      label='Show reconstructed volume after filtration?',
                      help='Volume after filtration')
        group.addParam('showBFactorCorrectedVolume', BooleanParam, default=False,
                       label='Show a b_factor corrected volume?',
                       help=""" This utility boost up the high frequencies. Do not use the automated 
                            mode [default] for maps with resolutions lower than 12-15 Angstroms.
                            It does not make sense to apply the Bfactor to the firsts iterations
                            see http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Correct_bfactor.
                            NOTE: bfactor will be applied ONLY to the reconstructed volumes NOT to
                            the reference ones
                            """)
        group.addParam('maxRes', FloatParam, default=12, 
                      condition='showBFactorCorrectedVolume',
                      label='Maximum resolution to apply B-factor (in Angstroms)')
        group.addParam('correctBfactorExtraCommand', StringParam, default='--auto',
                       condition='showBFactorCorrectedVolume', 
                       label='User defined flags for the correct_bfactor program',
                       help=""" See http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Correct_bfactor
                            for details. DEFAULT behaviour is --auto
                            """)
        
        group = form.addGroup('Library, classes and images')
        group.addParam('showProjectionMatchingLibrary', BooleanParam, default=False,
                      label='Display library?')
        group.addParam('showProjectionMatchingClasses', BooleanParam, default=False,
                      label='Display classes?')
        group.addParam('showProjectionMatchingLibraryAndClasses', BooleanParam, default=False,
                      label='Display library and classes in a single image?')
        group.addParam('showProjectionMatchingLibraryAndImages', BooleanParam, default=False,
                      label='Display library and experimental images in a single image?')
        group.addParam('showDiscardedImages', BooleanParam, default=False,
                      label='Display input discarded images?')
        group.addParam('showExperimentalImages', BooleanParam, default=False,
                      label='Display experimental images?',
                      help="""Display Input images with alignment and classification information
                           WARNING: the angles and shifts are the adecuate for reconstruction
                           but not for 2D aligment.
                           """)
        group = form.addGroup('Convergence')
        group.addParam('plotHistogramAngularMovement', BooleanParam, default=False,
                      label='Plot histogram with angular/shift changes?',
                      help=""" Plot histogram with angular changes from one iteration to next. 
                           Iteration 0 -> initial values
                           """)
        group.addParam('numberOfBins', IntParam, default=50, 
                      condition='plotHistogramAngularMovement',
                      label='Number of bins (for histogram)',
                      help='Number of bins in histograms')
        group.addParam('usePsi', BooleanParam, default=False,
                      label='Use Psi to compute angular distances?',
                      help='Use Psi')
        group.addParam('angleSort', BooleanParam, default=False,
                      label='Sort assigned angles?',
                      help='Sort by angles the experimental images.')
        group.addParam('shiftSort', BooleanParam, default=False,
                      label='Sort shift?',
                      help='Sort by shift the experimental images.')
        group = form.addGroup('Angular distribution and resolution plots')
        group.addParam('showAngDist', EnumParam, choices=['2D plot', 'chimera'], 
                      display=EnumParam.DISPLAY_LIST, default=ANGDIST_2DPLOT,
                      label='Display angular distribution',
                      help='*2D plot*: display angular distribution as interative 2D in matplotlib.\n'
                           '*chimera*: display angular distribution using Chimera with red spheres.') 
        group.addParam('showResolutionPlots', BooleanParam, default=True,
                      label='Display resolution plots (FSC)?',
                      help='')
        group.addParam('resolutionThreshold', FloatParam, default=0.5, 
                      expertLevel=LEVEL_ADVANCED,
                      label='Threshold in resolution plots',
                      help='')                                      

    def _getVisualizeDict(self):
        self._load()
        return {'showReference': self._showRefs,
                'showReconstruction': self._showRecons,
                'showFilteredReconstruction' : self._showFilVols,
                'showBFactorCorrectedVolume' : self._showBfactorVols,
                'showProjectionMatchingLibrary' : self._showProjMatchLibrary,
                'showProjectionMatchingClasses' : self._showProjMatchClasses,
                'showProjectionMatchingLibraryAndClasses' : self._showProjMatchLibAndClasses,
                'showProjectionMatchingLibraryAndImages' : self._showProjMatchLibAndImages,
                'showDiscardedImages' : self._showDiscardedImages,
                'showExperimentalImages' : self._showExperimentalImages,
                'plotHistogramAngularMovement' : self._plotHistogramAngularMovement,
                'showAngDist': self._showAngularDistribution,
                'showResolutionPlots': self._showFSC
                }
    
    def _viewAll(self, *args):
        pass
    
    def _validate(self):
        if self.lastIter is None:
            return ['There are not iterations completed.'] 
    
    def createDataView(self, filename, extraParams=''):
        return DataView(filename)
        
#     def createScipionView(self, filename, extraParams=''):
#         inputParticlesId = self.protocol.inputParticles.get().strId()
#         return Classes3DView(self._project.getName(), 
#                           self.protocol.strId(), filename, other=inputParticlesId,
#                           env=self._env)

    def _load(self):
        """ Load selected iterations and classes 3D for visualization mode. """
        self._refsList = [1]
        self.protocol._initialize() # Load filename templates
        if self.showRef3DNo == REF_ALL:
            self._refsList = range(1, self.protocol.numberOfReferences+1)
        else:
            self._refsList = self._getListFromRangeString(self.ref3DSelection.get())
        # ToDo: enhance this
        self.firstIter = 1
        self.lastIter = self.protocol.numberOfIterations.get()
        
        if self.viewIter.get() == ITER_LAST:
            self._iterations = [self.lastIter]
        else:
            self._iterations = self._getListFromRangeString(self.iterSelection.get())
            
        from matplotlib.ticker import FuncFormatter
        self._plotFormatter = FuncFormatter(self._formatFreq) 
    
    def _formatFreq(self, value, pos):
        """ Format function for Matplotlib formatter. """
        inv = 999
        if value:
            inv = int(1/value)
        return "1/%d" % inv

    def _getGridSize(self, n=None):
        """ Figure out the layout of the plots given the number of references. """
        if n is None:
            n = len(self._refsList)
        
        if n == 1:
            gridsize = [1, 1]
        elif n == 2:
            gridsize = [2, 1]
        else:
            gridsize = [(n+1)/2, 2]
            
        return gridsize
    
#     def _getPrefixes(self):
#         prefixes = self.protocol.PREFIXES
#         halves = getattr(self, 'showHalves', None)
#         if halves:
#             if halves == 0:
#                 prefixes = ['half1_']
#             elif halves == 1:
#                 prefixes = ['half2_']
#         return prefixes
    
#===============================================================================
# Show Volumes
#===============================================================================
    
    def _createVolumesMd(self, volumes):
        """ Write a metadata with all volumes selected for visualization. """
        mdPath = self.protocol._getTmpPath('viewer_volumes.xmd')
        cleanPath(mdPath)
        md = xmipp.MetaData()
        
        for volFn in volumes:
            md.clear()
            md.setValue(xmipp.MDL_IMAGE, volFn, md.addObject())
            blockName = volFn.split("/")[3]
            print "Volume: ", volFn, blockName
            md.write("%s@%s"% (blockName, mdPath), xmipp.MD_APPEND)
        return [self.createDataView(mdPath)]
    
    def _showVolumesChimera(self, volumes):
        """ Create a chimera script to visualize selected volumes. """
        
        if len(volumes) > 1:
            cmdFile = self.protocol._getTmpPath('chimera_volumes.cmd')
            cleanPath(cmdFile)
            f = open(cmdFile, 'w+')
            f.write('windowsize 800 600\n')
            for volFn in volumes:
                vol = os.path.relpath(volFn, self.protocol._getTmpPath())
                f.write("open %s\n" % vol)
            f.write('tile\n')
            f.close()
            view = CommandView('chimera %s &' % cmdFile)                    
        else:
            view = CommandView('xmipp_chimera_client --input "%s" --mode projector 256 &' % volumes[0])
        
        return [view]
    
    def _showVolumes(self, volumes, paramName):
        if self.showVolume == VOLUME_CHIMERA:
            return self._showVolumesChimera(volumes)
        elif self.showVolume == VOLUME_SLICES:
            return self._createVolumesMd(volumes)
    
    def _volFileNames(self, volTemplate):
        volumes = []
        for it in self._iterations:
            for ref3d in self._refsList:
                volFn = self.protocol._getFileName(volTemplate, iter=it, ref=ref3d)
                volumes.append(volFn)
        return volumes
    
    def _showRefs(self, paramName=None):
        volumes = self._volFileNames('maskedFileNamesIters')
        return self._showVolumes(volumes, paramName)
    
    def _showRecons(self, paramName=None):
        volumes = self._volFileNames('reconstructedFileNamesIters')
        return self._showVolumes(volumes, paramName)
    
    def _showFilVols(self, paramName=None):
        volumes = self._volFileNames('reconstructedFilteredFileNamesIters')
        return self._showVolumes(volumes, paramName)
    
    def _showBfactorVols(self, paramName=None):
        volsBfactor = []
        volumes = self._volFileNames('reconstructedFileNamesIters')
        for vol in volumes:
            volBfactor = vol + '.bfactor'
            volsBfactor.append(volBfactor)
            args = '-i %(vol)s --sampling %(samplingRate)s --maxres %(maxRes)s -o %(volBfactor)s ' + self.correctBfactorExtraCommand.get()
            maxRes = self.maxRes.get()
            samplingRate = self.protocol.resolSam
            args = args % locals()
            
            hostConfig = self.protocol.getHostConfig()
            # Create the steps executor
            executor = StepExecutor(hostConfig)
            self.protocol.setStepsExecutor(executor)
            # Finally run the protocol
            
            self.protocol.runJob("xmipp_volume_correct_bfactor", args, numberOfMpi=1)
        return self._showVolumes(volsBfactor, paramName)
#         pass
    
#===============================================================================
# Show Library, Classes and Images
#===============================================================================

    def _showProjMatchLibrary(self, paramName=None):
        #map stack position with ref number
        list = []
        mdIn  = xmipp.MetaData()
        mdOut = xmipp.MetaData()
        cleanPattern(self.protocol._getTmpPath('references_library*'))
        
        for ref3d in self._refsList:
            for it in self._iterations:
                convert_refno_to_stack_position = {}
                file_nameReferences = 'projectionDirections@' + self.protocol._getFileName('projectLibrarySampling', iter=it, ref=ref3d)
                #last reference name
                mdReferences     = xmipp.MetaData(file_nameReferences)
                mdReferencesSize = mdReferences.size()
                for id in mdReferences:
                    convert_refno_to_stack_position[mdReferences.getValue(xmipp.MDL_NEIGHBOR,id)]=id
                file_nameAverages   = self.protocol._getFileName('outClassesXmd', iter=it, ref=ref3d)
                if exists(file_nameAverages):
                    #print "OutClassesXmd", OutClassesXmd
                    mdIn.read(file_nameAverages)
                    mdOut.clear()
                    for i in mdIn:
                        #id1=mdOut.addObject()
                        #mdOut.setValue(MDL_IMAGE,mdIn.getValue(MDL_IMAGE,i),id1)
                        ref2D = mdIn.getValue(xmipp.MDL_REF,i)
                        file_references = self.protocol._getFileName('projectLibraryStk', iter=it, ref=ref3d)
                        file_reference = xmipp.FileName()
                        file_reference.compose(convert_refno_to_stack_position[ref2D],file_references)
                        id2=mdOut.addObject()
                        mdOut.setValue(xmipp.MDL_IMAGE, file_reference, id2)
                        
                    if mdOut.size() == 0:
                        print "Empty metadata: ", file_nameReferences
                    else:
                        file_nameReferences = self.protocol._getTmpPath('references_library.xmd')
                        sfn = createUniqueFileName(file_nameReferences)
                        file_nameReferences = 'projectionDirections@' + sfn
                        mdOut.write(file_nameReferences)
                        list.append(self.createDataView(file_nameReferences))
        return list
    
    def _showProjMatchClasses(self, paramName=None):
        classes = []
        for ref3d in self._refsList:
            for it in self._iterations:
                classesFn = self.protocol._getFileName('outClassesXmd', iter=it, ref=ref3d)
                if exists(classesFn):
                    _, _, _, _, size = xmipp.MetaDataInfo(classesFn)
                    
                    if size == 0:
                        print "Empty metadata: ", classesFn
                    else:
                        classes.append(self.createDataView(classesFn))
        return classes
    
    def _showProjMatchLibAndClasses(self, paramName=None):
        #map stack position with ref number
        list = []
        mdIn  = xmipp.MetaData()
        mdOut = xmipp.MetaData()
        
        for ref3d in self._refsList:
            for it in self._iterations:
                convert_refno_to_stack_position = {}
                file_nameReferences = 'projectionDirections@' + self.protocol._getFileName('projectLibrarySampling', iter=it, ref=ref3d)
                mdReferences     = xmipp.MetaData(file_nameReferences)
                mdReferencesSize = mdReferences.size()
                for id in mdReferences:
                    convert_refno_to_stack_position[mdReferences.getValue(xmipp.MDL_NEIGHBOR,id)]=id
                file_nameAverages   = self.protocol._getFileName('outClassesXmd', iter=it, ref=ref3d)
                file_references = self.protocol._getFileName('projectLibraryStk', iter=it, ref=ref3d)
                if exists(file_nameAverages):
                    mdIn.read(file_nameAverages)
                    mdOut.clear()
                    for i in mdIn:
                        ref2D = mdIn.getValue(xmipp.MDL_REF, i)
                        file_reference = xmipp.FileName()
                        file_reference.compose(convert_refno_to_stack_position[ref2D],file_references)
                        id1 = mdOut.addObject()
                        mdOut.setValue(xmipp.MDL_IMAGE, mdIn.getValue(xmipp.MDL_IMAGE,i), id1)
                        mdOut.setValue(xmipp.MDL_IMAGE2, file_reference, id1)
                        
                    if mdOut.size() == 0:
                        print "Empty metadata: ", file_nameReferences
                    else:
                        file_nameReferences = self.protocol._getFileName('projectLibrarySampling', iter=it, ref=ref3d)
                        sfn = createUniqueFileName(file_nameReferences)
                        file_nameReferences = 'projectionDirections@' + sfn
                        mdOut.merge(mdIn)
                        mdOut.write(file_nameReferences)
                        # ToDo: show the metadata in "metadata" form.
                        list.append(self.createDataView(file_nameReferences))
        return list
    
    def _showProjMatchLibAndImages(self, paramName=None):
        from numpy  import array, dot
        #map stack position with ref number
        imgAndClasses = []
        mdIn  = xmipp.MetaData()
        mdOut = xmipp.MetaData()
        mdTmp = xmipp.MetaData()
        for ref3d in self._refsList:
            for it in self._iterations:
                convert_refno_to_stack_position = {}
                file_nameReferences = 'projectionDirections@' + self.protocol._getFileName('projectLibrarySampling', iter=it, ref=ref3d)
                #last reference name
                mdReferences = xmipp.MetaData(file_nameReferences)
                mdReferencesSize = mdReferences.size()
                for id in mdReferences:
                    convert_refno_to_stack_position[mdReferences.getValue(xmipp.MDL_NEIGHBOR,id)]=id
                file_nameImages = "ctfGroup[0-9][0-9][0-9][0-9][0-9][0-9]@" + self.protocol._getFileName('docfileInputAnglesIters', iter=it)
                mdTmp.read(file_nameImages)#query with ref3D
                mdIn.importObjects(mdTmp, xmipp.MDValueEQ(xmipp.MDL_REF3D, ref3d))
                mdOut.clear()
                for i in mdIn:
                    id1 = mdOut.addObject()
                    mdOut.setValue(xmipp.MDL_IMAGE, mdIn.getValue(xmipp.MDL_IMAGE,i), id1)
                    psi = -1. * mdIn.getValue(xmipp.MDL_ANGLE_PSI, i)
                    flip = mdIn.getValue(xmipp.MDL_FLIP, i)
                    if(flip):
                        psi = -psi
                    eulerMatrix = xmipp.Euler_angles2matrix(0., 0., psi)
                    x = mdIn.getValue(xmipp.MDL_SHIFT_X, i)
                    y = mdIn.getValue(xmipp.MDL_SHIFT_Y, i)
                    shift = array([x, y, 0])
                    shiftOut = dot(eulerMatrix, shift)
                    [x,y,z] = shiftOut
                    if flip:
                        x = -x
                    mdOut.setValue(xmipp.MDL_ANGLE_PSI, psi, id1)
                    mdOut.setValue(xmipp.MDL_SHIFT_X, x, id1)
                    mdOut.setValue(xmipp.MDL_SHIFT_Y, y, id1)
                    mdOut.setValue(xmipp.MDL_FLIP, flip, id1)
                    ref2D = mdIn.getValue(xmipp.MDL_REF,i)
                    file_references = self.protocol._getFileName('projectLibraryStk', iter=it, ref=ref3d)
                    file_reference = xmipp.FileName()
                    file_reference.compose(convert_refno_to_stack_position[ref2D], file_references)
                    id2 = mdOut.addObject()
                    mdOut.setValue(xmipp.MDL_IMAGE, file_reference, id2)
                    mdOut.setValue(xmipp.MDL_ANGLE_PSI, 0., id2)
                if mdOut.size() == 0:
                    print "Empty metadata"
                else:
                    file_nameReferences = self.protocol._getFileName('projectLibrarySampling', iter=it, ref=ref3d)
                    sfn   = createUniqueFileName(file_nameReferences)
                    file_nameReferences = 'projectionDirections@' + sfn
                    mdOut.write(sfn)
                    imgAndClasses.append(self.createDataView(file_nameReferences))
        return imgAndClasses
    
    def _showDiscardedImages(self, paramName=None):
        md = xmipp.MetaData()
        for it in self._iterations:
            file_name = self.protocol._getFileName('outClassesDiscarded', iter=it)
            md.read(file_name)
            if exists(file_name):
                if md.size() == 0:
                    print "Empty metadata: ", file_name
                else:
                    return [self.createDataView(file_name)]
    
    def _showExperimentalImages(self, paramName=None):
        md = xmipp.MetaData()
        for ref3d in self._refsList:
            for it in self._iterations:
                file_name = "ctfGroup[0-9][0-9][0-9][0-9][0-9][0-9]@" + self.protocol._getFileName('docfileInputAnglesIters', iter=it)
                md.read(file_name)
                return [self.createDataView(file_name)]
    
#===============================================================================
# Convergence
#===============================================================================

    def _plotHistogramAngularMovement(self, paramName=None):
        from numpy import arange
        from matplotlib.ticker import FormatStrFormatter
        
        plots = []
        colors = ['g', 'b', 'r', 'y', 'c', 'm', 'k']
        lenColors=len(colors)
        
        numberOfBins = self.numberOfBins.get()
        md = xmipp.MetaData()
        for it in self._iterations:
            md.read(self.protocol._mdDevitationsFn(it))
            if not self.usePsi: 
                md.fillConstant(xmipp.MDL_ANGLE_PSI,0.)
            
            nrefs = len(self._refsList)
            gridsize = self._getGridSize(nrefs)
            xplotterShift = XmippPlotter(*gridsize, mainTitle='Iteration_%d\n' % it, windowTitle="ShiftDistribution")
            xplotter = XmippPlotter(*gridsize, mainTitle='Iteration_%d' % it, windowTitle="AngularDistribution")
            
            for ref3d in self._refsList:
                mDoutRef3D = xmipp.MetaData()
                mDoutRef3D.importObjects(md, xmipp.MDValueEQ(xmipp.MDL_REF3D, ref3d))
                _frequency = "Frequency (%d)" % mDoutRef3D.size()
                
                xplotterShift.createSubPlot("%s_ref3D_%d"%(xmipp.label2Str(xmipp.MDL_SHIFT_DIFF),ref3d), "pixels", _frequency) 
                xplotter.createSubPlot("%s_ref3D_%d"%(xmipp.label2Str(xmipp.MDL_ANGLE_DIFF),ref3d), "degrees", _frequency) 
                #mDoutRef3D.write("afterimportObject@mdIter.sqlite",MD_APPEND)
                xplotter.plotMd(mDoutRef3D, 
                                xmipp.MDL_ANGLE_DIFF, 
                                xmipp.MDL_ANGLE_DIFF, 
                                color=colors[ref3d%lenColors], 
                                #nbins=50
                                nbins=int(numberOfBins)
                                )#if nbins is present do an histogram
                xplotterShift.plotMd(mDoutRef3D, 
                                xmipp.MDL_SHIFT_DIFF, 
                                xmipp.MDL_SHIFT_DIFF, 
                                color=colors[ref3d%lenColors], 
                                nbins=int(numberOfBins)
                                )#if nbins is present do an histogram
                 
                if self.angleSort:
                    mDoutRef3D.sort(xmipp.MDL_ANGLE_DIFF)
                    fn = xmipp.FileName()
                    baseFileName   = self.protocol._getTmpPath("angle_sort.xmd")
                    fn = self.protocol._getRefBlockFileName("angle_iter", it, "ref3D", ref3d, baseFileName)
                    mDoutRef3D.write(fn, xmipp.MD_APPEND)
                    print "File with sorted angles saved in:", fn
                     
                if self.shiftSort:
                    mDoutRef3D.sort(xmipp.MDL_SHIFT_DIFF)
                    fn = xmipp.FileName()
                    baseFileName   = self.protocol._getTmpPath("angle_sort.xmd")
                    fn = self.protocol._getRefBlockFileName("shift_iter", it, "ref3D", ref3d, baseFileName)
                    mDoutRef3D.write(fn, xmipp.MD_APPEND)
                    print "File with sorted shifts saved in:", fn
                
                plots.append(xplotterShift)
                plots.append(xplotter)
        return plots
    
#===============================================================================
# Angular distribution and resolution plots
#===============================================================================
    
    def _createAngDistChimera(self, it):
        arguments = []
        outerRadius = self.protocol._outerRadius[it]
        radius = float(outerRadius) * 1.1

        for ref3d in self._refsList:
            file_name = self.protocol._getFileName('outClassesXmd', iter=it, ref=ref3d)
            file_name_rec_filt = self.protocol._getFileName('reconstructedFilteredFileNamesIters', iter=it, ref=ref3d)
            args = "--input '%s' --mode projector 256 -a %s red %f" %(file_name_rec_filt, file_name, radius)
            arguments.append(args)
                   
        return CommandView('xmipp_chimera_client %s &' % args)
    
    def _createAngDist2D(self, it):
        # Common variables to use
        nrefs = len(self._refsList)
        gridsize = self._getGridSize(nrefs)
        xplotter = XmippPlotter(*gridsize, mainTitle='Iteration %d' % it, windowTitle="Angular Distribution")
        for ref3d in self._refsList:
            md = xmipp.MetaData(self.protocol._getFileName('outClassesXmd', iter=it, ref=ref3d))
            plot_title = 'Ref3D_%d' % ref3d
            xplotter.plotMdAngularDistribution(plot_title, md)
        
        return xplotter
    
    def _plotFSC(self, a, mdFn):
        md = xmipp.MetaData(mdFn)
        resolution_inv = [md.getValue(xmipp.MDL_RESOLUTION_FREQ, id) for id in md]
        frc = [md.getValue(xmipp.MDL_RESOLUTION_FRC, id) for id in md]
        self.maxFrc = max(frc)
        self.minInv = min(resolution_inv)
        self.maxInv = max(resolution_inv)
        a.plot(resolution_inv, frc)
        a.xaxis.set_major_formatter(self._plotFormatter)
        a.set_ylim([-0.1, 1.1])
    
    def _showAngularDistribution(self, paramName=None):
        views = []
        
        if self.showAngDist == ANGDIST_CHIMERA:
            for it in self._iterations:
                views.append(self._createAngDistChimera(it))
                        
        elif self.showAngDist == ANGDIST_2DPLOT:
            for it in self._iterations:
                views.append(self._createAngDist2D(it))
                
        return views
    
    def _showFSC(self, paramName=None):
        threshold = self.resolutionThreshold.get()
        nrefs = len(self._refsList)
        gridsize = self._getGridSize(nrefs)
        
        xmipp.activateMathExtensions()
        
        xplotter = XmippPlotter(*gridsize, windowTitle='Resolution FSC')
        
        for ref3d in self._refsList:
            plot_title = 'Ref3D_%s' % ref3d
            a = xplotter.createSubPlot(plot_title, 'Armstrongs^-1', 'FSC', yformat=False)
            legends = []
            for it in self._iterations:
                file_name = self.protocol._getFileName('resolutionXmdFile', iter=it, ref=ref3d)
                if exists(file_name):
                    self._plotFSC(a, file_name)
                    legends.append('iter %d' % it)
            xplotter.showLegend(legends)
            if threshold < self.maxFrc:
                a.plot([self.minInv, self.maxInv],[threshold, threshold], color='black', linestyle='--')
            a.grid(True)
            
        return [xplotter]
