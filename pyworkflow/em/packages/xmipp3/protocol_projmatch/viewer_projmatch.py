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


from pyworkflow.protocol.executor import StepExecutor
from pyworkflow.viewer import CommandView, Viewer, ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO
from pyworkflow.em.viewer import DataView, ClassesView, Classes3DView
from pyworkflow.utils import createUniqueFileName, cleanPattern
from protocol_projmatch import XmippProtProjMatch
# from projmatch_initialize import createFilenameTemplates
from pyworkflow.em.packages.xmipp3.convert import * # change this
from pyworkflow.em.viewer import ChimeraDataView
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.protocol.params import (LabelParam, IntParam, FloatParam,
                                        StringParam, EnumParam, NumericRangeParam)
import xmipp
from pyworkflow.em.plotter import EmPlotter
from pyworkflow.em.packages.xmipp3.plotter import XmippPlotter
import numpy as np


ITER_LAST = 0
ITER_SELECTION = 1

ANGDIST_2DPLOT = 0
ANGDIST_CHIMERA = 1

VOLUME_SLICES = 0
VOLUME_CHIMERA = 1

DISPLAY_LIBRARY = 0

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
                      display=EnumParam.DISPLAY_HLIST,
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
                      label="Iteration list",
                      help="Write the iteration list to visualize.")

        group = form.addGroup('Particles')
        group.addParam('displayLibraryOrClasses', EnumParam, choices=['projections', 'classes', 'projections and classes'],
                          default=DISPLAY_LIBRARY, display=EnumParam.DISPLAY_COMBO,
                          label='Display',
                          help='Displays images with angular assignment')
        group.addParam('showProjectionMatchingLibraryAndImages', LabelParam, default=False,
                      label='Display projections and particles',
                      help="Display projections and particles")
        group.addParam('showExperimentalImages', LabelParam, default=False,
                      label='Display particles',
                      help="""Display particles with alignment and classification information
                           WARNING: the angles and shifts are the adequate for reconstruction
                           but not for 2D aligment.
                           """)
        group.addParam('showDiscardedImages', LabelParam, default=False,
                      label='Display discarded particles',
                      help='Display discarded particles.')
        
        group = form.addGroup('Volumes')
        group.addParam('showRef3DNo', EnumParam, choices=['all', 'selection'], default=REF_ALL,
                      display=EnumParam.DISPLAY_HLIST,
                      label='3D Class to visualize',
                      help='All: Display all 3D classes for each iteration'
                           'that you selected.\n'
                           'Selection: You may specify which 3D class (or classes)'
                           ' to visualize')
        group.addParam('ref3DSelection', NumericRangeParam, default='1',
                      condition='showRef3DNo == %d' % REF_SEL,
                      label='Classes list',
                      help='Write the 3d classes list to visualize.')
        group.addParam('matrixWidth', FloatParam, default=-1, 
                      expertLevel=LEVEL_ADVANCED,
                      label='Width of projection galleries',
                      help='Usually a multiple of 2 is the right value. -1 => authomatic')

        group.addParam('showAngDist', EnumParam, choices=['2D plot', 'chimera'],
                      display=EnumParam.DISPLAY_HLIST, default=ANGDIST_2DPLOT,
                      label='Display angular distribution',
                      help='*2D plot*: display angular distribution as interative 2D in matplotlib.\n'
                           '*chimera*: display angular distribution using Chimera with red spheres.')
        group.addParam('displayVolWith', EnumParam, choices=['slices', 'chimera'],
                      display=EnumParam.DISPLAY_HLIST, default=VOLUME_SLICES,
                      label='Display volume with',
                      help='*slices*: display volumes as 2D slices along z axis.\n'
                           '*chimera*: display volumes as surface with Chimera.')
        group.addParam('displayVolume', EnumParam, choices=['Reference', 'Reconstructed', 'Filtered', 'bfactor corrected'],
                          default=1, display=EnumParam.DISPLAY_COMBO,
                          label='Display volume',
                          help='Displays selected volume')
        group.addParam('spheresScale', IntParam, default=100, 
              expertLevel=LEVEL_ADVANCED,
              label='Spheres size',
              help='')
        group.addParam('maxRes', FloatParam, default=5, 
                      condition='displayVolume==3',
                      label='Maximum resolution to apply B-factor (in Angstroms)')
        group.addParam('correctBfactorExtraCommand', StringParam, default='--auto',
                       condition='displayVolume==3', 
                       label='User defined flags for the correct_bfactor program',
                       help=""" See http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/Correct_bfactor
                            for details. DEFAULT behaviour is --auto
                            """)

        group = form.addGroup('Resolution')

        group.addParam('showResolutionPlots', LabelParam, default=True,
                      label='Display resolution plots (FSC)',
                      help='')
        group.addParam('resolutionThreshold', FloatParam, default=0.5,
                      expertLevel=LEVEL_ADVANCED,
                      label='Threshold in resolution plots',
                      help='')

        form.addSection(label='Convergence')
        group = form.addGroup('Convergence')
        group.addParam('plotHistogramAngularMovement', LabelParam, default=False,
                      label='Plot histogram with angular/shift changes',
                      help=""" Plot histogram with angular changes from one iteration to next. 
                           Iteration 0 -> initial values
                           """)
        group.addParam('numberOfBins', IntParam, default=50, 
                      condition='plotHistogramAngularMovement',
                      label='Number of bins (for histogram)',
                      help='Number of bins in histograms')
        group.addParam('usePsi', BooleanParam, default=False,
                      label='Use Psi to compute angular distances',
                      help='Use Psi')
        group.addParam('angleSort', BooleanParam, default=False,
                      label='Sort assigned angles',
                      help='Sort by angles the experimental images.')
        group.addParam('shiftSort', BooleanParam, default=False,
                      label='Sort shift',
                      help='Sort by shift the experimental images.')
    
    def _getVisualizeDict(self):
        self._load()
        return {
                'displayVolume' : self._showVolume,
                'displayLibraryOrClasses' : self._showLibraryOrClasses,
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
    
    def createDataView(self, filename, viewParams={}):
        return DataView(filename, viewParams)
        
#     def createScipionView(self, filename, extraParams=''):
#         inputParticlesId = self.protocol.inputParticles.get().strId()
#         return Classes3DView(self._project,
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
        #self.lastIter = self.protocol.numberOfIterations.get()
        self.lastIter = self.protocol.getLastIter()
        
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
            inv = 1/value
        return "1/%0.2f" % inv

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
    
#===============================================================================
# Show Volumes
#===============================================================================

    def _showVolume(self, paramName=None):
        choice = self.displayVolume.get()
        if choice == 0:
            return self._showRefs(paramName)
        elif choice == 1:
            return self._showRecons(paramName)
        elif choice == 2:
            return self._showFilVols(paramName)
        else:
            return self._showBfactorVols(paramName)
            
    
    def _createVolumesMd(self, volumes):
        """ Write a metadata with all volumes selected for visualization. """
        mdPath = self.protocol._getTmpPath('viewer_volumes.xmd')
        cleanPath(mdPath)
        md = xmipp.MetaData()
        
        for volFn in volumes:
            md.clear()
            md.setValue(xmipp.MDL_IMAGE, volFn, md.addObject())
            blockName = volFn.split("/")[3]
            #print "Volume: ", volFn, blockName
            md.write("%s@%s"% (blockName, mdPath), xmipp.MD_APPEND)
        return [self.createDataView(mdPath)]

    def viewVolumesSqlite(self, volumes):
        path = self.protocol._getExtraPath('viewer_volumes.sqlite')
        samplingRate = self.protocol.inputParticles.get().getSamplingRate()
        self.createVolumesSqlite(volumes, path, samplingRate)
        return [ObjectView(self._project, self.protocol.strId(), path)]
    
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
            view = ChimeraView(cmdFile)
        else:
            
            #view = CommandView('xmipp_chimera_client --input "%s" --mode projector 256 &' % volumes[0])
            view = ChimeraClientView(volumes[0], showProjection=True)
        
        return [view]
    
    def _showVolumes(self, volumes):
        if self.displayVolWith == VOLUME_CHIMERA:
            return self._showVolumesChimera(volumes)
        elif self.displayVolWith == VOLUME_SLICES:
            #return self._createVolumesMd(volumes)
            return self.viewVolumesSqlite(volumes)
    
    def _volFileNames(self, volTemplate):
        volumes = []
        for it in self._iterations:
            for ref3d in self._refsList:
                volFn = self.protocol._getFileName(volTemplate, iter=it, ref=ref3d)
                if exists(volFn):
                    volumes.append(volFn)
                else:
                    print "Volume %s does not exist" % volFn
        return volumes
    
    def _showRefs(self, paramName=None):
        volumes = self._volFileNames('maskedFileNamesIters')
        return self._showVolumes(volumes)
    
    def _showRecons(self, paramName=None):
        volumes = self._volFileNames('reconstructedFileNamesIters')
        return self._showVolumes(volumes)
    
    def _showFilVols(self, paramName=None):
        volumes = self._volFileNames('reconstructedFilteredFileNamesIters')
        return self._showVolumes(volumes)
    
    def _showBfactorVols(self, paramName=None):
        volsBfactor = []
        guinierPlots = []
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
            guinierPlots.append(self._showGuinier(volBfactor))
            
        return self._showVolumes(volsBfactor) + guinierPlots
    
#===============================================================================
# Show Library, Classes and Images
#===============================================================================
    def _showLibraryOrClasses(self, paramName=None):
        choice = self.displayLibraryOrClasses.get()
        if choice == DISPLAY_LIBRARY:
            return self._showProjMatchLibrary(paramName)
        elif choice == 1:
            return self._showProjMatchClasses(paramName)
        else:
            return self._showProjMatchLibAndClasses(paramName)

    def _showProjMatchLibrary(self, paramName=None):
        #map stack position with ref number
        list = []
        mdIn  = xmipp.MetaData()
        mdOut = xmipp.MetaData()
        cleanPattern(self.protocol._getTmpPath('references_library*'))
        
        for ref3d in self._refsList:
            for it in self._iterations:
                convert_refno_to_stack_position = {}
                file_name = self.protocol._getFileName('projectLibrarySampling', iter=it, ref=ref3d)
                file_nameReferences = 'projectionDirections@' + file_name
                #last reference name
                if exists(file_name):
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
                else:
                    print "File %s does not exist" % file_name
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
                else:
                    print "File %s does not exist" % classesFn
        return classes
    
    def _showProjMatchLibAndClasses(self, paramName=None):
        #map stack position with ref number
        list = []
        mdIn  = xmipp.MetaData()
        mdOut = xmipp.MetaData()
        
        for ref3d in self._refsList:
            for it in self._iterations:
                convert_refno_to_stack_position = {}
                file_name = self.protocol._getFileName('projectLibrarySampling', iter=it, ref=ref3d)
                file_nameReferences = 'projectionDirections@' + file_name
                if exists(file_name):
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
                    else:
                        print "File %s does not exist" % file_name
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
                file_name = self.protocol._getFileName('projectLibrarySampling', iter=it, ref=ref3d)
                file_nameReferences = 'projectionDirections@' + file_name
                if exists(file_name):
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
                        mdOut.write(file_nameReferences)
                        imgAndClasses.append(self.createDataView(file_nameReferences))
                else:
                        print "File %s does not exist" % file_name
        return imgAndClasses
    
    def _showDiscardedImages(self, paramName=None):
        md = xmipp.MetaData()
        for it in self._iterations:
            file_name = self.protocol._getFileName('outClassesDiscarded', iter=it)
            if exists(file_name):
                md.read(file_name)
                if md.size() == 0:
                    print "Empty metadata: ", file_name
                else:
                    return [self.createDataView(file_name)]
            else:
                print "File %s does not exist" % file_name
                return []
    
    def _showExperimentalImages(self, paramName=None):
        views = []
        for ref3d in self._refsList:
            for it in self._iterations:
                partSet = self.protocol._getIterParticles(it)
                v = self.createScipionPartView(partSet)
                views.append(v)
        return views
    
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
            mdFn = self.protocol._mdDevitationsFn(it)
            if xmipp.existsBlockInMetaDataFile(mdFn):
                md.read(mdFn)
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
            else:
                print "File %s does not exist" % mdFn
        return plots
    
#===============================================================================
# showAngularDistribution
#===============================================================================
    def _showAngularDistribution(self, paramName=None):
        views = []
        if self.showAngDist == ANGDIST_CHIMERA:
            for it in self._iterations:
                angDist = self._createAngDistChimera(it)
                if angDist is not None:
                    views.append(angDist)
                        
        elif self.showAngDist == ANGDIST_2DPLOT:
            for it in self._iterations:
                angDist = self._createAngDist2D(it)
                if angDist is not None:
                    views.append(angDist)
        return views
    
    def _createAngDistChimera(self, it):
        radius = self.spheresScale.get()

        if len(self._refsList) == 1:
            ref3d = self._refsList[0]
            classesFn = self.protocol._getFileName('outClassesXmd', iter=it, ref=ref3d)
            vol = self.protocol._getFileName('reconstructedFilteredFileNamesIters', iter=it, ref=ref3d)
            if exists(classesFn):
                return ChimeraClientView(vol, showProjection=True, angularDistFile=classesFn, spheresDistance=radius)
            else:
                print "File %s does not exist" % classesFn
                return None
        else:
            return self.infoMessage('Please select only one class to display angular distribution')
    
    def _createAngDist2D(self, it):
        # Common variables to use
        nrefs = len(self._refsList)
        gridsize = self._getGridSize(nrefs)
        xplotter = XmippPlotter(*gridsize, mainTitle='Iteration %d' % it, windowTitle="Angular Distribution")
        for ref3d in self._refsList:
            classesFn = self.protocol._getFileName('outClassesXmd', iter=it, ref=ref3d)
            if exists(classesFn):
                md = xmipp.MetaData(classesFn)
                title = 'Ref3D_%d' % ref3d
                xplotter.plotMdAngularDistribution(title, md)
            else:
                print "File %s does not exist" % classesFn
                return None
        
        return xplotter
    
#===============================================================================
# plotFSC
#===============================================================================
    def _showFSC(self, paramName=None):
        threshold = self.resolutionThreshold.get()
        nrefs = len(self._refsList)
        gridsize = self._getGridSize(nrefs)
        xmipp.activateMathExtensions()
        
        for ref3d in self._refsList:
            xplotter = XmippPlotter(*gridsize, windowTitle='Resolution FSC')
            legends = []
            show = False
            plot_title = 'Ref3D_%s' % ref3d
            a = xplotter.createSubPlot(plot_title, 'frequency(1/A)', 'FSC', yformat=False)
            legends = []
            for it in self._iterations:
                file_name = self.protocol._getFileName('resolutionXmdFile', iter=it, ref=ref3d)
                if exists(file_name):
                    show = True
                    legends.append('iter %d' % it)
                    self._plotFSC(a, file_name)
                    xplotter.showLegend(legends)
            if show:
                if threshold < self.maxFrc:
                    a.plot([self.minInv, self.maxInv],[threshold, threshold], color='black', linestyle='--')
                a.grid(True)
            else:
                raise Exception("Set a valid iteration to show its FSC")
            
            return [xplotter]
    
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
    
#===============================================================================
# plotBfactor
#===============================================================================
    def _showGuinier(self, volume):
        nrefs = len(self._refsList)
        gridsize = self._getGridSize(nrefs)
        guinierFn = volume + ".guinier"
        
        d2 = self._getGuinierValue(guinierFn, 0)
        
        legends = ["lnFweighted ln(F)", "corrected ln(F)", "model"]
        xplotter = EmPlotter(*gridsize, windowTitle='Guinier Plots')
        subPlot = xplotter.createSubPlot(basename(volume), 'd^-2(A^-2)', 'ln(F)', yformat=False)
        for i, legend in enumerate(legends):
            y = self._getGuinierValue(guinierFn, i+2)
            subPlot.plot(d2, y)
            xplotter.showLegend(legends)
        subPlot.grid(True)
        return xplotter
    
    def _getGuinierValue(self, guinierFn, col):
        f1 = open(guinierFn)
        value = []
        for l in f1:
            if not "#" in l:
                valList = l.split()
                val = float(valList[col])
                value.append(val)
        f1.close()
        return value

#===============================================================================
# Utils Functions
#===============================================================================
    def createScipionPartView(self, partSet, viewParams={}):
        from pyworkflow.em import ObjectView
        inputParticlesId = self.protocol.inputParticles.get().strId()
        filename = partSet.getFileName()
        
        labels =  'enabled id _size _filename _transform._matrix'
        viewParams = {showj.ORDER:labels,
                      showj.VISIBLE: labels, showj.RENDER:'_filename',
                      'labels': 'id',
                      }
        
        return ObjectView(self._project, 
                          self.protocol.strId(), filename, other=inputParticlesId,
#                           env=self._env,
                          viewParams=viewParams)
