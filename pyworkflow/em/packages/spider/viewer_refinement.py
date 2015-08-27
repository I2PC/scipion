# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
from pyworkflow.em.packages.eman2.viewer import ANGDIST_CHIMERA
"""
This module implement the wrappers around xmipp_showj
visualization program.
"""

import os
from glob import glob

import pyworkflow.em as em
import pyworkflow.protocol.params as params
from pyworkflow.viewer import DESKTOP_TKINTER, WEB_DJANGO, ProtocolViewer
from pyworkflow.em.plotter import EmPlotter

from protocol import SpiderProtRefinement
from spider import SpiderDocFile


ITER_LAST = 0
ITER_SELECTION = 1

ANGDIST_2DPLOT = 0
ANGDIST_CHIMERA = 1

VOLUME_SLICES = 0
VOLUME_CHIMERA = 1

VOL = 0
VOL_HALF1 = 1
VOL_HALF2 = 2
VOL_FILTERED = 3

# Template volume names depending on the iteration
VOLNAMES = {
            VOL: 'vol%02d',
            VOL_HALF1: 'vol%02d_sub1',
            VOL_HALF2: 'vol%02d_sub2',
            VOL_FILTERED: 'vol%02d_filtered'
            }



    
class SpiderViewerRefinement(ProtocolViewer):
    """ Wrapper to visualize different type of objects
    with the Xmipp program xmipp_showj. """
    
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [SpiderProtRefinement]
    _label = 'viewer refinement'

    def _defineParams(self, form):
        form.addSection(label='Visualization')
        form.addParam('viewIter', params.EnumParam, 
                      choices=['last', 'selection'], default=ITER_LAST, 
                      display=params.EnumParam.DISPLAY_HLIST,
                      label="Iteration to visualize", important=True,
                      help="""
*last*: only the last iteration will be visualized.
*selection*: you may specify a range of iterations.
Examples:
"1,5-8,10" -> [1,5,6,7,8,10]
"2,6,9-11" -> [2,6,9,10,11]
"2 5, 6-8" -> [2,5,6,7,8]                      
                           """)
        form.addParam('iterSelection', params.NumericRangeParam, 
                      condition='viewIter==%d' % ITER_SELECTION, 
                      label="Iterations list", 
                      help="Write the iteration list to visualize.")

        group = form.addGroup('Angular assignment')
        group.addParam('showImagesAngularAssignment', params.LabelParam, default=True,
                       label='Particles angular assignment')
        group.addParam('displayAngDist', params.EnumParam, choices=['2D plot', 'chimera'], 
                      default=ANGDIST_2DPLOT, display=params.EnumParam.DISPLAY_HLIST, 
                      label='Display angular distribution',
                      help='*2D plot*: display angular distribution as interative 2D in matplotlib.\n'
                           '*chimera*: display angular distribution using Chimera with red spheres.') 
        group.addParam('spheresScale', params.IntParam, default=-1, 
                       condition="displayAngDist == %d" % ANGDIST_CHIMERA,
                       expertLevel=params.LEVEL_ADVANCED,
                       label='Spheres size',
                       help='')

        group = form.addGroup('Volumes')
        group.addParam('displayVol', params.EnumParam, 
                       choices=['slices', 'chimera'], 
                       default=VOLUME_SLICES, display=params.EnumParam.DISPLAY_HLIST, 
                       label='Display volume with',
                       help='*slices*: display volumes as 2D slices along z axis.\n'
                            '*chimera*: display volumes as surface with Chimera.')
        group.addParam('showVolumes', params.EnumParam, default=VOL,
                       choices=['reconstructed', 'half1', 'half2', 'filtered'],
                       label='Volume to visualize',
                       help='Select the volume to visualize')
        
        group = form.addGroup('Resolution')

        group.addParam('showFSC', params.LabelParam, default=True,
                       important=True,
                       label='Display resolution plots (FSC)',
                       help='')
        group.addParam('groupFSC', params.EnumParam, default=1,
                       choices=['iterations', 'defocus groups'],
                       display=params.EnumParam.DISPLAY_HLIST,
                       label='Group FSC plots by',
                       help='Select which FSC curve you want to '
                            'show together in the same plot.')
        group.addParam('groupSelection', params.NumericRangeParam, 
                      condition='groupFSC==%d' % 1, 
                      label="Groups list", 
                      help="Write the group list to visualize. See examples in iteration list")        
        group.addParam('resolutionThresholdFSC', params.FloatParam, default=0.5, 
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Threshold in resolution plots',
                      help='')
                                              
        
    def _getVisualizeDict(self):
        #self._load()
        return {
                'showVolumes': self._showVolumes,
                'showFSC': self._showFSC,
                'displayAngDist': self._displayAngDist
                }
        
    def _validate(self):
        return []
    
    def _formatFreq(self, value, pos):
        """ Format function for Matplotlib formatter. """
        inv = 999.
        if value:
            inv = 1/value
        return "1/%0.2f" % inv
    
    def _getIterations(self):
        if self.viewIter == ITER_LAST:
            return [self.protocol.numberOfIterations.get()] # FIXME: return the last completed iteration
        else:
            return self._getListFromRangeString(self.iterSelection.get(''))
        
    def _getGroups(self):
        return self._getListFromRangeString(self.groupSelection.get(''))
    
    def _getFinalPath(self, *paths):
        return self.protocol._getExtraPath('Refinement', 'final', *paths)
        
        
#===============================================================================
# ShowVolumes
#===============================================================================
    def _createVolumesSqlite(self):
        """ Write an sqlite with all volumes selected for visualization. """

        
        volSqlite = self.protocol._getExtraPath('viewer_volumes.sqlite')
        samplingRate = self.protocol.inputParticles.get().getSamplingRate()
        self.createVolumesSqlite(self.getVolumeNames(), 
                                 volSqlite, samplingRate)
        
        return [self.getObjectView(volSqlite)]
        
    def getVolumeNames(self, it=None):
        """ If it is not none, return the volume of this iteration only. """
        if it is None:
            iterations = self._getIterations()
        else:
            iterations = [it]
            
        volTemplate = VOLNAMES[self.showVolumes.get()]
        volumes = [self._getFinalPath(volTemplate % i) + '.stk'
                   for i in iterations]
        
        return volumes

    def _showVolumesChimera(self):
        """ Create a chimera script to visualize selected volumes. """
        volumes = self.getVolumeNames()
        
        if len(volumes) > 1:
            cmdFile = self._getFinalPath('chimera_volumes.cmd')
            f = open(cmdFile, 'w+')
            for vol in volumes:
                # We assume that the chimera script will be generated
                # at the same folder than spider volumes
                localVol = os.path.basename(vol)
                if os.path.exists(vol):
                    f.write("open spider:%s\n" % localVol)
            f.write('tile\n')
            f.close()
            view = em.ChimeraView(cmdFile)
        else:
            #view = CommandView('xmipp_chimera_client --input "%s" --mode projector 256 &' % volumes[0])
            #view = em.ChimeraClientView(volumes[0])
            view = em.ChimeraClientView(volumes[0], showProjection=True)#, angularDistFile=sqliteFn, spheresDistance=radius)
            
        return [view]
            
    def _showVolumes(self, paramName=None):
        if self.displayVol == VOLUME_CHIMERA:
            return self._showVolumesChimera()
        
        elif self.displayVol == VOLUME_SLICES:
            return self._createVolumesSqlite()#self._createVolumesMd()
    
#===============================================================================
# plotFSC            
#===============================================================================
    def _plotFSC(self, a, fscFile):
        resolution = []
        fsc = []
        
        fscDoc = SpiderDocFile(fscFile)
        for values in fscDoc:
            resolution.append(1/values[1])
            fsc.append(values[2])

        self.maxfsc = max(fsc)
        self.minInv = min(resolution)
        self.maxInv = max(resolution)
        a.plot(resolution, fsc)
        from matplotlib.ticker import FuncFormatter
        a.xaxis.set_major_formatter(FuncFormatter(self._formatFreq))
        a.set_ylim([-0.1, 1.1])
        fscDoc.close()
            
    def _showFSC(self, paramName=None):
        threshold = self.resolutionThresholdFSC.get()
        iterations = self._getIterations()
        groups = self._getGroups()
        
        if self.groupFSC == 0: # group by iterations           
            files = [(it, self._getFinalPath('fscdoc_%02d.stk' % it)) for it in iterations]
            legendPrefix = 'iter'
        else:
            it = iterations[-1]
            legendPrefix = 'group'
            def group(f): # retrieve the group number
                return int(f.split('_')[-1].split('.')[0])
            groupFiles = glob(self._getFinalPath('fscdoc_%02d_???.stk' % it))
            groupFiles.sort()
            files = [(group(f), f) for f in groupFiles if group(f) in groups]
            if not files: #empty files
                return [self.errorMessage("Please select valid groups to display", 
                                          title="Wrong groups selection")]
                
        plotter = EmPlotter(x=1, y=1, windowTitle='Resolution FSC')
        a = plotter.createSubPlot("FSC", 'Angstroms^-1', 'FSC', yformat=False)
        #fscFile = self._getFinalPath('fscdoc_%02d.stk' % iterations[0])
        legends = []
        for it, fscFile in files:
            if os.path.exists(fscFile):
                self._plotFSC(a, fscFile)
                legends.append('%s %d' % (legendPrefix, it))
            else:
                print "Missing file: ", fscFile
            
        if threshold < self.maxfsc:
            a.plot([self.minInv, self.maxInv],[threshold, threshold], 
                   color='black', linestyle='--')
        
        plotter.showLegend(legends)
        a.grid(True)
        
        return [plotter]
    
    def _iterAngles(self, it):
        """ Iterate over the angular distribution for a given iteration. """
        # Get the alignment files of each group for this iteration
        files = glob(self._getFinalPath('align_%02d_???.stk' % it))
        for anglesFile in files:
            fscDoc = SpiderDocFile(anglesFile)
            for values in fscDoc:
                theta = values[1]
                phi = values[2]
                
                if theta > 90:
                    theta = abs(180. - theta)
                    phi += 180
                yield phi, theta
            fscDoc.close()
        
    def _displayAngDist(self, *args):
        #print "_displayAngDist...."
        iterations = self._getIterations()
        nparts = self.protocol.inputParticles.get().getSize()
        views = []
        
        if self.displayAngDist == ANGDIST_2DPLOT:
	    #print " self.displayAngDist == ANGDIST_2DPLO "
            for it in iterations:
                anglesSqlite = self._getFinalPath('angular_dist_%03d.sqlite' % it)
                title = 'Angular distribution iter %03d' % it
                plotter = EmPlotter(x=1, y=1, windowTitle=title)
                self.createAngDistributionSqlite(anglesSqlite, nparts, 
                                                 itemDataIterator=self._iterAngles(it))
                plotter.plotAngularDistributionFromMd(anglesSqlite, title)
                views.append(plotter)
        else:
            it = iterations[-1]
            print "Using last iteration: ", it
            anglesSqlite = self._getFinalPath('angular_dist_%03d.sqlite' % it)
            self.createAngDistributionSqlite(anglesSqlite, nparts, 
                                                 itemDataIterator=self._iterAngles(it))
            volumes = self.getVolumeNames(it)
            views.append(em.ChimeraClientView(volumes[0], 
                                              showProjection=True, 
                                              angularDistFile=anglesSqlite, 
                                              spheresDistance=2))#self.spheresScale.get()))
            
        return views
    
