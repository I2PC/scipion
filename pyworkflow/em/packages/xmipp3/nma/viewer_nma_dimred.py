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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
from pyworkflow.em.packages.xmipp3.nma.data import PathData
"""
This module implement the wrappers aroung Xmipp CL2D protocol
visualization program.
"""

from os.path import basename, join, exists
import numpy as np

from pyworkflow.utils.path import cleanPath, makePath, cleanPattern
from pyworkflow.viewer import (ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO)
from pyworkflow.protocol.params import StringParam, LabelParam
from pyworkflow.em.data import SetOfParticles
from pyworkflow.utils.process import runJob
from pyworkflow.em.viewer import VmdView
from pyworkflow.gui.browser import FileBrowserWindow

from protocol_nma_dimred import XmippProtDimredNMA, DIMRED_MAPPINGS
from data import Point, Data
from plotter import XmippNmaPlotter
from gui import ClusteringWindow, TrajectoriesWindow

        
class XmippDimredNMAViewer(ProtocolViewer):
    """ Visualization of results from the NMA protocol
    """
    _label = 'viewer nma dimred'
    _targets = [XmippProtDimredNMA]
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    
    def __init__(self, **kwargs):
        ProtocolViewer.__init__(self, **kwargs)
        self._data = None

    def getData(self):
        if self._data is None:
            self._data = self.loadData()
        return self._data

    def _defineParams(self, form):
        form.addSection(label='Visualization')
        form.addParam('displayRawDeformation', StringParam, default='1',
                      label='Display raw deformation',
                      help='Type 1 to see the histogram of raw deformation number 1; \n'
                           'type 2 to see the histogram of raw deformation number 2, etc.\n'
                           'Type 1 2 to see the 2D plot of raw deformations number 1 vs 2.\n'
                           'Type 1 2 3 to see the 3D plot of raw deformations 1, 2 and 3; etc.'
                           )
        
        form.addParam('displayClustering', LabelParam,
                      label='Open clustering tool?',
                      help='Open a GUI to visualize the images as points'
                           'and select some of them to create new clusters.')
         
        form.addParam('displayTrajectories', LabelParam,
                      label='Open trajectories tool?',
                      help='Open a GUI to visualize the images as points'
                           'to draw and ajust trajectories.')       
        
    def _getVisualizeDict(self):
        return {'displayRawDeformation': self._viewRawDeformation,
                'displayClustering': self._displayClustering,
                'displayTrajectories': self._displayTrajectories,
                } 
                        
    def _viewRawDeformation(self, paramName):
        components = self.displayRawDeformation.get()
        return self._doViewRawDeformation(components)
        
    def _doViewRawDeformation(self, components):
#        components = map(int, self.displayRawDeformation.get().split())
        components = map(int, components.split())
        dim = len(components)
        views = []
        
        if dim > 0:
            modeList = [m - 1 for m in components]
            modeNameList = ['Mode %d' % m for m in components]
            missingList = []
                    
            if missingList:
                return [self.errorMessage("Invalid mode(s) *%s*\n." % (', '.join(missingList)), 
                              title="Invalid input")]
            
            # Actually plot
            plotter = XmippNmaPlotter(data=self.getData())
            baseList = [basename(n) for n in modeNameList]
            
            if dim == 1:
                self.getData().XIND = modeList[0]
                plotter.plotArray1D("Histogram for %s" % baseList[0], 
                                    "Deformation value", "Number of images")
            else:
                self.getData().YIND = modeList[1]
                if dim == 2:
                    plotter.plotArray2D("%s vs %s" % tuple(baseList), *baseList)
                elif dim == 3:
                    self.getData().ZIND = modeList[2]
                    plotter.plotArray3D("%s %s %s" % tuple(baseList), *baseList)
            views.append(plotter)
            
        return views
    
    def _displayClustering(self, paramName):
        self.clusterWindow = self.tkWindow(ClusteringWindow, 
                                           title='Clustering Tool',
                                           dim=self.protocol.reducedDim.get(),
                                           data=self.getData(),
                                           callback=self._createCluster
                                           )
        return [self.clusterWindow]
    
    
    def _displayTrajectories(self, paramName):
        self.trajectoriesWindow = self.tkWindow(TrajectoriesWindow, 
                              title='Trajectories Tool',
                              dim=self.protocol.reducedDim.get(),
                              data=self.getData(),
                              callback=self._generateAnimation,
                              loadCallback=self._loadAnimation,
                              numberOfPoints=10
                              )
        return [self.trajectoriesWindow]
        
    def _createCluster(self):
        """ Create the cluster with the selected particles
        from the cluster. This method will be called when
        the button 'Create Cluster' is pressed.
        """
        # Write the particles
        prot = self.protocol
        project = prot.getProject()
        inputSet = prot.getInputParticles()
        fnSqlite = prot._getTmpPath('cluster_particles.sqlite')
        cleanPath(fnSqlite)
        partSet = SetOfParticles(filename=fnSqlite)
        partSet.copyInfo(inputSet)
        for point in self.getData():
            if point.getState() == Point.SELECTED:
                particle = inputSet[point.getId()]
                partSet.append(particle)
        partSet.write()
        partSet.close()
                
        from protocol_batch_cluster import BatchProtNMACluster
        newProt = project.newProtocol(BatchProtNMACluster)
        clusterName = self.clusterWindow.getClusterName()
        if clusterName:
            newProt.setObjLabel(clusterName)
        newProt.inputNmaDimred.set(prot)
        newProt.sqliteFile.set(fnSqlite)
        
        project.launchProtocol(newProt)
        
    def _loadAnimationData(self, obj):
        prot = self.protocol
        animationName = obj.getFileName() # assumes that obj.getFileName is the folder of animation
        animationPath = prot._getExtraPath(animationName)
        #animationName = animationPath.split('animation_')[-1]
        animationRoot = join(animationPath, animationName)
        
        animationSuffixes = ['.vmd', '.pdb', 'trajectory.txt']
        for s in animationSuffixes:
            f = animationRoot + s
            if not exists(f):
                self.errorMessage('Animation file "%s" not found. ' % f)
                return
        
        # Load animation trajectory points
        trajectoryPoints = np.loadtxt(animationRoot + 'trajectory.txt')
        data = PathData(dim=trajectoryPoints.shape[1])
        
        for i, row in enumerate(trajectoryPoints):
            data.addPoint(Point(pointId=i+1, data=list(row), weight=1))
            
        self.trajectoriesWindow.setPathData(data)
        self.trajectoriesWindow.setAnimationName(animationName)
        self.trajectoriesWindow._onUpdateClick()
        
        def _showVmd():
            vmdFn = animationRoot + '.vmd'
            VmdView(' -e %s' % vmdFn).show()        
            
        self.getTkRoot().after(500, _showVmd)
        
        
    def _loadAnimation(self):
        prot = self.protocol
        browser = FileBrowserWindow("Select the animation folder (animation_NAME)", 
                                    self.getWindow(), prot._getExtraPath(), 
                                    onSelect=self._loadAnimationData)
        browser.show()
        
    def _generateAnimation(self):
        prot = self.protocol
        projectorFile = prot.getProjectorFile()
        
        animation = self.trajectoriesWindow.getAnimationName()
        animationPath = prot._getExtraPath('animation_%s' % animation)
        
        cleanPath(animationPath)
        makePath(animationPath)
        animationRoot = join(animationPath, 'animation_%s' % animation)
        
        trajectoryPoints = np.array([p.getData() for p in self.trajectoriesWindow.pathData])
        
        if projectorFile:
            M = np.loadtxt(projectorFile)
            deformations = np.dot(trajectoryPoints, np.linalg.pinv(M))
            np.savetxt(animationRoot + 'trajectory.txt', trajectoryPoints)
        else:
            Y = np.loadtxt(prot.getOutputMatrixFile())
            X = np.loadtxt(prot.getDeformationFile())
            # Find closest points in deformations
            deformations = [X[np.argmin(np.sum((Y - p)**2, axis=1))] for p in trajectoryPoints]
            
        pdb = prot.getInputPdb()
        pdbFile = pdb.getFileName()
        
        modesFn = prot.inputNMA.get()._getExtraPath('modes.xmd')
        
        for i, d in enumerate(deformations):
            atomsFn = animationRoot + 'atomsDeformed_%02d.pdb' % (i+1)
            cmd = '-o %s --pdb %s --nma %s --deformations %s' % (atomsFn, pdbFile, modesFn, str(d)[1:-1])
            runJob(None, 'xmipp_pdb_nma_deform', cmd, env=prot._getEnviron())
            
        # Join all deformations in a single pdb
        # iterating going up and down through all points
        # 1 2 3 ... n-2 n-1 n n-1 n-2 ... 3, 2
        n = len(deformations)
        r1 = range(1, n+1)
        r2 = range(2, n) # Skip 1 at the end
        r2.reverse()
        loop = r1 + r2
        
        trajFn = animationRoot + '.pdb'
        trajFile = open(trajFn, 'w')
        
        for i in loop:
            atomsFn = animationRoot + 'atomsDeformed_%02d.pdb' % i
            atomsFile = open(atomsFn)
            for line in atomsFile:
                trajFile.write(line)
            trajFile.write('TER\nENDMDL\n')
            atomsFile.close()

        trajFile.close()
        # Delete temporary atom files
        cleanPattern(animationRoot + 'atomsDeformed_??.pdb')    
        
        # Generate the vmd script
        vmdFn = animationRoot + '.vmd'
        vmdFile = open(vmdFn, 'w')
        vmdFile.write("""
        mol new %s
        animate style Loop
        display projection Orthographic
        mol modcolor 0 0 Index
        mol modstyle 0 0 Beads 1.000000 8.000000
        animate speed 0.5
        animate forward
        """ % trajFn)
        vmdFile.close() 
        
        VmdView(' -e ' + vmdFn).show()
            
        
    def loadData(self):
        """ Iterate over the images and the output matrix txt file
        and create a Data object with theirs Points.
        """
        matrix = np.loadtxt(self.protocol.getOutputMatrixFile())
        particles = self.protocol.getInputParticles()
        
        data = Data()
        for i, particle in enumerate(particles):
            data.addPoint(Point(pointId=particle.getObjId(),
                                data=matrix[i, :],
                                weight=particle._xmipp_cost.get()))
            
        return data
    
