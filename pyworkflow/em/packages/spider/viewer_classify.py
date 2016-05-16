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
"""
This module implements the visualization program
for Spider classify protocols.
"""

from os.path import join

import Tkinter as tk

from pyworkflow.em import ProtUserSubSet, SetOfClasses2D
from pyworkflow.protocol.params import IntParam, FloatParam, LabelParam
from pyworkflow.protocol.constants import STATUS_FINISHED
from pyworkflow.utils.properties import Icon
from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO
from pyworkflow.utils.graph import Graph
from pyworkflow.utils.path import cleanPath
from pyworkflow.gui import Window
from pyworkflow.gui.widgets import HotButton
from pyworkflow.gui.graph import LevelTree
from pyworkflow.gui.canvas import Canvas, ImageBox
from pyworkflow.em.viewer import ClassesView
from protocol import SpiderProtClassifyWard, SpiderProtClassifyDiday
from pyworkflow.gui.dialog import askString

from spider import SpiderDocFile


class SpiderViewerClassify(ProtocolViewer):
    """ Visualization of some Classification protocols in Spider. """

    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    
    def _defineParams(self, form):
        form.addSection(label='Visualization')
        group1 = form.addGroup('Dendrogram')
        group1.addParam('doShowDendrogram', LabelParam,
                        label="Show dendrogram", default=True,
                        help='In a dendrogram larger vertical bars signify a '
                             'greater difference between classes. Many small '
                             'differences at the bottom can be eliminated with '
                             'an increase of the "Minimum height" setting.')
        group1.addParam('minHeight', FloatParam, default=0.5,
                      label='Minimum height',
                      help='The dendrogram will be cut at this level')
        self.groupClass = form.addGroup('Classes')
        self.groupClass.addParam('doShowClasses', LabelParam,
                                 label="Visualize class averages", default=True,
                                 help='Display class averages')

    def _getVisualizeDict(self):
        return {'doShowDendrogram': self._plotDendrogram,
                'doShowClasses': self.visualizeClasses
                }

    def _plotDendrogram(self, e=None):
        from pyworkflow.em.packages.xmipp3.plotter import XmippPlotter
        xplotter = XmippPlotter()
        self.plt = xplotter.createSubPlot("Dendrogram", "", "")
        self.step = 0.25
        self.rightMost = 0.0 # Used to arrange leaf nodes at the bottom
        
        node = self.protocol.buildDendrogram()
        self.plotNode(node, self.minHeight.get())    
        self.plt.set_xlim(0., self.rightMost + self.step)
        self.plt.set_ylim(-10, 105)
        
        return [xplotter]
    
    def plotNode(self, node, minHeight=-1):
        childs = node.getChilds()
        h = node.height
        if h > minHeight and len(childs) > 1:
            x1, y1 = self.plotNode(childs[0], minHeight)
            x2, y2 = self.plotNode(childs[1], minHeight)
            xm = (x1 + x2)/2
            x = [x1, x1, x2, x2]
            y = [y1, h, h, y2]
            self.plt.plot(x, y, color='b')
            point = (xm, h)
        else:
            self.rightMost += self.step
            point = (self.rightMost, 0.)
        
        length = node.length
        index = node.index
        self.plt.annotate("%d(%d)" % (index, length), point, xytext=(0, -5),
                          textcoords='offset points', va='top', ha='center',
                          size='x-small')
        self.plt.plot(point[0], point[1], 'ro')
        
        return point
    
    def _createNode(self, canvas, node, y):
        node.selected = False
        node.box = SpiderImageBox(canvas, node, y)
        return node.box
            
            
class SpiderViewerWard(SpiderViewerClassify):
    """ Visualization of Spider - classify Ward protocol results. """
    
    _targets = [SpiderProtClassifyWard]
    _label = "viewer ward"
    
    def _defineParams(self, form):
        SpiderViewerClassify._defineParams(self, form)

        self.groupClass.addParam('maxLevel', IntParam, default=4,
                      label='Maximum level',
                      help='Maximum level of classes to show')

    def visualizeClasses(self, e=None):
        classTemplate = "class_%03d"
        averages = '%03d@' + self.protocol._getFileName('averages')
        
        def getInfo2(level, classNo):
            return classTemplate % classNo, averages % classNo
        
        node = self.protocol.buildDendrogram(writeAverages=False,
                                             stripSingleChild=True)
        
        g = Graph(root=node)
        self.graph = g
               
        self.win = Window("Select classes", self.formWindow, minsize=(1000, 600))
        root = self.win.root
        canvas = Canvas(root)
        canvas.grid(row=0, column=0, sticky='nsew')
        root.grid_columnconfigure(0, weight=1)
        root.grid_rowconfigure(0, weight=1) 
        
        self.buttonframe = tk.Frame(root)
        self.buttonframe.grid(row=2, column=0, columnspan=2)  
        self.win.createCloseButton(self.buttonframe).grid(row=0, column=0,
                                                          sticky='n',
                                                          padx=5, pady=5)
        saveparticlesbtn = HotButton(self.buttonframe, "Particles",
                                     Icon.PLUS_CIRCLE,
                                     command=self._askCreateParticles)
        saveparticlesbtn.grid(row=0, column=1, sticky='n', padx=5, pady=5)  
        btn = HotButton(self.buttonframe, "Classes", Icon.PLUS_CIRCLE,
                        command=self._askCreateClasses)
        btn.grid(row=0, column=2, sticky='n', padx=5, pady=5)
            
        lt = LevelTree(g)
        lt.DY = 135 # TODO: change in percent of the image size
        lt.setCanvas(canvas)
        lt.paint(self._createNode, maxLevel=self.maxLevel.get()-1)
        canvas.updateScrollRegion()
        
        return [self.win]

    def _askCreateParticles(self):
        self._askCreateSubset('Particles', self.getSelectedNodesCount(2))

    def _askCreateClasses(self):
        self._askCreateSubset('Classes', self.getSelectedNodesCount(1))

    def _askCreateSubset(self, output, size):

        if self._selectionOverlap():
            self.win.showError("Classes could not overlap in the tree.")
            return

        s = '' if size == 1 else 's'
        headerLabel = 'Are you sure you want to create a new set of ' \
                      ' %s with %s element%s?' % (output, size, s)
        runname =  askString('Question','Run name:', self.win.getRoot(), 30,
                             defaultValue='ProtUserSubSet',
                             headerLabel=headerLabel)
        if runname:
            createFunc = getattr(self, 'save' + output)
            createFunc(runname)

    def _createSubsetProtocol(self, createOutputFunc, label=None):
        """ Create a subset of classes or particles. """
        try:
            project = self.getProject()
            prot = project.newProtocol(ProtUserSubSet)
            prot.setObjLabel(label)
            prot.inputObject.set(self.protocol)
            project._setupProtocol(prot)
            prot.makePathsAndClean()
            createOutputFunc(prot)
            prot.setStatus(STATUS_FINISHED)
            project._storeProtocol(prot)
            #self.project.launchProtocol(prot, wait=True)

        except Exception, ex:
            import traceback
            traceback.print_exc()    
            self.win.showError(str(ex))
        
    def getSelectedNodesCount(self, depth):
        if depth == 1:
            return len([node for node in self.graph.getNodes() if node.selected])
        else:
            count = 0
            for node in self.graph.getNodes():
                if node.selected:
                    count += len(node.imageList)
            return count

    def _selectionOverlap(self):
        """ Check if selected classes do overlap. """
        allImages = set()
        for n in self._selectedNodes():
            for i in n.imageList:
                if i in allImages:
                    return True
                allImages.add(i)
        return False

    def _selectedNodes(self):
        return [node for node in self.graph.getNodes() if node.selected]

    def saveClasses(self, runname=None):
        """ Store selected classes. """
        def createClasses(prot):
            classes = prot._createSetOfClasses2D(self.protocol.inputParticles.get(),
                                                 suffix='Selection')
            self.protocol._fillClassesFromNodes(classes, self._selectedNodes())
            prot._defineOutputs(outputClasses=classes)

        self._createSubsetProtocol(createClasses, runname)
            
    def saveParticles(self, runname=None):
        """ Store particles from selected classes. """
        def createParticles(prot):
            inputParticles = self.protocol.inputParticles.get()
            particles = prot._createSetOfParticles(suffix='Selection')
            particles.copyInfo(inputParticles)
            self.protocol._fillParticlesFromNodes(inputParticles,
                                                  particles,
                                                  self._selectedNodes())
            prot._defineOutputs(outputParticles=particles)

        self._createSubsetProtocol(createParticles, runname)


class SpiderImageBox(ImageBox):
    def __init__(self, canvas, node, y):
        ImageBox.__init__(self, canvas, node.path, text=node.getName(), y=y)
        
    def _onClick(self, e=None):
        if self.node.path is None:
            return
        # On click change the selection state
        self.node.selected = not self.node.selected
        if self.node.selected:
            self.label.config(bd=2, bg='green')
        else:
            self.label.config(bd=0, bg='grey')
            
            
class SpiderViewerDiday(SpiderViewerClassify):
    """ Visualization of Spider - classify Diday protocol results. """
    
    _targets = [SpiderProtClassifyDiday]
    _label = "viewer diday"
    
    def _defineParams(self, form):
        SpiderViewerClassify._defineParams(self, form)

        self.groupClass.addParam('numberOfClasses', IntParam, default=4,
                      label='Number of classes',
                      help='Desired number of classes.')

    def visualizeClasses(self, e=None):
        prot = self.protocol
        classDir = prot.getClassDir()
        classAvg = 'classavg'
        classVar = 'classvar'
        classDoc = 'docclass'
        
        params = {'[class_dir]': classDir,
                  '[desired-classes]': self.numberOfClasses.get(),
                  '[particles]': prot._params['particles'] + '@******',
                  '[class_doc]': join(classDir, classDoc + '***'), 
                  '[class_avg]': join(classDir, classAvg + '***'),
                  '[class_var]': join(classDir, classVar + '***'),        
                  }
        
        prot.runTemplate('mda/classavg.msa', prot.getExt(), params)

        particles = prot.inputParticles.get()
        particles.load()
        sampling = particles.getSamplingRate()
        
        setFn = prot._getTmpPath('classes2D.sqlite')
        cleanPath(setFn)
        classes2D = SetOfClasses2D(filename=setFn)
        classes2D.setImages(particles)

        # We need to first create a map between the particles index and
        # the assigned class number
        classDict = {}
        for classId in range(1, self.numberOfClasses.get()+1):
            docClass = prot._getPath(classDir, classDoc + '%03d.stk' % classId)
            doc = SpiderDocFile(docClass)
            for values in doc.iterValues():
                imgIndex = int(values[0])
                classDict[imgIndex] = classId
            doc.close()

        updateItem = lambda p, i: p.setClassId(classDict[i])

        def updateClass(cls):
            rep = cls.getRepresentative()
            rep.setSamplingRate(particles.getSamplingRate())
            avgFn = prot._getPath(classDir,
                                  classAvg + '%03d.stk' % cls.getObjId())
            rep.setLocation(1, avgFn)

        particlesRange = range(1, particles.getSize()+1)
        classes2D.classifyItems(updateItemCallback=updateItem,
                                updateClassCallback=updateClass,
                                itemDataIterator=iter(particlesRange))

        classes2D.write()
        classes2D.close()

        return [ClassesView(self.getProject(), prot.strId(),
                            classes2D.getFileName(), particles.strId())]
 
