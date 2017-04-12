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
Some Spider protocol base classes.
"""
from os.path import basename

from pyworkflow.em import ProtClassify2D, Particle, Class2D
from pyworkflow.em.convert import ImageHandler
from pyworkflow.protocol.params import PointerParam, IntParam
from pyworkflow.utils import removeExt, copyFile
import pyworkflow.utils.graph as graph

from ..spider import SpiderDocFile
from protocol_base import SpiderProtocol



class SpiderProtClassify(ProtClassify2D, SpiderProtocol):
    """ Base for Spider classification protocols
    """
    def __init__(self, script, classDir, **kwargs):
        ProtClassify2D.__init__(self, **kwargs)
        SpiderProtocol.__init__(self, **kwargs)
        # To avoid showing MPI box due to duplicated init
        self.allowMpi = False
        self._script = script
        self._classDir = classDir
        
        self._params = {'ext': 'stk',
                        '[class_dir]': self._classDir,
                        'particles': 'input_particles',
                        'particlesSel': 'input_particles_sel',
                        'dendroPs': 'dendrogram',
                        'dendroDoc': '%s/docdendro' % self._classDir,
                        'averages': 'averages'
                        }  
        
    def getClassDir(self):
        return self._classDir
    
    def getNumberOfClasses(self):
        return None
    
    #--------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        self._defineBasicParams(form)
        form.addParallelSection(threads=4, mpi=0)

    def _defineBasicParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputParticles', PointerParam,
                      label="Input particles", important=True,
                      pointerClass='SetOfParticles',
                      help='Input images to perform PCA')
        form.addParam('pcaFile', PointerParam, pointerClass='PcaFile',
                      label="IMC/SEQ file", 
                      help='The IMC file contains the coordinates of each image '
                           'in the reduced-dimension space. '
                           'The SEQ file contains, for all images, '
                           'the pixel values under the mask. ')
        form.addParam('numberOfFactors', IntParam, default=10,
                      label='Number of factors',
                      help='After running, examine the eigenimages and decide '
                           'which ones to use. \n'
                           'Typically all but the first few are noisy.')
        

    
    #--------------------------- INSERT steps functions -----------------------
    
    def _insertAllSteps(self):    
        
        pcaFile = self.pcaFile.get().filename.get()
        
        self._insertFunctionStep('convertInput', 'inputParticles',
                                 self._getFileName('particles'),
                                 self._getFileName('particlesSel'))
        
        self._insertFunctionStep('classifyStep', pcaFile, 
                                 self.numberOfFactors.get(),
                                 self.getNumberOfClasses())
        
        self._insertFunctionStep('createOutputStep')
        
    #--------------------------- STEPS functions ------------------------------

    def _updateParams(self):   
        pass 
        
    def classifyStep(self, imcFile, numberOfFactors, numberOfClasses):
        """ Apply the selected filter to particles. 
        Create the set of particles.
        """
        # Copy file to working directory, it could be also a link
        imcLocalFile = basename(imcFile)
        copyFile(imcFile, self._getPath(imcLocalFile))
        self.info("Copied file '%s' to '%s' " % (imcFile, imcLocalFile))
        # Spider automatically add _IMC to the ca-pca result file
        # JMRT: I have modify the kmeans.msa script to not add _IMC
        # automatically, it can be also used with _SEQ suffix,
        # so we will pass the whole cas_file
        imcBase = removeExt(imcLocalFile)#.replace('_IMC', '')
        imcPrefix = imcBase.replace('_IMC', '').replace('_SEQ', '')
        
        self._params.update({'x27': numberOfFactors,
                             'x30': self.numberOfThreads.get(),
                             '[cas_prefix]': imcPrefix,
                             '[cas_file]': imcBase,
                             })
        self._updateParams()

        self.runTemplate(self.getScript(), self.getExt(), self._params)


class SpiderProtClassifyCluster(SpiderProtClassify):
    """ Base for Clustering Spider classification protocols.
    """
    def __init__(self, script, classDir, **kwargs):
        SpiderProtClassify.__init__(self, script, classDir, **kwargs)

    #--------------------------- STEPS functions ------------------------------
       
    def createOutputStep(self):
        self.buildDendrogram(True)
         
    #--------------------------- UTILS functions ------------------------------
    
    def _fillClassesFromNodes(self, classes2D, nodeList):
        """ Create the SetOfClasses2D from the images of each node
        in the dendrogram.
        """
        particles = classes2D.getImages()
        sampling = classes2D.getSamplingRate()

        # We need to first create a map between the particles index and
        # the assigned class number
        classDict = {}
        nodeDict = {}
        classCount = 0
        for node in nodeList:
            if node.path:
                classCount += 1
                node.classId = classCount
                nodeDict[classCount] = node
                for i in node.imageList:
                    classDict[int(i)] = classCount

        def updateItem(p, i):
            classId = classDict.get(i, None)
            if classId is None:
                p._appendItem = False
            else:
                p.setClassId(classId)

        def updateClass(cls):
            node = nodeDict[cls.getObjId()]
            rep = cls.getRepresentative()
            rep.setSamplingRate(sampling)
            rep.setLocation(node.avgCount, self.dendroAverages)

        particlesRange = range(1, particles.getSize()+1)

        classes2D.classifyItems(updateItemCallback=updateItem,
                                updateClassCallback=updateClass,
                                itemDataIterator=iter(particlesRange))
                
    def _fillParticlesFromNodes(self, inputParts, outputParts, nodeList):
        """ Create the SetOfClasses2D from the images of each node
        in the dendrogram.
        """
        allImages = set()
        for node in nodeList:
            if node.path:
                for i in node.imageList:
                    allImages.add(i)

        def updateItem(item, index):
            item._appendItem = index in allImages

        particlesRange = range(1, inputParts.getSize()+1)

        outputParts.copyItems(inputParts,
                              updateItemCallback=updateItem,
                              itemDataIterator=iter(particlesRange))

    def buildDendrogram(self, writeAverages=False):
        """ Parse Spider docfile with the information to build the dendrogram.
        Params:
            writeAverages: whether to write class averages or not.
        """
        dendroFile = self._getFileName('dendroDoc')
        # Dendrofile is a docfile with at least 3 data colums (class, height, id)
        doc = SpiderDocFile(dendroFile)
        values = []
        indexes = []
        for _, h, i in doc.iterValues():
            indexes.append(i)
            values.append(h)
        doc.close()
        
        self.dendroValues = values
        self.dendroIndexes = indexes
        self.dendroImages = self._getFileName('particles')
        self.dendroAverages = self._getFileName('averages')
        self.dendroAverageCount = 0 # Write only the number of needed averages
        self.dendroMaxLevel = 10 # FIXME: remove hard coding if working the levels
        self.ih = ImageHandler()

        return self._buildDendrogram(0, len(values)-1, 1, writeAverages)

    def getImage(self, particleNumber):
        return self.ih.read((particleNumber, self.dendroImages))
        
    def addChildNode(self, node, leftIndex, rightIndex, index,
                     writeAverages, level, searchStop):
        child = self._buildDendrogram(leftIndex, rightIndex, index,
                                      writeAverages, level+1, searchStop)
        node.addChild(child)
        node.extendImageList(child.imageList)

        if writeAverages:
            node.addImage(child.image)
            del child.image # Allow to free child image memory
                
    def _buildDendrogram(self, leftIndex, rightIndex, index,
                         writeAverages=False, level=0, searchStop=0):
        """ This function is recursively called to create the dendrogram
        graph (binary tree) and also to write the average image files.
        Params:
            leftIndex, rightIndex: the indexes within the list where to search.
            index: the index of the class average.
            writeAverages: flag to select when to write averages
            searchStop: this could be 1, means that we will search until the
                last element (used for right childs of the dendrogram or,
                can be 0, meaning that the last element was already the max
                (used for left childs )
        From self:
            self.dendroValues: the list with the heights of each node
            self.dendroImages: image stack filename to read particles
            self.dendroAverages: stack name where to write averages
        It will search for the max in values list (between minIndex and maxIndex).
        Nodes to the left of the max are left childs and the other right childs.
        """
        if level < self.dendroMaxLevel:
            avgCount = self.dendroAverageCount + 1
            self.dendroAverageCount += 1
                    
        if rightIndex == leftIndex: # Just only one element
            height = self.dendroValues[leftIndex]
            node = DendroNode(index, height)
            node.extendImageList([self.dendroIndexes[leftIndex]])
            node.addImage(self.getImage(node.imageList[0]))

        elif rightIndex == leftIndex + 1: # Two elements
            height = max(self.dendroValues[leftIndex], 
                         self.dendroValues[rightIndex])
            node = DendroNode(index, height)
            node.extendImageList([self.dendroIndexes[leftIndex],
                                  self.dendroIndexes[rightIndex]])
            node.addImage(self.getImage(node.imageList[0]),
                          self.getImage(node.imageList[1]))
        else: # 3 or more elements
            # Find the max value (or height) of the elements
            maxValue = self.dendroValues[leftIndex]
            maxIndex = 0
            # searchStop could be 0 (do not consider last element, coming from
            # left child, or 1 (consider also the last one, coming from right)
            values = self.dendroValues[leftIndex+1:rightIndex+searchStop]
            for i, v in enumerate(values):
                if v > maxValue:
                    maxValue = v
                    maxIndex = i+1
            m = maxIndex + leftIndex
            node = DendroNode(index, maxValue)
            hasRightChild = m < rightIndex

            if maxValue > 0:
                nextIndex = 2 * index if hasRightChild else index
                self.addChildNode(node, leftIndex, m, nextIndex,
                                  writeAverages, level, 0)

                if hasRightChild:
                    self.addChildNode(node, m+1, rightIndex, 2 * index + 1,
                                      writeAverages, level, 1)
                else:
                    # If the node has a single child, we will remove a node
                    # just to advance in the level of the tree to get more
                    # different class averages
                    if node.getChilds():
                        child = node.getChilds()[0]
                        child.image = node.image
                        child.parents = []
                        node = child
            else:
                node.extendImageList(self.dendroIndexes[leftIndex:rightIndex+1])
                node.addImage(*[self.getImage(img) for img in node.imageList])

        if level < self.dendroMaxLevel:
            node.avgCount = avgCount
            node.path = '%d@%s' % (node.avgCount, self.dendroAverages)
            
            if writeAverages:
                # normalize the sum of images depending on the number of particles
                # assigned to this classes
                avgImage = node.image / float(node.getSize()) 
                self.ih.write(avgImage, (node.avgCount, self.dendroAverages))
                fn = self._getTmpPath('doc_class%03d.stk' % index)
                doc = SpiderDocFile(fn, 'w+')
                for i in node.imageList:
                    doc.writeValues(i)
                doc.close()
                
        return node
    

class DendroNode(graph.Node):
    """ Special type of Node to store dendrogram values. """
    def __init__(self, index, height):
        graph.Node.__init__(self, 'class_%03d' % index)
        self.index = index
        self.height = height
        self.length = 0
        self.path = None
        self.selected = False
        self.imageList = []
        self.image = None
        
    def getChilds(self):
        return [c for c in self._childs if c.path]
    
    def getSize(self):
        """ Return the number of images assigned to this class. """
        return len(self.imageList)

    def addImage(self, *images):
        """ Add some images to this node. """
        for img in images:
            if self.image is None:
                self.image = img
            else:
                self.image.inplaceAdd(img)

    def extendImageList(self, imageList):
        self.imageList.extend(imageList)
        self.length = self.getSize()