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
""" 
This module implements a simple algorithm to display a graph(mainly a tree)
level by level, using only Tkinter.
"""
import Tkinter as tk

#TODO: all LevelTree code is DEPRECATED...remove it after cleanning
# and include all code from graph_layout

from graph_layout import *


class LevelTree(object):
    """ Class to render the Graph in a Canvas. """
    def __init__(self, graph):
        self.DY = 65
        self.DX = 15
        self.FONT = "sans-serif"
        self.FONTSIZE = 9
        self.graph = graph
        self.canvas = None
        
    def setCanvas(self, canvas):
        self.canvas = canvas
        
    def paint(self, createNode=None, createEdge=None, maxLevel=9999, usePositions=False):
        """ Paint the Graph, nodes will be positioned by levels.
        Params:
            canvas: the canvas object to paint the graph.
            createNode: function to build the Item that represents a Node.
                    the Item created should have the following methods:
                        getDimensions: return the width and height
                        moveTo: change the position of the Item
            createEdge: function to build an Edge connection two Nodes
            usePositions: if this is True, use the nodes positions without
                recomputing them.
            
            If createNode and createEdge are None, the default ones will be used,
            that requires the setCanvas method had to be called first.
        """
        self.createNode = createNode or self._defaultCreateNode
        self.createEdge = createEdge or self._defaultCreateEdge
        self.maxLevel = maxLevel
        rootNode = self.graph.getRoot()
        
        if usePositions:
            self._paintNodeWithPosition(rootNode)
            self._paintEdges(rootNode)
        else:            
            self._setLevel(rootNode, 0, None)
            self._paintNodeWithChilds(rootNode, 1)
            m = 9999
            for left, right in rootNode.hLimits:
                m = min(m, left)
            self._createEdges(rootNode, -m + self.DY)
        
    def _setLevel(self, node, level, parent):
        """ Set the level of the nodes. """
        node.level = level
        node.parent = parent
        nextLevel = level + 1
        if nextLevel > self.maxLevel:
            return
        for child in node.getChilds():
            if nextLevel > getattr(child, 'level', 0):
                self._setLevel(child, nextLevel, node)
                 
    def _paintNodeWithChilds(self, node, level):
        y = level * self.DY
        
        self._paintNode(node, y)
        
        if level > self.maxLevel:
            return
        
        childs = [c for c in node.getChilds() if c.parent is node]
        n = len(childs)

        if n > 0:
            #width = (xmax - xmin) / n
            for c in childs:
                self._paintNodeWithChilds(c, level + 1)
                
            if n > 1:
                offset = 0
                for i in range(n-1):
                    sep = self._getChildsSeparation(childs[i], childs[i+1])
                    offset += sep
                    c = childs[i+1]
                    c.offset = offset
                
                total = childs[0].half + offset + childs[-1].half
                half = total/2
                for c in childs:
                    c.offset -= half - childs[0].half
                
            else:
                childs[0].offset = 0
            self._getHLimits(node)
        
    def _defaultCreateNode(self, canvas, node, y):
        """ If not createNode is specified, this one will be used
        by default. 
        """
        if canvas is None:
            raise Exception("method setCanvas should be called before using _defaultCreateNode")
        nodeText = node.getLabel()
        textColor = 'black'
        if nodeText.startswith('Project'):
            textColor='white'
        
        return canvas.createTextbox(nodeText, 100, y, bgColor='light blue', textColor=textColor, margin=0)
        
    def _defaultCreateEdge(self, srcItem, dstItem):
        if self.canvas is None:
            raise Exception("method setCanvas should be called before using _defaultCreateEdge")
        self.canvas.createEdge(srcItem, dstItem)
        
    def _paintNode(self, node, y):
        """ Paint a node of the graph.
        Params:
           canvas: the canvas in which to paint.
           node: the node of the graph to be painted.
           y: level in the tree where to paint.
        Returns:
           the create item in the canvas.
        """
        item = self.createNode(self.canvas, node, y)
        node.width, node.height = item.getDimensions()
        node.half = node.width / 2
        node.hLimits = [[-node.half, node.half]]
        node.y = item.y
        node.offset = 0
        # Create link from both sides to reach
        node.item = item 
        item.node = node
        
        return item
    
    def _printHLimits(self, node, msg):
        print "\n=====%s========" % msg
        print " dd: %s" % node.t.text.replace('\n', '_')
        print "  offset: %d, width: %d" % (node.offset, node.width)
        print "  hlimits:"
        for l, r in node.hLimits:
            print "   [%d, %d]" % (l, r)
            
    def _getHLimits(self, node):
        '''This function will traverse the tree
        from node to build the left and right profiles(hLimits)
        for each level of the tree'''
        node.hLimits = [[-node.half, node.half]]
        #print "getHLimits, parent: ", node.t.text
        childs = [c for c in node.getChilds() if c.parent is node]
        for child in childs:
            count = 1
            #printHLimits(c, " child")
            if not hasattr(child, 'hLimits'):
                print "node %s has no hLimits" % child.label
                raise Exception()
            
            for l, r in child.hLimits:
                l += child.offset
                r += child.offset
                #print "  level ", count
                #print "   l, r", l, r
                if count < len(node.hLimits):
                    #print "   node.hLimits(A)", node.hLimits[count]
                    if l < node.hLimits[count][0]:
                        node.hLimits[count][0] = l
                    if r > node.hLimits[count][1]:
                        node.hLimits[count][1] = r
                else:
                    node.hLimits.append([l, r])
                    #print "   node.hLimits(A)", node.hLimits[count]
                #print "   node.hLimits(B)", node.hLimits[count]
                count += 1
                
    def _getChildsSeparation(self, child1, child2):
        '''Calcualte separation between siblings
        at each height level'''
        sep = 0 
        hL1 = child1.hLimits
        hL2 = child2.hLimits
        n1 = len(hL1)
        n2 = len(hL2)
        h = min(n1, n2)
            
        for i in range(h):
            right = hL1[i][1]
            left = hL2[i][0]            
            if left + sep < right:
                sep = right - left                
  
        return sep + self.DX
    
    def _createEdges(self, node, x):
        """ Adjust the position of the nodes
        and create the edges between them.
        """
        nx = x + node.offset
        node.item.moveTo(nx, node.y)
        
        if node.level == self.maxLevel:
            return 
        
        for c in node.getChilds():
            if c.parent is node:
                self._createEdges(c, nx)
            self.createEdge(node.item, c.item)
            
    def _paintNodeWithPosition(self, node):
        """ Paint nodes using its position. """
        self._paintNode(node, None)
        
        for child in node.getChilds():
            # parent = None for nodes that have been not traversed
            parent = getattr(child, 'parent', None)
            if parent is None:
                child.parent = node
                self._paintNodeWithPosition(child)
                
    def _paintEdges(self, node):
        """ Paint only the edges between nodes, asumming they are 
        already well positioned. 
        """
        for child in node.getChilds():
            if child.parent is node:
                self._paintEdges(child)        
            self.createEdge(node.item, child.item)
            


def createNode(canvas, node, y):
    from canvas import ImageBox
    return ImageBox(canvas, node.path, text=node.getName(), y=y)
    
def createEdge(srcItem, dstItem):
    pass
    
    
if __name__ == '__main__':
    from canvas import Canvas
    from pyworkflow.utils.graph import Graph
    
    root = tk.Tk()
    canvas = Canvas(root, width=600, height=500)
    canvas.grid(row=0, column=0, sticky='nsew')
    root.grid_columnconfigure(0, weight=1)
    root.grid_rowconfigure(0, weight=1)
    
    import sys
    
    maxLevel = 2
    if len(sys.argv) > 1:
        maxLevel = int(sys.argv[1])
    
    def example1():
        """ Display a very simple graph: A, B and C. """
        g = Graph(rootName='A')
        a = g.getRoot()
        b = g.createNode('B')
        c = g.createNode('C')
        a.addChild(b, c)
        lt = LevelTree(g)
        lt.setCanvas(canvas)
        #lt.paint(createNode, createEdge)
        lt.paint()

    def example2():
        avg = '%03d@/home/josem/Scipion/projects/analu_MCM/Runs/SpiderProtClassifyWard4236/averages.stk'
        name2 = "class_%03d"
        path = '%d@/home/josem/Scipion/projects/analu_MCM/Runs/XmippProtCL2D4394/extra/level_%02d/classes_sorted.stk'
        name = 'level_%02d_%02d'
        
        def getInfo(level, classNo):
            return name % (level, classNo), path % (classNo, level)
        
        def getInfo2(level, classNo):
            return name2 % classNo, avg % classNo
        
        n, p = getInfo2(0, 1)
        g = Graph(rootName=n)
        a = g.getRoot()
        a.path = p
        
        def addChilds(node, nodeNumber, level):
            imgNo = nodeNumber * 2 
            for off in [0, 1]:
                n, p = getInfo2(level+1, imgNo + off)
                b = g.createNode(n)
                b.path = p
                node.addChild(b)
                if level < maxLevel - 2:
                    addChilds(b, imgNo + off, level + 1)
                        
            
        addChilds(a, 1, 0)
            
        lt = LevelTree(g)
        lt.DY = 135
        lt.setCanvas(canvas)
        lt.paint(createNode)
        
    example2()
    canvas.updateScrollRegion()
    
    root.mainloop()             
