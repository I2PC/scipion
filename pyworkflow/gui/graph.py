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

class LevelTree(object):
    """ Class to render the Graph in a Canvas. """
    def __init__(self, graph):
        self.DY = 50
        self.DX = 15
        self.FONT = "sans-serif"
        self.FONTSIZE = 9
        self.graph = graph
        
    def paint(self, canvas, createNodeItem=None):
        """ Paint the Graph.
        Params:
            canvas: the canvas object to paint the graph.
            createNodeItem: function to build the item in the canvas to represent the node.
                if not function is passed, a default one will be used.  
        """
        self.canvas = canvas
        self.createNodeItem = createNodeItem or self._defaultNodeItem
        rootNode = self.graph.getRoot()
        self._paintNodeWithChilds(rootNode, 1)
        m = 9999
        for left, right in rootNode.hLimits:
            m = min(m, left)
        self._createEdges(rootNode, -m + self.DY)
        
    def _paintNodeWithChilds(self, node, level):
        y = level * self.DY
        childs = node.getChilds()
        n = len(childs)
        
        self._paintNode(node, y)
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
                
#                print "\n\n----A------"
#                print "Parent: %s" % node.t.text.replace('\n', '_')
#                for c in childs:
#                    print "  child: %s, width: %d, offset: %d" % (c.t.text.replace('\n', '_'), c.width, c.offset)
                
                total = childs[0].half + offset + childs[-1].half
                half = total/2
                for c in childs:
                    c.offset -= half - childs[0].half
                
#                print "\n----B------"
#                print "childs[0].half: ", childs[0].half
#                print "Parent: %s" % node.t.text.replace('\n', '_')
#                for c in childs:
#                    print "  child: %s, width: %d, offset: %d" % (c.t.text.replace('\n', '_'), c.width, c.offset)
            else:
                childs[0].offset = 0
            self._getHLimits(node)
#            print "\n=====C========"
#            print "Parent: %s" % node.t.text.replace('\n', '_')
#            print " offset: %d, width: %d" % (node.offset, node.width)
#            for l, r in node.hLimits:
#                print "[%d, %d]" % (l, r)
        
    def _defaultNodeItem(self, node, y):
        """ If not createNodeItem is specified, this one will be used
        by default. 
        """
        nodeText = node.getName()
        textColor = 'black'
        if nodeText.startswith('Project'):
            textColor='white'
        
        return self.canvas.createTextbox(nodeText, 100, y, bgColor='light blue', textColor=textColor)
        
        
    def _paintNode(self, node, y):
        """ Paint a node of the graph.
        Params:
           canvas: the canvas in which to paint.
           node: the node of the graph to be painted.
           y: level in the tree where to paint.
        Returns:
           the create item in the canvas.
        """
        item = self.createNodeItem(node, y)
        node.width, node.height = item.getDimensions()
        node.half = node.width / 2
        node.hLimits = [[-node.half, node.half]]
        node.y = y
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
        for child in node.getChilds():
            count = 1
            #printHLimits(c, " child")
            
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
        #print "node: ", node.t.text, " x:", nx
        for c in node.getChilds():
            self._createEdges(c, nx)
            self.canvas.createEdge(node.item, c.item)


if __name__ == '__main__':
    from canvas import Canvas
    from pyworkflow.utils.graph import Graph
    
    root = tk.Tk()
    canvas = Canvas(root, width=400, height=400)
    canvas.grid(row=0, column=0, sticky='nsew')
    root.grid_columnconfigure(0, weight=1)
    root.grid_rowconfigure(0, weight=1)
    
    g = Graph(rootName='A')
    a = g.getRoot()
    b = g.createNode('B')
    c = g.createNode('C')
    a.addChild(b, c)
    
    lt = LevelTree(g)
    lt.paint(canvas)
    
    root.mainloop()             
