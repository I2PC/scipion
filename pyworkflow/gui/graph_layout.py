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


class GraphLayout(object):
    """ Base class for all algorithm that implement
    functions to organize a graph in a plane.
    """
    
    def draw(self, graph, **kwargs):
        """ Setup the nodes position in the plane. """
        pass


class LevelTreeLayout(GraphLayout):
    """ Organize the nodes of the graph by levels.
    It will recursively organize childs and then
    fit two sibling trees. """
    
    def __init__(self):
        GraphLayout.__init__(self)
        self.DY = 65
        self.DX = 15
        self.maxLevel = 9999
        
    def draw(self, graph, **kwargs):
        """ Organize nodes of the graph in the plane.
        Nodes should have: x, y, width and height attributes
        x and y will be modified.
        """
        rootNode = graph.getRoot()
        
        # Do some level initialization on each node
        self._setLayoutLevel(rootNode,  1, None)
        self._computeNodeOffsets(rootNode, 1)
        # Compute extreme left limit
        m = 9999
        for left, _ in rootNode._layout['hLimits']:
            m = min(m, left)
        
        self._applyNodeOffsets(rootNode, -m + self.DY)
        
    
    def _setLayoutLevel(self, node, level, parent):
        """ Iterate over all nodes and set _layout dict.
        Also set the node level, which is defined
        as the max level of a parent + 1
        """
        if level > self.maxLevel:
            return 
        
        if not hasattr(node, '_layout'):
            node._layout = {}
            
        # Calculate the y-position depending on the level
        # and the deltha-Y (DY)
        node.y = level * self.DY
        
        #print "DEBUG: setLayoutLevel: node: ", node.getName()
        #print "DEBUG: setLayoutLevel: parent: ", parent.getName() if parent else "None"
        
        
        layout = node._layout
        if level > layout.get('level', 0):
            #print "DEBUG: setLayoutLevel:   updating level to: ", level
            layout['level'] = level
            layout['parent'] = parent
            half = node.width / 2
            layout['half'] = half
            layout['hLimits'] = [[-half, half]]
            layout['offset'] = 0
    
            for child in node.getChilds():
                self._setLayoutLevel(child, level+1, node)
    
    def __setNodeOffset(self, node, offset):
        node._layout['offset'] = offset
        
    def __getNodeHalf(self, node):
        return node._layout['half']
    
    def __getNodeChilds(self, node):
        """ Return the node's childs that have been 
        visited by this node first (its 'parent')
        """
        return [c for c in node.getChilds() if c._layout['parent'] is node]
        
    def _computeNodeOffsets(self, node, level):
        """ Position a parent node and its childs.
        Only this sub-tree will be considered at this point.
        Then it will be adjusted with node siblings.
        """
        if level > self.maxLevel:
            return
        
        #print "DEBUG: _computeNodeOffsets: node", node.getName()
        
        childs = self.__getNodeChilds(node)
        n = len(childs)

        if n > 0:
            #print "DEBUG: _computeNodeOffsets: childs: ", n
            for c in childs:
                self._computeNodeOffsets(c, level + 1)
                
            if n > 1:
                offset = 0
                for i in range(n-1):
                    sep = self._getChildsSeparation(childs[i], childs[i+1])
                    offset += sep
                    c = childs[i+1]
                    self.__setNodeOffset(c, offset)
                
                half0 = self.__getNodeHalf(childs[0])
                total = half0 + offset + self.__getNodeHalf(childs[1])
                half = total / 2
                for c in childs:
                    self.__setNodeOffset(c, c._layout['offset'] - half + half0)
            else:
                self.__setNodeOffset(childs[0], 0)
            
            self._computeHLimits(node)      
    
    def _computeHLimits(self, node):
        """ This function will traverse the tree
        from node to build the left and right profiles(hLimits)
        for each level of the tree
        """
        layout = node._layout
        hLimits = layout['hLimits']

        childs = self.__getNodeChilds(node)
        
        for child in childs:
            count = 1
            layout = child._layout
            for l, r in layout['hLimits']:
                l += layout['offset']
                r += layout['offset']

                if count < len(hLimits):
                    if l < hLimits[count][0]:
                        hLimits[count][0] = l
                    if r > hLimits[count][1]:
                        hLimits[count][1] = r
                else:
                    hLimits.append([l, r])
                count += 1
                
    def _getChildsSeparation(self, child1, child2):
        '''Calcualte separation between siblings
        at each height level'''
        sep = 0 
        hL1 = child1._layout['hLimits']
        hL2 = child2._layout['hLimits']
        n1 = len(hL1)
        n2 = len(hL2)
        h = min(n1, n2)
            
        for i in range(h):
            right = hL1[i][1]
            left = hL2[i][0]            
            if left + sep < right:
                sep = right - left                
  
        return sep + self.DX
    
    def _applyNodeOffsets(self, node, x):
        """ Adjust the x-position of the nodes by applying the offsets.
        """
        if node._layout['level'] == self.maxLevel:
            return 
        
        #print "DEBUG: _applyNodeOffsets: node", node.getName()
        
        layout = node._layout
        #print "DEBUG: _applyNodeOffsets:     bofore node.x: ", node.x
        node.x = x + layout['offset']
        
        childs = self.__getNodeChilds(node)
        
        #print "DEBUG: _applyNodeOffsets:     after node.x: ", node.x
        #print "DEBUG: _applyNodeOffsets:     childs: ", len(childs)
        
        for child in childs:
            self._applyNodeOffsets(child, node.x)

