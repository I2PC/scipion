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


class GraphLayout(object):
    """ Base class for all algorithm that implement
    functions to organize a graph in a plane.
    """
    
    def draw(self, graph, **kwargs):
        """ Setup the nodes position in the plane. """
        pass


class BasicLayout(GraphLayout):
    """ This layout will keep node position as much as possible.
    It will try to allocate the nodes with x=0 and y=0.
    """
    
    def __init__(self):
        GraphLayout.__init__(self)
        self.DY = 65
        self.DX = 15
        
    def draw(self, graph, **kwargs):
        """ Organize nodes of the graph in the plane.
        Nodes should have: x, y, width and height attributes
        x and y will be modified.
        """
        for node in graph.getNodes():
            if hasattr(node, 'x') and hasattr(node, 'y'):
                if getattr(node, 'x', 0) == 0 or node.y == 0:
                    self._drawNode(node)
                
    def _drawNode(self, node):
        """ Allocate node with x=0 and y=0. """

        parents = node.getParents()
        if not parents:
            print "EMPTY NODE ask JM"
            return
        maxParent = parents[0]
        
        for p in parents[1:]:
            if p.y > maxParent.y:
                maxParent = p
                
        siblings = maxParent.getChilds()
        
        if len(siblings) == 1:
            node.x = maxParent.x
            node.y = maxParent.y + self.DY
        else:
            rightSibling = siblings[0]
            for s in siblings:
                if s.x > rightSibling.x:
                    rightSibling = s
            node.x = rightSibling.x + rightSibling.width/2 + self.DX + node.width/2
            node.y = rightSibling.y
        

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
        
        # Setup empty layout for each node
        for node in graph.getNodes():
            node._layout = {}
            
        # Do some level initialization on each node
        self._setLayoutLevel(rootNode,  1, None)
        self._computeNodeOffsets(rootNode, 1)
        # Compute extreme left limit
        m = 9999
        for left, _ in rootNode._layout['hLimits']:
            m = min(m, left)
        
        self._applyNodeOffsets(rootNode, -m + self.DY)
        
        # Clean temporarly _layout attributes
        for node in graph.getNodes():
            del node._layout
        
    def _setLayoutLevel(self, node, level, parent):
        """ Iterate over all nodes and set _layout dict.
        Also set the node level, which is defined
        as the max level of a parent + 1
        """
        if level > self.maxLevel:
            return 
        
        layout = node._layout
        
        if level > layout.get('level', 0):
            # Calculate the y-position depending on the level
            # and the deltha-Y (DY)
            node.y = level * self.DY
            layout['level'] = level
            layout['parent'] = parent
            half = node.width / 2
            layout['half'] = half
            layout['hLimits'] = [[-half, half]]
            layout['offset'] = 0
    
            if self.__isNodeExpanded(node):
                for child in node.getChilds():
                    self._setLayoutLevel(child, level+1, node)
    
    def __isNodeExpanded(self, node):
        """ Check if the status of the node is expanded or collapsed. """
        return getattr(node, 'expanded', True)
        
    def __setNodeOffset(self, node, offset):
        node._layout['offset'] = offset
        
    def __getNodeHalf(self, node):
        return node._layout['half']
    
    def __getNodeChilds(self, node):
        """ Return the node's childs that have been 
        visited by this node first (its 'parent')
        """
        if self.__isNodeExpanded(node):
            return [c for c in node.getChilds() if c._layout['parent'] is node]
        else:
            return [] # treat collapsed nodes as if they have no childs
        
    def _computeNodeOffsets(self, node, level):
        """ Position a parent node and its childs.
        Only this sub-tree will be considered at this point.
        Then it will be adjusted with node siblings.
        """
        if level > self.maxLevel:
            return
        
        childs = self.__getNodeChilds(node)
        n = len(childs)

        if n > 0:
            for c in childs:
                self._computeNodeOffsets(c, level + 1)
                
            if n > 1:
                offset = 0
                # Keep right limits to compute the separation between siblings
                # some times it not enough to compare with the left sibling
                # for some child levels of the node
                rightLimits = [r for l, r in childs[0]._layout['hLimits']]
                
                for i in range(n-1):
                    sep = self._getChildsSeparation(childs[i], childs[i+1], rightLimits)
                    offset += sep
                    c = childs[i+1]
                    self.__setNodeOffset(c, offset)
                
                half0 = self.__getNodeHalf(childs[0])
                total = half0 + offset + self.__getNodeHalf(childs[-1])
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
                
    def _getChildsSeparation(self, child1, child2, rightLimits):
        '''Calcualte separation between siblings
        at each height level'''
        sep = 0 
        hL2 = child2._layout['hLimits']
        n1 = len(rightLimits)
        n2 = len(hL2)
        h = min(n1, n2)
            
        for i in range(h):
            right = rightLimits[i]
            left = hL2[i][0]            
            if left + sep < right:
                sep = right - left
            rightLimits[i] = hL2[i][1]
                
        if n1 > n2:
            # If there are more levels in the rightLimits
            # updated the last ones like if they belong 
            # to next sibling is is now (sep + self.DX) away
            for i in range(h, n1):
                rightLimits[i] -= sep + self.DX
        else:
            # If the current right sibling has more levels
            # just add them to the current rightLimits
            for i in range(h, n2):
                rightLimits.append(hL2[i][1])
                
        return sep + self.DX
    
    def _applyNodeOffsets(self, node, x):
        """ Adjust the x-position of the nodes by applying the offsets.
        """
        if node._layout['level'] == self.maxLevel:
            return 
        
        layout = node._layout
        node.x = x + layout['offset']
        
        childs = self.__getNodeChilds(node)
        
        for child in childs:
            self._applyNodeOffsets(child, node.x)

