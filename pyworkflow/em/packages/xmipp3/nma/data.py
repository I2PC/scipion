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
Define some classes to store Data points for clustering.
"""


class Point():
    """ Return x, y 2d coordinates and some other properties
    such as weight and state.
    """
    # Selection states
    DISCARDED = -1
    NORMAL = 0
    SELECTED = 1
    
    def __init__(self, pointId, data, weight, state=0):
        self._id = pointId
        self._data = data
        self._weight = weight
        self._state = state
        self._container = None
        
    def getId(self):
        return self._id
    
    def getX(self):
        return self._data[self._container.XIND]
    
    def setX(self, value):
        self._data[self._container.XIND] = value
    
    def getY(self):
        return self._data[self._container.YIND]
    
    def setY(self, value):
        self._data[self._container.YIND] = value
        
    def getZ(self):
        return self._data[self._container.ZIND]
        
    def setZ(self, value):
        self._data[self._container.ZIND] = value    
    
    def getWeight(self):
        return self._weight
    
    def getState(self):
        return self._state
    
    def setState(self, newState):
        self._state = newState
        
    def eval(self, expression):
        localDict = {}
        for i, x in enumerate(self._data):
            localDict['x%d' % (i+1)] = x
        return eval(expression, {"__builtins__":None}, localDict)

    def setSelected(self):
        self.setState(Point.SELECTED)
        
    def isSelected(self):
        return self.getState()==Point.SELECTED    
    
    def setDiscarded(self):
        self.setState(Point.DISCARDED)
           
    def isDiscarded(self):
        return self.getState()==Point.DISCARDED
    
    def getData(self):
        return self._data 

    
class Data():
    """ Store data points. """
    def __init__(self, **kwargs):
        # Indexes of data
        self._dim = kwargs.get('dim') # The points dimensions
        self.clear()
        
    def addPoint(self, point, position=None):
        point._container = self
        if position is None:
            self._points.append(point)
        else:
            self._points.insert(position, point)
            
    def getPoint(self, index):
        return self._points[index]
        
    def __iter__(self):
        for point in self._points:
            if not point.isDiscarded():
                yield point
                
    def iterAll(self):
        """ Iterate over all points, including the discarded ones."""
        return iter(self._points)
            
    def getXData(self):
        return [p.getX() for p in self]
    
    def getYData(self):
        return [p.getY() for p in self]
    
    def getZData(self):
        return [p.getZ() for p in self]
    
    def getWeights(self):
        return [p.getWeight() for p in self]
    
    def getSize(self):
        return len(self._points)
    
    def getSelectedSize(self):
        return len([p for p in self if p.isSelected()])
    
    def getDiscardedSize(self):
        return len([p for p in self.iterAll() if p.isDiscarded()])
    
    def clear(self):
        self.XIND = 0
        self.YIND = 1
        self.ZIND = 2
        self._points = []


class PathData(Data):
    """ Just contains two list of x and y coordinates. """
    
    def __init__(self, **kwargs):
        Data.__init__(self, **kwargs)
    
    def splitLongestSegment(self):
        """ Split the longest segment by adding the midpoint. """
        maxDist = 0
        n = self.getSize()
        # Find the longest segment and its index
        for i in range(n-1):
            p1 = self.getPoint(i)
            x1, y1 = p1.getX(), p1.getY()
            p2 = self.getPoint(i+1)
            x2, y2 = p2.getX(), p2.getY()
            dist = (x1-x2)**2 + (y1-y2)**2
            if dist > maxDist:
                maxDist = dist
                maxIndex = i+1
                midX = (x1+x2)/2
                midY = (y1+y2)/2
        # Add a midpoint to it
        point = self.createEmptyPoint()
        point.setX(midX)
        point.setY(midY)
        self.addPoint(point, position=maxIndex)
        
    def createEmptyPoint(self):
        data = [0.] * self._dim # create 0, 0...0 point
        point = Point(0, data, 0)
        point._container = self
        
        return point
    
    def removeLastPoint(self):
        del self._points[-1]
