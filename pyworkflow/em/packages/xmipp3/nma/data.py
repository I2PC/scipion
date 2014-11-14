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
Define some classes to store Data points for clustering.
"""


class Point():
    """ Return x, y 2d coordinates and some other properties
    such as weight and state.
    """
    # Indexes of data
    XIND = 0
    YIND = 1
    ZIND = 2
    # Selection states
    DISCARDED = -1
    NORMAL = 0
    SELECTED = 1
    
    def __init__(self, pointId, data, weight, state=0):
        self._id = pointId
        self._data = data
        self._weight = weight
        self._state = state
        
    def getId(self):
        return self._id
    
    def getX(self):
        return self._data[Point.XIND]
    
    def getY(self):
        return self._data[Point.YIND]
    
    def getZ(self):
        return self._data[Point.ZIND]
    
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
    
    
class Data():
    """ Store data points. """
    def __init__(self):
        self._points = []
        
    def addPoint(self, point):
        self._points.append(point)
        
    def __iter__(self):
        for point in self._points:
            yield point
            
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
        return len([p for p in self if p.getState()==Point.SELECTED])
    
    
    

            
