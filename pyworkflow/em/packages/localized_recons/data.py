# **************************************************************************
# *
# * Authors:     Josue Gomez Blanco (jgomez@cnb.csic.es)
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
from math import sqrt
from pyworkflow.object import OrderedObject, Float, CsvList, Set

class Vector(OrderedObject):
    """Define Vector class and method to obtain individual vectors from lists
       with 3 values"""
    """Base object for all EM classes"""
    
    def __init__(self, vector=[0,0,1], dist=0, **args):
        OrderedObject.__init__(self, **args)
        self._distance = Float(dist)
        self._vector = CsvList(pType=float)
        self._vector.set(vector)
    
    def setVector(self, v):
        self._vector.set(v)
    
    def getVector(self):
        return self._vector
    
    def setDistance(self, d):
        if d <= 0:
            d = self.length()
        self._distance.set(d)
    
    def getDistance(self):
        return self._distance.get()
    
    def xAxis(self):
        return self.getVector()[0]
    
    def yAxis(self):
        return self.getVector()[1]
    
    def zAxis(self):
        return self.getVector()[2]
    
    def length(self):
        [x, y, z] = self.getVector()
        return sqrt(x**2 + y**2 + z**2)
    
    def normalize(self):
        if self.length() > 0:
            length = self.length()
            v = [x / length for x in self.getVector()]
        else:
            v = [0,0,0]
        self.setVector(v)


class SetOfVectors(Set):
    """ Represents a set of Images """
    ITEM_TYPE = Vector
    
    def __init__(self, **kwargs):
        Set.__init__(self, **kwargs)

