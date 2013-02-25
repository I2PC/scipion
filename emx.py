# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Roberto Marabini       (roberto@cnb.csic.es)
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

from object import Object, Float, Integer, String, Pointer
   
class EmxScalar(Object):
    '''Store a value and units'''
    def __init__(self, ElemType=Float, unit=None, value=None, **args):
        self.value = ElemType(store=False)
        Object.__init__(self, value, **args)        
        self.unit = String(unit, tag='attribute')
        
    def setUnit(self, unit):
        self.unit.set(unit)
        
    def getUnit(self):
        return self.unit.get()
        
    def hasValue(self):    
        return self.value.hasValue()
    
    def set(self, value):
        self.value.set(value)
    
class EmxVector(Object):
    '''This class will implement a list of named elements.
    It will serve for vectors or matrices.
    The elements should have a name like: (X, Y, Z)'''
    def __init__(self, ElemType=Float, **args):
        Object.__init__(self, **args)
        # Get the named list if exists, by default X,Y and Z
        self.keys = getattr(self, 'keys', ['X', 'Y', 'Z'])
        # Build items
        for k in self.keys:
            setattr(self, k, EmxScalar(ElemType, args.get(k, None)))
        # Set unit for all items if comes as argument
        if 'unit' in args:
            self.setUnits(args.get('unit'))
            
    def setUnit(self, varname, unit):
        getattr(self, varname).setUnit(unit)
        
    def setUnits(self, unit):
        for k in self.keys:
            self.setUnit(k, unit)
        
    def __str__(self):
        slots = [''] * len(self.keys)
        for i, k in enumerate(self.keys):
            item = getattr(self, k)
            if item.hasValue():
                slots[i] = '%s' % str(item.get())
        return '(%s)' % ','.join(slots)

    def hasValue(self):
        for k in self.keys:
            if getattr(self, k).hasValue():
                return True
        return False
    
    
class EmxMatrix(EmxVector):
    '''Class based on vector but arranges items in a Matrix'''
    def __init__(self, ElemType=Float, xdim=3, **args):
        prefix = args.get('prefix', 't')
        self.xdim = xdim
        self.ydim = args.get('ydim', xdim)
        self.keys = []
        for i in range(1, self.ydim+1):
            for j in range(1, self.xdim+1):
                self.keys.append('%(prefix)s%(i)d%(j)d' % locals())
        EmxVector.__init__(self, ElemType, **args)  
                    
    def __str__(self):
        vectorStr = EmxVector.__str__(self)
        matrixStr = ""
        counter = 0
        for c in vectorStr:
            if c == ',':
                counter += 1
                if counter % self.xdim == self.xdim - 1:
                    counter = 0
                    c = '\n'
            matrixStr += c
        return matrixStr 
            
                    
class TransformationMatrix(EmxMatrix):
    def __init__(self, ElemType=Float, **args):
        EmxMatrix.__init__(self, ElemType, prefix='t', xdim=4, **args)
        self.setUnit('t14', 'A')
        self.setUnit('t24', 'A')
        self.setUnit('t34', 'A')


class PixelSpacing(EmxVector):
    def __init__(self, **args):
        EmxVector.__init__(self, Float, unit='A/px', **args)
        
        
class BoxSize(EmxVector):
    def __init__(self, **args):
        EmxVector.__init__(self, Integer, unit='px', **args)
        
        
class CenterCoord(EmxVector):
    def __init__(self, **args):
        EmxVector.__init__(self, Float, unit='px', **args)
        
    
class Micrograph(Object):
    def __init__(self, **args):
        Object.__init__(self, **args)
        self.acceleratingVoltage = EmxScalar(Float, 'kV')
        self.activeFlag = Integer(1)
        self.cs = EmxScalar(Float, 'mm')
        self.pixelSpacing = PixelSpacing()
        self.defocusU = EmxScalar(Float,'nn')
        self.defocusV = EmxScalar(Float,'nn')
        self.defocusUAngle = EmxScalar(Float,'deg')
        self.amplitudeContrast = Float()
        self.fom = Float()

    def hasValue(self):
        ''' A micrograph must have always PK'''
        return True

#    def set(self,micrographFK):
#        '''assign Primary Key to Foreign Key'''
#        self.id=micrographFK.id
        
class Particle(Object):
    def __init__(self, **args):
        Object.__init__(self, **args)
        self.activeFlag = Integer(1)
        self.boxSize = BoxSize()
        self.centerCoord = CenterCoord()
        self.defocusU = EmxScalar(Float,'nn')
        self.defocusV = EmxScalar(Float,'nn')
        self.defocusUAngle = EmxScalar(Float,'deg')
        self.fom = Float()
        self.pixelSpacing = PixelSpacing()
        self.transformationMatrix=TransformationMatrix()
        #define foreign keys
        self.micrograph = Pointer()
            
#    def __str__(self):
#        return "File %s\n Particles: %d" % (self.Path.get(), self.N.get())             
