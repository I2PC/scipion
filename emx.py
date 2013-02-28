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
from pydoc import classname
   
class EmxScalar(Object):
    '''Store a value and units.
    This class will behave like a simple type(Integer,Float,String)
    but can set units and hold any of these types as value'''
    def __init__(self, ElemType=Float, unit=None, value=None, **args):
        self.converter = ElemType(store=False)        
        Object.__init__(self, value, **args)        
        self.unit = String(unit, tag='attribute')
        
    def setUnit(self, unit):
        self.unit.set(unit)
        
    def getUnit(self):
        return self.unit.get()
    
    def set(self, value):
        self.converter.set(value)
        self.value = self.converter.value
        
    def convert(self, value):
        return self.converter.convert(value)
    
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
            setattr(self, k, EmxScalar(ElemType, value=args.get(k, None)))
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
    
    def getAttributesToStore(self):
        '''Return the list of attributes than are
        subclasses of Object and will be stored'''
        for key in self.keys:
            yield (key, getattr(self, key))
    
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
    
    def setValue(self, **args):
        '''Set values of different items'''
        for k, v in args.iteritems():
            if k in self.keys:
                getattr(self, k).set(v)
                
    def setAllValues(self, value):
        '''Set all items to the same value'''
        for k in self.keys:
            getattr(self, k).set(value)
            
                    
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
        
class Entity(Object):
    def __init__(self, **args):
        Object.__init__(self, **args)
        self.activeFlag = Integer()
        self.fom = Float()

    def hasValue(self):
        ''' An Entity must have always PK'''
        return True
    def __str__(self):
        className                = self.getClassName()
        partStr = "%s id = %s" % (className,str(self.id))
        for key, attr in self.getAttributesToStore():
            if attr.hasValue():
                partStr += '\n %s: %s' % (key, attr)
            if attr.pointer:
#                partStr += '\n %s %s' % (key, str(attr.id))
                partStr += '\n %s %s' % (key, str(attr.get().id))
        return partStr

class micrograph(Entity):
    def __init__(self, **args):
        Entity.__init__(self, **args)
        self.acceleratingVoltage = EmxScalar(Float, 'kV')
        self.cs = EmxScalar(Float, 'mm')
        self.pixelSpacing = PixelSpacing()
        self.defocusU = EmxScalar(Float,'nn')
        self.defocusV = EmxScalar(Float,'nn')
        self.defocusUAngle = EmxScalar(Float,'deg')
        self.amplitudeContrast = Float()
        
class particle(Entity):
    def __init__(self, **args):
        Entity.__init__(self, **args)
        self.boxSize = BoxSize()
        self.centerCoord = CenterCoord()
        self.defocusU = EmxScalar(Float,'nn')
        self.defocusV = EmxScalar(Float,'nn')
        self.defocusUAngle = EmxScalar(Float,'deg')
        self.pixelSpacing = PixelSpacing()
        self.transformationMatrix=TransformationMatrix()
        #define foreign keys
        self.micrograph = Pointer()

    def getMicrograph(self):
        return self.micrograph.get()

    def setMicrograph(self, mic):
        self.micrograph.set(mic)

class EmxData():
    ''' Class to group EMX objects'''
    def __init__(self):
        self.objLists = {'micrograph' : [], 
                          'particle'  : []}

    def addObject(self,obj):
        className                = obj.getClassName()
        self.objLists[className].append(obj)
        
    def getObjectwithID(self, className,id):
        for obj in self.objLists[className]:
            if (obj.id == id):
                return obj

    def clean(self):
        for list in self.objLists:
            objLists[list]=[]

    def __str__(self):
        partStr=""
        for k,v in self.objLists.iteritems():
            if len(v):
                partStr += "\n****\n%sS\n****\n"% k.upper()
            for obj in v:
                partStr += obj.__str__()+"\n"
        return partStr