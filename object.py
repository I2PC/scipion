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
# *  e-mail address 'xmipp@cnb.csic.es'
# *
# **************************************************************************

class Object():
    ''' All objects in our Domain should inherit from this class
    that will contains all base properties'''
    def __init__(self, value=None, **args):
        self.set(value)
        self.id = args.get('id', None)
        self.parent_id = args.get('parent_id', None)
        self.name = args.get('name', '')
        self.tag = args.get('tag', None) # True if the object serves as input to his parent
        self.store = args.get('store', True) # True if this object will be stored from his parent
        self.reference = args.get('reference', False) # True if will be treated as a reference for storage
        
    def getClassName(self):
        return self.__class__.__name__
    
    def getAttributesToStore(self):
        '''Return the list of attributes than are
        subclasses of Object and will be stored'''
        for key, attr in self.__dict__.iteritems():
            if issubclass(attr.__class__, Object) and attr.store:
                yield (key, attr)
    
    def convert(self, value):
        '''Convert a value to desired scalar type'''
        return value
    
    def setId(self,obj):
        print "kkkk", obj.id,
        self.id=obj.id
        print self.id
        
        
    def set(self, value):
        '''Set the internal value, if it is different from None
        call the convert function in subclasses'''
        if not value is None:
            value = self.convert(value)            
        self.value = value
    
    def get(self):
        '''Return internal value'''
        return self.value
    
    def __str__(self):
        '''String representation of the scalar value'''
        return str(self.value)
        
    def hasValue(self):
        return not self.value is None
    
        
class Integer(Object):
    '''Integer object'''
    def convert(self, value):
        return int(value)
    
        
class String(Object):
    '''String object'''
    def convert(self, value):
        return str(value)
    
        
class Float(Object):
    '''Float object'''
    def convert(self, value):
        return float(value)
    
    
class Boolean(Object):
    '''Boolean object'''
    def convert(self, value):
        return bool(value)    
       

class Array(Object):
    '''Class for holding fixed len array'''
    def __init__(self, size=10, **args):
        Object.__init__(self, size, **args)
        
    def set(self, size):
        '''Set the array size'''
        self.value = int(size)  
        for i in range(int(size)):
            self.__setitem__(i, None)                 
        
    def strIndex(self, i):
        return 'item_%04d' % i
    
    def __setitem__(self, index, value):
        self.__dict__[self.strIndex(index)] = value
        
    def __getitem__(self, index):
        return self.__dict__[self.strIndex(index)]
    
    def __len__(self):
        return self.value

    
class Vector3D(Object):
    def __init__(self, ElemType=Float, x=None, y=None, z=None, **args):
        Object.__init__(self, **args)
        self.X = ElemType(x)
        self.Y = ElemType(y)     
        self.Z = ElemType(z)
            
    def setUnit(self, varname, unit):
        getattr(self, varname).id = {'unit': unit}     
        
    def __str__(self):
        slots = ['']*3
        for i, v in enumerate([self.X, self.Y, self.Z]):
            if v.hasValue():
                slots[i] = '%f' % v.get()
        return '(%s)' % ','.join(slots)

    def hasValue(self):
        return \
               self.X.hasValue() \
            or self.Y.hasValue() \
            or self.Z.hasValue()
            
class TransformationMatrix(Object):
    def __init__(self, ElemType=Float, t11=None, t12=None, t13=None, t14=None,\
                                       t21=None, t22=None, t23=None, t24=None,\
                                       t31=None, t32=None, t33=None, t34=None, **args):
        Object.__init__(self, **args)
        self.t11 = ElemType(t11)
        self.t12 = ElemType(t12)     
        self.t13 = ElemType(t13)
        self.t14 = ElemType(t14)
        self.setUnit('t14', 'A')

        self.t21 = ElemType(t21)
        self.t22 = ElemType(t22)     
        self.t23 = ElemType(t23)
        self.t24 = ElemType(t24)
        self.setUnit('t24', 'A')

        self.t31 = ElemType(t31)
        self.t32 = ElemType(t32)     
        self.t33 = ElemType(t33)
        self.t34 = ElemType(t34)
        self.setUnit('t34', 'A')

    def setUnit(self, varname, unit):
        getattr(self, varname).id = {'unit': unit}     
        
    def __str__(self):
        return (('%f %f %f %f\n%f %f %f %f\n%f %f %f %f') % \
        (self.t11.get(), self.t12.get(), self.t13.get(), self.t14.get(),\
         self.t21.get(), self.t22.get(), self.t23.get(), self.t24.get(),\
         self.t31.get(), self.t32.get(), self.t33.get(), self.t34.get()))

    def hasValue(self):
        return \
               self.t11.hasValue() \
            or self.t12.hasValue() \
            or self.t13.hasValue() \
            or self.t14.hasValue() \
            or self.t21.hasValue() \
            or self.t22.hasValue() \
            or self.t23.hasValue() \
            or self.t24.hasValue() \
            or self.t31.hasValue() \
            or self.t32.hasValue() \
            or self.t33.hasValue() \
            or self.t34.hasValue() 

class PixelSpacing(Vector3D):
    def __init__(self, x=None, y=None, z=None, **args):
        Vector3D.__init__(self, Float, x, y, z, **args)
        self.setUnit('X', 'A/px')
        self.setUnit('Y', 'A/px')
        self.setUnit('Z', 'A/px')
        
class BoxSize(Vector3D):
    def __init__(self, x=None, y=None, z=None, **args):
        Vector3D.__init__(self, Integer, x, y, z, **args)
        self.setUnit('X', 'px')
        self.setUnit('Y', 'px')
        self.setUnit('Z', 'px')
        
class CenterCoord(Vector3D):
    def __init__(self, x=None, y=None, z=None, **args):
        Vector3D.__init__(self, Float, x, y, z, **args)
        self.setUnit('X', 'px')
        self.setUnit('Y', 'px')
        self.setUnit('Z', 'px')
        
    
class Micrograph(Object):
    def __init__(self, **args):
        Object.__init__(self, **args)
        self.acceleratingVoltage = Float(id={'unit':'kV'})
        self.activeFlag = Integer(1)
        self.cs = Float(id={'unit': 'mm'})
        self.pixelSpacing = PixelSpacing()
        self.defocusU = Float(id={'unit':'nn'})
        self.defocusV = Float(id={'unit':'nn'})
        self.defocusUAngle = Float(id={'unit':'deg'})
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
        self.defocusU = Float(id={'unit':'nn'})
        self.defocusV = Float(id={'unit':'nn'})
        self.defocusUAngle = Float(id={'unit':'deg'})
        self.fom = Float()
        self.pixelSpacing = PixelSpacing()
        self.transformationMatrix=TransformationMatrix()
        #define foreign keys
        self.micrographFK = Micrograph()
            
#    def __str__(self):
#        return "File %s\n Particles: %d" % (self.Path.get(), self.N.get())             
