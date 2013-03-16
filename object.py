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
This modules holds the base classes for the ORM implementation.
The Object class is the root in the hierarchy and some other
basic classes.
"""

class Object(object):
    """ All objects in our Domain should inherit from this class
    that will contains all base properties"""
    def __init__(self, value=None, **args):
        self.set(value)
        self.id        = args.get('id', None)
        self.parent_id = args.get('parent_id', None)
        self.name      = args.get('name', '')
        self.tag       = args.get('tag', None) # True if the object serves as input to his parent
        self.store     = args.get('store', True) # True if this object will be stored from his parent
        self.pointer   = args.get('pointer', False) # True if will be treated as a reference for storage
        
    def getClassName(self):
        return self.__class__.__name__
    
    def getAttributesToStore(self):
        """Return the list of attributes than are
        subclasses of Object and will be stored"""
        for key, attr in self.__dict__.iteritems():
            if issubclass(attr.__class__, Object) and attr.store:
                yield (key, attr)
                
    def isPointer(self):
        """If this is true, the value field is a pointer 
        to anothor object"""
        return self.pointer
    
    def convert(self, value):
        """Convert a value to desired scalar type"""
        return value
    
    def set(self, value):
        """Set the internal value, if it is different from None
        call the convert function in subclasses"""
        if not value is None:
            value = self.convert(value)            
        self.value = value
    
    def get(self):
        """Return internal value"""
        return self.value
    
    def __str__(self):
        """String representation of the scalar value"""
        return str(self.value)
        
    def hasValue(self):        
        return not self.value is None
    
    def __eq__(self, other):
        """Comparison for scalars should be by value
        and for other objects by reference"""
        if self.value is None:
            return object.__eq__(other)
        return self.value == other.value
    
        
class Integer(Object):
    """Integer object"""
    def convert(self, value):
        return int(value)
    
        
class String(Object):
    """String object"""
    def convert(self, value):
        return str(value)
    
        
class Float(Object):
    """Float object"""
    def convert(self, value):
        return float(value)
    
    
class Boolean(Object):
    """Boolean object"""
    def convert(self, value):
        return bool(value)    
    
    
class Pointer(Object):
    """Reference object to other one"""
    def __init__(self, value=None, **args):
        Object.__init__(self, value, pointer=True, **args)    
       

class Array(Object):
    """Class for holding fixed len array"""
    def __init__(self, size=10, **args):
        Object.__init__(self, size, **args)
        
    def set(self, size):
        """Set the array size"""
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
           
