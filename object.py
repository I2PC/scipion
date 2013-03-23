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
        object.__init__(self)
        self.set(value)
        self.id =  args.get('id', None)
        self.parent_id =  args.get('parent_id', None)
        self.name =  args.get('name', '')
        self.tag =  args.get('tag', None) # True if the object serves as input to his parent
        self.store =  args.get('store', True) # True if this object will be stored from his parent
        self.pointer =  args.get('pointer', False) # True if will be treated as a reference for storage
        
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
    
    def getId(self):
        """Return object id"""
        return self.id
    
    def strId(self):
        """String representation of id"""
        return str(self.id)
    
    def __str__(self):
        """String representation of the scalar value"""
        return str(self.value)
        
    def hasValue(self):        
        return True
    
    def __eq__(self, other):
        """Comparison for scalars should be by value
        and for other objects by reference"""
        if self.value is None:
            return object.__eq__(other)
        return self.value == other.value
    
    def equalAttributes(self, other):
        """Compare that all attributes are equal"""
        for k, v in self.getAttributesToStore():
            v1 = getattr(self, k) # This is necessary because of FakedObject simulation of getattr
            v2 = getattr(other, k)
            if issubclass(type(v1), Object):
                comp = v1.equalAttributes(v2)
            else:
                comp = v1 == v2
            if not comp:
                return False
        return True
    

class OrderedObject(Object):
    """This is based on Object, but keep the list
    of the attributes to store in the same order
    of insertion, this can be useful where order matters"""
    def __init__(self, value=None, **args):
        object.__setattr__(self, '_attributes', [])
        Object.__init__(self, value, **args)
        
    def __setattr__(self, name, value):
        if not name in self._attributes and issubclass(value.__class__, Object) and value.store:
            self._attributes.append(name)
        Object.__setattr__(self, name, value)
    
    def getAttributesToStore(self):
        """Return the list of attributes than are
        subclasses of Object and will be stored"""
        for key in self._attributes:
            yield (key, getattr(self, key))
            
class FakedObject(Object):
    """This is based on Object, but will hide the set and get
    access to the attributes, they need to be defined with addAttribute"""
    def __init__(self, value=None, **args):
        object.__setattr__(self, '_attributes', {})
        Object.__init__(self, value, **args)
        
    def addAttribute(self, name, attrClass, **args):
        self._attributes[name] = attrClass(**args)
           
    def __setattr__(self, name, value):
        if name in self._attributes:
            if issubclass(type(value), Object):
                self._attributes[name] = value
            else:
                self._attributes[name].set(value)
        else:
            object.__setattr__(self, name, value)
    
    def __getattr__(self, name):
        if name in self._attributes:
            return self._attributes[name].get()
        return None
    
    def getAttributesToStore(self):
        """Return the list of attributes than are
        subclasses of Object and will be stored"""
        return self._attributes.iteritems()

                
class Scalar(Object):
    """Base class for basic types"""
    def hasValue(self):        
        return not self.value is None
    
    def equalAttributes(self, other):
        """Compare that all attributes are equal"""
        return self.value == other.value
    
    
class Integer(Scalar):
    """Integer object"""
    def convert(self, value):
        return int(value)
    
        
class String(Scalar):
    """String object"""
    def convert(self, value):
        return str(value)
    
        
class Float(Scalar):
    """Float object"""
    def convert(self, value):
        return float(value)
    
    
class Boolean(Scalar):
    """Boolean object"""
    def convert(self, value):
        t = type(value)
        if t is bool:
            return value
        if t is str:
            v = value.strip().lower()
            return v == 'true' or v == '1'
        return bool(value)    
    
    
class Pointer(Object):
    """Reference object to other one"""
    def __init__(self, value=None, **args):
        Object.__init__(self, value, pointer=True, **args)   
       

class List(Object, list):
    """Class to store a list of objects"""
    def __init__(self, **args):
        Object.__init__(self, **args)
        list.__init__(self)
        
    def __setattr__(self, name, value):
        if name.startswith('__item__') or len(name)==0:
            self.append(value)
        else:
            object.__setattr__(self, name, value)
            
    def __getattr__(self, name):
        if name.startswith('__item__'):
            index = int(name.split('__item__')[1]) - 1
            if index < len(self):
                return self[index]
        return None

    def getAttributesToStore(self):
        for key, attr in self.__dict__.iteritems():
            if issubclass(attr.__class__, Object) and attr.store:
                yield (key, attr)
        
        for i, item in enumerate(self):
            yield ("__item__%06d" % (i+1), item)
            
    def __str__(self):
        return list.__str__(self)
        
        
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
           
