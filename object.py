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
        self._objId =  args.get('objId', None)
        self._objParentId =  args.get('objParentId', None)
        self._objName =  args.get('objName', '')
        self._objTag =  args.get('objTag', None) # True if the object serves as input to his parent
        self._objDoStore =  args.get('objDoStore', True) # True if this object will be stored from his parent
        self._objIsPointer =  args.get('objIsPointer', False) # True if will be treated as a reference for storage
        self._objCreationTime = None
        
    def getClassName(self):
        return self.__class__.__name__
    
    def getAttributesToStore(self):
        """Return the list of attributes than are
        subclasses of Object and will be stored"""
        for key, attr in self.__dict__.iteritems():
            if issubclass(attr.__class__, Object) and attr._objDoStore:
                yield (key, attr)
                
    def isPointer(self):
        """If this is true, the value field is a pointer 
        to anothor object"""
        return self._objIsPointer
    
    def _convertValue(self, value):
        """Convert a value to desired scalar type"""
        return value
    
    def set(self, value):
        """Set the internal value, if it is different from None
        call the convert function in subclasses"""
        if not value is None:
            value = self._convertValue(value)            
        self._objValue = value
    
    def get(self):
        """Return internal value"""
        return self._objValue
    
    def getInternalValue(self):
        """Return the internal value for storage.
        This is a good place to do some update of the
        internal value before been stored"""
        return self._objValue
    
    def getObjId(self):
        """Return object id"""
        return self._objId
    
    def setObjId(self, newId):
        """Set the object id"""
        self._objId = newId
        
    def hasObjId(self):
        return not self._objId is None
    
    def strId(self):
        """String representation of id"""
        return str(self._objId)
    
    def getName(self):
        return self._objName
    
    def getLastName(self):
        """ If the name contains parent path, remove it
        and only return the attribute name in its parent. 
        """
        if '.' in self._objName:
            return self._objName.split('.')[-1]
        return self._objName 
    
    def setName(self, name):
        self._objName = name
        
    def hasValue(self):        
        return True
    
    def getStore(self):
        """Return True if the object will be stored by the mapper"""
        return self._objDoStore
    
    def setStore(self, value):
        """set the store flag"""
        self._objDoStore = value
    
    def __eq__(self, other):
        """Comparison for scalars should be by value
        and for other objects by reference"""
        if self._objValue is None:
            return object.__eq__(other)
        return self._objValue == other._objValue
    
    def equalAttributes(self, other):
        """Compare that all attributes are equal"""
        for k, _ in self.getAttributesToStore():
            v1 = getattr(self, k) # This is necessary because of FakedObject simulation of getattr
            v2 = getattr(other, k)
            if issubclass(type(v1), Object):
                comp = v1.equalAttributes(v2)
            else:
                comp = v1 == v2
            if not comp:
                return False
        return True
    
    def printAll(self, name=None, level=0):
        """Print object and all its attributes.
        Main for debugging"""
        tab = ' ' * (level*3)
        if name is None:
            print tab, self.getClassName()
        else:
            print tab, '%s = %s' % (name, self._objValue)
        for k, v in self.getAttributesToStore():
            v.printAll(k, level + 1)
            
    def copyAttributes(self, other, *attrNames):
        """ Copy attributes in attrNames from other to self. 
        If the name X is in attrNames, it would be equivalent to:
        self.X.set(other.X.get())
        """
        for name in attrNames:
            getattr(self, name).set(getattr(other, name).get())
    

class OrderedObject(Object):
    """This is based on Object, but keep the list
    of the attributes to store in the same order
    of insertion, this can be useful where order matters"""
    def __init__(self, value=None, **args):
        object.__setattr__(self, '_attributes', [])
        Object.__init__(self, value, **args)
        
    def __setattr__(self, name, value):
        if not name in self._attributes and issubclass(value.__class__, Object) and value._objDoStore:
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
            attr = self._attributes[name]
            if issubclass(type(attr), Scalar):
                return attr.get()
            else:
                return attr
        return None
    
    def getAttributesToStore(self):
        """Return the list of attributes than are
        subclasses of Object and will be stored"""
        return self._attributes.iteritems()

                
class Scalar(Object):
    """Base class for basic types"""
    def hasValue(self):        
        return self._objValue is not None
    
    def equalAttributes(self, other):
        """Compare that all attributes are equal"""
        return self._objValue == other._objValue
    
    def __str__(self):
        """String representation of the scalar value"""
        return str(self._objValue)
    
    def __eq__(self, other):
        """Comparison for scalars should be by value
        and for other objects by reference"""
        if isinstance(other, Object):
            return self._objValue == other._objValue
        return self._objValue == other

    def __ne__(self, other):
        return not self.__eq__(other)
    
    def __cmp__(self, other):
        """ Comparison implementation for scalars. """
        if isinstance(other, Object):
            return cmp(self._objValue, other._objValue)
        return cmp(self._objValue, other)        
       
    def get(self, default=None):
        """Get the value, if internal value is None
        the default argument passed is returned"""
        if self.hasValue():
            return self._objValue
        return default
        
    
class Integer(Scalar):
    """Integer object"""
    def _convertValue(self, value):
        return int(value)
    
        
class String(Scalar):
    """String object"""
    def _convertValue(self, value):
        return str(value)
    
        
class Float(Scalar):
    """Float object"""
    def _convertValue(self, value):
        return float(value)
    
    
class Boolean(Scalar):
    """Boolean object"""
    def _convertValue(self, value):
        t = type(value)
        if t is bool:
            return value
        if t is str or t is unicode:
            v = value.strip().lower()
            return v == 'true' or v == '1'
        return bool(value) 
    
    def __nonzero__(self):
        return self.get()
    
    def __bool__(self):
        return self.get()  
    
    
class Pointer(Scalar):
    """Reference object to other one"""
    def __init__(self, value=None, **args):
        Object.__init__(self, value, objIsPointer=True, **args)   
       
    def __str__(self):
        """String representation of the scalar value"""
        return '-> %s (%s)' % (self.get().getClassName(), self.get().strId())
    
    def _convertValue(self, value):
        """Avoid storing _objValue and future objects
        obtained from .get()"""
        value.setStore(False)
        return value
    
    
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

    def getAttributesToStore(self):
        for key, attr in self.__dict__.iteritems():
            if issubclass(attr.__class__, Object) and attr._objDoStore:
                yield (key, attr)
        
        for i, item in enumerate(self):
            yield (self.getIndexStr(i+1), item)
            
    def getIndexStr(self, i):
        """Return the way the string index is generated"""
        return "__item__%06d" % i
            
    
    #TODO: check if needed
    def __len__(self):
        return list.__len__(self)
    
    def isEmpty(self):
        return len(self) > 0
    
            
class CsvList(Scalar, list):
    """This class will store a list of objects
    in a single DB row separated by comma.
    pType: the type of the list elememnts, int, bool, str"""
    def __init__(self, pType=str, **args):
        Scalar.__init__(self, **args)
        list.__init__(self)
        self._pType = pType
        
    def _convertValue(self, value):
        """Value should be a str with comman separated values"""

        del self[:] # Clear the list before setting new value
        for s in value.split(','):
            self.append(self._pType(s))
            
    def getInternalValue(self):
        return ','.join(map(str, self))
    
    def get(self):
        return self
    
    def __str__(self):
        return list.__str__(self)
    
        
class Array(Object):
    """Class for holding fixed len array"""
    def __init__(self, size=10, **args):
        Object.__init__(self, size, **args)
        
    def set(self, size):
        """Set the array size"""
        self._objValue = int(size)  
        for i in range(int(size)):
            self.__setitem__(i, None)                 
        
    def strIndex(self, i):
        return 'item_%04d' % i
    
    def __setitem__(self, index, value):
        self.__dict__[self.strIndex(index)] = value
        
    def __getitem__(self, index):
        return self.__dict__[self.strIndex(index)]
    
    def __len__(self):
        return self._objValue
    
    
def ObjectWrap(value):
    """This function will act as a simple Factory
    to create objects from Python basic types"""
    t = type(value)
    if issubclass(t, Object):
        return value
    if t is int:
        return Integer(value)
    if t is bool:
        return Boolean(value)
    if t is float:
        return Float(value)
    if t is None:
        return None
    #If it is str, unicode or unknown type, convert to string
    return String(value)
         
           
