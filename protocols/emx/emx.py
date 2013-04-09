'''
/***************************************************************************
 * Authors:     Roberto Marabini (roberto@cnb.csic.es)
 *              Jose Miguel de la Rosa
 *
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/
MODIFICATION ADVICE:

Please,  do not  generate or  distribute 
a modified version of this file under its original name. 
 '''
 
VERSION='EMX_1.0'
import sys
try:
    import collections
except ImportError:
    sys.stderr.write('Could not import OrderedDict. For Python versions '
                     'earlier than 2.7 this module may be missing. '
                     )

'''naming conventions:
1) 
 xxxxx__yy will be written in XML as
<xxxx>
   <yy> value </yy>
</xxxx>

anything starting by the name of a class (i.e. micrographXXXX) 
 is a pointer to that class
'''
'''
GLOSARY:

 primary key: A primary key is a set of labels/attributes
 that uniquely identify an object (i.e: a micrograph)
 
 foreign key: Given two objects (i.e one micrograph an one particle) 
 a foreign key is a set of labels/attributes in the first object 
 that uniquely identify the second object 
'''
EMX_SEP = '__'
#classes 
MICROGRAPH = 'micrograph'
PARTICLE   = 'particle'
#order in which items should be read
CLASSLIST  = [MICROGRAPH,PARTICLE]
#primary keys
FILENAME   = 'fileName'
INDEX      = 'index'

class EmxLabel:
    '''auxiliary class to assign data type (i.e.: int, str, etc)
    and unit to each attribute/label
    '''
    def __init__(self, type, unit=None):
        self.type = type
        self.unit = unit
        
    def hasUnit(self):
        return not self.unit is None

    def getUnit(self):
        return self.unit

    def getType(self):
        return self.type

#Dictionary with attribute names, data types and units
#By default an attribute does not has unit assign to it
emxDataTypes={
              FILENAME:EmxLabel(str)
              ,INDEX:EmxLabel(int)
              ,'acceleratingVoltage':EmxLabel(float,'kV')
              ,'activeFlag':EmxLabel(int)
              ,'amplitudeContrast':EmxLabel(float)
              ,'boxSize__X':EmxLabel(int,'px')
              ,'boxSize__Y':EmxLabel(int,'px')
              ,'boxSize__Z':EmxLabel(int,'px')
              ,'centerCoord__X':EmxLabel(float,'px')
              ,'centerCoord__Y':EmxLabel(float,'px')
              ,'centerCoord__Z':EmxLabel(float,'px')
              ,'cs':EmxLabel(float,'mm')
              ,'defocusU':EmxLabel(float,'nm')
              ,'defocusV':EmxLabel(float,'nm')
              ,'defocusUAngle':EmxLabel(float,'deg')
              ,'fom':EmxLabel(float)
              ,'pixelSpacing__X':EmxLabel(float,'A/px')
              ,'pixelSpacing__Y':EmxLabel(float,'A/px')
              ,'pixelSpacing__Z':EmxLabel(float,'A/px')
              ,'transformationMatrix__t11':EmxLabel(float)
              ,'transformationMatrix__t12':EmxLabel(float)
              ,'transformationMatrix__t13':EmxLabel(float)
              ,'transformationMatrix__t14':EmxLabel(float,'A')
              ,'transformationMatrix__t21':EmxLabel(float)
              ,'transformationMatrix__t22':EmxLabel(float)
              ,'transformationMatrix__t23':EmxLabel(float)
              ,'transformationMatrix__t24':EmxLabel(float,'A')
              ,'transformationMatrix__t31':EmxLabel(float)
              ,'transformationMatrix__t32':EmxLabel(float)
              ,'transformationMatrix__t33':EmxLabel(float)
              ,'transformationMatrix__t34':EmxLabel(float,'A')
}

class EmxObject:
    '''Base class for all emx objects/classes
       name is the class type so far micrograph or particles
    '''    
    _foreignKeys=[]
    _primaryKey=[]
    _attributes=[]
    def __init__(self,name):
        #try:
            #dictionaries with labels used as 1) primary keys
            #2)foreing keys and 3) plain attributes
        self.dictPrimaryKeys = collections.OrderedDict()
        self.dictForeignKeys = collections.OrderedDict()
        self.dictAttributes  = collections.OrderedDict()
#        except ImportError:
#            #ordereddict was introduced in python 2.7
#            #use plain dictionaries if no available
#            self.dictPrimaryKeys = {}
#            self.dictForeignKeys = {}
#            self.dictAttributes  = {}
        #chlid class name
        self.name=name


    def clear(self):
        """ Generic clean. PK cannot be modified or clean
        """
        for key in self.dictAttributes:
            self.dictAttributes[key] = None
        for key, value in self.iterForeignKeys():
            self.dictForeignKeys[key] = None

    def pprint_pk(self, printNone=False):
        '''Print primary keys
        '''
        out = "fileName: %(fileName)s"
        if (self.get(INDEX) != None) or printNone:
            out += ", index:%(index)s"
        return (out % self.dictPrimaryKeys) +'\n'
        
    def pprint_od(self, printNone=False):
        '''print ordered dictionaries, default routine is ugly.
        '''
        #primary key
        out = "\nObject type: %s\n"% self.name
        out += "Primary key:\n  "
        out += self.pprint_pk(printNone)
        #foreign key
        if len(self.dictForeignKeys) and\
               self.dictForeignKeys[self.dictForeignKeys.keys()[0]]:
            out += "Foreign keys:\n  "
            for key, value in self.dictForeignKeys.iteritems():
                if (value != None) or printNone:
                    out += "%s -> " % key
                    if value is None:
                        out += "None"
                    else:
                        out += value.pprint_pk(printNone)
                        #out += str(value) +"\n"
            
        #other attributes
        out += "Other Attributes:\n"
        for key, value in self.iterAttributes():
            _unit = emxDataTypes[key].getUnit()
            if _unit is None:
                _unit=""
            else:
                _unit= "(%s)"%_unit
            if (value != None) or printNone:
                out += "  %s %s:%s,\n" % (key, _unit, str(value) )
        return out
    
    def _initAttribute_ (self,key, value=None):
        '''Private function do not use outside this file
        '''
        self.dictAttributes[key] = value
        
    def _initPrimaryKey_ (self,key, value=None):
        '''Private function do not use outside this file
        '''
        if value is None:
            pass
        else:
            self.dictPrimaryKeys[key] = emxDataTypes[key].getType()(value)
        
    def get(self, key):
        '''given a key (attribute name) returns 
           the value assigned to it'''
        if key in self.dictPrimaryKeys:
            return self.dictPrimaryKeys[key]
        if key in self.dictAttributes:
            return self.dictAttributes[key]
        return None
        #raise Exception("Key %s not found in: %s" % (key, self.name))
    
    def set(self, key, value=None):
        '''given a key (attribute name) assigns a value to it'''
        if value is None:
            if key in self._primaryKey:
                self.dictPrimaryKeys[key] = None
            elif key in self._attributes:
                self.dictAttributes[key] = None
        else:
            if key in self._primaryKey:
                self.dictPrimaryKeys[key] = emxDataTypes[key].getType()(value)
            elif key in self._attributes:
                self.dictAttributes[key] = emxDataTypes[key].getType()(value)
            else:
                raise Exception("Key %s not allowed in: %s" % (key, self.name))

    def iterAttributes(self):
        '''Returns list with valid keys (attribute names) 
           and values for this class. Primary keys are ignored'''
        return self.dictAttributes.iteritems()

    def iterPrimaryKeys(self):
        '''Returns list with valid primary keys (attribute names) 
        and values for this class'''
        return self.dictPrimaryKeys.iteritems()

    def iterForeignKeys(self):
        '''Returns list with valid primary keys (attribute names) 
        and values for this class'''
        return self.dictForeignKeys.iteritems()

    def __str__(self):
        
        #strPrint  = self.pprint_pk()
        return self.pprint_od()

    def __eq__(self, other):
        """ If the primary keys of two objects are identical 
        then both objects are identical"""
        return (self.dictPrimaryKeys == other.dictPrimaryKeys)

    def strongEq(self, other):
        """ true if both objects are truly identical"""
        return (self.dictPrimaryKeys == other.dictPrimaryKeys) and\
               (self.dictForeignKeys == other.dictForeignKeys) and\
               (self.dictAttributes == other.dictAttributes)
        
    def getForeignKey(self,object):
        """ return empty object with foreign key
        """
        if object.name not in self._foreignKeys:
            raise Exception("class %s does not have FK of type: %s"%(self.name, object.name))
        object.clear(self.dictForeignKeys[object.name])


    def setForeignKey(self, object):
        """Set object as foreign key
        """
        self._setForeignKey(object.name, object)

    def _setForeignKey(self, className, object):
        """Set object as foreign key
        """
        if className not in self._foreignKeys:
            raise Exception("class %s does not have FK of type: %s"%(self.name, object.name))
        self.dictForeignKeys[className] = object

    def addForeignKey(self, className, object):
        """Set object as foreign key
        """
        if className not in self._foreignKeys:
            raise Exception("class %s does not have FK of type: %s"%(self.name, object.name))
        self.dictForeignKeys[className] = object

class Emxmicrograph(EmxObject):
    '''Class for Micrographs
    '''    
    _primaryKey=[FILENAME,INDEX]
    _foreignKey=None
    _attributes=[
        'acceleratingVoltage'
        ,'activeFlag'
        ,'amplitudeContrast'
        ,'cs'
        ,'defocusU'
        ,'defocusV'
        ,'defocusUAngle'
        ,'fom'
        ,'pixelSpacing__X'
        ,'pixelSpacing__Y'
        ,'pixelSpacing__Z'
        ]
    def __init__(self,fileName=None,index=None):
        #init emx object
        EmxObject.__init__(self,MICROGRAPH)
        if fileName is None and index is None:
            raise Exception(MICROGRAPH + "cannot be created with fileName=None and index=None")
        #set primary keys. At least one of this must be different from None
        self._initPrimaryKey_(FILENAME,fileName)
        self._initPrimaryKey_(INDEX, index)
            
        
class Emxparticle(EmxObject):
    '''Class for Particles
    '''
    _primaryKey=[FILENAME,INDEX]
    _foreignKeys=[MICROGRAPH]
    _foreignKeysMap={MICROGRAPH:Emxmicrograph('a',1)}
    _attributes=['activeFlag'
                ,'boxSize__X'
                ,'boxSize__Y'
                ,'boxSize__Z'
                ,'centerCoord__X'
                ,'centerCoord__Y'
                ,'centerCoord__Z'
                ,'defocusU'
                ,'defocusV'
                ,'defocusUAngle'
                ,'pixelSpacing__X'
                ,'pixelSpacing__Y'
                ,'pixelSpacing__Z'
                ,'transformationMatrix__t11'
                ,'transformationMatrix__t12'
                ,'transformationMatrix__t13'
                ,'transformationMatrix__t14'
                ,'transformationMatrix__t21'
                ,'transformationMatrix__t22'
                ,'transformationMatrix__t23'
                ,'transformationMatrix__t24'
                ,'transformationMatrix__t31'
                ,'transformationMatrix__t32'
                ,'transformationMatrix__t33'
                ,'transformationMatrix__t34'
                ]
    def __init__(self,fileName,index=None):
        #init emx object
        EmxObject.__init__(self,PARTICLE)
        if fileName is None and index is None:
            raise Exception(PARTICLE + "cannot be created with fileName=None and index=None")
        #set primary keys
        self._initPrimaryKey_(FILENAME, fileName)
        self._initPrimaryKey_(INDEX, index) # Index of an image in a multi-image 
        #set foreign keys
#        if micrograph is not None:
#            self._setForeignKey(MICROGRAPH,micrograph)
    def _setXMLForeignKey(self, className, mapFK):
        """Never ever use this function unless you know what you are doing
        """
        if className not in self._foreignKeys:
            raise Exception("class %s does not have FK of type: %s"%(self.name, object.name))
        self.dictForeignKeys[className] = mapFK

    def _getXMLForeignKey(self, className):
        """Never ever use this function unless you know what you are doing
        """
        if className not in self._foreignKeys:
            raise Exception("class %s does not have FK of type: %s"%(self.name, object.name))
        return self.dictForeignKeys[className]

class EmxData():
    ''' Class to group EMX objects'''
    def __init__(self):
        self.objLists = {MICROGRAPH : [], 
                          PARTICLE  : []}
        self.mapObjectPK={}

    def addObject(self,obj):
        self.objLists[obj.name].append(obj)
        self.mapObjectPK[str(obj.dictPrimaryKeys)]=obj
        
    def getObjectwithPK(self, mapPK):
        return self.mapObjectPK[str(mapPK)]

    def clear(self):
        for listName in self.objLists:
            self.objLists[listName]=[]
        self.mapObjectPK={}

    def size(self):
        return len(self.mapObjectPK)

    def __iter__(self):
        for node, nodelist in self.objLists.iteritems():
            for subnode in nodelist:
                yield subnode

    def iterClasses(self,className):
        """Iterate through a particular class
        """
        return self.objLists[className]

    def __str__(self):
        partStr=""
        for k,v in self.objLists.iteritems():
            if len(v):
                partStr += "\n****\n%sS\n****\n"% k.upper()
            for obj in v:
                partStr += obj.__str__()+"\n"
        return partStr