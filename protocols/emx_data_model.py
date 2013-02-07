'''
/***************************************************************************
 * Authors:     Roberto Marabini (roberto@cnb.csic.es)
 *              Jose Miguel 
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
                     'Program will work but attributes '
                     'are not guarantee to be sorted by alphabetical order' 
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
              ,INDEX:EmxLabel(int)
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
    def __init__(self,name):
        try:
            #dictionaries with labels used as 1) primary keys
            #2)foreing keys and 3) plain attributes
            self.dictPrimaryKeys = collections.OrderedDict()
            self.dictForeignKeys = collections.OrderedDict()
            self.dictAttributes  = collections.OrderedDict()
        except ImportError:
            #ordereddict was introduced in python 2.7
            #use plain dictionaries if no available
            self.dictPrimaryKeys = {}
            self.dictForeignKeys = {}
            self.dictAttributes  = {}
        #chlid class name
        self.name=name
    
    def pprint_pk(self, printNone=False):
        '''Print primary keys
        '''
        out = "fileName: %(fileName)s"
        if (self.get(INDEX) != None) or printNone:
            out += " (index=%(index)s)"
        return (out % self.dictPrimaryKeys) +'\n'
        
    def pprint_od(self, printNone=False):
        '''print ordered dictionaries, default routine is ugly.
        INPUT: ordered dictionary
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
            
        #other attributes
        out += "Other Attributes:\n"
        for key, value in self.iterAttributes():
            if (value != None) or printNone:
                out += "  %s:%s,\n" % (key, str(value) )
        return out
    
    def _initAttribute_ (self,key, value=None):
        '''Private function do not use outside this file
        '''
        self.dictAttributes[key] = value
        
    def _initPrimaryKey_ (self,key, value=None):
        '''Private function do not use outside this file
        '''
        self.dictPrimaryKeys[key] = value
        
    def get(self, key):
        '''given a key (attribute name) returns 
           the value assigned to it'''
        if key in self.dictPrimaryKeys:
            return self.dictPrimaryKeys[key]
        if key in self.dictAttributes:
            return self.dictAttributes[key]
        raise Exception("Key %s not found in: %s" % (key, self.name))
    
    def set(self, key, value=None):
        '''given a key (attribute name) assigns a value to it'''
        if key in self.dictPrimaryKeys:
            self.dictPrimaryKeys[key] = value
        elif key in self.dictAttributes:
            self.dictAttributes[key] = value
        else:
            raise Exception("Key %s not allowed in: %s" % (key, self.name))

    def iterAttributes(self):
        '''Returns list with valid keys (attribute names) 
           for this class. Primary keys are ignored'''
        return self.dictAttributes.iteritems()

    def iterPrimaryKeys(self):
        '''Returns list with valid primary keys (attribute names) 
        for this class'''
        return self.dictPrimaryKeys.iteritems()

    def __eq__(self, other):
        '''equality operator'''
        return self.dictAttributes == other.dictAttributes\
               and self.dictPrimaryKeys == other.dictPrimaryKeys\
               and self.dictForeignKeys == other.dictForeignKeys
    
    def __str__(self):
        '''print operator'''
        return self.pprint_od()
    
    def comparePK(self, **args):
        return self.dictPrimaryKeys == args
        

class EmxMicrograph(EmxObject):
    '''Class for Micrographs
    '''    
    def __init__(self,fileName,index=None,activeFlag=1):
        #init emx object
        EmxObject.__init__(self,MICROGRAPH)
        #define primary keys
        self._initPrimaryKey_(FILENAME,fileName)
        self._initPrimaryKey_(INDEX, index)
        #define rest of attributes
        self._initAttribute_('acceleratingVoltage')
        self._initAttribute_('activeFlag',activeFlag)
        self._initAttribute_('amplitudeContrast')
        self._initAttribute_('cs')
        self._initAttribute_('defocusU')
        self._initAttribute_('defocusV')
        self._initAttribute_('defocusUAngle')
        self._initAttribute_('fom')
        for i in ['X','Y','Z']:
            self._initAttribute_('pixelSpacing__'+i)

class EmxParticle(EmxObject):
    '''Class for Particles
    '''    
    def __init__(self,fileName,index=1,micrograph=None,activeFlag=1):
        #init emx object
        EmxObject.__init__(self,PARTICLE)
        #define primary keys
        self._initPrimaryKey_(FILENAME, fileName)
        self._initPrimaryKey_(INDEX, index) # Index of an image in a multi-image 
                                                     # file. Defaults to 1 if omitted.'''
        #define foreign keys
        self.setMicrograph(micrograph)
        #define rest of attributes
        self._initAttribute_('activeFlag', activeFlag)
        for i in ['X','Y','Z']:
            self._initAttribute_('boxSize__'+i)
            self._initAttribute_('centerCoord__'+i)
        self._initAttribute_('defocusU')
        self._initAttribute_('defocusV')
        self._initAttribute_('defocusUAngle')
        for i in ['X','Y','Z']:
            self._initAttribute_('pixelSpacing__'+i)
        for i in ['11','12','13','14','21','22','23','24','31','32','33','34']:
            self._initAttribute_('transformationMatrix__t'+i)
            
    def getMicrograph(self):
        return self.dictForeignKeys[MICROGRAPH]
    
    def setMicrograph(self, micrograph):
        self.dictForeignKeys[MICROGRAPH] = micrograph
        
class EmxData():
    '''EMX data model. version 1.0
       No file format information here
    '''    
    def __init__(self):
        self.version = 1.0
        self.listParticles = []
        self.listMicrographs = []
        self.dictLists = {MICROGRAPH: self.listMicrographs, 
                          PARTICLE: self.listParticles}
    
    def __eq__(self,other):
        ''' equality operator'''
        micBool = all (
                       self.listMicrographs[m] == 
                       other.listMicrographs[m] 
                       for m in range(len(self.listMicrographs))
                       )
        parBool= all(
                    self.listParticles[p] == 
                    other.listParticles[p] 
                    for p in range(len(self.listParticles)) 
                     )
        return parBool and micBool

    def appendMicrograph(self,object):
        self.listMicrographs.append(object)

    def appendParticle(self,object):
        self.listParticles.append(object)
    
    def findObject(self, objList, **objPK):
        for obj in objList:
            if obj.comparePK(**objPK):
                return obj
        return None

    def __str__(self):
        ''' print operator'''
        str = 'MICROGRAPH\n '
        for micrograph in self.listMicrographs:
            str += micrograph.__str__()
        str += '\n\nPARTICLES\n '
        for particle in self.listParticles:
            str += particle.__str__()
        return str

