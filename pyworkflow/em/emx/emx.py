"""
/***************************************************************************
 * Authors:     Roberto Marabini (roberto@cnb.csic.es)
 *              J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
"""
 

import os
import sys
try:
    import collections
except ImportError:
    sys.stderr.write('Could not import OrderedDict. For Python versions '
                     'earlier than 2.7 this module may be missing. '
                     )

"""
NAMING CONVENTIONS:
1) 
 xxxxx__yy will be written in XML as
<xxxx>
   <yy> value </yy>
</xxxx>

anything starting by the name of a class (i.e. micrographXXXX) 
 is a pointer to that class

GLOSARY:

 primary key: A primary key is a set of labels/attributes
 that uniquely identify an object (i.e: a micrograph)
 
 foreign key: Given two objects (i.e one micrograph an one particle) 
 a foreign key is a set of labels/attributes in the first object 
 that uniquely identify the second object 
"""
EMX_SEP = '__'
#classes 
MICROGRAPH = 'micrograph'
PARTICLE   = 'particle'
#order in which items should be read
CLASSLIST  = [MICROGRAPH, PARTICLE]
#primary keys
FILENAME   = 'fileName'
INDEX      = 'index'


class EmxLabel:
    """
    Auxiliary class to assign data type (i.e.: int, str, etc)
    and unit to each attribute/label pair.
    """
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
              ,'transformationMatrix__t14':EmxLabel(float,'px')
              ,'transformationMatrix__t21':EmxLabel(float)
              ,'transformationMatrix__t22':EmxLabel(float)
              ,'transformationMatrix__t23':EmxLabel(float)
              ,'transformationMatrix__t24':EmxLabel(float,'px')
              ,'transformationMatrix__t31':EmxLabel(float)
              ,'transformationMatrix__t32':EmxLabel(float)
              ,'transformationMatrix__t33':EmxLabel(float)
              ,'transformationMatrix__t34':EmxLabel(float,'px')
}


class EmxObject:
    """
    Base class for all EMX objects/classes
    name is the class type so far micrograph or particles
    """    
    _foreignKeys = []
    _primaryKey = []
    _attributes = []
    _name = None
    
    def __init__(self):
        self.dictPrimaryKeys = collections.OrderedDict()
        self.dictForeignKeys = collections.OrderedDict()
        self.dictAttributes  = collections.OrderedDict()

    #---------- Public object methods ------------------------------------
    def __str__(self):
        return self._pprint()
    
    def clear(self):
        """ Generic clean. PK cannot be modified or clean. """
        for key in self.dictAttributes:
            self.dictAttributes[key] = None
        for key, value in self.iterForeignKeys():
            self.dictForeignKeys[key] = None
            
    def get(self, key, default=None):
        """ Given a key (attribute name) returns the value assigned to it.
        If not present, the default will be returned.
        """
        if key in self.dictPrimaryKeys:
            return self.dictPrimaryKeys[key]
        if key in self.dictAttributes:
            return self.dictAttributes[key]
        return default
        
    def has(self, key):
        return self.get(key) is not None
    
    def set(self, key, value=None):
        """ Given a key (attribute name) assigns a value to it. """
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
                raise Exception("Key %s not allowed in: %s" % (key, self._name))

    def iterAttributes(self):
        """Returns list with valid keys (attribute names) 
           and values for this class. Primary keys are ignored"""
        return self.dictAttributes.iteritems()

    def iterPrimaryKeys(self):
        """Returns list with valid primary keys (attribute names) 
        and values for this class"""
        return self.dictPrimaryKeys.iteritems()

    def iterForeignKeys(self):
        """Returns list with valid primary keys (attribute names) 
        and values for this class"""
        return self.dictForeignKeys.iteritems()

    def __eq__(self, other):
        """ If the primary keys of two objects are identical 
        then both objects are identical"""
        return (self.dictPrimaryKeys == other.dictPrimaryKeys)

    def strongEq(self, other):
        """ true if both objects are truly identical"""
        return (self.dictPrimaryKeys == other.dictPrimaryKeys and
                self.dictForeignKeys == other.dictForeignKeys and
                self.dictAttributes == other.dictAttributes)
                   
    #---------- Internal object methods----------------------------------
    def _pprintPK(self, printNone=False):
        """ Print primary keys ."""
        out = "fileName: %(fileName)s"
        if (self.get(INDEX) != None) or printNone:
            out += ", index:%(index)s"
        return (out % self.dictPrimaryKeys) +'\n'
        
    def _pprint(self, printNone=False):
        """print ordered dictionaries, default routine is ugly.
        """
        #primary key
        out = "\nObject type: %s\n"% self._name
        out += "Primary key:\n  "
        out += self._pprintPK(printNone)
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
                        out += value._pprintPK(printNone)
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
    
    def _initAttribute(self,key, value=None):
        """Private function do not use outside this file
        """
        self.dictAttributes[key] = value
        
    def _initPrimaryKey(self,key, value=None):
        """Private function do not use outside this file
        """
        if value is None:
            pass
        else:
            self.dictPrimaryKeys[key] = emxDataTypes[key].getType()(value)
        
    def _validateForeignKey(self, className):
        """ Raise an Exception if className is not in foreign keys dict. """
        if className not in self._foreignKeys:
            raise Exception("class %s does not have FK of type: %s" % (self._name, className))
        
    def _setForeignKey(self, className, object):
        """ Set another object as foreign key of self for a given className. """
        self._validateForeignKey(className)
        self.dictForeignKeys[className] = object
        
    def _getForeignKey(self, className, validate=False):
        """ Return the object assigned as foreign key for a given className. """
        if validate:
            self._validateForeignKey(className)
        return self.dictForeignKeys.get(className, None)


class EmxImage(EmxObject):
    """ Base class for EmxMicrograph and EmxParticle.
    Both share that have FILENAME and INDEX as primary keys.
    """
    _primaryKey = [FILENAME, INDEX]
    
    def __init__(self, fileName=None, index=None):
        #init emx object
        EmxObject.__init__(self)
        if fileName is None and index is None:
            raise Exception(self._name + "cannot be created with fileName=None and index=None")
        #set primary keys. At least one of this must be different from None
        self._initPrimaryKey(FILENAME, fileName)
        self._initPrimaryKey(INDEX, index)
    
    
class EmxMicrograph(EmxImage):
    """Class for Micrographs
    """    
    _foreignKey = None
    _name = MICROGRAPH 
    _attributes = [
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
            
        
class EmxParticle(EmxImage):
    """Class for Particles
    """
    _foreignKeys = [MICROGRAPH]
    _foreignKeysMap = {MICROGRAPH:EmxMicrograph('a',1)}
    _name = PARTICLE
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
                ,'fom'
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
    
    def setMicrograph(self, micrograph):
        """ Set the micrograph associated with this particle. """
        self._setForeignKey(MICROGRAPH, micrograph)
        
    def getMicrograph(self):
        """ Return the micrograph associated with this particles. """
        return self._getForeignKey(MICROGRAPH)


class EmxData():
    """ Class to group EMX objects"""
    def __init__(self):
        self.objLists = {MICROGRAPH : [], 
                          PARTICLE  : []}
        self.mapObjectPK = {}
        self._mapper = EmxXmlMapper(self)

    def addObject(self, obj):
        self.objLists[obj._name].append(obj)
        self.mapObjectPK[str(obj.dictPrimaryKeys)] = obj
        
    def getObject(self, mapPK):
        """ Return an object given its primary key. """
        return self.mapObjectPK[str(mapPK)]

    def clear(self):
        for listName in self.objLists:
            self.objLists[listName] = []
        self.mapObjectPK = {}

    def size(self):
        return len(self.mapObjectPK)

    def __iter__(self):
        for node, nodelist in self.objLists.iteritems():
            for subnode in nodelist:
                yield subnode

    def iterClasses(self, className):
        """ Iterate through objects of a particular class. """
        return self.objLists[className]

    def __str__(self):
        partStr=""
        for k, v in self.objLists.iteritems():
            if len(v):
                partStr += "\n****\n%sS\n****\n"% k.upper()
            for obj in v:
                partStr += obj.__str__()+"\n"
        return partStr

    def read(self, emxFile):
        """ Read data from an emxFile. """
        self._mapper.readEMXFile(emxFile)
        
    def readFirstObject(self, className, emxFile):
        """ Read only the first object of a given className from file."""
        return self._mapper.firstObject(className, emxFile)
        
    def write(self, emxFile):
        """ Write data to an emxFile. """
        self._mapper.writeEMXFile(emxFile)

#------------------- XmlMapper implementation -----------------------------------
"""
Following from here there is the implementation of the XmlMapper to store
EmxData objects in XML files as described in the EMX format.
This mapper is the default one used by EmxData (and no other make sense at this moment)
"""
try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET
    
ERR_VALIDATION_WRONG=1
VERSION = 1.0
ROOTNAME = 'EMX'
HEADER = """<%(ROOTNAME)s version="%(VERSION)s">
  <!--
  ##########################################################################
  #               EMX Exchange file 
  #               Produced using the emx library 
  #               (http://i2pc.cnb.csic.es/emx/LoadTools.htm?type=Library)
  # 
  #  Information on this file format is available at 
  #  http://i2pc.cnb.csic.es/emx
  ##########################################################################
  #  One of the best ways you can help us to improve this software
  #  is to let us know about any problems you find with it.
  #  Please report bugs to: emx@cnb.csic.es
  ##########################################################################
  -->
""" % globals()
  
EMXSCHEMA10 = 'http://sourceforge.net/p/emexchange/code/ci/master/tree/trunk/resourcesEmx/schemas/emx.xsd?format=raw'
EMXSCHEMA11 = 'http://sourceforge.net/p/emexchange/code/ci/master/tree/trunk/resourcesEmx/schemas/emx_11.xsd?format=raw'


class ValidateError(Exception):
    def __init__(self, code, message):
        self.errorMessage = message
        self.errorCode = code
    def __str__(self):
        return "Error Code: %d. Message: %s" % (self.errorCode, self.errorMessage)
    def getCode(self):
        return self.errorCode
    def getMessage(self):
        return self.errorMessage
    
    
class EmxXmlMapper():
    """Mapper for XML"""
    def __init__(self, emxData):
        self.emxData = emxData
        self.classObject = {MICROGRAPH: EmxMicrograph, PARTICLE: EmxParticle}

    def __del__(self):
        pass
    
    def objectToXML(self, object):
        """ Given an object persist it in XML dataBase.
        Each object goes to a different element. Much much faster...
        """
        
        # write primary key
        xmlString = r"  <%s" % object._name
        for key, value in object.dictPrimaryKeys.iteritems():
            if value is None:
                continue
            xmlString += ' %(key)s="%(value)s"' % ({'key':key, 'value':str(value)})
        xmlString += ">\n"
        
        # write attributes
        oldParent = ""
        for key, value in object.iterAttributes():
            if value is None:
                continue
            unit = emxDataTypes[key].getUnit()
            # is this an special case, that is, 
            # does the label contains '__'?
            # I asumme there is no grandchild
            if EMX_SEP in key:
                (parent, child) = key.split(EMX_SEP)
                # take care of cases like:
                # <pixelSpacing>
                # <X>5.6</X>
                # <Y>5.7</Y>
                # </pixelSpacing>
                # second entry
                if oldParent == parent:
                    xmlString = xmlString.replace("  </%s>\n" % parent, "")
                # first entry
                else:
                    xmlString += "    <%s>\n  " % parent
                if unit is None:
                    xmlString += "    <%(child)s>%(value)s"\
                              "</%(child)s>\n  </%(parent)s>\n" % ({'parent':parent,
                                                           'child':child,
                                                           'value':str(value)})
                else:
                    xmlString += '    <%(child)s unit="%(unit)s">%(value)s'\
                              "</%(child)s>\n    </%(parent)s>\n" % ({'parent':parent,
                                                           'child':child,
                                                           'value':str(value),
                                                           'unit':unit})
                oldParent = parent
            # simple attributes with no child
            else:
                if unit is None:
                    xmlString += "    <%(key)s>%(value)s</%(key)s>\n" % ({'key':key, 'value':str(value)})
                else:
                    xmlString += '    <%(key)s unit="%(unit)s">%(value)s</%(key)s>\n' % ({'key':key, 'value':str(value), 'unit':unit})
                    
        # write foreign key
        if len(object.dictForeignKeys) and\
               object.dictForeignKeys[object.dictForeignKeys.keys()[0]]:
            pointedObject = object.dictForeignKeys[object.dictForeignKeys.keys()[0]]
            xmlString += "    <%s" % pointedObject._name
            for key, value in pointedObject.dictPrimaryKeys.iteritems():
                if value is None:
                    continue
                xmlString += ' %(key)s="%(value)s"' % ({'key':key, 'value':str(value)})
            xmlString += "/>\n"
        xmlString += "  </%s>\n" % object._name
        # print xmlString
        return xmlString

    _attributes = [
        'acceleratingVoltage'
        , 'activeFlag'
        , 'amplitudeContrast'
        , 'cs'
        , 'defocusU'
        , 'defocusV'
        , 'defocusUAngle'
        , 'fom'
        , 'pixelSpacing__X'
        , 'pixelSpacing__Y'
        , 'pixelSpacing__Z'
        ] 
    
    def readEMXFile(self, fileName, classElement=None):
        """ create tree from xml file 
        If classElement is not None, the first element of this class
        will be returned
        """
        # get context
        context = ET.iterparse(fileName, events=('start', 'end'))
        # turn it into an iterator
        context = iter(context)
        # get the root element
        event, root = context.next()
        
        # self.classObject = globals()['Emx'+element]
        doItPK = True
        skipLabelPK = 'kk'
        
        mergeParent = False
        parentLabel = 'kk'
        
        lastStartTagA = 'kk'
        lastEventStartA = False

        listObjectWithForeignKey = []

        for event, elem in context:
            tag = elem.tag
            if tag == 'EMX':
                continue
            if event == 'start':
                # primary key and FK
                if tag in CLASSLIST:
                    # only primary key
                    if(doItPK):
                        self.createObject(elem)
                        doItPK = False
                        skipLabelPK = tag
                    # foreign key
                    else:
                        # get PF and save the map for the first pass 
                        # since the actual pointed object may not exists
                        FK = self.readObjectPK(elem)
                        self._object._setForeignKey(tag, FK)
                        listObjectWithForeignKey.append(self._object)
                else:
                    if lastEventStartA == True:
                        mergeParent = True
                        parentLabel = lastStartTagA
                    lastStartTagA = tag
                    lastEventStartA = True

            elif event == 'end':
                # PK or FG
                if tag in CLASSLIST and skipLabelPK == tag:
                    doItPK = True
                    if tag == classElement:
                        return
                # other attributes
                else:
                    # simple element
                    if lastStartTagA == tag:
                        if elem.text is None:
                            raise Exception ("Element: " + tag + " is empty")
                        else:
                            text = elem.text.strip(' \n\t')
                        if(len(text) < 1):
                            raise Exception ("ZERO for tag=%s, value=%s" % (tag, text))
                        if  mergeParent: 
                            self._object.set(parentLabel + EMX_SEP + tag, text)
                        else:
                            self._object.set(tag, text)
                    elif parentLabel == tag:
                        mergeParent = False
                        parentLabel = 'kk'
                    lastEventStartA = False
            else:
                raise Exception ("Unknown event type %s" % event)
            root.clear()
        # Now loop Trough all objects and fix the FK
        for object in listObjectWithForeignKey:
            for key in object._foreignKeys:
                fk = object._getForeignKey(key)
                object._setForeignKey(key, self.emxData.getObject(fk))

    def createObject(self, elem):
        self.myClass = self.classObject[elem.tag]
        # primary key 
        # get PK
        self.dict = self.readObjectPK(elem)
        # create object
        self._object = self.classObject[elem.tag](**(self.dict))
        # add it to emxData
        self.emxData.addObject(self._object)
        
    def readObjectPK(self, elem):
        """ read primary key. So far all entities has the same PK. 
        We may need to specialize or use dictPrimaryKeys in the future
        """
        mapPK = collections.OrderedDict()
        for attribute in elem.attrib:
            mapPK[attribute] = emxDataTypes[attribute].getType()(elem.get(attribute))
        if mapPK:
            return collections.OrderedDict(sorted(mapPK.items(), key=lambda t: t[0]))
        else:
            raise Exception("readObjectPK: No fileName or index provided")

    def firstObject(self, classname, fileName):
        """ Iterate over the tags elements and find the 
        first one of type 'classname', build the object
        and return it. The foreing keys will be not updated.
        """
        self._object = None
#        context = ET.iterparse(fileName, events=('start', 'end'))
#        for event, elem in iter(context):
#            tag = elem.tag
#            if event == 'start':
#                # print "tag: '%s'" % tag, "class: '%s'" % classname
#                if tag == classname:
#                    # print "tag==class"
#                    self.createObject(elem)
#                    # print "self._object: ", self._object
#                    return self._object
        self.readEMXFile(fileName, classElement=classname)
        
        return self._object
        
    def writeEMXFile(self, fileName):
        """read xml file and store it in a document
        """
        xmlFile = open(fileName, "w")
        xmlFile.write("<?xml version='1.0' encoding='utf-8'?>\n")
        xmlFile.write(HEADER)

        for object in self.emxData:
            text = self.objectToXML(object)  # 
#            #implement this with a regular expression
#            #format matrices properly
            for i, j in {
                          '</t11>\n    ':'</t11> '
                         , '</t12>\n    ':'</t12> '
                         , '</t13>\n    ':'</t13> '
                         , '</t21>\n    ':'</t21> '
                         , '</t22>\n    ':'</t22> '
                         , '</t23>\n    ':'</t23> '
                         , '</t31>\n    ':'</t31> '
                         , '</t32>\n    ':'</t32> '
                         , '</t33>\n    ':'</t33> '
                         }.iteritems():
                text = text.replace(i, j)
            xmlFile.write(text)
        xmlFile.write("</EMX>")
        xmlFile.close()

def validateSchema(filename, schema_file=None):
    """
    Code from astropy project released under BSD licence
    Validates an XML file against a schema or DTD.

    Functions to do XML schema and DTD validation.  At the moment, this
    makes a subprocess call first to xerces then to  xmllint. 
    This could use a Python-based
    library at some point in the future, if something appropriate could be
    found. lxml is a possibility but has too many dependences if anyone
    knows about a pure python validator let my know

    Parameters
    ----------
    filename : str
        The path to the XML file to validate

    schema : str
        The path to the XML schema or DTD

    Returns
    -------
    returncode, stdout, stderr : int, str, str
        Returns the returncode from validator and the stdout and stderr
        as strings
    """
    import subprocess, os
    ###########
    # try xerces
    ###########
    answerSize = 1024  # avoid overflow in web
    endding = ''
    if schema_file is None:
        _schema = EMXSCHEMA11
    else:
        _schema = schema_file
    # print "java jaxp.SourceValidator -a %s -i %s -xsd11"% (_schema, filename)
    p = subprocess.Popen("java jaxp.SourceValidator -a %s -i %s -xsd11" 
                         % (_schema, filename),
            shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    # xerces exists but is error
    if p.returncode == 0 and (stderr != ""):
        if len(stderr) > answerSize:
           endding = '... (too many errors, displayed first %d characters)' % (answerSize)
        raise ValidateError(ERR_VALIDATION_WRONG, """Error: when validating file %s with schema %s.
        \nError:%s""" % (filename, _schema, stderr[:answerSize] + endding))
    #######
    # no xerces available, let us try xmlint
    ######
    if p.returncode != 0:
        print "validating with xmllint"
        if schema_file is None:
            _schema = EMXSCHEMA10
        else:
            _schema = schema_file
        schema_part = '--schema ' + _schema
        p = subprocess.Popen(
            "xmllint --noout %s %s" % (schema_part, filename),
            shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = p.communicate()
        if p.returncode == 127:
            raise ValidateError(127,
                """Error: neither xerces-f nor xmllint could be found,  I cannot validate schema. 
    Schema validation is based either on the xmllint program that belongs to the libxml2-tools package.
    or on the xerces-f project""")
            
        if p.returncode != 0:
            if len(stderr) > answerSize:
                 endding = '... (too many errors, displayed first %d characters)' % (answerSize)
            message = """Error: when validating file %s with schema %s.
            \nError:%s""" % (filename, _schema, stderr[:answerSize] + endding)
            # print "message", message
            raise ValidateError(ERR_VALIDATION_WRONG, message)
    return p.returncode, stdout[:answerSize], stderr[:answerSize]

