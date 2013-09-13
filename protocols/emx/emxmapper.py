# **************************************************************************
# *
# * Authors:     Roberto Marabini       (roberto@cnb.csic.es)
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
# *  e-mail address 'roberto@cnb.csic.es'
# *
# **************************************************************************
try:
    import collections
except ImportError:
    sys.stderr.write('Could not import OrderedDict. For Python versions '
                     'earlier than 2.7 this module may be missing. '
                     )

VERSION=1.0
ROOTNAME = 'EMX'
HEADER = '''
  ##########################################################################
  #               EMX Exchange file 
  #               Produced using the emxLibrary 
  #               (http://i2pc.cnb.csic.es/emx/LoadTools.htm?type=Library)
  # 
  #  This is a EMX file.
  #
  #  Information on this file format is available at 
  #  http://i2pc.cnb.csic.es/emx
  ##########################################################################
  #  One of the best ways you can help us to improve this software
  #  is to let us know about any problems you find with it.
  #  Please report bugs to: emx@cnb.csic.es
  ##########################################################################
  '''
EMXSCHEMA10 ='http://sourceforge.net/p/emexchange/code/ci/master/tree/trunk/resourcesEmx/schemas/emx.xsd?format=raw'
EMXSCHEMA11 ='http://sourceforge.net/p/emexchange/code/ci/master/tree/trunk/resourcesEmx/schemas/emx_11.xsd?format=raw'
from emx import MICROGRAPH, EMX_SEP, CLASSLIST, emxDataTypes
from emx import Emxmicrograph, Emxparticle

from os.path import exists
from os      import remove
try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET

#def convert2Dictionary(mydict):
#    s = ''
#    for value in mydict.values(): 
#        s += value
#    return s
class ValidateError(Exception):
    def __init__(self, code, message):
        self.errorMessage = message
        self.errorCode = code
    def __str__(self):
        return "Error Code: %d. Message: %s" % (self.errorCode, self.errorMessage)
    def getCode(self):
        return self.code
    def getMessage(self):
        return self.message
    
class XmlMapper():
    '''Mapper for XML'''
    def __init__(self, emxData):
        self.emxData = emxData
        self.classObject={}
        for element in CLASSLIST:
            self.classObject[element] = globals()['Emx'+element]

    def __del__(self):
        pass
    
    def objectToXML(self, object):
        """ given an object persist it in XML dataBase.
        Each object goes to a different element. Much much faster...
        """
        
        #write primary key
        xmlString = r"  <%s" % object.name
        for key, value in object.dictPrimaryKeys.iteritems():
            if value is None:
                continue
            xmlString += ' %(key)s="%(value)s"'%({'key':key,'value':str(value)})
        xmlString += ">\n"
        
        #write attributes
        oldParent=""
        for key, value in object.iterAttributes():
            if value is None:
                continue
            unit = emxDataTypes[key].getUnit()
            #is this an special case, that is, 
            # does the label contains '__'?
            # I asumme there is no grandchild
            if EMX_SEP in key:
                (parent,child) = key.split(EMX_SEP)
                #take care of cases like:
                #<pixelSpacing>
                # <X>5.6</X>
                # <Y>5.7</Y>
                #</pixelSpacing>
                #second entry
                if oldParent == parent:
                    xmlString = xmlString.replace("  </%s>\n"%parent,"")
                #first entry
                else:
                    xmlString += "    <%s>\n  "%parent
                if unit is None:
                    xmlString += "    <%(child)s>%(value)s"\
                              "</%(child)s>\n  </%(parent)s>\n"%({'parent':parent,
                                                           'child':child,
                                                           'value':str(value)})
                else:
                    xmlString += '    <%(child)s unit="%(unit)s">%(value)s'\
                              "</%(child)s>\n    </%(parent)s>\n"%({'parent':parent,
                                                           'child':child,
                                                           'value':str(value),
                                                           'unit':unit})
                oldParent = parent
            #simple attributes with no child
            else:
                if unit is None:
                    xmlString += "    <%(key)s>%(value)s</%(key)s>\n"%({'key':key,'value':str(value)})
                else:
                    xmlString += '    <%(key)s unit="%(unit)s">%(value)s</%(key)s>\n'%({'key':key,'value':str(value),'unit':unit})
                    
        #write foreign key
        if len(object.dictForeignKeys) and\
               object.dictForeignKeys[object.dictForeignKeys.keys()[0]]:
            pointedObject = object.dictForeignKeys[object.dictForeignKeys.keys()[0]]
            xmlString += "    <%s"%pointedObject.name
            for key, value in pointedObject.dictPrimaryKeys.iteritems():
                if value is None:
                    continue
                xmlString += ' %(key)s="%(value)s"'%({'key':key,'value':str(value)})
            xmlString += "/>\n"
        xmlString += "  </%s>\n"%object.name
        #print xmlString
        return xmlString

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
    
    def readEMXFile(self, fileName):
        """ create tree from xml file 
        """
        #get context
        context     = ET.iterparse(fileName, events=('start','end'))
        # turn it into an iterator
        context     = iter(context)
        # get the root element
        event, root = context.next()
        
        #self.classObject = globals()['Emx'+element]
        doItPK = True
        skipLabelPK='kk'
        
        mergeParent = False
        parentLabel='kk'
        
        lastStartTagA='kk'
        lastEventStartA = False

        listObjectWithForeignKey = []

        for event, elem in context:
            tag = elem.tag
            if tag == 'EMX':
                continue
            if event == 'start':
                #primary key and FK
                if tag in CLASSLIST:
                    #only primary key
                    if(doItPK):
                        self.createObject(elem)
                        doItPK = False
                        skipLabelPK = tag
                    #foreign key
                    else:
                        #get PF and save the map for the first pass 
                        #since the actual pointed object may not exists
                        FK=self.readObjectPK(elem)
                        self._object._setXMLForeignKey(tag,FK )
                        listObjectWithForeignKey.append(self._object)
                else:
                    if lastEventStartA == True:
                        mergeParent   = True
                        parentLabel   = lastStartTagA
                    lastStartTagA=tag
                    lastEventStartA = True

            elif event == 'end':
                #PK or FG
                if tag in CLASSLIST and skipLabelPK == tag:
                    doItPK = True
                #other attributes
                else:
                    #simple element
                    if lastStartTagA == tag:
                        if elem.text is None:
                            raise Exception ("Element: "+tag+" is empty")
                        else:
                            text = elem.text.strip(' \n\t')
                        if(len(text)<1):
                            raise Exception ("ZERO for tag=%s, value=%s"%(tag,text))
                        if  mergeParent: 
                            self._object.set(parentLabel+EMX_SEP+tag, text)
                        else:
                            self._object.set(tag, text)
                    elif parentLabel == tag:
                        mergeParent   = False
                        parentLabel   = 'kk'
                    lastEventStartA = False
            else:
                raise Exception ("Unknown event type %s"%event)
            root.clear()
        #Now loop Trough all objects and fix the FK
        for object in listObjectWithForeignKey:
            for key in object._foreignKeys:
                object.setForeignKey(self.emxData.getObjectwithPK(object._getXMLForeignKey(key)))

    def createObject(self,elem):
        self.myClass = self.classObject[elem.tag]
        #primary key 
        #get PK
        self.dict         = self.readObjectPK(elem)
        #create object
        self._object = self.classObject[elem.tag](**(self.dict))
        #add it to emxData
        self.emxData.addObject(self._object)
        
    def readObjectPK(self, elem):
        ''' read primary key. So far all entities has the
        same PK. We may need to specialize or use dictPrimaryKeys
        in the future
        '''
        mapPK=collections.OrderedDict()
        for attribute in elem.attrib:
            mapPK[attribute] = emxDataTypes[attribute].getType()(elem.get(attribute))
        if mapPK:
            return collections.OrderedDict(sorted(mapPK.items(), key=lambda t: t[0]))
        else:
            raise Exception("readObjectPK: No fileName or index provided" )

    def firstObject(self, classname, fileName):
        """ Iterate over the tags elements and find the 
        first one of type 'classname', build the object
        and return it. The foreing keys will be not updated.
        """
        context = ET.iterparse(fileName, events=('start','end'))
        for event, elem in iter(context):
            tag = elem.tag
            if event == 'start':
                print "tag: '%s'" % tag, "class: '%s'" % classname
                if tag == classname:
                    print "tag==class"
                    self.createObject(elem)
                    print "self._object: ", self._object
                    return self._object
        return None
        
    def writeEMXFile(self, fileName):
        """read xml file and store it in a document
        """
        xmlFile = open(fileName, "w")
        xmlFile.write("<?xml version='1.0' encoding='utf-8'?>\n")
        xmlFile.write('''<EMX version="1.0">
  <!--
  ##########################################################################
  #               EMX Exchange file 
  #               Produced by the emxLibrary module
  #
  #
  #  Information on this file format is available at 
  #  http://i2pc.cnb.csic.es/emx
  ##########################################################################
  #  One of the best ways you can help us to improve this software
  #  is to let us know about any problems you find with it.
  #  Please report bugs to: emx@cnb.csic.es
  ##########################################################################
  -->
''')
        for object in self.emxData:
            text = self.objectToXML(object)# 
#            #implement this with a regular expression
#            #format matrices properly
            for i, j in {
                          '</t11>\n    ':'</t11> '
                         ,'</t12>\n    ':'</t12> '
                         ,'</t13>\n    ':'</t13> '
                         ,'</t21>\n    ':'</t21> '
                         ,'</t22>\n    ':'</t22> '
                         ,'</t23>\n    ':'</t23> '
                         ,'</t31>\n    ':'</t31> '
                         ,'</t32>\n    ':'</t32> '
                         ,'</t33>\n    ':'</t33> '
                         }.iteritems():
                text = text.replace(i, j)
            xmlFile.write( text)
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
    #try xerces
    ###########
    answerSize=1024 # avoid overflow in web
    endding=''
    if schema_file is None:
        _schema = EMXSCHEMA11
    else:
        _schema = schema_file
    #print "java jaxp.SourceValidator -a %s -i %s -xsd11"% (_schema, filename)
    p = subprocess.Popen("java jaxp.SourceValidator -a %s -i %s -xsd11" 
                         % (_schema, filename),
            shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    #xerces exists but is error
    if p.returncode == 0 and (stderr!=""):
        if len(stderr) > answerSize:
           endding='... (too many errors, displayed first %d characters)'%(answerSize)
        print ">>>>>>>>>>>>>>>>>>>>>>>>1"
        raise ValidateError(p.returncode, """Error: when validating file %s with schema %s.
        \nError:%s"""%(filename,_schema,stderr[:answerSize]+endding))
    #######
    #no xerces available, let us try xmlint
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
            print ">>>>>>>>>>>>>>>>>>>>>>>>2"
            raise ValidateError(127,
                """Error: neither xerces-f nor xmllint could be found,  I cannot validate schema. 
    Schema validation is based either on the xmllint program that belongs to the libxml2-tools package.
    or on the xerces-f project""")
            
        if p.returncode != 0:
            if len(stderr) > answerSize:
                 endding='... (too many errors, displayed first %d characters)'%(answerSize)
            print ">>>>>>>>>>>>>>>>>>>>>>>>3"
            message = """Error: when validating file %s with schema %s.
            \nError:%s"""%(filename,_schema,stderr[:answerSize]+endding)
            print "message", message
            raise ValidateError(p.returncode, message)
    return p.returncode, stdout[:answerSize], stderr[:answerSize]

