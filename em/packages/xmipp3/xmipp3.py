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
This sub-package will contains Xmipp3.0 specific protocols
"""

import os
import xmipp


def getXmippPath(*paths):
    '''Return the path the the Xmipp installation folder
    if a subfolder is provided, will be concatenated to the path'''
    if os.environ.has_key('XMIPP_HOME'):
        return os.path.join(os.environ['XMIPP_HOME'], *paths)  
    else:
        raise Exception('XMIPP_HOME environment variable not set')
    
    
class XmippProtocol():
    """ This class groups some common functionalities that
    share some Xmipp protocols, like converting steps.
    """
    
    def _insertConvertStep(self, inputName, xmippClass, resultFn):
        """ Insert the convertInputToXmipp if the inputName attribute
        is not an instance of xmippClass.
        It will return the result filename, if the 
        convertion is needed, this will be input resultFn.
        If not, it will be inputAttr.getFileName()
        """
        inputAttr = getattr(self, inputName)
        if not isinstance(inputAttr, xmippClass):
            self._insertFunctionStep('convertInputToXmipp', inputName, resultFn)
            return resultFn
        return inputAttr.getFileName()
        
    def convertInputToXmipp(self, inputName, xmippClass, resultFn):
        """ This step can be used whenever a convertion is needed.
        It will receive the inputName and get this attribute from self,
        invoke the convert function and check the result files if
        convertion was done (otherwise the input was already in Xmipp format).
        """
        inputAttr = getattr(self, inputName)
        inputXmipp = xmippClass.convert(inputAttr, resultFn)
        
        if inputXmipp != inputAttr:
            self._insertChild(inputName + 'Xmipp', inputXmipp)
            return [resultFn] # validate resultFn was produced if converted
        
    def getConvertedInput(self, inputName):
        """ Retrieve the converted input, it can be the case that
        it is the same as input, when not convertion was done. 
        """
        return getattr(self, inputName + 'Xmipp', getattr(self, inputName))
        

class XmippMdRow():
    """ Support Xmipp class to tore label and value pairs 
    corresponding to a Metadata row. It can be used as base
    for classes that maps to a MetaData row like XmippImage, XmippMicrograph..etc. 
    """
    def __init__(self):
        self._labelDict = {} # Dictionary containing labels and values
    
    def setValue(self, label, value):
        """args: this list should contains tuples with 
        MetaData Label and the desired value"""
        self._labelDict[label] = value
            
    def getValue(self, label):
        return self._labelDict[label]
    
    def getFromMd(self, md, objId):
        """ Get all row values from a given id of a metadata. """
        self._labelDict.clear()
        for label in md.getActiveLabels():
            self._labelDict[label] = md.getValue(label, objId)
            
    def setToMd(self, md, objId):
        """ Set back row values to a metadata row. """
        for label, value in self._labelDict.iteritems():
            # TODO: Check how to handle correctly unicode type
            # in Xmipp and Scipion
            if type(value) is unicode:
                value = str(value)
            md.setValue(label, value, objId)
            
  
class XmippSet():
    """ Support class to store sets in Xmipp base on a MetaData. """
    def __init__(self):
        """ Create new set, base on a Metadata.
        itemLabel: main key xmipp label of the items
        itemClass: Class that represent the items.
        A method .getFileName should be available to store the md.
        Items contained in XmippSet are suposed to inherit from XmippMdRow.
        """
        self._itemLabel = self._getItemLabel()
        self._itemClass = self._getItemClass()
        self._md = xmipp.MetaData()
        
    def _getItemLabel(self):
        """ Should be implemented to return items MetaData label. """
        pass
    
    def _getItemClass(self):
        """ Should be implemented to return items Class. """
        pass
        
    def __iter__(self):
        """Iterate over the set of images in the MetaData"""
        self._md.read(self.getFileName())
        
        for objId in self._md:  
            item = self._itemClass()
            item.readFromMd(self._md, objId)  
            #m = Image(md.getValue(xmipp.MDL_IMAGE, objId))
            #if self.hasCTF():
            #    m.ctfModel = XmippCTFModel(md.getValue(xmipp.MDL_CTF_MODEL, objId)) 
            yield item
            
    def setMd(self, md):
        self._md = md
        
    def sort(self):
        """Sort the set according to MDL_IMAGE"""
        self._md.sort(self._itemLabel)
        
    def write(self):
        self._md.write(self.getFileName())
        
    def append(self, item):
        """Add a new item to the set"""
        objId = self._md.addObject()
        # Convert to xmipp micrograph if necessary
        itemXmipp = self._itemClass.convert(item)
        itemXmipp.setToMd(self._md, objId)
        
    @staticmethod
    def convert(inputSet, xmippSetClass, filename):
        """ Convert from a generic set to a XmippSet subclass(xmippSetClass).
        In particular a filename is requiered to store the result MetaData.
        It is also asummed that this class have a .copyInfo method.
        """
        if isinstance(inputSet, xmippSetClass):
            return inputSet
        
        setOut = xmippSetClass(filename)
        setOut.copyInfo(inputSet)
        
        for item in inputSet:
            setOut.append(item)
        setOut.write()
        
        return setOut
                  
