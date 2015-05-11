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
Add functions related to metadata
"""

import sys
from collections import OrderedDict

from pyworkflow.object import ObjectWrap
from xmipp import MetaData, label2Str, str2Label, MD_APPEND



class Row():
    """ Support Xmipp class to store label and value pairs 
    corresponding to a Metadata row. 
    """
    def __init__(self):
        self._labelDict = OrderedDict() # Dictionary containing labels and values
        self._objId = None # Set this id when reading from a metadata
        
    def getObjId(self):
        return self._objId
    
    def hasLabel(self, label):
        return self.containsLabel(label)
    
    def containsLabel(self, label):
        # Allow getValue using the label string
        if isinstance(label, basestring):
            label = str2Label(label)
        return label in self._labelDict
    
    def removeLabel(self, label):
        if self.hasLabel(label):
            del self._labelDict[label]
    
    def setValue(self, label, value):
        """args: this list should contains tuples with 
        MetaData Label and the desired value"""
        # Allow setValue using the label string
        if isinstance(label, basestring):
            label = str2Label(label)
        self._labelDict[label] = value
            
    def getValue(self, label, default=None):
        """ Return the value of the row for a given label. """
        # Allow getValue using the label string
        if isinstance(label, basestring):
            label = str2Label(label)
        return self._labelDict.get(label, default)
    
    def getValueAsObject(self, label, default=None):
        """ Same as getValue, but making an Object wrapping. """
        return ObjectWrap(self.getValue(label, default))
    
    def readFromMd(self, md, objId):
        """ Get all row values from a given id of a metadata. """
        self._labelDict.clear()
        self._objId = objId
        
        for label in md.getActiveLabels():
            self._labelDict[label] = md.getValue(label, objId)
            
    def writeToMd(self, md, objId):
        """ Set back row values to a metadata row. """
        for label, value in self._labelDict.iteritems():
            # TODO: Check how to handle correctly unicode type
            # in Xmipp and Scipion
            if type(value) is unicode:
                value = str(value)
            try:
                md.setValue(label, value, objId)
            except Exception, ex:
                print >> sys.stderr, "XmippMdRow.writeToMd: Error writting value to metadata."
                print >> sys.stderr, "                     label: %s, value: %s, type(value): %s" % (label2Str(label), value, type(value))
                raise ex
            
    def readFromFile(self, fn):
        md = MetaData(fn)
        self.readFromMd(md, md.firstObject())
        
    def writeToFile(self, fn):
        md = MetaData()
        self.writeToMd(md, md.addObject())
        md.write(fn)
        
    def copyFromRow(self, other):
        for label, value in other._labelDict.iteritems():
            self.setValue(label, value)
            
    def clone(self):
        """ Return another Row that have exactly the same
        values as self.
        """
        row = Row()
        row.copyFromRow(self)
        row._objId = self._objId
    
        return row
        
    def __str__(self):
        s = '{'
        for k, v in self._labelDict.iteritems():
            s += '  %s = %s\n' % (label2Str(k), v)
        return s + '}'
    
    def __iter__(self):
        return self._labelDict.iteritems()
        
    def containsAll(self, labels):
        """ Check if all labels are present in the row.
        Params:
            row: the Row object.
            labels: either a dict or list object containing the labels
                (in the case of dicts, label are the dict.values())
        """
        values = labels.values() if isinstance(labels, dict) else labels
        return all(self.containsLabel(l) for l in values)
    
    def containsAny(self, labels):
        """ Check if at least one of labels is present in the row.
        Params:
            row: the Row object.
            labels: either a dict or list object containing the labels
                (in the case of dicts, label are the dict.values())
        """
        values = labels.values() if isinstance(labels, dict) else labels
        return any(self.containsLabel(l) for l in values)
            
    def printDict(self):
        """ Fancy printing of the row, mainly for debugging. """
        print str(self)
        
        
class RowMetaData():
    """ This class is a wrapper for MetaData in row mode.
    Where only one object is used.
    """
    def __init__(self, filename=None):
        self._md = MetaData()
        self._md.setColumnFormat(False)
        self._id = self._md.addObject()
        
        if filename:
            self.read(filename)
        
    def setValue(self, label, value):
        self._md.setValue(label, value, self._id)
        
    def getValue(self, label):
        return self._md.getValue(label, self._id)
        
    def write(self, filename, mode=MD_APPEND):
        self._md.write(filename, mode)
        
    def read(self, filename):
        self._md.read(filename)
        self._md.setColumnFormat(False)
        self._id = self._md.firstObject()
        
    def containsLabel(self, label):
        return self._md.containsLabel(label)
     
    def __str__(self):
        return str(self._md)