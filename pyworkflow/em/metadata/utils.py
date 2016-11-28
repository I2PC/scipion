# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Laura del Cano (ldelcano@cnb.csic.es)
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
This module contains utils functions to operate over xmipp metadata files.
"""

from classes import MetaData, Row
from constants import LABEL_TYPES, MDL_ITEM_ID
from functions import labelType, getBlocksInMetaDataFile


def label2Python(label):
    """ Return the Python type (int, float, bool) for a given 
    metadata label (LABEL_INT, LABEL_DOUBLE..etc)
    """
    return LABEL_TYPES.get(labelType(label), str)


def getFirstRow(mdOrFn):
    """ Return the first object of a metadata.
    Params:
        mdOrFn: you can pass a metadata or a filename as argument.
    """

    if isinstance(mdOrFn, basestring):
        md = MetaData()
        md.read(mdOrFn, 1)
    else: # mdOrFn is MetaData
        md = mdOrFn
        
    if md.getParsedLines():
        row = Row()
        row.readFromMd(md, md.firstObject())
    else:
        row = None
    
    return row


def getSize(filename):
    """ Return the metadata size without parsing entirely. """
    md = MetaData()
    md.read(filename, 1)
    return md.getParsedLines()


def isEmpty(filename):
    """ Use getMdSize to check if metadata is empty. """
    return getSize(filename) == 0


def iterRows(md, sortByLabel=None):
    """ Iterate over the rows of the given metadata.
    Params:
        md: a MetaData object or a filename (MetaData will be read)
        sortByLabel: a label to sort the metadata before iterate.
    """
    # If md is string, take as filename and create the metadata

    if isinstance(md, basestring):
        md = MetaData(md)

    if sortByLabel is not None:
        md.sort(sortByLabel)

    row = Row()
    
    for objId in md:
        row.readFromMd(md, objId)
        yield row


def joinBlocks(inputMd, blockPrefix=None):
    mdImages = MetaData()
    mdAll = MetaData()
    mdBlocks = getBlocksInMetaDataFile(inputMd)
    
    for mdBlock in mdBlocks:
        if blockPrefix is not None:
            if mdBlock.startswith(blockPrefix):
                mdImages.read(mdBlock + "@" + inputMd)
                mdAll.unionAll(mdImages)
        else:
            mdImages.read(mdBlock + "@" + inputMd)
            mdAll.unionAll(mdImages)
    return mdAll


class SetMdIterator():
    """ Class to iterate over an input set and skip
    elements not present in metadata.
    This class can be used in copyItems when the number
    of elements in the set is higher that in metadata.
    """
    def __init__(self, md, sortByLabel=None, 
                 keyLabel=MDL_ITEM_ID,
                 updateItemCallback=None):
        
        if updateItemCallback is None:
            raise Exception('Set an updateItemCallback')
        
        self.iterMd = iterRows(md, sortByLabel) 
        self.keyLabel = keyLabel
        self.updateItemCallback = updateItemCallback   
        self.__nextRow()
        
    def __nextRow(self):
        try:
            self.lastRow = next(self.iterMd)
        except StopIteration:
            self.lastRow = None
            
    def updateItem(self, item, row):
        """ This function should be passed to copyItems
        as callback and it will filter the items
        not present in the metadata.
        """
        row = self.lastRow
        if (row is None or
            item.getObjId() != row.getValue(self.keyLabel)):
            item._appendItem = False
        
        else:
            item._appendItem = True
            self.updateItemCallback(item, row)
            self.__nextRow()
                
             
    
