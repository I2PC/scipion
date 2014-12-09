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
from constants import LABEL_TYPES
from functions import labelType


def label2Python(label):
    """ Return the Python type (int, float, bool) for a given 
    metadata label (LABEL_INT, LABEL_DOUBLE..etc)
    """
    labelType = labelType(label)
    return LABEL_TYPES.get(labelType, str)


def getFirstRow(filename):
    """ Create a MetaData but only read the first row.
    This method should be used for validations of labels
    or metadata size, but the full metadata is not needed.
    """
    md = MetaData()
    md.read(filename, 1)
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


def iterRows(md):
    """ Iterate over the rows of the given metadata. """
    # If md is string, take as filename and create the metadata
    if isinstance(md, basestring):
        md = MetaData(md)
        
    row = Row()
    
    for objId in md:
        row.readFromMd(md, objId)
        yield row
