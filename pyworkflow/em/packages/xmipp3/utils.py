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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
"""
This module contains utils functions to operate over xmipp metadata files.
"""

import xmipp
from xmipp3 import XmippMdRow


def getMdFirstRow(filename):
    """ Create a MetaData but only read the first row.
    This method should be used for validations of labels
    or metadata size, but the full metadata is not needed.
    """
    md = xmipp.MetaData()
    md.read(filename, 1)
    if md.getParsedLines():
        row = XmippMdRow()
        row.readFromMd(md, md.firstObject())
    else:
        row = None
    
    return row


def getMdSize(filename):
    """ Return the metadata size without parsing entirely. """
    md = xmipp.MetaData()
    md.read(filename, 1)
    return md.getParsedLines()


def isMdEmpty(filename):
    """ Use getMdSize to check if metadata is empty. """
    return getMdSize(filename) == 0


def iterMdRows(md):
    """ Iterate over the rows of the given metadata. """
    # If md is string, take as filename and create the metadata
    if isinstance(md, basestring):
        md = xmipp.MetaData(md)
        
    row = XmippMdRow()
    
    for objId in md:
        row.readFromMd(md, objId)
        yield row
