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
This module contains converter functions 
from base classes to xmipp classes
"""

from data import XmippMicrograph, XmippSetOfMicrographs
import xmipp
    
def convertMicrograph(mic):
    """Convert from Micrograph to XmippMicrograph"""
    if type(mic) is XmippMicrograph:
        return mic
    
    micXmipp = XmippMicrograph(mic.getFileName())
    # TODO: copyInfo??
    # from mic to micXmipp??  
    return micXmipp
       
def convertSetOfMicrographs(setOfMics, filename):
    """Method to convert from a general SetOfMicrographs to XmippSetOfMicrographs"""
    if type(setOfMics) is XmippSetOfMicrographs:
        return setOfMics
        
    micsOut = XmippSetOfMicrographs(filename)
    micsOut.copyInfo(setOfMics)

    for mic in setOfMics:
        micsOut.append(mic)

    micsOut.write()
        
    return micsOut
    
    
def convertCTFModel(setOfMics, filename):
    """Method to convert from a general SetOfMicrographs to XmippSetOfMicrographs"""
    if type(setOfMics) is XmippSetOfMicrographs:
        return setOfMics
        
    return None
