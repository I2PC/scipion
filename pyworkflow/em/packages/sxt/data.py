# **************************************************************************
# *
# * Authors:     Mohsen Kazemi  (mkazemi@cnb.csic.es)
# *              Joaquin Oton   (joton@cnb.csic.es)
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
This modules contains basic hierarchy
for ET data objects like: Tilt Series, Focal Series and others
"""

import os
from pyworkflow.object import *
#import pyworkflow.utils as pwutils

from constants import *
from pyworkflow.em.convert import *

from xmipp import *

"""
    
class TiltSeries(SetOfImages):
    """Represents an ET Tilt-series object"""
    
    def __init__(self, **kwargs):
        SetOfImages.__init__(self, **kwargs)
        
    
    
    def getFileName(self):
        """ Use the _objValue attribute to store filename. """
        return self._filename.get()
    
    def setFileName(self, filename):
        """ Use the _objValue attribute to store filename. """
        self._filename.set(filename) 
        
        
    def getAngles(self):
        """ Return ET angles """
        return self._angles.get()
    
    def getNormalized(self):
        """ Return true if data are normalized """
        return self._normalized.get()
    
    

  
class FocalSeries(TiltSeries):
    """ Represents a set of Tomograms """
    ITEM_TYPE = TiltSeries
    
    def __init__(self, **kwargs):
        TiltSeries.__init__(self, **kwargs)
        
"""




