# **************************************************************************
# *
# * Authors:     Mohsen Kazemi  (mkazemi@cnb.csic.es)
# *              
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
 
import pyworkflow.object as pwobj
import pyworkflow.em as em
from h5py import File


    
class TiltSeries(em.Image):    
    
    def __init__(self, **kwargs):        
        em.Image.__init__(self, **kwargs)
        self._angles = pwobj.CsvList()
        self._normalized = pwobj.Boolean()
        self._focalSeries = None    
      
        
    def getAngles(self):
        return self._angles.get()
    
    def setNormalized(self,fileName):
        fhHdf5 = File(fileName, 'r')
        if "TomoNormalized" in fhHdf5:
            self._normalized.set(True)
        else:    
            self._normalized.set(False)
            
    def getNormalized(self):
        return self._normalized.get()
    
    

  
#class FocalSeries(em.EmObject):
#    ITEM_TYPE = TiltSeries
   
#    def __init__(self, **kwargs):
#        TiltSeries.__init__(self, **kwargs)
        





