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

import pyworkflow.object as pwobj
import pyworkflow.em as em

    
class TiltSeries(em.Image):
    
    def __init__(self, **kwargs):        
        em.Image.__init__(self, **kwargs)
        self._angles = pwobj.CsvList(pType=float)
        self._focalSeries = None
        self._acquisition = XrayAcquisition()
        self._size = pwobj.Integer(0)      
        
    def getAngles(self):
        return self._angles.get()
    
    def setAngles(self, angles):
        """ Set the angles of the tilt series.
        Params:
            angles: Python list with angles. """
        self._angles.set(angles)
    
    def hasXrayAcquisition(self):
        return (self._acquisition is not None and
                self._acquisition.getLensLabel() is not None and
                self._acquisition.getEnergy() is not None and
                self._acquisition.getDate() is not None        
                )
        
    def getXrayAcquisition(self):
        return self._acquisition

    def setXrayAcquisition(self, acquisition):
        self._acquisition = acquisition
        
    def getSize(self):
        return  self._size.get()   
    
    def setSize(self, value):
        self._size.set(value)             
    
class XrayAcquisition(em.EMObject):
    """Soft Xray acquisition information"""
    def __init__(self, **kwargs):
        em.EMObject.__init__(self, **kwargs)
        self._lensLabel = pwobj.String(kwargs.get('lensLabel', None))
        self._energy = pwobj.Float(kwargs.get('energy', None))
        self._date = pwobj.String(kwargs.get('date', None))
        
    def copyInfo(self, other):
        self.copyAttributes(other, '_lensLabel', '_energy', '_date')
        
    def getLensLabel(self):
        return self._lensLabel.get()
        
    def setLensLabel(self, value):
        self._lensLabel.set(value)
        
    def getEnergy(self):
        return self._energy.get()
        
    def setEnergy(self, value):
        self._energy.set(value)
        
    def getDate(self):
        return self._date.get()
    
    def setDate(self, value):
        self._date.set(value)
        
    def __str__(self):
        return "\n    lensLabel=%s\n    energy= %f\n    date=%d\n\n" % (
                self._lensLabel.get(),
                self._energy.get(),
                self._date.get()
                )
    
#class FocalSeries(em.EmObject):
#    ITEM_TYPE = TiltSeries
   
#    def __init__(self, **kwargs):
#        TiltSeries.__init__(self, **kwargs)
        

    #def append(self, image):
    #    """ Add a image to the set. """
    #    if self.hasXrayAcquisition(): 
    #        self.setXrayAcquisition(self.getXrayAcquisition())
    #    if self.getSize() == 0: 
    #        if self._firstDim.isEmpty():
    #            self._firstDim.set(image.getDim())
    #    em.EMSet.append(self, image)   
        


