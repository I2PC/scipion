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
        self._size = pwobj.Integer(0)
        self._acquisition = XrayAcquisition()
        self._focalSeries = FocalSeries()
        
    def getAngles(self):
        return self._angles.get()
    
    def setAngles(self, angles):
        """ Set the angles of the tilt series.
        Params:
            angles: Python list with angles. """
        self._angles.set(angles)
    
    def getSize(self):
        return  self._size.get()   
    
    def setSize(self, value):
        self._size.set(value)             
    
    def getXrayAcquisition(self):
        return self._acquisition

    def setXrayAcquisition(self, acquisition):
        self._acquisition = acquisition
        
    def getFocalSeries(self):
        return self._focalSeries

    def setFocalSeries(self, focalSeries):
        self._focalSeries = focalSeries
        
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
        return "\n    lensLabel=%s\n    energy= %f\n    date=%s\n\n" % (
                self._lensLabel.get(),
                self._energy.get(),
                self._date.get()
                )
    
class FocalSeries(em.EMObject):
    """Info related to the setOfTiltSeries"""   
    def __init__(self, **kwargs):
        em.EMObject.__init__(self, **kwargs)
        self._tiltSeriesGroup = pwobj.Integer(kwargs.get('tiltSeriesGroup', None))
        self._index = pwobj.Integer(kwargs.get('index', None))
        self._defocus = pwobj.Float(kwargs.get('defocus', None))
        self._reference = pwobj.Integer(kwargs.get('reference', None))
        
    def copyInfo(self, other):
        self.copyAttributes(other, '_tiltSeriesGroup', '_index', '_defocus', '_reference')
    
    def gettiltSeriesGroup(self):
        return self._tiltSeriesGroup.get()
    
    def settiltSeriesGroup(self, tiltSeriesGroup):
        self._tiltSeriesGroup.set(tiltSeriesGroup)
        
    def getIndex(self):
        return self._index.get()
    
    def setIndex(self, index):
        self._index.set(index)
        
    def getDefocus(self):
        return self._defocus.get()
    
    def setDefocus(self, defocus):
        self._defocus.set(defocus)
        
    def getReference(self):
        return  self._reference.get()   
    
    def setReference(self, reference):
        self._reference.set(reference)
        
    def __str__(self):
        return "\n    tiltSeriesGroup = %d\n    index = %d\n    defocus = %f\n    reference = %d\n\n" % (
                self._tiltSeriesGroup.get(),
                self._index.get(),
                self._defocus.get(),
                self._reference.get()
                )

class SetOfTiltSeries(em.SetOfImages):
    
    ITEM_TYPE = TiltSeries
    
    def __init__(self, **kwargs):        
        em.SetOfImages.__init__(self, **kwargs)
        
    def append(self, tiltSeries):
        """ Add a tilt series to the set. """
        if self.getSize() == 0: # only check this for first time append is called
            if self._firstDim.isEmpty():
                self._firstDim.set(tiltSeries.getDim())
        em.EMSet.append(self, tiltSeries)
        
    
    def iterItems(self, orderBy='id', direction='ASC'):
        
        for tiltSeries in em.Set.iterItems(self, orderBy=orderBy, direction=direction):
            yield tiltSeries

