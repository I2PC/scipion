# **************************************************************************
# *
# * Authors:     J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
# *  e-mail address 'jgomez@cnb.csic.es'
# *
# **************************************************************************
"""
This module implement the import/export of Micrographs and Particles to EMX
"""
from pyworkflow.utils import *
from pyworkflow.em.data import *
import emx
from collections import OrderedDict
    
        
def exportData(emxDir, inputSet, ctfSet=None):
    """ Export micrographs, coordinates or particles to  EMX format. """
    cleanPath(emxDir)
    makePath(emxDir) 
    emxData = emx.EmxData()

    if isinstance(inputSet, SetOfMicrographs):
        emxSetOfMicrographs(emxData, inputSet, emxDir, ctfSet)
        
    elif isinstance(inputSet, SetOfCoordinates):
        micSet = inputSet.getMicrographs()
        emxSetOfMicrographs(emxData, micSet, emxDir, ctfSet)
        emxSetOfParticles(emxData, inputSet, None, micSet)
        
    elif isinstance(inputSet, SetOfParticles):
        if inputSet.hasCoordinates():
            micSet = inputSet.getCoordinates().getMicrographs()
            emxSetOfMicrographs(emxData, micSet, emxDir, writeData=False)
        fnMrcs = join(emxDir, 'data.mrc')
        emxSetOfParticles(emxData, inputSet, fnMrcs, micSet)
        
    fnXml = join(emxDir, 'data.emx')
    emxData.write(fnXml)
    
    
def importData(protocol, emxFile):
    """ Import objects into Scipion from a give EMX file. 
    Returns:
        a dictionary with key and values as outputs sets
        (Micrographs, Coordinates or Particles)
    """
    emxData = emx.EmxData()
    emxData.read(emxFile)
    
    emxMic = emxData.getFirstObject(emx.MICROGRAPH)
    outputDict = {}
    
    if emxMic is not None:
        hasMicrographs = True
        micSet = protocol._createSetOfMicrographs()
        mic = Micrograph()
        mic.setAcquisition(Acquisition())
        if emxMic.has('acceleratingVoltage'):
            mic.setCTF(CTFModel())
        _micrographFromEmx(emxMic, mic)
        acq = mic.getAcquisition().clone()
        if acq.getMagnification() is None:
            acq.setMagnification(6000)
            acq.setScannedPixel
        micSet.setAcquisition(acq)        
        micSet.setSamplingRate(mic.getSamplingRate())
    
        for emxMic in emxData.iterClasses(emx.MICROGRAPH):
            _micrographFromEmx(emxMic, mic)
            micSet.append(mic)
            
        outputDict['outputMicrographs'] = micSet
        
    return outputDict


#---------------- Export related functions -------------------------------

def _writeDict(emxObj, dictValues):
    """ Set values of an EMX object. 
    Key-values pairs should be provided in dictValues.
    """
    for k, v in dictValues.iteritems():
        emxObj.set(k, v)
            
def _writeCTF(emxObj, ctf):
    """ Write the CTF values in the EMX object. """
    _writeDict(emxObj, {'defocusU': ctf.getDefocusU(),
                        'defocusV': ctf.getDefocusV(),
                        'defocusUAngle': ctf.getDefocusAngle()
                        })
    
def _writeSamplingRate(emxObj, sampling):
    """ Write sampling rate info as expected by EMX. """
    emxObj.set('pixelSpacing__X', sampling)
    emxObj.set('pixelSpacing__Y', sampling)
    
def _writeAcquisition(emxObj, acq):
    _writeDict(emxObj, {'acceleratingVoltage': acq.getVoltage(),
                        'amplitudeContrast': acq.getAmplitudeContrast(),
                        'cs': acq.getSphericalAberration()
                        })
    
def emxMicrograph(mic):
    """ Create an EMX micrograph and fill all the values. """
    index, fn = mic.getLocation()
    i = index or 1
    emxMic = emx.EmxMicrograph(fileName=fn, index=i)
    _writeAcquisition(emxMic, mic.getAcquisition())
    _writeSamplingRate(emxMic, mic.getSamplingRate())
    if mic.hasCTF():
        _writeCTF(emxMic, mic.getCTF())
    
    return emxMic
    
def _writeCenter(emxObj, coord):
    _writeDict(emxObj, {'X': coord.getX(), 'Y': coord.getY()})
    
def _emxParticle(emxData, coordinate, index, filename, micSet):
    """ This helper function will serve to write common xml string
    of Coordinate and Particle, actually, in EMX both are 'particle'
    """
    i = index or 1
    emxParticle = emx.EmxParticle(fileName=filename, index=i)
    if coordinate:
        mic = micSet[coordinate.getMicId()] # get micrograph
        index, filename = mic.getLocation()
        i = index or 1
        filename = basename(filename)
        mapKey = OrderedDict([('fileName', filename), ('index', i)])
        emxMic = emxData.getObject(mapKey)
        emxParticle.setMicrograph(emxMic)
        #TODO: ADD foreign key
        _writeDict(emxParticle, {'centerCoord__X': coordinate.getX(), 
                                 'centerCoord__Y': coordinate.getY()})

    return emxParticle
           
def emxParticle(emxData, particle, micSet):
    index, filename = particle.getLocation()
    emxParticle = _emxParticle(emxData, particle.getCoordinate(), index, filename, micSet)
    if particle.hasCTF():
        _writeCTF(emxParticle, particle.getCTF())
    _writeSamplingRate(emxParticle, particle.getSamplingRate())
    
    return emxParticle

def emxCoordinate(emxData, coordinate, micSet):
    return  _emxParticle(emxData, coordinate, coordinate.getObjId(), 
                         "coordinate", micSet)  
    
    
def emxSetOfMicrographs(emxData, micSet, emxDir, ctfSet=None, writeData=True):
    """ Write a SetOfMicrograph as expected in EMX format (xml file)
    Params:
        micSet: input set of micrographs
        filename: the EMX file where to store the micrographs information.
    """
    from pyworkflow.em.convert import ImageHandler
    from pyworkflow.em.constants import NO_INDEX
    
    ih = ImageHandler()

    for mic in micSet:
        loc = mic.getLocation()
        fnMicBase = replaceBaseExt(loc[1], 'mrc')
        if writeData:
            newLoc = join(emxDir, fnMicBase)
            ih.convert(loc, newLoc)
        mic.setLocation(NO_INDEX, fnMicBase)
        if ctfSet:
            mic.setCTF(ctfSet[mic.getObjId()])
        emxMic = emxMicrograph(mic)
        emxData.addObject(emxMic)
    
    
def emxSetOfParticles(emxData, partSet, stackFn=None, micSet=None):
    """ Write a SetOfMicrograph as expected in EMX format (xml file)
    Params:
        micSet: input set of micrographs
        filename: the EMX file where to store the micrographs information.
    """
    from pyworkflow.em.convert import ImageHandler
    
    ih = ImageHandler()
    
    for i, particle in enumerate(partSet):
        if stackFn:
            loc = particle.getLocation()
            newLoc = (i+1, stackFn)
            ih.convert(loc, newLoc)
            newFn = basename(stackFn)
            particle.setLocation(i+1, newFn)
            emxObj = emxParticle(emxData, particle, micSet)
        else:
            emxObj = emxCoordinate(emxData, particle, micSet)
        emxData.addObject(emxObj)
        


 #---------------- Export related functions ------------------------------- 
    
   
def _setLocationFromEmx(emxObj, img):
    """ Set image location from attributes "index" and "fileName". """
    img.setLocation(emxObj.get(emx.INDEX, 1), 
                    emxObj.get(emx.FILENAME))
    
def _setSamplingFromEmx(emxObj, img):
    """ Set the sampling rate from the pixelSpacing element. """
    img.setSamplingRate(emxObj.get('pixelSpacing'))
        
def _setCoordinatesFromEmx(emxObj, coordinate):
    if emxObj.has('centerCoord__X'):
        coordinate.setX(emxObj.get('centerCoord__X'))
        coordinate.setY(emxObj.get('centerCoord__Y'))
    
def _acquisitionFromEmx(emxObj, acquisition):
    """ Create an acquistion from elem. """
    acquisition.setVoltage(emxObj.get('acceleratingVoltage'))
    acquisition.setAmplitudeContrast(emxObj.get('amplitudeContrast'))
    acquisition.setSphericalAberration(emxObj.get('cs'))    
    
def _ctfFromEmx(emxObj, ctf):
    """ Create a CTF model from the values in elem. """
    for tag in ['defocusU', 'defocusV', 'defocusUAngle']:
        if not emxObj.has(tag):
            return None
    ctf.setDefocusU(emxObj.get('defocusU'))
    ctf.setDefocusV(emxObj.get('defocusV'))
    ctf.setDefocusAngle(emxObj.get('defocusUAngle'))
    
def _imageFromEmx(emxObj, img):
    """ Create either Micrograph or Particle from xml elem. """
    _setLocationFromEmx(emxObj, img)
    _setSamplingFromEmx(emxObj, img)
    _ctfFromEmx(emxObj, img.getCTF())
    
def _micrographFromEmx(emxMic, mic):
    """ Create a micrograph and set the values properly
    from the corresponding xml tree element.
    """
    _imageFromEmx(emxMic, mic)
    _acquisitionFromEmx(emxMic, mic.getAcquisition())
    
    #mic.printAll()
    
def _particleFromEmx(emxObj, particle):
    _imageFromEmx(emxObj, particle)
    _setCoordinatesFromEmx(emxObj, particle.getCoordinate())
    
def _coordinateFromEmx(emxObj, coordinate):
    #_imageFromEmx(emxObj, coordinate)
    _setCoordinatesFromEmx(emxObj, coordinate)
  
   