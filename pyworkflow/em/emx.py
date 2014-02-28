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
from data import *
from pyworkflow import emx
from collections import OrderedDict
    

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
    from convert import ImageHandler
    from constants import NO_INDEX
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
    
    
    
    
def _getFloat(elem, key):
    """ Parse a float value from elem child text. """
    return float(elem.find(key).text.strip())
    
def _setLocationFromElem(elem, img):
    """ Set image location from attributes "index" and "fileName". """
    index = elem.get('index', 1)
    img.setLocation(index, elem.get('fileName'))
    
def _setSamplingFromElem(elem, img):
    """ Set the sampling rate from the pixelSpacing element. """
    s = elem.find('pixelSpacing')
    if s is None:
        print "NO pixelSpacing:", elem
    else:
        img.setSamplingRate(_getFloat(s, 'X'))
        
def _setCoordinatesFromElem(elem, coordinate):
    c = elem.find('centerCoord')
    if c is not None:
        coordinate.setX(_getFloat(c, 'X'))
        coordinate.setY(_getFloat(c, 'Y'))
    
def _acquisitionFromElem(elem, acquisition):
    """ Create an acquistion from elem. """
    acquisition.setVoltage(_getFloat(elem, 'acceleratingVoltage'))
    acquisition.setAmplitudeContrast(_getFloat(elem, 'amplitudeContrast'))
    acquisition.setSphericalAberration(_getFloat(elem, 'cs'))    
    
def _ctfFromElem(elem, ctf):
    """ Create a CTF model from the values in elem. """
    for tag in ['defocusU', 'defocusV', 'defocusUAngle']:
        if elem.find(tag) is None:
            return None
    ctf.setDefocusU(_getFloat(elem, 'defocusU'))
    ctf.setDefocusV(_getFloat(elem, 'defocusV'))
    ctf.setDefocusAngle(_getFloat(elem, 'defocusUAngle'))
    
def _imageFromElem(elem, img):
    """ Create either Micrograph or Particle from xml elem. """
    _setLocationFromElem(elem, img)
    _setSamplingFromElem(elem, img)
    _ctfFromElem(elem, img.getCTF())
    
def _micrographFromElem(elem, mic):
    """ Create a micrograph and set the values properly
    from the corresponding xml tree element.
    """
    _imageFromElem(elem, mic)
    _acquisitionFromElem(elem, mic.getAcquisition())
    
    #mic.printAll()
    
def _particleFromElem(elem, particle):
    _imageFromElem(elem, particle)
    _setCoordinatesFromElem(elem, particle.getCoordinate())
    
def _coordinateFromElem(elem, coordinate):
    #_imageFromElem(elem, coordinate)
    _setCoordinatesFromElem(elem, coordinate)
  
def _iterXml(xmlFn, tagsCallback):
    """ Iterate incrementally a give XML file. 
    This implementation is used having in mind huge xml files.
    Recipe taken from:
    http://effbot.org/zone/element-iterparse.htm
    Params:
        xmlFn: filename of the xml file.
    """
    # get an iterable
    context = ET.iterparse(xmlFn, events=("start", "end"))
    # turn it into an iterator
    context = iter(context)
    # get the root element
    _, root = context.next()
    
    # This flag will serve to check whether we are parsing micrographs
    # or we are parsing particles or coordinates
    # After the first <particle> tag is found, the flag should be set to False
    onlyMics = True
    # Check if pixelSpacing is present to create Particle or Coordinate
    foundPixelSpacing = False
    # Create single instances of objects that will be populated with 
    # the values parsed from the EMX file
    mic = Micrograph()
    mic.setAcquisition(Acquisition())    
    coordinate = Coordinate()
    particle = Particle() # particle should be of type Particle or Coordinate
    particle.setCoordinate(coordinate)
    
    for event, elem in context:
        if event == 'start':
            # After found the first particle, set onlyMics flag to False
            # after this, micrographs tags are just particle references
            if elem.tag == 'particle':
                onlyMics = False
            # If pixelSpacing is found(and we require it for particles)
            # we asume that element are Particle, if not, there are Coordinates
            if elem.tag == 'pixelSpacing':
                if not onlyMics:
                    foundPixelSpacing = True
            # If found CTF information, set the CTFModel
            # only once for either micrograph or particle
            if elem.tag == 'defocusU':
                if onlyMics:
                    if not mic.hasCTF():
                        mic.setCTF(CTFModel())
                else:
                    if not particle.hasCTF():
                        particle.setCTF(CTFModel())
        else: # event == 'end'
            if elem.tag == 'micrograph':
                if onlyMics:
                    _micrographFromElem(elem, mic)
                    root.clear()
                else:
                    pass #TODO: set particle micrograph reference
            elif elem.tag == 'particle':
                if foundPixelSpacing:
                    _particleFromElem(elem, particle)
                    #kparticle.printAll()
                else:
                    _coordinateFromElem(elem, coordinate)
                    #coordinate.printAll()
                root.clear()
    