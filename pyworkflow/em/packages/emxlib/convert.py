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

from __future__ import print_function
from os.path import join, dirname, basename, exists

from pyworkflow.utils.path import createLink, makePath, cleanPath, replaceBaseExt
from pyworkflow.em.convert import ImageHandler, NO_INDEX
from pyworkflow.em.data import (Micrograph, CTFModel, Particle, 
                                Coordinate, Transform,
                                SetOfMicrographs, SetOfCoordinates, 
                                SetOfParticles)
import emxlib
from collections import OrderedDict

        
def exportData(emxDir, inputSet, ctfSet=None, xmlFile='data.emx', binaryFile=None):
    """ Export micrographs, coordinates or particles to  EMX format. """
    cleanPath(emxDir)
    makePath(emxDir) 
    emxData = emxlib.EmxData()
    micSet=None
    
    if binaryFile is None:
        binaryFile = xmlFile.replace('.emx', '.mrc')

    if isinstance(inputSet, SetOfMicrographs):
        _micrographsToEmx(emxData, inputSet, emxDir, ctfSet)
        
    elif isinstance(inputSet, SetOfCoordinates):
        micSet = inputSet.getMicrographs()
        _micrographsToEmx(emxData, micSet, emxDir, ctfSet)
        _particlesToEmx(emxData, inputSet, None, micSet)
        
    elif isinstance(inputSet, SetOfParticles):
        if inputSet.hasCoordinates():
            micSet = inputSet.getCoordinates().getMicrographs()
            _micrographsToEmx(emxData, micSet, emxDir, writeData=False)
        fnMrcs = join(emxDir, binaryFile)
        _particlesToEmx(emxData, inputSet, fnMrcs, micSet)
        
    fnXml = join(emxDir, xmlFile)
    emxData.write(fnXml)
    
    
def importData(protocol, emxFile, outputDir, acquisition, 
               samplingRate=None, copyOrLink=createLink):
    """ Import objects into Scipion from a given EMX file. 
    Returns:
        a dictionary with key and values as outputs sets
        (Micrographs, Coordinates or Particles)
    """
    emxData = emxlib.EmxData()
    emxData.read(emxFile)

    _micrographsFromEmx(protocol, emxData, emxFile, outputDir, acquisition, 
                        samplingRate, copyOrLink)
    _particlesFromEmx(protocol, emxData, emxFile, outputDir, acquisition,
                      samplingRate, copyOrLink)


#---------------- Export related functions -------------------------------

def _dictToEmx(emxObj, dictValues):
    """ Set values of an EMX object. 
    Key-values pairs should be provided in dictValues.
    """
    for k, v in dictValues.iteritems():
        emxObj.set(k, v)
    
            
def _ctfToEmx(emxObj, ctf):
    """ Write the CTF values in the EMX object. """
    # Divide by 10 the defocus, since we have them in A
    # and EMX wants it in nm
    _dictToEmx(emxObj, {'defocusU': ctf.getDefocusU()/10.0,
                        'defocusV': ctf.getDefocusV()/10.0,
                        'defocusUAngle': ctf.getDefocusAngle()
                        })
    
    
def _samplingToEmx(emxObj, sampling):
    """ Write sampling rate info as expected by EMX. """
    emxObj.set('pixelSpacing__X', sampling)
    emxObj.set('pixelSpacing__Y', sampling)
    

def _acquisitionToEmx(emxObj, acq):
    _dictToEmx(emxObj, {'acceleratingVoltage': acq.getVoltage(),
                        'amplitudeContrast': acq.getAmplitudeContrast(),
                        'cs': acq.getSphericalAberration()
                        })
    
    
def _micrographToEmx(mic):
    """ Create an EMX micrograph and fill all the values. """
    index, fn = mic.getLocation()
    i = index or 1
    #TODO: this basename(fn) if potentially dangerous if two micrographs
    #have the same name but it is not easy to do a more general approach

    emxMic = emxlib.EmxMicrograph(fileName=basename(fn), index=i)
    _acquisitionToEmx(emxMic, mic.getAcquisition())
    _samplingToEmx(emxMic, mic.getSamplingRate())
    if mic.hasCTF():
        _ctfToEmx(emxMic, mic.getCTF())
    
    return emxMic
    
    
def _centerToEmx(emxObj, coord):
    _dictToEmx(emxObj, {'X': coord.getX(), 'Y': coord.getY()})
    
    
def _setupEmxParticle(emxData, coordinate, index, filename, micSet):
    """ This helper function will serve to setup common attributes
    of Coordinate and Particle, actually, in EMX both are 'particles'
    """
    i = index or 1
    emxParticle = emxlib.EmxParticle(fileName=filename, index=i)
    if coordinate:
        if micSet is not None:
            mic = micSet[coordinate.getMicId()] # get micrograph
            index, filename = mic.getLocation()
            i = index or 1
            filename = basename(filename)  ## Careful here if two micrographs have the same name...
            mapKey = OrderedDict([('fileName', filename), ('index', i)])

            emxMic = emxData.getObject(mapKey)
            emxParticle.setMicrograph(emxMic)
        #TODO: ADD foreign key
        _dictToEmx(emxParticle, {'centerCoord__X': coordinate.getX(), 
                                 'centerCoord__Y': coordinate.getY()})

    return emxParticle
           
           
def _particleToEmx(emxData, particle, micSet):
    #import pdb
    #pdb.set_trace()

    index, filename = particle.getLocation()
    emxParticle = _setupEmxParticle(emxData, particle.getCoordinate(), index, filename, micSet)
    if particle.hasCTF():
        _ctfToEmx(emxParticle, particle.getCTF())
    _samplingToEmx(emxParticle, particle.getSamplingRate())
    
    return emxParticle


def _coordinateToEmx(emxData, coordinate, micSet):
    return  _setupEmxParticle(emxData, coordinate, coordinate.getObjId(), 
                         "coordinate", micSet)  
    
    
def _micrographsToEmx(emxData, micSet, emxDir, ctfSet=None, writeData=True):
    """ Write a SetOfMicrograph as expected in EMX format (xml file)
    Params:
        micSet: input set of micrographs
        filename: the EMX file where to store the micrographs information.
    """
    ih = ImageHandler()

    for mic in micSet:
        if writeData:
            loc = mic.getLocation()
            fnMicBase = replaceBaseExt(loc[1], 'mrc')
            newLoc = join(emxDir, fnMicBase)
            ih.convert(loc, newLoc)
            mic.setLocation(NO_INDEX, fnMicBase)
        if ctfSet:
            mic.setCTF(ctfSet[mic.getObjId()])
        emxMic = _micrographToEmx(mic)
        emxData.addObject(emxMic)
    
    
def _particlesToEmx(emxData, partSet, stackFn=None, micSet=None):
    """ Write a SetOfMicrograph as expected in EMX format 
    Params:
        micSet: input set of micrographs
        filename: the EMX file where to store the micrographs information.
    """

    ih = ImageHandler()

    for i, particle in enumerate(partSet):
        if stackFn:
            loc = particle.getLocation()
            newLoc = (i+1, stackFn)
            ih.convert(loc, newLoc)
            newFn = basename(stackFn)
            particle.setLocation(i+1, newFn)
            emxObj = _particleToEmx(emxData, particle, micSet)
        else:
            emxObj = _coordinateToEmx(emxData, particle, micSet)
        emxData.addObject(emxObj)
        
        
#---------------Function related to import from EMX ---------------

def _setLocationFromEmx(emxObj, img):
    """ Set image location from attributes "index" and "fileName". """
    index = emxObj.get(emxlib.INDEX, 1)
#     if index == 1:
#         index = NO_INDEX
        
    img.setLocation(index, emxObj.get(emxlib.FILENAME))
    
    
def _setSamplingFromEmx(emxObj, img):
    """ Set the sampling rate from the pixelSpacing element. """
    img.setSamplingRate(emxObj.get('pixelSpacing__X'))
        
        
def _setCoordinatesFromEmx(emxObj, coordinate):
    if emxObj.has('centerCoord__X'):
        coordinate.setX(emxObj.get('centerCoord__X'))
        coordinate.setY(emxObj.get('centerCoord__Y'))
    
def _hasAcquisitionLabels(emxObj):
    """ Check that needed labels for CTF are present in object dict. """
    return (emxObj is not None and 
            all([emxObj.has(tag) for tag in ['acceleratingVoltage', 'amplitudeContrast', 'cs']]))
    
def _acquisitionFromEmx(emxObj, acquisition):
    """ Create an acquistion from elem. """
    if _hasAcquisitionLabels(emxObj):
        acquisition.setVoltage(emxObj.get('acceleratingVoltage'))
        acquisition.setAmplitudeContrast(emxObj.get('amplitudeContrast'))
        acquisition.setSphericalAberration(emxObj.get('cs'))    
    
    
def _hasCtfLabels(emxObj):
    """ Check that needed labels for CTF are present in object dict. """
    return (emxObj is not None and 
            all([emxObj.has(tag) for tag in ['defocusU', 'defocusV', 'defocusUAngle']]))
    
    
def _ctfFromEmx(emxObj, ctf):
    """ Create a CTF model from the values in elem. """
    if _hasCtfLabels(emxObj):
        # Multiply by 10 to convert from nm to A
        ctf.setStandardDefocus(emxObj.get('defocusU')*10.0,
                               emxObj.get('defocusV')*10.0,
                               emxObj.get('defocusUAngle'))
    else: # Consider the case of been a Particle and take CTF from Micrograph
        if isinstance(emxObj, emxlib.EmxParticle):
            _ctfFromEmx(emxObj.getMicrograph(), ctf)

    
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
    
    
def _particleFromEmx(emxObj, particle):
    _imageFromEmx(emxObj, particle)
    _setCoordinatesFromEmx(emxObj, particle.getCoordinate())
    
    
def _transformFromEmx(emxParticle, part, transform):
    """ Read the transformation matrix values from EMX tags. """
    m = transform.getMatrix()
    
    for i in range(3):
        for j in range(4):
            m[i, j] = emxParticle.get('transformationMatrix__t%d%d' % (i+1, j+1))
    
    transform.setObjId(part.getObjId())
            
def _coordinateFromEmx(emxObj, coordinate):
    #_imageFromEmx(emxObj, coordinate)
    _setCoordinatesFromEmx(emxObj, coordinate)
    emxMic = emxObj.getMicrograph()
    coordinate.setMicId(emxMic._micId)
    
    
def _micrographsFromEmx(protocol, emxData, emxFile, outputDir, 
                        acquisition, samplingRate, copyOrLink):
    """ Create the output SetOfMicrographs given an EMXData object.
    If there is information of the CTF, also the SetOfCTF will
    be registered as output of the protocol.
    """
    emxMic = emxData.getFirstObject(emxlib.MICROGRAPH)
    if emxMic is not None:
        micSet = protocol._createSetOfMicrographs()
        mic = Micrograph()
        mic.setAcquisition(acquisition)

        if _hasCtfLabels(emxMic):
            mic.setCTF(CTFModel())
            ctfSet = protocol._createSetOfCTF()
        else:
            ctfSet = None

        _micrographFromEmx(emxMic, mic)
        acq = mic.getAcquisition().clone()
        micSet.setAcquisition(acq)
        if not samplingRate:
            samplingRate = mic.getSamplingRate()
        micSet.setSamplingRate(samplingRate)
        micDir = dirname(emxFile)

    for emxMic in emxData.iterClasses(emxlib.MICROGRAPH):
        _micrographFromEmx(emxMic, mic)
        _, fn = mic.getLocation()
        micFn = join(micDir, fn)
        if copyOrLink is not None:
            micBase = basename(micFn)
            newFn = join(outputDir, micBase)
            copyOrLink(micFn, newFn)
            mic.setLocation(newFn)
        else:
            mic.setLocation(micFn)
        micSet.append(mic)
        emxMic._micId = mic.getObjId()
        mic.cleanObjId()
        if ctfSet is not None:
            ctf = mic.getCTF().clone()
            ctf.setMicrograph(mic)
            #TODO I do not think next line is needed
            #ctf.setMicFile(newFn)
            ctfSet.append(ctf)

    if emxMic is not None:
        protocol._defineOutputs(outputMicrographs=micSet)

        if ctfSet is not None:
            protocol._defineOutputs(outputCTF=ctfSet)
            ctfSet.setMicrographs(micSet)
            protocol._defineCtfRelation(micSet, ctfSet)


def _particlesFromEmx(protocol
                      , emxData
                      , emxFile
                      , outputDir
                      , acquisition
                      , samplingRate
                      , copyOrLink):
    """ Create the output SetOfCoordinates or SetOfParticles given an EMXData object.
    Add CTF information to the particles if present.
    """    
    emxParticle = emxData.getFirstObject(emxlib.PARTICLE)
    partDir = dirname(emxFile)
    micSet = getattr(protocol, 'outputMicrographs', None)

    if emxParticle is not None:     
        # Check if there are particles or coordinates
        fn = emxParticle.get(emxlib.FILENAME)
        if exists(join(partDir, fn)): # if the particles has binary data, means particles case
            partSet = protocol._createSetOfParticles()
            partSet.setAcquisition(acquisition)
            
            part = Particle()
            if _hasCtfLabels(emxParticle) or _hasCtfLabels(emxParticle.getMicrograph()):
                part.setCTF(CTFModel())
                partSet.setHasCTF(True)
            if emxParticle.has('centerCoord__X'):
                part.setCoordinate(Coordinate())
            _particleFromEmx(emxParticle, part)
            if emxParticle.has('transformationMatrix__t11'):
                #FIXME: Detect if the alignment is 2D or 3D
                partSet.setAlignment3D()
            if not samplingRate:
                samplingRate = part.getSamplingRate()
            partSet.setSamplingRate(samplingRate) 
            particles = True
        else: # if not binary data, the coordinate case
            if micSet is None:
                raise Exception('Could not import Coordinates from EMX, micrographs not imported.')
            partSet = protocol._createSetOfCoordinates(micSet)
            part = Coordinate()
            particles = False
            
        copiedFiles = {} # copied or linked
        
        for emxParticle in emxData.iterClasses(emxlib.PARTICLE):
            if particles:
                _particleFromEmx(emxParticle, part)
                i, fn = part.getLocation()
                partFn = join(partDir, fn)
                newFn = join(outputDir, basename(partFn))
                newLoc = (i, newFn)
                
                if not partFn in copiedFiles:
                    copyOrLink(partFn, newFn)
                    copiedFiles[partFn] = newFn
                    
                part.setLocation(newLoc)
                
                if partSet.hasAlignment():
                    transform = Transform()
                    _transformFromEmx(emxParticle, part, transform)
                    part.setTransform(transform)
            else:
                _coordinateFromEmx(emxParticle, part)
                
            partSet.append(part)
            part.cleanObjId()
            
        if particles:
            protocol._defineOutputs(outputParticles=partSet)
        else:
            protocol._defineOutputs(outputCoordinates=partSet)
        
        if micSet is not None:
            protocol._defineSourceRelation(micSet, partSet)
            
