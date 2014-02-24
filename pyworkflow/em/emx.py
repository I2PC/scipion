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
import os
from datetime import datetime as dt
try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET

from pyworkflow.utils import *
from data import *
    

def writeHeader(f):
    """ Write the header of the EMX file. """
    f.write("""<?xml version='1.0' encoding='utf-8'?>
<EMX version="1.0">
 <!--
  ##########################################################################
  #               EMX Exchange file 
  #               Produced by Scipion (%s)
  # 
  #  Information about EMX file format is available at 
  #  http://i2pc.cnb.csic.es/emx
  ##########################################################################
  --> """ % dateStr(dt.now()))
    
    
def writeFooter(f):
    f.write("""
    
</EMX>""")
    
    
def _writeCTF(f, ctf):
    """ Write the CTF values in the EMX file. """
    ctfDict = {'defocusU': ctf.getDefocusU(),
               'defocusV': ctf.getDefocusV(),
               'defocusUAngle': ctf.getDefocusAngle()
               }
    f.write("""
  <defocusU unit="nm">%(defocusU)0.2f</defocusU>
  <defocusV unit="nm">%(defocusV)0.2f</defocusV>
  <defocusUAngle unit="deg">%(defocusUAngle)0.2f</defocusUAngle>""" % ctfDict)
    
def _writeSamplingRate(f, sampling):
    """ Write sampling rate info as expected by EMX. """
    f.write("""
  <pixelSpacing>
    <X unit="A/px">%(sampling)0.2f</X> <Y unit="A/px">%(sampling)0.2f</Y>
  </pixelSpacing>""" % locals())
    
def writeMicrograph(f, mic):
    """ Write the micrograph xml to the file """
    index, filename = mic.getLocation()
    acq = mic.getAcquisition()
    micDict = {'index': index or 1,
               'fileName': filename,
               'voltage': acq.getVoltage(),
               'amplitudeContrast': acq.getAmplitudeContrast(),
               'cs': acq.getSphericalAberration(),
               }
    f.write("""
<micrograph index="%(index)d" fileName="%(fileName)s">    
  <acceleratingVoltage unit="kV">%(voltage)0.2f</acceleratingVoltage>
  <amplitudeContrast>%(amplitudeContrast)0.2f</amplitudeContrast>
  <cs unit="mm">%(cs)0.2f</cs>""" % micDict)
    _writeSamplingRate(f, mic.getSamplingRate())
    if mic.hasCTF():
        _writeCTF(f, mic.getCTF())
    f.write("""
</micrograph> """)
    
    
def _writeCenter(f, coord):
    coordDict = {'X': coord.getX(), 'Y': coord.getY()}
    f.write("""
    <centerCoord>
      <X unit="px">%(X)d</X> <Y unit="px">%(Y)d</Y>
    </centerCoord> """ % coordDict)
    
    
def _writeCommon(f, coordinate, index, filename, micSet):
    """ This helper function will serve to write common xml string
    of Coordinate and Particle, actually, in EMX both are 'particle'
    """
    partDict = {'index': index or 1,
               'fileName': filename,
               }
    """ Write the partice xml to the file. """
    f.write("""
<particle fileName="%(fileName)s" index="%(index)d">""" % partDict)
    if coordinate:
        mic = micSet[coordinate.getMicId()] # get micrograph
        index, filename = mic.getLocation()
        partDict.update({'index': index or 1, 'fileName': basename(filename)})
        f.write("""
    <micrograph fileName="%(fileName)s" index="%(index)d"/>""" % partDict)
        _writeCenter(f, coordinate)
        
        
def _writeCommonFooter(f):
    f.write("""
</particle> """ )
           
           
def writeParticle(f, particle, micSet):
    index, filename = particle.getLocation()
    _writeCommon(f, particle.getCoordinate(), index, filename, micSet)
    if particle.hasCTF():
        _writeCTF(f, particle.getCTF())
    _writeSamplingRate(f, particle.getSamplingRate())
    _writeCommonFooter(f)


def writeCoordinate(f, coordinate, micSet):
    _writeCommon(f, coordinate, 1, "coordinate", micSet)  
    _writeCommonFooter(f)
    
    
def writeSetOfMicrographs(f, micSet, emxDir, ctfSet=None, writeData=True):
    """ Write a SetOfMicrograph as expected in EMX format (xml file)
    Params:
        micSet: input set of micrographs
        filename: the EMX file where to store the micrographs information.
    """
    from convert import ImageHandler
    from constants import NO_INDEX
    ih = ImageHandler()
    f.write(""" 
 <!--
  ##########################################################################
  #               Micrographs
  ##########################################################################
  --> """)
    for mic in micSet:
        loc = mic.getLocation()
        fnMicBase = replaceBaseExt(loc[1], 'mrc')
        if writeData:
            newLoc = join(emxDir, fnMicBase)
            ih.convert(loc, newLoc)
        mic.setLocation(NO_INDEX, fnMicBase)
        if ctfSet:
            mic.setCTF(ctfSet[mic.getObjId()])
        writeMicrograph(f, mic)
    
    
def writeSetOfParticles(f, partSet, stackFn=None, micSet=None):
    """ Write a SetOfMicrograph as expected in EMX format (xml file)
    Params:
        micSet: input set of micrographs
        filename: the EMX file where to store the micrographs information.
    """
    f.write("""
 <!--
  ##########################################################################
  #               Particles
  ##########################################################################
  --> """)
    from pyworkflow.em.convert import ImageHandler
    ih = ImageHandler()
    
    for i, particle in enumerate(partSet):
        if stackFn:
            loc = particle.getLocation()
            newLoc = (i+1, stackFn)
            ih.convert(loc, newLoc)
            newFn = basename(stackFn)
            particle.setLocation(i+1, newFn)
            writeParticle(f, particle, micSet)
        else:
            writeCoordinate(f, particle, micSet)
        
        
def exportData(emxDir, inputSet, ctfSet=None):
    """ Export micrographs, coordinates or particles to  EMX format. """
    cleanPath(emxDir)
    makePath(emxDir) 
    fnXml = join(emxDir, 'data.emx')
    f = open(fnXml, 'w+')
    writeHeader(f)

    if isinstance(inputSet, SetOfMicrographs):
        writeSetOfMicrographs(f, inputSet, emxDir, ctfSet)
        
    elif isinstance(inputSet, SetOfCoordinates):
        micSet = inputSet.getMicrographs()
        writeSetOfMicrographs(f, micSet, emxDir, ctfSet)
        writeSetOfParticles(f, inputSet, None, micSet)
        
    elif isinstance(inputSet, SetOfParticles):
        if inputSet.hasCoordinates():
            micSet = inputSet.getCoordinates().getMicrographs()
            writeSetOfMicrographs(f, micSet, emxDir, writeData=False)
        fnMrcs = join(emxDir, 'data.mrc')
        writeSetOfParticles(f, inputSet, fnMrcs, micSet)
        
    writeFooter(f)
    f.close()
    
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
    