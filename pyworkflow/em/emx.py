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
try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET

from pyworkflow.utils import *
from data import SetOfMicrographs, SetOfCoordinates, SetOfParticles
    

def writeHeader(f):
    """ Write the header of the EMX file. """
    f.write("""
<?xml version='1.0' encoding='utf-8'?>
<EMX version="1.0">
 <!--
  ##########################################################################
  #               EMX Exchange file 
  #               Produced by Scipion
  # 
  #  This is a EMX file.
  #
  #  Information on this file format is available at 
  #  http://i2pc.cnb.csic.es/emx
  ##########################################################################
  -->
     """)
    
    
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
    
    
def writeMicrograph(f, mic):
    """ Write the micrograph xml to the file """
    index, filename = mic.getLocation()
    acq = mic.getAcquisition()
    micDict = {'index': index or 1,
               'fileName': filename,
               'voltage': acq.getVoltage(),
               'amplitudeContrast': acq.getAmplitudeContrast(),
               'cs': acq.getSphericalAberration(),
               'samplingRate': mic.getSamplingRate(),
               }
    f.write("""
<micrograph index="%(index)d" fileName="%(fileName)s">    
  <acceleratingVoltage unit="kV">%(voltage)0.2f</acceleratingVoltage>
  <amplitudeContrast>%(amplitudeContrast)0.2f</amplitudeContrast>
  <cs unit="mm">%(cs)0.2f</cs>
  <pixelSpacing>
    <X unit="A/px">%(samplingRate)0.2f</X> <Y unit="A/px">%(samplingRate)0.2f</Y>
  </pixelSpacing>""" % micDict)
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
    <micrograph fileName="%(fileName)s" index="%(index)d">""" % partDict)
        _writeCenter(f, coordinate)
        
        
def _writeCommonFooter(f):
    f.write("""
</particle> """ )
           
           
def writeParticle(f, particle, micSet):
    index, filename = particle.getLocation()
    _writeCommon(f, particle.getCoordinate(), index, filename, micSet)
    if particle.hasCTF():
        particle.getCTF().printAll()
        _writeCTF(f, particle.getCTF())
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
    
    for mic in micSet:
        loc = mic.getLocation()
        fnMicBase = basename(loc[1])
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
        fnMrcs = join(emxDir, 'data.mrcs')
        writeSetOfParticles(f, inputSet, fnMrcs, micSet)
        
    writeFooter(f)
    f.close()
    