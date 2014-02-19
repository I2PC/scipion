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
    

def writeCTF(ctf, f):
    """ Write the CTF values in the EMX file. """
    ctfDict = {'defocusU': ctf.getDefocusU(),
               'defocusV': ctf.getDefocusV(),
               'defocusUAngle': ctf.getDefocusAngle()
               }
    f.write("""
  <defocusU unit="nm">%(defocusU)0.2f</defocusU>
  <defocusV unit="nm">%(defocusV)0.2f</defocusV>
  <defocusUAngle unit="deg">%(defocusUAngle)0.2f</defocusUAngle>""" % ctfDict)
    
    
def writeMicrograph(mic, f):
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
    <X unit="A/px">%(samplingRate)0.2f</X>
    <Y unit="A/px">%(samplingRate)0.2f</Y>
  </pixelSpacing>""" % micDict)
    if mic.hasCTF():
        writeCTF(mic.getCTF(), f)
    f.write("""
</micrograph>
    """)
    
    
def writeParticle(particle, f):
    """ Write the partice xml to the file. """
    particleXmlStr = """
<particle fileName="particles.mrc" index="3">
    <micrograph fileName="image002.mrc"/>
    <centerCoord>
      <X unit="px">522</X>
      <Y unit="px">194</Y>
   </centerCoord>
</particle> """
    
    
def writeSetOfMicrographs(micSet, f):
    """ Write a SetOfMicrograph as expected in EMX format (xml file)
    Params:
        micSet: input set of micrographs
        filename: the EMX file where to store the micrographs information.
    """
    for mic in micSet:
        writeMicrograph(mic, f)
    

def exportSetOfMicrographs(emxDir, micSet=None, ctfSet=None, partSet=None):
    """ Export micrographs as EMX format. """
    fnXml = os.path.join(emxDir, 'data.emx')
    f = open(fnXml, 'w+')
    
    from pyworkflow.em.convert import ImageHandler
    ih = ImageHandler()
    
    for i, mic in enumerate(micSet):
        loc = mic.getLocation()
        fnMicBase = basename(loc[1])
        fnMicNew = join(emxDir, fnMicBase)
        newLoc = (i+1, fnMicNew)
        ih.convert(loc, newLoc)
        mic.setLocation(None, fnMicBase)
        if ctfSet:
            mic.setCTF(ctfSet[mic.getObjId()])
        writeMicrograph(mic, f)
        
    f.close()
    