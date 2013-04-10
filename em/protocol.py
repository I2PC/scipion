# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
Protocols related to EM
"""
import os
import shutil
from pyworkflow.object import String, Float
from pyworkflow.protocol import Protocol
from pyworkflow.protocol.params import *
from pyworkflow.em import SetOfMicrographs


def defineImportMicrograph():
    """Create the definition of parameters for
    the ImportMicrographs protocol"""
    f = Form()
    
    f.addSection(label='Input')
    f.addParam('pattern', StringParam, label="Pattern")
    f.addParam('tiltPairs', BooleanParam, default=False, important=True,
               label='Are micrographs tilt pairs?')
    
    f.addSection(label='Microscope description')
    f.addParam('voltage', FloatParam, default=200,
               label='Microscope voltage (in kV)')
    f.addParam('sphericalAberration', FloatParam, default=2.26,
               label='Spherical aberration (in mm)')
    f.addParam('samplingRateMode', EnumParam, default=0,
               label='Sampling rate mode',
               choices=['From image', 'From scanner'])
    f.addParam('samplingRate', FloatParam,
               label='Sampling rate (A/px)', condition='samplingRateMode==0')
    f.addParam('magnification', IntParam, default=60000,
               label='Magnification rate', condition='samplingRateMode==1')
    f.addParam('ScannedPixelSize', FloatParam, 
               label='Scanned pixel size', condition='samplingRateMode==1')
    
    return f


class ProtImportMicrographs(Protocol):
    """Protocol to import a set of micrographs in the project"""
    _definition = defineImportMicrograph()
    
    def __init__(self, **args):
        Protocol.__init__(self, **args)
        
#        self.pattern = String(args.get('pattern', None))
#        self.voltage = Float(args.get('voltage', 300))
#        self.sphericalAberration = Float(args.get('sphericalAberration', 1.2))
#        self.samplingRate = Float(args.get('samplingRate', 1.2))
#        self.magnification = Float(args.get('magnification', 50000))             
        
    def defineSteps(self):
        self.insertFunctionStep('importMicrographs', self.pattern.get(),
                                self.voltage.get(), self.sphericalAberration.get(), 
                                self.samplingRate.get())
        
    def importMicrographs(self, pattern, voltage, sphericalAberration, samplingRate):
        """Copy micrographs matching the filename pattern
        Register other parameters"""
        from glob import glob
        files = glob(pattern)
        if len(files) == 0:
            raise Exception('importMicrographs:There is not files matching pattern')
        path = self.getPath('micrographs.txt')
        micFile = open(path, 'w+')
        for f in files:
            dst = self.getPath(os.path.basename(f))
            print >> micFile, dst
            shutil.copyfile(f, dst)
        micFile.close()
        
        micSet = SetOfMicrographs(value=path)
        micSet.microscope.voltage.set(voltage)
        micSet.microscope.sphericalAberration.set(sphericalAberration)
        micSet.samplingRate.set(samplingRate)
        self.defineOutputs(micrograph=micSet)
        
        return path
        

class ProtCTFMicrographs(Protocol):
    pass


class ProtDownsampleMicrographs(Protocol):
    pass


class ProtParticlePicking(Protocol):
    pass


class ProtAlign(Protocol):
    pass


class ProtClassify(Protocol):
    pass


class ProtAlignClassify(Protocol):
    pass
