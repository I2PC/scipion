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
from pyworkflow.object import String
from pyworkflow.protocol import Protocol
from pyworkflow.em import SetOfMicrographs

class ProtImportMicrographs(Protocol):
    def __init__(self, **args):
        Protocol.__init__(self, **args)
        self.pattern = String(args.get('pattern', None))
              
        
    def defineSteps(self):
        self.insertFunctionStep('copyMicrographs', self.pattern.get())
        
    def copyMicrographs(self, pattern):
        """Copy micrographs matching the filename pattern"""
        from glob import glob
        files = glob(pattern)
        path = self.getPath('micrographs.txt')
        micFile = open(path, 'w+')
        for f in files:
            dst = self.getPath(os.path.basename(f))
            print >> micFile, dst
            shutil.copyfile(f, dst)
        micFile.close()
        
        self.defineOutputs(micrograph=SetOfMicrographs(value=path),
                           micrographFiltered=SetOfMicrographs(value='KKK'))
        
        return path
        

class ProtScreenMicrographs(Protocol):
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
