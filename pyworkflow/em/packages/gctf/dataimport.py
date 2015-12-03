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

import os

import pyworkflow.utils as pwutils 
from pyworkflow.em.data import CTFModel
from convert import readCtfModel


class ZhangLabImportCTF():
    """ Import CTF estimated with GCTF. """
    def __init__(self, protocol):
        self.protocol = protocol
        self.copyOrLink = self.protocol.getCopyOrLink()

    def importCTF(self, mic, fileName):
        ctf = CTFModel()
        ctf.setMicrograph(mic)
        readCtfModel(ctf, fileName, ctf4=False)
        
        # Try to find the given PSD file associated with the cttfind log file
        # we handle special cases of .ctf extension and _ctffindX prefix for Relion runs
        fnBase = pwutils.removeExt(fileName)
        for suffix in ['_psd.mrc', '.mrc', '.ctf']:
            psdPrefixes = [fnBase, 
                           fnBase.replace('_ctffind3', '')]
            for prefix in psdPrefixes:
                psdFile =  prefix + suffix
                if os.path.exists(psdFile):
                    if psdFile.endswith('.ctf'):
                        psdFile += ':mrc'
                    ctf.setPsdFile(psdFile)
        return ctf

