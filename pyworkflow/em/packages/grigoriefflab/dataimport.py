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
from pyworkflow.utils.path import removeExt
from pyworkflow.em.data import CTFModel
from convert import readCtfModel, ctffindOutputVersion



class BrandeisImport():
    """ Class used to import different kind of objects
    from Brandeis projects into Scipion.
    """
    def __init__(self, protocol):
        self.protocol = protocol
        self.copyOrLink = self.protocol.getCopyOrLink()

    def importCTF(self, mic, fileName):
        ctf = CTFModel()
        ctf.setMicrograph(mic)
        readCtfModel(ctf, fileName, ctf4=ctffindOutputVersion(fileName)==4)
        
        for suffix in ['_psd.mrc', '.mrc']:
            if os.path.exists(removeExt(fileName) + suffix):
                ctf.setPsdFile(removeExt(fileName) + suffix)
        return ctf


    

                