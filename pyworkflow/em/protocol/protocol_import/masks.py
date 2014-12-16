# **************************************************************************
# *
# * Authors:     Laura del Caño (ldelcano@cnb.csic.es)
# *              Ignacio Foche (ifoche@cnb.csic.es)
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
In this module are protocol related to EM imports of Masks...
"""

from os.path import join

from pyworkflow.protocol.params import (PathParam, BooleanParam, 
                                        EnumParam, StringParam, LEVEL_ADVANCED)
from pyworkflow.utils.path import expandPattern, createLink, copyFile
from pyworkflow.em.protocol import EMProtocol



class ProtImportMask(EMProtocol):
    """ Class for import masks from existing files. """
    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        
        form.addSection(label='Import')

        form.addParam('maskPath', PathParam, 
                      label="Mask path",
                      help="Select the file path of the mask\n")
        form.addParam('samplingRate', FloatParam, default=1.,
                   label=Message.LABEL_SAMP_RATE)
        
    #--------------------------- INFO functions ----------------------------------------------------
    def _validate(self):
        return []
    
    #--------------------------- UTILS functions ---------------------------------------------------


