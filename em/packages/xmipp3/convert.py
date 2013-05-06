# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *              Laura del Cano (ldelcano@cnb.csic.es)
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
This module contains converter functions 
from base classes to xmipp classes
"""

from data import *

    
    
def convertCTFModel(ctfModel, filename):
    """ Method to convert from a general ctfModel to XmippCtfModel """
    if type(ctfModel) is XmippCTFModel:
        return ctfModel
        
    xmippCTFModel = XmippCTFModel()
    xmippCTFModel.copyInfo(ctfModel)
    
    xmippCTFModel.write(filename)
    
    return xmippCTFModel


    
    
