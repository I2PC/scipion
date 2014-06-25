# **************************************************************************
# *
# * Authors:     Airen
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
# *  e-mail address 'azaldivar@cnb.csic.es'
# *
# **************************************************************************
"""
This sub-package will contains Eman3.0 specific protocols
"""

import os

from pyworkflow.utils import environAdd



def loadEnvironment():
    BSOFT_HOME = os.environ['BSOFT_HOME']
    BSOFT_BIN = os.path.join(BSOFT_HOME, 'bin')
    environAdd('PATH', BSOFT_BIN)
    os.environ['BSOFT'] = BSOFT_HOME
    
    
    
