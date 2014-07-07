# **************************************************************************
# *
# * Authors:     Airen Zaldivar Peraza (azaldivar@cnb.csic.es)
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
In this module are protocol base classes related to 2D processing

"""
from pyworkflow.em.protocol import *



class Prot2D(EMProtocol):
    pass


class ProtAlign2D(Prot2D):
    """ This class will serve as a base for all protocols that align a set of 2D images.
    All Align protocols receive as input:
        A set of partices
    and will allow the option to generate the aligned particles.
    """
    def _defineParams(self, form):
        form.addSection(label='Input')
        
        form.addParam('inputParticles', PointerParam, important=True,
                      label=Message.LABEL_INPUT_PART, pointerClass='SetOfParticles')
        form.addParam('writeAlignedParticles', BooleanParam, default=True,
                      label=Message.LABEL_ALIG_PART, 
                      help=Message.TEXT_ALIG_PART)
        # Hook that should be implemented in subclasses
        self._defineAlignParams(form)
        
    def _defineAlignParams(self, form):
        """ This method should be implemented by subclasses
        to add other parameter relatives to the specific align protocol."""
        pass 
    

class ProtCreateMask2D(Prot2D):
    """ For those protocols who create mask as output. """
    pass


class ProtClassify2D(Prot2D):
    pass

class ProtAnalysis2D(Prot2D):
    pass