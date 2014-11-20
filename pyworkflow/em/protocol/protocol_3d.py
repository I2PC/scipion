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


class Prot3D(EMProtocol):
    pass

class ProtPreprocessVolumes(Prot3D):
    """ This class will serve as a base for all protocol
    that performs some operation on Volumes (i.e. filters, mask, resize, etc)
    It is mainly defined by an inputVolumes and outputVolumes.
    """
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        
        form.addParam('inputVolumes', PointerParam, important=True,
                      label=Message.LABEL_INPUT_VOLS, pointerClass='Volume, SetOfVolumes',
                      help='Can be a density volume or a SetOfVolumes.')
        # Hook that should be implemented in subclasses
        self._defineProcessParams(form)
        form.addParallelSection(threads=2, mpi=1)
        
    def _defineProcessParams(self, form):
        """ This method should be implemented by subclasses
        to add other parameter relatives to the specific operation."""
        pass


class ProtFilterVolumes(ProtPreprocessVolumes):
    """ This is the base for the branch of filters, 
    between the ProtPreprocessVolumes """
    pass

class ProtOperateVolumes(ProtPreprocessVolumes):
    """ This is the base for the branch of filters,
    between the ProtPreprocessParticles """
    pass

class ProtMaskVolumes(ProtPreprocessVolumes):
    """ This is the base for the branch of mask, 
    between the ProtPreprocessVolumes """
    pass


class ProtCreateMask3D(ProtPreprocessVolumes):
    pass


class ProtAlignVolume(ProtPreprocessVolumes):
    """Protocol base for Align volumes protocols"""
    pass


class ProtReconstruct3D(Prot3D):
    pass

class ProtRefine3D(Prot3D):
    pass


class ProtClassify3D(Prot3D):
    pass


class ProtInitialVolume(Prot3D):
    """Protocol base for Initial volumes protocols"""
    pass


class ProtAnalysis3D(Prot3D):
    pass