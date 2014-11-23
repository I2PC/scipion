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
In this module are protocol base classes related to EM imports of Micrographs, Particles, Volumes...
"""

from pyworkflow.utils.properties import Message
from pyworkflow.protocol.params import FloatParam, FileParam

from images import ProtImportImages



class ProtImportParticles(ProtImportImages):
    """Protocol to import a set of particles to the project"""
    _label = 'import particles'
    _outputClassName = 'SetOfParticles'

    IMPORT_FROM_EMX = 1
    IMPORT_FROM_XMIPP3 = 2
    IMPORT_FROM_RELION = 3

    def _getImportChoices(self):
        """ Return a list of possible choices
        from which the import can be done.
        (usually packages formas such as: xmipp3, eman2, relion...etc.
        """
        return ['emx', 'xmipp3', 'relion']
    
    def _defineBasicParams(self, form):
        ProtImportImages._defineBasicParams(self, form)
        
        form.addParam('micrographsEMX', FileParam,
              condition = '(importFrom == %d)' % self.IMPORT_FROM_EMX,
              label='Input EMX file',
              help="Select the EMX file containing particles information.\n"
                   "See more about [[http://i2pc.cnb.csic.es/emx][EMX format]]")
        
        form.addParam('particlesMd', FileParam,
                      condition = '(importFrom == %d)' % self.IMPORT_FROM_XMIPP3,
                      label='Particles metadata file',
                      help="Select the particles Xmipp metadata file.\n"
                           "It is usually a images.xmd_ file result\n"
                           "from Xmipp protocols execution.")
        
        form.addParam('dataStar', FileParam,
                      condition = '(importFrom == %d)' % self.IMPORT_FROM_RELION,
                      label='Particles metadata file',
                      help="Select the particles Xmipp metadata file.\n"
                           "It is usually a images.xmd_ file result\n"
                           "from Xmipp protocols execution.")
         
    def _defineParams(self, form):
        ProtImportImages._defineParams(self, form)
        group = self.acquisitionGroup
        group.addParam('samplingRate', FloatParam,
                   label=Message.LABEL_SAMP_RATE)
        
    def setSamplingRate(self, imgSet):
        imgSet.setSamplingRate(self.samplingRate.get())
    

class ProtImportAverages(ProtImportParticles):
    """Protocol to import a set of averages to the project"""
    _label = 'import averages'
    _outputClassName = 'SetOfAverages'    
    

