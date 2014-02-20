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
This module contains the protocol base class for Relion.
"""

from pyworkflow.em import *  
from pyworkflow.utils.which import which
from pyworkflow.utils.path import makePath, replaceBaseExt, join, basename
from convert import createRelionInputImages, createRelionInputVolume


class ProtRelionBase(EMProtocol):
    """ Base class for all relion protocols"""
    def __init__(self, **args):
        EMProtocol.__init__(self, **args)
        self.ParamsStr = ''
        self.program = 'relion_refine'
        if self.numberOfMpi > 1:
            self.program += '_mpi'

#TODO: move this from here
#    def loadEnvironment(self):
#        """ Load the environment variables needed for use RELION tools. """
#        RELION_DIR = os.environ['RELION_HOME']
#        os.environ['PATH'] = os.environ['PATH'] + os.pathsep + join(RELION_DIR, 'bin', '')
#        os.environ['LD_LIBRARY_PATH'] = os.environ['LD_LIBRARY_PATH'] + os.pathsep + join(RELION_DIR, 'lib64', '')
#        #os.environ['RELION_QSUB_TEMPLATE'] = join(RELION_DIR, 'bin', 'qsub.csh')

            
    def defineSteps(self): 
        # Convert input reference and images to Spider/MRC format as needed by Relion
        self._insertFunctionStep('convertInputStep')
        # Create a STAR files with the input images
        self._insertFunctionStep('createStarStep')
        
        if self.doContinue:
            self._insertContinueSteps()
        else:
            self._insertNormalSteps()
        self._insertFunctionStep('createOutputStep')
        
        self.loadEnvironment()
        print "doCTF: ", self.doCtf.get()
        self.imgStar    = createRelionInputImages(self, self.inputParticles.get(), doCTF=self.doCtf.get())        
        self.volumeName = createRelionInputVolume(self, self.input3DReferences.get())

    
    def convertInputStep(self):
        pass
    
    def createStarStep(self):
        pass


                        
