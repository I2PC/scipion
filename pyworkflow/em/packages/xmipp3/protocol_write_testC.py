# **************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
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

import pyworkflow.protocol.params as params

from pyworkflow.em.protocol import ProtProcessParticles
from pyworkflow.protocol.params import IntParam

class XmippProtWriteTestC(ProtProcessParticles):
    """    
    using mpi write data to a large file (C++ level)
    """
    _label = 'write_test_C'
    
    def __init__(self, *args, **kwargs):
        ProtProcessParticles.__init__(self, *args, **kwargs)

    #--------------------------- DEFINE param functions --------------------------------------------   
    def _defineParams(self, form):
        form.addSection(label='Size is important')
        form.addParam('xDim', IntParam, default=128,label='X dimension')
        form.addParam('yDim', IntParam, default=128,label='Y dimension')
        form.addParam('nDim', IntParam, default=1024,label='Number of frames')
        form.addParallelSection(threads=1, mpi=4)


    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('testWriteStep')


    def testWriteStep(self):
        params =  ' -i kk.mrcs '
        params += ' --xdim %d '%self.xDim.get()
        params += ' --ydim %d '%self.yDim.get()
        params += ' --ndim %d '%self.nDim.get()

        #change to directory
        nproc = self.numberOfMpi.get()
        #nT=self.numberOfThreads.get()

        self.runJob('xmipp_write_test',
                    params, numberOfMpi=nproc,cwd=self._getExtraPath())
    
    #--------------------------- INFO functions --------------------------------------------
    def _validate(self):
        pass
    
    def _summary(self):
        pass
    
    def _methods(self):
        messages = []
        return messages
    
    def _citations(self):
        return ['Vargas2014a']
    
    #--------------------------- UTILS functions -------------------------------------------- 
    def _updateLocation(self, item, row):
        index, filename = xmippToLocation(row.getValue(md.MDL_IMAGE))
        item.setLocation(index, filename)

