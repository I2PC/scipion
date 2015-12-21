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
from xmipp import *
import os
from pyworkflow.em.protocol import ProtProcessParticles
from pyworkflow.protocol.constants import STEPS_PARALLEL
from pyworkflow.protocol.params import IntParam

class XmippProtWriteTestP(ProtProcessParticles):
    """    
    using mpi write data to a large file (python level)
    """
    _label = 'write_test_P'
    
    def __init__(self, *args, **kwargs):
        ProtProcessParticles.__init__(self, *args, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL

    #--------------------------- DEFINE param functions --------------------------------------------   
    def _defineParams(self, form):
        form.addSection(label='Size is important')
        form.addParam('xDim', IntParam, default=128,label='X dimension')
        form.addParam('yDim', IntParam, default=128,label='Y dimension')
        form.addParam('nframes', IntParam, default=1024,label='Number of frames')
        form.addParallelSection(threads=1, mpi=4)

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
       stepId = self._insertFunctionStep('testWriteStep', self.nframes.get())
       pre=[]
       for counter in range(self.nframes.get(),0,-1):
            id = self._insertFunctionStep('testWriteStep', counter, prerequisites=[stepId] )
            pre.append(id)
       stepId = self._insertFunctionStep('testWriteValidate',  prerequisites=pre)
       #check step...

    def testWriteStep(self, counter):
        #change directory
        #cwd=self._getExtraPath()
        #os.chdir(cwd)
	img = Image()
	img.setDataType(DT_FLOAT)
	img.resize(self.xDim.get(), self.yDim.get())
	img.initConstant(counter)
	img.write("%d@%s"%(counter,os.path.join(self._getExtraPath(),"kk.mrcs")))

    def testWriteValidate(self):
	img = Image()
        for counter in range(self.nframes.get(),0,-1):
	    img.read("%d@%s"%(counter,os.path.join(self._getExtraPath(),"kk.mrcs")))
            mean, dev, min, max = img.computeStats()
            if (min == max) and (min == counter) and (dev == 0):
                print("frame %d ok"%counter)
            else:
                print("frame %d no OK mean dev, min max"%counter, mean,dev,min,max)
                break
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

