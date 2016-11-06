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

from protocol import Protocol
from pyworkflow.protocol.params import IntParam, STEPS_PARALLEL

        
                
class ProtTestParallel(Protocol):
    """ A parallel test protocol.
    """    
    _label = "parallel test"
    
    def __init__(self, **args):
        Protocol.__init__(self, **args)        
        self.stepsExecutionMode = STEPS_PARALLEL
        
    #--------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('numberOfIterations', IntParam, default=2, 
                      label="Number of iterations", 
                      help='Repeat the insertion of steps N times.')
        form.addParam('numberOfParallelSleeps', IntParam, default=2, 
                      label="Number of parallel sleeps", 
                      help='How many sleep steps can be done at the same time.')
        form.addParam('failAfter', IntParam, default=0,
                      label="Fail after", 
                      help='If you set an id, the next step should fail')
        form.addParam('sleepSecs', IntParam, default=2,
                      label='Seconds to sleep')
        
        form.addParallelSection(threads=4, mpi=1)
            
         
    #--------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        n = self.numberOfIterations.get()
        m = self.numberOfParallelSleeps.get()
        secs = self.sleepSecs.get()
        failAfter = self.failAfter.get()
        
        deps = []
        
        for i in range(n):
            initId = self._insertFunctionStep('initStep', i, prerequisites=deps)
            deps = []
            for j in range(m):
                fail = deps[-1] == failAfter if deps else False
                tag = 'iter %d, n: %d' % (i, j)
                sleepId = self._insertFunctionStep('sleepStep', secs, fail, tag,
                                                   prerequisites=[initId])
                deps.append(sleepId)
            endId = self._insertFunctionStep('endStep', i, prerequisites=deps)
            deps = [endId]
            
            
    #--------------------------- STEPS functions -------------------------------

    def initStep(self, iterN):
        """ All subsequent sleep steps should depend on this. """
        self._log.info("Starting iteration: %d" % iterN)
    
    def sleepStep(self, secs=5, forceFail=False, tag=''):
        if forceFail:
            self.runJob('echo', " 'Failing for testing purposes...'; exit 1")
        else:
            from pyworkflow.em.packages.xmipp3 import getEnviron
            self.runJob('xmipp_work_test',
                        "--time %d --tag '%s'" % (secs, tag), env=getEnviron())
            
    def awakeStep(self):
        print "Awaked after an sleep step"
        
    def endStep(self, iterN):
        self._log.info("Ending iteration: %d" % iterN)
        
