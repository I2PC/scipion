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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

from pyworkflow.object import *
from pyworkflow.em import *
from tests import *
from pyworkflow.mapper import SqliteMapper
from pyworkflow.utils import dateStr
from pyworkflow.protocol.constants import MODE_RESUME, STATUS_FINISHED
from pyworkflow.protocol.executor import StepExecutor

    
#Protocol for tests, runs in resume mode, and sleeps for??
class MyProtocol(Protocol):
    def __init__(self, **args):
        Protocol.__init__(self, **args)
        self.name = String(args.get('name', None))
        self.numberOfSleeps = Integer(args.get('n', 1))
        self.runMode = Integer(MODE_RESUME)
        
    def sleepStep(self, t, s):
        log = self._getPath("step_%02d.txt" % t)
        import time, datetime
        f = open(log, 'w+')
        f.write("Going to sleep at %s\n" % dateStr(datetime.datetime.now(), True))
        time.sleep(t)
        f.write("  Slept: %d seconds\n" % t)
        f.write("Awaked at %s\n" % dateStr(datetime.datetime.now(), True))
        f.close()
        return [log]
        
    def _insertAllSteps(self):
        for i in range(self.numberOfSleeps.get()):
            self._insertFunctionStep('sleepStep', i+1, 'sleeping %d'%i)
            
            
class MyParallelProtocol(MyProtocol):
    def _insertAllSteps(self):
        step1 = self._insertFunctionStep('sleepStep', 1, '1')
        n = 2
        deps = [step1]
        for i in range(n):
            self._insertFunctionStep('sleepStep')
    
            
# TODO: this test seems not to be finished.
class TestProtocolExecution(BaseTest):
    
    @classmethod
    def setUpClass(cls):
        setupTestOutput(cls)
    
    def test_StepExecutor(self):
        """Test the list with several Complex"""
        fn = self.getOutputPath("protocol.sqlite")   
        mapper = SqliteMapper(fn, globals())
        prot = MyProtocol(mapper=mapper, n=2, workingDir=self.getOutputPath(''))
        prot._stepsExecutor = StepExecutor(hostConfig=None)
        prot.run()
        
        self.assertEqual(prot._steps[0].getStatus(), STATUS_FINISHED)
        
        mapper2 = SqliteMapper(fn, globals())
        prot2 = mapper2.selectById(prot.getObjId())
        
        self.assertEqual(prot.endTime.get(), prot2.endTime.get())
