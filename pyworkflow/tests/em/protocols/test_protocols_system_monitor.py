# ***************************************************************************
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *
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
# ***************************************************************************/
from pyworkflow.em.protocol.monitors.protocol_monitor_system import SYSTEM_LOG_SQLITE
from pyworkflow.em.protocol.protocol_stress import STRESS_NG
from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pyworkflow.em.protocol import ProtStress, ProtMonitorSystem
from pyworkflow.utils import commandExists


class TestStress(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

    def test_pattern(self):
        """ Import several Particles from a given pattern.
        """

        # Check stress-ng is installed
        print "Command exists?&&&&&&&&&&&&&&&&&&& %s" % commandExists(STRESS_NG)
        self.assertTrue(commandExists(STRESS_NG), "This test requires to have %s installed" % STRESS_NG)

        kwargs = {'noCpu': 2,
                'noMem': 0,
                'amountMem': 256,
                'timeout':10,
                'delay':1
                }
        
        #put some stress on the system
        prot1 = self.newProtocol(ProtStress, **kwargs)
        prot1.setObjLabel('stress')
        self.proj.launchProtocol(prot1,wait=False)

        #TODO fill protol pointer
        kwargs = {'samplingInterval': 1,
                'interval': 20
                }
        prot2 = self.newProtocol(ProtMonitorSystem, **kwargs)
        prot2.inputProtocols.append(prot1)
        self.launchProtocol(prot2)

        baseFn = prot2._getPath(SYSTEM_LOG_SQLITE)
        import os.path
        #not sure what to test here
        print baseFn
        self.assertTrue(os.path.isfile(baseFn))
