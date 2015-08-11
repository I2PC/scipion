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
# *  e-mail address 'xmipp@cnb.csic.es'
# ***************************************************************************/


import os
from base import ProgramTest

import pyworkflow.em.packages.xmipp3 as xmipp3



class XmippTestProgAngularDiscreteAssign(ProgramTest):
    @classmethod
    def setUpClass(cls):
        #cls.setTestDir(os.path.join(os.environ['SCIPION_TESTS', 'testXmipp']))
        cls.program = 'xmipp_angular_discrete_assign'
        cls.env = xmipp3.getEnviron()
        
        
        """
                <CASE
            arguments="-i input/aFewProjections.sel -o %o/assigned_angles.txt --ref %o/reference.doc"
            changeDir="FALSE">
            <PRERUN
                command="xmipp_angular_project_library -i input/phantomBacteriorhodopsin.vol -o %o/reference.stk --sampling_rate 10" />
                        <FILE filename="assigned_angles.txt"/>
        </CASE>
        """
        
    def test_case1(self):
        """ Check that an ouput was generated and
        the condition is valid. 
        """
        self.runCase("-i input/aFewProjections.sel -o %o/assigned_angles.txt --ref %o/reference.doc",
                     preruns=["xmipp_angular_project_library -i input/phantomBacteriorhodopsin.vol -o %o/reference.stk --sampling_rate 10"],
                     outputs=['assigned_angles.txt'])
        
