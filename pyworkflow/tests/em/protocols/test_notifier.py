# ***************************************************************************
# * Authors:    Roberto Marabini (roberto@cnb.csic.es)
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
import  time
import urllib2
import json

from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pyworkflow.em.protocol import ProtStress, ProtMonitorSystem
from pyworkflow.gui.project.notifier import ProjectNotifier

class TestNotifier(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

    def test_projectNotifier(self):
        """ Execute a protocol and then report on it
        """
        #run any protocol
        kwargs = {'noCpu': 2,
                  'noMem': 0,
                  'amountMem': 8,
                  'timeout':3,
                  'delay':1
                  }

        prot1 = self.newProtocol(ProtStress, **kwargs)
        prot1.setObjLabel('stress')
        self.proj.launchProtocol(prot1,wait=True)

        # remove uuid file
        uuidFn=self.proj.getLogPath("uuid.log")
        os.remove(uuidFn) if os.path.exists(uuidFn) else None

        #we want to test this class: ProjectNotifier
        projectNotifier = ProjectNotifier(self.proj)
        #stores uuid enry in database
        projectNotifier.notifyWorkflow()
        #wait 15 second so heroku can wake up
        time.sleep(30)
        #get uuid from local file
        uuid = projectNotifier._getUuid()
        #get protocol list from database
        url = os.environ.get('SCIPION_NOTIFY_URL', '').strip()
        #url = "http://secret-reaches-65198.herokuapp.com/report_protocols/api/workflow/workflow/?project_uuid="
        url += "?project_uuid=" + uuid
        results = json.loads(urllib2.urlopen(url).read())
        print "url",url
        print "res", results
        #print "resObj", results["objects"]
        project_workflowRemote = results["objects"][0]['project_workflow']
        #get protocol list local
        project_workflowLocal  = self.proj.getProtocolsNameJson()
        #print "project_workflowRemote", project_workflowRemote
        #print "project_workflowLocal", project_workflowLocal

        #not sure what to test here
        self.assertEqual(project_workflowLocal, project_workflowRemote)
