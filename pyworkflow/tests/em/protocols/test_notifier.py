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

import pyworkflow.utils as pwutils
from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pyworkflow.em.protocol import ProtStress, ProtMonitorSystem
from pyworkflow.gui.project.notifier import ProjectNotifier

class TestNotifier(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

    def _getUrl(self):
        if not pwutils.envVarOn('SCIPION_NOTIFY'):
            return ''

        return os.environ.get('SCIPION_NOTIFY_URL', '').strip()

    def test_projectNotifier(self):
        """ Execute a protocol and then report on it
        """
        # Connect to heroku and retrieve protocol list
        # then store number of times protstress has been executed
        url = self._getUrl()

        if not url:
            print "SCIPION_NOTIFY and SCIPION_NOTIFY_URL should be defined!"
            print "A test server is at Heroku, you can setup as:"
            print "export SCIPION_NOTIFY=1"
            print "export SCIPION_NOTIFY_URL=http://calm-shelf-73264.herokuapp.com/report_protocols/api/workflow/workflow/"
            return

        url = url.replace("workflow/workflow/",
                          "workflow/protocol/?name=ProtStress")

        results = json.loads(urllib2.urlopen(url).read())
        times_protocolRemote = results["objects"][0]["timesUsed"]

        #run a project that executes one time protStress
        kwargs = {'noCpu': 2,
                  'noMem': 0,
                  'amountMem': 8,
                  'timeout':3,
                  'delay':1
                  }
        #create and execute protocol stress
        prot1 = self.newProtocol(ProtStress, **kwargs)
        prot1.setObjLabel('stress')
        self.proj.launchProtocol(prot1, wait=True)

        # remove uuid file (so we start always in the same conditions)
        uuidFn=self.proj.getLogPath("uuid.log")
        os.remove(uuidFn) if os.path.exists(uuidFn) else None

        #we want to test this class: ProjectNotifier
        projectNotifier = ProjectNotifier(self.proj)
        #store workflow  in database
        projectNotifier.notifyWorkflow()
        #wait until protocol gui has connected to heroku and sent the information
        time.sleep(5)
        #get uuid from local file
        uuid = projectNotifier._getUuid()
        #get protocol list from database
        url = self._getUrl()
        #url = "http://secret-reaches-65198.herokuapp.com/report_protocols/api/workflow/workflow/?project_uuid="
        url += "?project_uuid=" + uuid
        results = json.loads(urllib2.urlopen(url).read())
        project_workflowRemote = results["objects"][0]['project_workflow']
        # get protocol list local
        project_workflowLocal  = self.proj.getProtocolsJson(namesOnly=True)
        #test that stored protocol and local one are identical
        self.assertEqual(project_workflowLocal, project_workflowRemote)

        #now check that the  number of times protstress has been executed
        #has increased in 1
        url = self._getUrl()
        url = url.replace("workflow/workflow/",
                          "workflow/protocol/?name=ProtStress")
        results = json.loads(urllib2.urlopen(url).read())
        times_protocolRemote_2 = results["objects"][0]["timesUsed"]
        self.assertEqual(times_protocolRemote_2, times_protocolRemote +1)
