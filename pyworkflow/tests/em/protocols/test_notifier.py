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
# *  e-mail address 'scipion@cnb.csic.es'
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
        return os.environ.get('SCIPION_NOTIFY_URL', 'http://calm-shelf-73264.herokuapp.com/report_protocols/api/workflow/workflow/').strip()

    def test_projectNotifier(self):
        """ Execute a protocol and then report on it
        """
        # Connect to heroku and retrieve protocol list
        # then store number of times protstress has been executed
        url = self._getUrl()

        #change periodicy in notification
        os.environ["SCIPION_NOTIFY_SECONDS"] = "30"
        os.environ["SCIPION_NOTIFY"] = "True"
        url = url.replace("workflow/workflow/",
                          "workflow/protocol/?name=ProtStress")
        results = json.loads(urllib2.urlopen(url).read())
        objects = results["objects"]
        if  (len(objects)!=0):
            objects = results["objects"][0]
            times_protocolRemote = objects["timesUsed"]
        else:
            times_protocolRemote = 0
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

        #we want to test this class: ProjectNotifier
        projectNotifier = ProjectNotifier(self.proj)
        # remove uuid and data files (so we start always in the same conditions)
        uuidFn=projectNotifier._getUuidFileName()
        os.remove(uuidFn) if os.path.exists(uuidFn) else None
        dataFn=projectNotifier._getDataFileName()
        os.remove(dataFn) if os.path.exists(dataFn) else None

        #store workflow  in database
        projectNotifier.notifyWorkflow()
        #wait until protocol gui has connected to heroku and sent the information
        time.sleep(5)
        #get uuid from local file
        uuid = projectNotifier._getUuid()
        #get protocol list from database
        urlWork = self._getUrl()
        #url = "http://secret-reaches-65198.herokuapp.com/report_protocols/api/workflow/workflow/?project_uuid="
        urlWork += "?project_uuid=" + uuid
        results = json.loads(urllib2.urlopen(urlWork).read())

        objects = results["objects"]
        if  (len(objects)!=0):
            objects = results["objects"][0]
            project_workflowRemote = objects["project_workflow"]
        else:
            project_workflowRemote = ""
        # get protocol list local
        project_workflowLocal  = self.proj.getProtocolsJson(namesOnly=True)
        #test that stored protocol and local one are identical
        self.assertEqual(project_workflowLocal, project_workflowRemote)

        #now check that the  number of times postgress has been executed
        #has increased in 1
        urlProt = self._getUrl()
        urlProt = urlProt.replace("workflow/workflow/",
                          "workflow/protocol/?name=ProtStress")
        time.sleep(5)# notifier runs in a thread so wait a bit
        results = json.loads(urllib2.urlopen(urlProt).read())
        objects = results["objects"]
        if  (len(objects)!=0):
            objects = results["objects"][0]
            times_protocolRemote_2 = objects["timesUsed"]
        else:
            times_protocolRemote_2 = 0
        self.assertEqual(times_protocolRemote_2, times_protocolRemote +1)

        #try to resend the project, as we have a 30 sec time this should
        #not go through number times should not change
        projectNotifier.notifyWorkflow()
        time.sleep(5)# notifier runs in a thread so wait a bit
        results = json.loads(urllib2.urlopen(urlProt).read())
        objects = results["objects"]
        if  (len(objects)!=0):
            objects = results["objects"][0]
            times_protocolRemote_2 = objects["timesUsed"]
        else:
            times_protocolRemote_2 = 0
        self.assertEqual(times_protocolRemote_2, times_protocolRemote +1)

        #try to resend project, as no modification has been made number of
        # times should not change
        time.sleep(25)
        projectNotifier.notifyWorkflow()
        time.sleep(5)# notifier runs in a thread so wait a bit
        results = json.loads(urllib2.urlopen(urlProt).read())
        objects = results["objects"]
        if  (len(objects)!=0):
            objects = results["objects"][0]
            times_protocolRemote_2 = objects["timesUsed"]
        else:
            times_protocolRemote_2 = 0
        self.assertEqual(times_protocolRemote_2, times_protocolRemote +1)

        #add new protocol and resend
        prot2 = self.newProtocol(ProtStress, **kwargs)
        prot2.setObjLabel('stress2')
        self.proj.launchProtocol(prot2, wait=True)
        projectNotifier.notifyWorkflow()
        time.sleep(5)# notifier runs in a thread so wait a bit
        results = json.loads(urllib2.urlopen(urlProt).read())
        objects = results["objects"]
        if  (len(objects)!=0):
            objects = results["objects"][0]
            times_protocolRemote_2 = objects["timesUsed"]
        else:
            times_protocolRemote_2 = 0
        self.assertEqual(times_protocolRemote_2, times_protocolRemote +2)
