# **************************************************************************
# *
# * Authors:    Roberto Marabini       (roberto@cnb.csic.es)
#               J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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

import threading
import uuid
import urllib, urllib2
import os
import time
from datetime import timedelta, datetime

import pyworkflow.utils as pwutils


class ProjectNotifier(object):
    """ Implement different types of notifications about a given
    project. Currently, the protocols in a workflow will be sent.
    """
    def __init__(self, project):
        self.project = project

    def _getUuidFileName(self):
        return self.project.getLogPath("uuid.log")

    def _getUuid(self):
        # Load (or create if not exits) a file
        # in the project Logs folder to store an unique
        # project identifier
        uuidFn = self._getUuidFileName()
        try:
            with open(uuidFn) as f:
                uuidValue = f.readline()
        except IOError:
            uuidValue = str(uuid.uuid4())
            with open(uuidFn,'w') as f:
                 f.write(uuidValue)

        return uuidValue

    def _modifiedBefore(self, seconds):
        """ Return True if the uuid.log file has been modified within a given
        number of seconds. """
        uuidFn = self._getUuidFileName()
        if not os.path.exists(uuidFn):
            return False
        mTime = datetime.fromtimestamp(os.path.getmtime(uuidFn))
        delta = datetime.now() - mTime

        return delta < timedelta(seconds=seconds)

    def _sendData(self, url, dataDict=None):
        #then connect to webserver a send json
        opener = urllib2.build_opener(urllib2.HTTPHandler(debuglevel=0))#no messages
        data = urllib.urlencode(dataDict)
        content = opener.open(url, data=data).read()
        now = time.time()
        #print "Notifying...."
        #pwutils.prettyDate(now)
        #print "dataDict: ", dataDict
        os.utime(self._getUuidFileName(), (now, now))

    def notifyWorkflow(self):
        #check if enviroment exists otherwise abort
        if not pwutils.envVarOn('SCIPION_NOTIFY'):
            return

        # Check the seconds range of the notify, by default one day
        seconds = int(os.environ.get('SCIPION_NOTIFY_SECONDS', '86400'))

        if self._modifiedBefore(seconds): # notify not more than once a day
            #print "sec, no notification", seconds
            return

        # INFO: now we are only sending the protocols names in the project.
        # We could pass namesOnly=False to get the full workflow template
        projectWorfklow = self.project.getProtocolsJson(namesOnly=True)

        dataDict = {'project_uuid': self._getUuid(),
                    'project_workflow': projectWorfklow}

        urlName = os.environ.get('SCIPION_NOTIFY_URL', '').strip()
        urlName += "addOrUpdateWorkflow/"
        t = threading.Thread(target=lambda: self._sendData(urlName, dataDict))
        t.start() # will execute function in a separate thread

    def getEntryFromWebservice(self,uuid):
        if not pwutils.envVarOn('SCIPION_NOTIFY'):
            return
        urlName = os.environ.get('SCIPION_NOTIFY_URL').strip()
        #remove last directory
        urlName = os.path.split(urlName)[0]
        url = urlName + "/?project_uuid=" + uuid
        resultDict = self._sendData(url)
