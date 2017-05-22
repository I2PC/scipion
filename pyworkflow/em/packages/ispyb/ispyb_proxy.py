# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es) [1]
# *              Kevin Savage (kevin.savage@diamond.ac.uk) [2]
# *
# * [1] Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# * [2] Diamond Light Source, Ltd
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

import os
import sys
from collections import OrderedDict
from datetime import datetime


SCRIPT = os.path.realpath(__file__)
INPUT_END = 'END'
OUTPUT_OK = 'OK'


def timeStamp(dt=None, format='%Y-%m-%d_%H%M%S'):
    if dt is None:
        dt = datetime.now()
    return dt.strftime(format)

class ISPyBdb():
    """ This is a Facade to provide access to the ispyb_api to store movies.
    In case ispyb_api could not be loaded, we used a simple file log
    to still test things are working. In the normal case, it is supposed
    to be used from ISPyBProcess, which would be executed using the proper
    Python with the required modules for ispyb_api to work.
    """
    def __init__(self, db, experimentParams):
        self.f = open("/tmp/ispybdb.txt", 'a')
        self.experimentParams = experimentParams
        print >> self.f, "ENV : %s" % os.environ

        try:
            from ispybapi.ispyb.dbconnection import dbconnection
            from ispybapi.ispyb.core import core
            from ispybapi.ispyb.mxacquisition import mxacquisition
            self.dbconnection = dbconnection
            self.core = core
            self.mxacquisition = mxacquisition
        except:
            self.dbconnection = None
            self.core = None
            self.mxacquisition = None
            print('Could not import ispyb api')
            print >> self.f, "\n%s: >>> OPENING DB" % timeStamp()

        self._loadCursor(db)
        self._create_group()

    def _loadCursor(self, db):
        # db should be one of: 'prod', 'dev' or 'test'
        # let's try to connect to db and get a cursor
        self.cursor = None

        if self.dbconnection:
            # connect = getattr(self.dbconnection, 'connect' + db)
            self.cursor = self.dbconnection.connect(conf=db)
        else:
            print ('Couldnt connect to db')

    def _create_group(self):
        if self.mxacquisition:
            self.visit_id = self.core.retrieve_visit_id(self.cursor, self.experimentParams['visit'])
            params = self.mxacquisition.get_data_collection_group_params()
            params['parentid'] = self.visit_id
            self.group_id = self.mxacquisition.insert_data_collection_group(self.cursor, params.values())
        else:
            self.visit_id = 'no-visit'
            self.group_id = 'no-group'

    def get_data_collection_params(self):
        if self.mxacquisition:
            params = self.mxacquisition.get_data_collection_params()
            params['parentid'] = self.group_id
            params['visitid'] = self.visit_id
            return params
        else:
            return OrderedDict()

    def update_data_collection(self, params):
        if self.mxacquisition:
            return self.mxacquisition.insert_data_collection(self.cursor,
                                                             params.values())
        else:
            s = "params = {"
            for k, v in params.iteritems():
                s += "   %s:%s\n " % (k, v)
            s += "}"
            print >> self.f, s
            self.f.flush()
            return -1 # Fake id

    def disconnect(self):
        if self.dbconnection:
            self.dbconnection.disconnect()
        else:
            print >> self.f, "%s: <<< CLOSING DB\n" % timeStamp()
            self.f.flush()
            self.f.close()
