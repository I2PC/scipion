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

import sys
import os
from os.path import abspath, dirname
from collections import OrderedDict
import json
import subprocess
from datetime import datetime


SCRIPT = os.path.realpath(__file__)
INPUT_END = 'END'
OUTPUT_OK = 'OK'


def timeStamp(dt=None, format='%Y-%m-%d_%H%M%S'):
    if dt is None:
        dt = datetime.now()
    return dt.strftime(format)


class ISPyBProxy():
    """ This class could be used from Scipion environment. It would be used
    from the MonitorISPyB to spawn a new process that will communicate with
    the database. The communication to the process will be through input/output
    pipes. The data lines will be in json format. """
    def __init__(self, db, experimentParams):
        self._createISPyBProcess(db)
        self._sendDict(experimentParams)

    # <<<<<<< Updated upstream
    # class ISPyBdb():
    #     def __init__(self, db, parentid, visit, sampleid, detectorid):
    #         self.params =  {}
    #         self.params['parentid'] = parentid
    #         self.params['visit'] = visit
    #         self.params['sampleid'] = sampleid
    #         self.params['detectorid'] = detectorid
    #
    #     def put_movie(self, movie):
    #         movieFn = movie.getFileName()
    #         self.params['mfile'] = movieFn
    #
    #         # TODO: Use ispyb-api some day
    #         script = pw.join('em', 'packages', 'ispyb', 'em_put_movie.py')
    #         cmd = 'python ' + script
    #         for k, v in self.params.iteritems():
    #            cmd += ' --%s %s' % (k, v)
    #         pwutils.runJob(None,
    #                        'source /etc/profile.d/modules.sh; '
    #                        'module load python/ana; '
    #                        'module load ispyb-api/ana; ',
    #                        cmd)
    # =======
    def _createISPyBProcess(self, db):
        cmd = ''
        cmd += 'python %s %s' % (SCRIPT, db)

        #    gcmd = greenStr(cmd)
        print "** Running: '%s'" % cmd
        self.proc = subprocess.Popen(cmd, shell=True,
                                     stdin=subprocess.PIPE,
                                     stdout=subprocess.PIPE)

    def _sendLine(self, line):
        print "Sending line: ", line
        print >> self.proc.stdin, line
        self.proc.stdin.flush()
        r =  self.proc.stdout.readline()
        print "Recv: ", r
        return r

    def _sendDict(self, paramsDict):
        return self._sendLine(json.dumps(paramsDict))

    def sendMovieParams(self, movieParams):
        return self._sendDict(movieParams)

    def close(self):
        self._sendLine(INPUT_END)
        self.proc.kill()


class ISPyBProccess():
    def __init__(self, db):
        # self.db should be set after this
        self._loadDb(db, self._readExperimentParams())
        self.processInputMovies()

    def readLine(self):
        line = sys.stdin.readline()
        return line

    def readDict(self, line=None):
        """ Read a line with json format and convert to a dict. """
        try:
            line = line or self.readLine()
            lineDict = json.loads(line)
            if not isinstance(lineDict, dict):
                self.error("Expecting line with json dict. Got: " + line)
            return lineDict
        except:
            self.error("Could not parse line as json dict. Line: " + line)

    def writeLine(self, line):
        print >> sys.stdout, line
        sys.stdout.flush()

    def error(self, msg):
        self.writeLine("ERROR: %s" % msg)
        sys.exit(1)

    def ok(self):
        self.writeLine(OUTPUT_OK)

    def _loadDb(self, db, experimentParams):
        # db should be one of: 'prod', 'dev' or 'test'
        # let's try to connect to db and get a cursor
        try:
            self.db = ISPyBdb(db, experimentParams)
        except Exception as ex:
            self.error(str(ex))

    def _readExperimentParams(self):
        # We should read experiment params
        # from an input line in json format.
        # expected values: parentid, visit, sampleid, detectorid
        experimentParams = self.readDict()
        for key in ['parentid', 'visit', 'sampleid', 'detectorid']:
            if key not in experimentParams:
                self.error("Missing key '%s' in experiment params" % key)
        self.ok()
        return experimentParams

    def processInputMovies(self):
        """ Read lines with json dict (per movie) until we find the
        special line INPUT_END. """
        try:
            line = self.readLine()

            while line.strip() != INPUT_END:
                movieParams = self.readDict(line)
                params = self.db.get_data_collection_params()
                params.update(movieParams)

                self.db.update_data_collection(params)
                self.ok()

                line = self.readLine()

            self.db.disconnect()
            self.ok()

        except Exception as ex:
            self.db.disconnect()
            self.error(str(ex))


class ISPyBdb():
    """ This is a Facade to provide access to the ispyb_api to store movies.
    In case ispyb_api could not be loaded, we used a simple file log
    to still test things are working. In the normal case, it is supposed
    to be used from ISPyBProcess, which would be executed using the proper
    Python with the required modules for ispyb_api to work.
    """
    def __init__(self, db, experimentParams):
        try:
            from ispyb_api.dbconnection import dbconnection
            from ispyb_api.core import core
            from ispyb_api.mxacquisition import mxacquisition
            self.dbconnection = dbconnection
            self.core = core
            self.mxacquisition = mxacquisition
        except:
            self.dbconnection = None
            self.core = None
            self.mxacquisition = None
            self.f = open("/tmp/ispybdb.txt", 'a')
            print >> self.f, "\n%s: >>> OPENING DB" % timeStamp()

        self._loadCursor(db)

    def _loadCursor(self, db):
        # db should be one of: 'prod', 'dev' or 'test'
        # let's try to connect to db and get a cursor
        self.cursor = None

        if self.dbconnection:
            connect = getattr(self.dbconnection, 'connect_to_' + db)
            self.cursor = connect()

    def get_data_collection_params(self):
        if self.mxacquisition:
            return self.mxacquisition.get_data_collection_params()
        else:
            return OrderedDict()

    def update_data_collection(self, params):
        if self.mxacquisition:
            dc_id = self.mxacquisition.update_data_collection(self.cursor,
                                                              params.values())
        else:
            s = "params = {"
            for k, v in params.iteritems():
                s += "   %s:%s\n " % (k, v)
            s += "}"
            print >> self.f, s
            self.f.flush()

    def disconnect(self):
        if self.dbconnection:
            self.dbconnection.disconnect()
        else:
            print >> self.f, "%s: <<< CLOSING DB\n" % timeStamp()
            self.f.flush()
            self.f.close()


if __name__ == "__main__":
    db = sys.argv[1]
    ISPyBProccess(db)


