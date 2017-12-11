from pyworkflow.em import ProtProcessMovies
import json
import os
import pyworkflow.protocol.params as params


class ProtArchive(ProtProcessMovies):
    """
    Sends a message to the dials archiving queue
    """
    _label = 'Archive'

    def _defineParams(self, form):
        form.addSection(label='Params')
        group = form.addGroup('Message Queue')
        group.addParam('queueHost', params.StringParam,
                       label="Host", important=True,
                       help='The host name for a queue broker (e.g. activemq)')
        group.addParam('queuePort', params.IntParam,
                       label="Port", important=True, default=61613,
                       help='The stomp port for a queue broker (e.g. activemq)')
        group.addParam('queueName', params.StringParam,
                       label='Queue Name',
                       help='The name of the queue within the broker')
        group.addParam('dropFileDir', params.FolderParam,
                       label='Dropfile dir',
                       help='The location of the dropfile dir')

    def processMovieStep(self, movieDict, hasAlignment):
        SCIPION_HOME = os.environ['SCIPION_HOME']
        LOGFILE = os.path.join(SCIPION_HOME, 'pyworkflow', 'em', 'packages','archiving','recipe.json')
        with open(LOGFILE) as recipe:
            recipe_data = json.load(recipe)

        real_path = os.path.realpath(movieDict['_filename'])

        self._setValues(recipe_data, 'files', [real_path ])
        self._setValues(recipe_data, 'dropfile-dir', dropFileDir)

        json_data = json.dumps(recipe_data)
        self._sendMessage(json_data)

    def _setValues(self, data, key, value):
        for k in data:
            if k == key:
                data[k] = value
            elif type(data[k]) is dict:
                self._setValues(data[k], key, value)

    def _sendMessage(self, json_data):
        # type: (str) -> void
        try:
            try:
                from stomp import Connection
            except ImportError as ex:
                print "You need stomp to run the archiver, run 'scipion run pip stomp.py' first"
                raise ex
            self.connection.send(self.queue_name, json_data)
        except AttributeError:
            self.connection = Connection([(self.queueHost, self.queuePort)])
            self.connection.connect(wait=True)
            self._sendMessage(json_data)

    def createOutputStep(self):
        pass
